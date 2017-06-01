/*!*************************************************************************
 * \file field2d.cxx
 *
 * Class for 2D X-Y profiles
 *
 **************************************************************************
 * Copyright 2010 B.D.Dudson, S.Farley, M.V.Umansky, X.Q.Xu
 *
 * Contact: Ben Dudson, bd512@york.ac.uk
 * 
 * This file is part of BOUT++.
 *
 * BOUT++ is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * BOUT++ is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with BOUT++.  If not, see <http://www.gnu.org/licenses/>.
 * 
 **************************************************************************/

#include <boutcomm.hxx>
#include <bout/rvec.hxx>

#include <globals.hxx> // for mesh

#include <field2d.hxx>

#include <utils.hxx>

#include <boundary_op.hxx>
#include <boundary_factory.hxx>

#include <boutexception.hxx>
#include <msg_stack.hxx>

#include <cmath>
#include <output.hxx>

#include <bout/assert.hxx>

Field2D::Field2D(Mesh *msh) : fieldmesh(msh), deriv(nullptr) { 

  boundaryIsSet = false;

  if(fieldmesh) {
    nx = fieldmesh->LocalNx;
    ny = fieldmesh->LocalNy;
  }
  
#ifdef TRACK
  name = "<F2D>";
#endif
}

Field2D::Field2D(const Field2D& f) : fieldmesh(f.fieldmesh), // The mesh containing array sizes
                                     data(f.data), // This handles references to the data array
                                     deriv(nullptr) {
  if(fieldmesh) {
    nx = fieldmesh->LocalNx;
    ny = fieldmesh->LocalNy;
  }
  
  boundaryIsSet = false;
  *this = f;
}

Field2D::Field2D(BoutReal val) : fieldmesh(nullptr), deriv(nullptr) {
  boundaryIsSet = false;
  
  fieldmesh = mesh;
  nx = fieldmesh->LocalNx;
  ny = fieldmesh->LocalNy;
  
  *this = val;
}

Field2D::~Field2D() {
  if(deriv)
    delete deriv;
}

void Field2D::allocate() {
  if(data.empty()) {
    if(!fieldmesh) {
      /// If no mesh, use the global
      fieldmesh = mesh;
      nx = fieldmesh->LocalNx;
      ny = fieldmesh->LocalNy;
    }
    data = Array<BoutReal>(nx*ny);
  }else
    data.ensureUnique();
}

Field2D* Field2D::timeDeriv() {
  if(deriv == nullptr)
    deriv = new Field2D();
  return deriv;
}

///////////// OPERATORS ////////////////

Field2D & Field2D::operator=(const Field2D &rhs) {
  // Check for self-assignment
  if(this == &rhs)
    return(*this); // skip this assignment

  TRACE("Field2D: Assignment from Field2D");

  checkData(rhs);
  
#ifdef TRACK
  name = rhs.name;
#endif

  // Copy the data and data sizes
  fieldmesh = rhs.fieldmesh;
  nx = rhs.nx; ny = rhs.ny; 

  // Copy reference to data
  data = rhs.data;

  return *this;
}

Field2D & Field2D::operator=(const BoutReal rhs) {
  TRACE("Field2D = BoutReal");
#ifdef TRACK
  name = "<r2D>";
#endif
  
  allocate();
  for(auto i : (*this))
    (*this)[i] = rhs;
  
  return *this;
}

////////////// Indexing ///////////////////

const DataIterator Field2D::iterator() const {
  return DataIterator(0, mesh->LocalNx-1, 
                      0, mesh->LocalNy-1,
                      0, 0);
}

const DataIterator Field2D::begin() const {
  /*return DataIterator(0, 0, mesh->LocalNx-1,
                      0, 0, mesh->LocalNy-1,
                      0, 0, 0);*/
  return Field2D::iterator();
}

const DataIterator Field2D::end() const {
  return DataIterator(0, mesh->LocalNx-1, 
                      0, mesh->LocalNy-1,
                      0, 0, DI_GET_END);
}

const IndexRange Field2D::region(REGION rgn) const {
  switch(rgn) {
  case RGN_ALL: {
    return IndexRange{0, nx-1,
        0, ny-1,
        0, 0};
    break;
  }
  case RGN_NOBNDRY: {
    return IndexRange{fieldmesh->xstart, fieldmesh->xend,
        fieldmesh->ystart, fieldmesh->yend,
        0, 0};
    break;
  }
  case RGN_NOX: {
    return IndexRange{fieldmesh->xstart, fieldmesh->xend,
        0, ny-1,
        0, 0};
    break;
  }
  case RGN_NOY: {
    return IndexRange{0, nx-1,
        fieldmesh->ystart, fieldmesh->yend,
        0, 0};
    break;
  }
  default: {
    throw BoutException("Field2D::region() : Requested region not implemented");
  }
  };
}

///////// Operators

#define F2D_UPDATE_FIELD(op,bop,ftype)                       \
  Field2D & Field2D::operator op(const ftype &rhs) {         \
    msg_stack.push("Field2D: %s %s", #op, #ftype);           \
    checkData(rhs) ;                                         \
    checkData(*this);                                        \
    if(data.unique()) {                                      \
      /* This is the only reference to this data */          \
      for(auto i : (*this))                                  \
        (*this)[i] op rhs[i];                                \
    }else {                                                  \
      /* Shared data */                                      \
      (*this) = (*this) bop rhs;                             \
    }                                                        \
    msg_stack.pop();                                         \
    return *this;                                            \
  }

F2D_UPDATE_FIELD(+=, +, Field2D); // operator+=(const Field2D &rhs)
F2D_UPDATE_FIELD(-=, -, Field2D); // operator-=(const Field2D &rhs)
F2D_UPDATE_FIELD(*=, *, Field2D); // operator*=(const Field2D &rhs)
F2D_UPDATE_FIELD(/=, /, Field2D); // operator/=(const Field2D &rhs)

#define F2D_UPDATE_REAL(op,bop)                              \
  Field2D & Field2D::operator op(const BoutReal rhs) {       \
    msg_stack.push("Field2D: %s Field2D", #op);              \
    if(!finite(rhs))                                         \
      throw BoutException("Field2D: %s operator passed non-finite BoutReal number", #op); \
    checkData(*this);                                        \
                                                             \
    if(data.unique()) {                                      \
      /* This is the only reference to this data */          \
      for(auto i : (*this))                                  \
        (*this)[i] op rhs;                                   \
    }else {                                                  \
      /* Need to put result in a new block */                \
      (*this) = (*this) bop rhs;                             \
    }                                                        \
    msg_stack.pop();                                         \
    return *this;                                            \
  }

F2D_UPDATE_REAL(+=,+);    // operator+= BoutReal
F2D_UPDATE_REAL(-=,-);    // operator-= BoutReal
F2D_UPDATE_REAL(*=,*);    // operator*= BoutReal
F2D_UPDATE_REAL(/=,/);    // operator/= BoutReal

////////////////////// STENCILS //////////////////////////

void Field2D::getXArray(int y, int UNUSED(z), rvec &xv) const {
  ASSERT0(isAllocated());

  xv.resize(nx);
  
  for(int x=0;x<nx;x++)
    xv[x] = operator()(x,y);
}

void Field2D::getYArray(int x, int UNUSED(z), rvec &yv) const {
  ASSERT0(isAllocated());

  yv.resize(ny);
  
  for(int y=0;y<ny;y++)
    yv[y] = operator()(x,y);
}

void Field2D::getZArray(int x, int y, rvec &zv) const {
  ASSERT0(isAllocated());

  zv.resize(mesh->LocalNz);
  
  for(int z=0;z<mesh->LocalNz;z++)
    zv[z] = operator()(x,y);
}

void Field2D::setXArray(int y, int UNUSED(z), const rvec &xv) {
  allocate();

  ASSERT0(xv.capacity() == (unsigned int) nx);

  for(int x=0;x<nx;x++)
    operator()(x,y) = xv[x];
}

void Field2D::setYArray(int x, int UNUSED(z), const rvec &yv) {
  allocate();

  ASSERT0(yv.capacity() == (unsigned int) mesh->LocalNy);

  for(int y=0;y<mesh->LocalNy;y++)
    operator()(x,y) = yv[y];
}

void Field2D::setXStencil(stencil &fval, const bindex &bx, CELL_LOC UNUSED(loc)) const {
  fval.mm = operator()(bx.jx2m,bx.jy);
  fval.m  = operator()(bx.jxm,bx.jy);
  fval.c  = operator()(bx.jx,bx.jy);
  fval.p  = operator()(bx.jxp,bx.jy);
  fval.pp = operator()(bx.jx2p,bx.jy);
}

void Field2D::setXStencil(forward_stencil &fval, const bindex &bx, CELL_LOC UNUSED(loc)) const {
  fval.m  = operator()(bx.jxm,bx.jy);
  fval.c  = operator()(bx.jx,bx.jy);
  fval.p  = operator()(bx.jxp,bx.jy);
  fval.p2 = operator()(bx.jx2p,bx.jy);
  fval.p3 = operator()(bx.jx+3,bx.jy);
  fval.p4 = operator()(bx.jx+4,bx.jy);
}

void Field2D::setXStencil(backward_stencil &fval, const bindex &bx, CELL_LOC UNUSED(loc)) const {
  fval.m4 = operator()(bx.jx-4,bx.jy);
  fval.m3 = operator()(bx.jx-3,bx.jy);
  fval.m2 = operator()(bx.jx2m,bx.jy);
  fval.m  = operator()(bx.jxm,bx.jy);
  fval.c  = operator()(bx.jx,bx.jy);
  fval.p  = operator()(bx.jxp,bx.jy);
}

void Field2D::setYStencil(stencil &fval, const bindex &bx, CELL_LOC UNUSED(loc)) const {
  fval.mm = operator()(bx.jx,bx.jy2m);
  fval.m  = operator()(bx.jx,bx.jym);
  fval.c  = operator()(bx.jx,bx.jy);
  fval.p  = operator()(bx.jx,bx.jyp);
  fval.pp = operator()(bx.jx,bx.jy2p);
}

void Field2D::setYStencil(forward_stencil &fval, const bindex &bx, CELL_LOC UNUSED(loc)) const {
  fval.m  = operator()(bx.jx,bx.jym);
  fval.c  = operator()(bx.jx,bx.jy);
  fval.p  = operator()(bx.jx,bx.jyp);
  fval.p2 = operator()(bx.jx,bx.jy2p);
  fval.p3 = operator()(bx.jx,bx.jy+3);
  fval.p4 = operator()(bx.jx,bx.jy+4);
}

void Field2D::setYStencil(backward_stencil &fval, const bindex &bx, CELL_LOC UNUSED(loc)) const {
  fval.m4 = operator()(bx.jx,bx.jy-4);
  fval.m3 = operator()(bx.jx,bx.jy-3);
  fval.m2 = operator()(bx.jx,bx.jy2m);
  fval.m  = operator()(bx.jx,bx.jym);
  fval.c  = operator()(bx.jx,bx.jy);
  fval.p  = operator()(bx.jx,bx.jyp);
}

void Field2D::setZStencil(stencil &fval, const bindex &bx, CELL_LOC UNUSED(loc)) const {
  fval = operator()(bx.jx,bx.jy);
}

///////////////////// FieldData VIRTUAL FUNCTIONS //////////

int Field2D::getData(int x, int y, int z, void *vptr) const {
  ASSERT0(isAllocated()); // Check data set
  
#if CHECK > 2
  // check ranges
  if((x < 0) || (x >= nx) || (y < 0) || (y >= ny)) {
    throw BoutException("Field2D: getData (%d,%d,%d) out of bounds\n", x, y, z);
  }
#endif
  BoutReal *ptr = (BoutReal*) vptr;
  *ptr = operator()(x,y);
  
  return sizeof(BoutReal);
}

int Field2D::getData(int x, int y, int z, BoutReal *rptr) const {
  ASSERT0(isAllocated()); // Check data set
  
#if CHECK > 2
  // check ranges
  if((x < 0) || (x >= nx) || (y < 0) || (y >= ny)) {
    throw BoutException("Field2D: getData (%d,%d,%d) out of bounds\n", x, y, z);
  }
#endif

  *rptr = operator()(x,y);
  return 1;
}

int Field2D::setData(int x, int y, int z, void *vptr) {
  allocate();
  
#if CHECK > 2
  // check ranges
  if((x < 0) || (x >= nx) || (y < 0) || (y >= ny)) {
    throw BoutException("Field2D: setData (%d,%d,%d) out of bounds\n", x, y, z);
  }
#endif
  BoutReal *ptr = (BoutReal*) vptr;
  operator()(x,y) = *ptr;
  
  return sizeof(BoutReal);
}

int Field2D::setData(int x, int y, int UNUSED(z), BoutReal *rptr) {
  allocate();
#if CHECK > 2
  // check ranges
  if((x < 0) || (x >= nx) || (y < 0) || (y >= ny) ) {
    throw BoutException("Field2D: setData (%d,%d) out of bounds\n", x, y);
  }
#endif

  operator()(x,y) = *rptr;
  return 1;
}

///////////////////// BOUNDARY CONDITIONS //////////////////

void Field2D::applyBoundary(bool init) {
#ifdef CHECK
  if (init) {
    msg_stack.push("Field2D::applyBoundary()");
    if(!boundaryIsSet)
      output << "WARNING: Call to Field2D::applyBoundary(), but no boundary set" << endl;
  }
#endif
  for(const auto& bndry : bndry_op)
    if ( !bndry->apply_to_ddt || init) // Always apply to the values when initialising fields, otherwise apply only if wanted
      bndry->apply(*this);
  msg_stack.pop();
}

void Field2D::applyBoundary(const string &condition) {
  TRACE("Field2D::applyBoundary(condition)");
  
#ifdef CHECK
  if(!isAllocated())
    output << "WARNING: Empty data in Field2D::applyBoundary(condition)" << endl;
#endif
  
  if(!isAllocated())
    return;

  /// Get the boundary factory (singleton)
  BoundaryFactory *bfact = BoundaryFactory::getInstance();
  
  /// Loop over the mesh boundary regions
  for(const auto& reg : mesh->getBoundaries()) {
    BoundaryOp* op = static_cast<BoundaryOp*>(bfact->create(condition, reg));
    op->apply(*this);
    delete op;
  }
  
  // Set the corners to zero
  for(int jx=0;jx<mesh->xstart;jx++) {
    for(int jy=0;jy<mesh->ystart;jy++) {
      operator()(jx,jy) = 0.;
    }
    for(int jy=mesh->yend+1;jy<mesh->LocalNy;jy++) {
      operator()(jx,jy) = 0.;
    }
  }
  for(int jx=mesh->xend+1;jx<mesh->LocalNx;jx++) {
    for(int jy=0;jy<mesh->ystart;jy++) {
      operator()(jx,jy) = 0.;
    }
    for(int jy=mesh->yend+1;jy<mesh->LocalNy;jy++) {
      operator()(jx,jy) = 0.;
    }
  }
}

void Field2D::applyBoundary(const string &region, const string &condition) {
  if(!isAllocated())
    return;

  /// Get the boundary factory (singleton)
  BoundaryFactory *bfact = BoundaryFactory::getInstance();
  
  /// Loop over the mesh boundary regions
  for(const auto& reg : mesh->getBoundaries()) {
    if(reg->label.compare(region) == 0) {
      BoundaryOp* op = static_cast<BoundaryOp*>(bfact->create(condition, reg));
      op->apply(*this);
      delete op;
      break;
    }
  }
  
  // Set the corners to zero
  for(int jx=0;jx<mesh->xstart;jx++) {
    for(int jy=0;jy<mesh->ystart;jy++) {
      operator()(jx,jy) = 0.;
    }
    for(int jy=mesh->yend+1;jy<mesh->LocalNy;jy++) {
      operator()(jx,jy) = 0.;
    }
  }
  for(int jx=mesh->xend+1;jx<mesh->LocalNx;jx++) {
    for(int jy=0;jy<mesh->ystart;jy++) {
      operator()(jx,jy) = 0.;
    }
    for(int jy=mesh->yend+1;jy<mesh->LocalNy;jy++) {
      operator()(jx,jy) = 0.;
    }
  }
}

void Field2D::applyTDerivBoundary() {
  for(const auto& bndry : bndry_op)
    bndry->apply_ddt(*this);
}

void Field2D::setBoundaryTo(const Field2D &f2d) {
  TRACE("Field2D::setBoundary(const Field2D&)");
  allocate(); // Make sure data allocated
  
  ASSERT0(f2d.isAllocated());

  /// Loop over boundary regions
  for(const auto& reg : mesh->getBoundaries()) {
    /// Loop within each region
    for(reg->first(); !reg->isDone(); reg->next()) {
      // Get value half-way between cells
      BoutReal val = 0.5*(f2d(reg->x,reg->y) + f2d(reg->x-reg->bx, reg->y-reg->by));
      // Set to this value
      (*this)(reg->x,reg->y) = 2.*val - (*this)(reg->x-reg->bx, reg->y-reg->by);
    }
  }
}

////////////// NON-MEMBER OVERLOADED OPERATORS //////////////

#define F2D_OP_F2D(op)                                     \
  const Field2D operator op(const Field2D &lhs, const Field2D &rhs) { \
    Field2D result;                                                 \
    result.allocate();                                              \
    for(auto i : result)                                            \
      result[i] = lhs[i] op rhs[i];                                 \
    return result;                                                  \
  }

F2D_OP_F2D(+);  // Field2D + Field2D
F2D_OP_F2D(-);  // Field2D - Field2D
F2D_OP_F2D(*);  // Field2D * Field2D
F2D_OP_F2D(/);  // Field2D / Field2D

#define F2D_OP_F3D(op)                                     \
  const Field3D operator op(const Field2D &lhs, const Field3D &rhs) { \
    Field3D result;                                                 \
    result.allocate();                                              \
    for(auto i : result)                                            \
      result[i] = lhs[i] op rhs[i];                                 \
    return result;                                                  \
  }

F2D_OP_F3D(+);  // Field2D + Field3D
F2D_OP_F3D(-);  // Field2D - Field3D
F2D_OP_F3D(*);  // Field2D * Field3D
F2D_OP_F3D(/);  // Field2D / Field3D

#define F2D_OP_REAL(op)                                     \
  const Field2D operator op(const Field2D &lhs, BoutReal rhs) {     \
    Field2D result;                                                 \
    result.allocate();                                              \
    for(auto i : result)                                            \
      result[i] = lhs[i] op rhs;                                    \
    return result;                                                  \
  }

F2D_OP_REAL(+);  // Field2D + BoutReal
F2D_OP_REAL(-);  // Field2D - BoutReal
F2D_OP_REAL(*);  // Field2D * BoutReal
F2D_OP_REAL(/);  // Field2D / BoutReal

#define REAL_OP_F2D(op)                                     \
  const Field2D operator op(BoutReal lhs, const Field2D &rhs) {     \
    Field2D result;                                                 \
    result.allocate();                                              \
    for(auto i : result)                                            \
      result[i] = lhs op rhs[i];                                    \
    return result;                                                  \
  }

REAL_OP_F2D(+);  // BoutReal + Field2D
REAL_OP_F2D(-);  // BoutReal - Field2D
REAL_OP_F2D(*);  // BoutReal * Field2D
REAL_OP_F2D(/);  // BoutReal / Field2D

// Unary minus
const Field2D operator-(const Field2D &f) {
  return -1.0*f;
}

//////////////// NON-MEMBER FUNCTIONS //////////////////

BoutReal min(const Field2D &f, bool allpe) {
  TRACE("min(Field2D)");
  
  ASSERT2(f.isAllocated());

  BoutReal result = f(mesh->xstart,mesh->ystart);

  for(auto i : f.region(RGN_NOBNDRY))
    if(f[i] < result)
      result = f[i];
  
  if(allpe) {
    // MPI reduce
    BoutReal localresult = result;
    MPI_Allreduce(&localresult, &result, 1, MPI_DOUBLE, MPI_MIN, BoutComm::get());
  }
  
  return result;
}

BoutReal max(const Field2D &f, bool allpe) {
  TRACE("max(Field2D)");
  
  ASSERT2(f.isAllocated());

  BoutReal result = f(mesh->xstart,mesh->ystart);

  for(auto i : f.region(RGN_NOBNDRY))
    if(f[i] > result)
      result = f[i];
  
  if(allpe) {
    // MPI reduce
    BoutReal localresult = result;
    MPI_Allreduce(&localresult, &result, 1, MPI_DOUBLE, MPI_MAX, BoutComm::get());
  }
  
  return result;
}

bool finite(const Field2D &f) {
  TRACE("finite(Field2D)");
  ASSERT0(f.isAllocated());

  for(auto i : f)
    if(!::finite(f[i]))
      return false;
  
  return true;
}

/////////////////////////////////////////////////
// functions

/*!
 * This macro takes a function \p func, which is
 * assumed to operate on a single BoutReal and return
 * a single BoutReal, and wraps it up into a function
 * of a Field2D called \p name.
 *
 * @param name  The name of the function to define
 * @param func  The function to apply to each value
 *
 * If CHECK >= 1, checks if the Field2D is allocated
 *
 * Loops over the entire domain, applies function,
 * and if CHECK >= 3 then checks result for non-finite numbers
 *
 */
#define F2D_FUNC(name, func)                               \
  const Field2D name(const Field2D &f) {                   \
    msg_stack.push(#name "(Field2D)");                     \
    /* Check if the input is allocated */                  \
    ASSERT1(f.isAllocated());                              \
    /* Define and allocate the output result */            \
    Field2D result;                                        \
    result.allocate();                                     \
    /* Loop over domain */                                 \
    for(auto d : result) {                                 \
      result[d] = func(f[d]);                              \
      /* If checking is set to 3 or higher, test result */ \
      ASSERT3(finite(result[d]));                          \
    }                                                      \
    msg_stack.pop();                                       \
    return result;                                         \
  }

F2D_FUNC(abs, ::fabs);

F2D_FUNC(sqrt, ::sqrt);

F2D_FUNC(exp, ::exp);
F2D_FUNC(log, ::log);

F2D_FUNC(sin, ::sin);
F2D_FUNC(cos, ::cos);
F2D_FUNC(tan, ::tan);

F2D_FUNC(sinh, ::sinh);
F2D_FUNC(cosh, ::cosh);
F2D_FUNC(tanh, ::tanh);

const Field2D copy(const Field2D &f) {
  Field2D result = f;
  result.allocate();
  return result;
}

const Field2D floor(const Field2D &var, BoutReal f) {
  Field2D result = copy(var);

  for(auto d : result)
    if(result[d] < f)
      result[d] = f;
  
  return result;
}

Field2D pow(const Field2D &lhs, const Field2D &rhs) {
  TRACE("pow(Field2D, Field2D)");
  // Check if the inputs are allocated
  ASSERT1(lhs.isAllocated());
  ASSERT1(rhs.isAllocated());

  // Define and allocate the output result
  Field2D result;
  result.allocate();

  // Loop over domain
  for(auto i: result) {
    result[i] = ::pow(lhs[i], rhs[i]);
    ASSERT3(finite(result[i]));
  }
  return result;
}

Field2D pow(const Field2D &lhs, BoutReal rhs) {
  TRACE("pow(Field2D, BoutReal)");
  // Check if the inputs are allocated
  ASSERT1(lhs.isAllocated());

  // Define and allocate the output result
  Field2D result;
  result.allocate();

  // Loop over domain
  for(auto i: result) {
    result[i] = ::pow(lhs[i], rhs);
    ASSERT3(finite(result[i]));
  }
  return result;
}

Field2D pow(BoutReal lhs, const Field2D &rhs) {
  TRACE("pow(lhs, Field2D)");
  // Check if the inputs are allocated
  ASSERT1(rhs.isAllocated());

  // Define and allocate the output result
  Field2D result;
  result.allocate();

  // Loop over domain
  for(auto i: result) {
    result[i] = ::pow(lhs, rhs[i]);
    ASSERT3(finite(result[i]));
  }
  return result;
}

#ifdef CHECK
/// Check if the data is valid
void checkData(const Field2D &f) {
  if(!f.isAllocated()) {
    throw BoutException("Field2D: Operation on empty data\n");
  }
  
#if CHECK > 2
  // Do full checks
  for(auto i : f.region(RGN_NOBNDRY)){
    if(!::finite(f[i])) {
      throw BoutException("Field2D: Operation on non-finite data at [%d][%d]\n", i.x, i.y);
    }
  }
#endif
}
#endif
