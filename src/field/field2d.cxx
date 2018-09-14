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

Field2D::Field2D(Mesh *localmesh) : Field(localmesh), deriv(nullptr) {

  boundaryIsSet = false;

  if(fieldmesh) {
    nx = fieldmesh->LocalNx;
    ny = fieldmesh->LocalNy;
  }
#if CHECK > 0
  else {
    nx=-1;
    ny=-1;
  }
#endif

#ifdef TRACK
  name = "<F2D>";
#endif
}

Field2D::Field2D(const Field2D& f) : Field(f.fieldmesh), // The mesh containing array sizes
                                     data(f.data), // This handles references to the data array
                                     deriv(nullptr) {
  TRACE("Field2D(Field2D&)");

#ifdef TRACK
  name = f.name;
#endif

#if CHECK > 2
  checkData(f);
#endif

  setLocation(f.location);

  if(fieldmesh) {
    nx = fieldmesh->LocalNx;
    ny = fieldmesh->LocalNy;
  }
#if CHECK > 0
  else {
    nx=-1;
    ny=-1;
  }
#endif

  boundaryIsSet = false;
}

Field2D::Field2D(BoutReal val, Mesh *localmesh) : Field(localmesh), deriv(nullptr) {
  boundaryIsSet = false;

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
#if CHECK > 2
    invalidateGuards(*this);
#endif
  }else
    data.ensureUnique();
}

Field2D* Field2D::timeDeriv() {
  if(deriv == nullptr)
    deriv = new Field2D(fieldmesh);
  return deriv;
}

////////////// Indexing ///////////////////

const DataIterator Field2D::iterator() const {
  return DataIterator(0, nx-1,
                      0, ny-1,
                      0, 0);
}

const DataIterator Field2D::begin() const {
  return Field2D::iterator();
}

const DataIterator Field2D::end() const {
  return DataIterator(0, nx-1,
                      0, ny-1,
                      0, 0, DI_GET_END);
}

const IndexRange Field2D::region(REGION rgn) const {
  switch(rgn) {
  case RGN_ALL: {
    return IndexRange{0, nx-1,
        0, ny-1,
        0, 0};
  }
  case RGN_NOBNDRY: {
    return IndexRange{fieldmesh->xstart, fieldmesh->xend,
        fieldmesh->ystart, fieldmesh->yend,
        0, 0};
  }
  case RGN_NOX: {
    return IndexRange{fieldmesh->xstart, fieldmesh->xend,
        0, ny-1,
        0, 0};
  }
  case RGN_NOY: {
    return IndexRange{0, nx-1,
        fieldmesh->ystart, fieldmesh->yend,
        0, 0};
  }
  default: {
    throw BoutException("Field2D::region() : Requested region not implemented");
  }
  };
}

void Field2D::setLocation(CELL_LOC new_location) {
  if (getMesh()->StaggerGrids) {
    if (new_location == CELL_VSHIFT) {
      throw BoutException(
          "Field2D: CELL_VSHIFT cell location only makes sense for vectors");
    }
    if (new_location == CELL_DEFAULT) {
      new_location = CELL_CENTRE;
    }
    location = new_location;
    if (new_location != location)
      fieldCoordinates = getMesh()->coordinates(location);
  } else {
#if CHECK > 0
    if (new_location != CELL_CENTRE && new_location != CELL_DEFAULT) {
      throw BoutException("Field2D: Trying to set off-centre location on "
                          "non-staggered grid\n"
                          "         Did you mean to enable staggerGrids?");
    }
#endif
    location = CELL_CENTRE;
  }

  /// Would like to do something like the following but can lead to problems
  /// at least in unit testing.
  /// if (fieldCoordinates == nullptr) fieldCoordinates = getMesh()->coordinates(location);  
}

CELL_LOC Field2D::getLocation() const {
  return location;
}

// Not in header because we need to access fieldmesh
BoutReal& Field2D::operator[](const Ind3D &d) {
  return operator[](fieldmesh->map3Dto2D(d));
}

const BoutReal& Field2D::operator[](const Ind3D &d) const {
  return operator[](fieldmesh->map3Dto2D(d));
}

///////////// OPERATORS ////////////////

Field2D &Field2D::operator=(const Field2D &rhs) {
  // Check for self-assignment
  if (this == &rhs)
    return (*this); // skip this assignment

  TRACE("Field2D: Assignment from Field2D");

  checkData(rhs);

#ifdef TRACK
  name = rhs.name;
#endif

  // Copy the data and data sizes
  fieldmesh = rhs.fieldmesh;
  nx = rhs.nx;
  ny = rhs.ny;

  // Copy reference to data
  data = rhs.data;

  // Copy location
  setLocation(rhs.location);

  return *this;
}

Field2D &Field2D::operator=(const BoutReal rhs) {
#ifdef TRACK
  name = "<r2D>";
#endif

  TRACE("Field2D = BoutReal");
  allocate();

#if CHECK > 0
  if (!finite(rhs))
    throw BoutException("Field2D: Assignment from non-finite BoutReal\n");
#endif
  for (const auto &i : (*this))
    (*this)[i] = rhs;

  return *this;
}

///////////////////// BOUNDARY CONDITIONS //////////////////

void Field2D::applyBoundary(bool init) {
  TRACE("Field2D::applyBoundary()");

#if CHECK > 0
  if (init) {

    if(!boundaryIsSet)
      output_warn << "WARNING: Call to Field2D::applyBoundary(), but no boundary set" << endl;
  }
#endif

  checkData(*this);

  for(const auto& bndry : bndry_op)
    if ( !bndry->apply_to_ddt || init) // Always apply to the values when initialising fields, otherwise apply only if wanted
      bndry->apply(*this);
}

void Field2D::applyBoundary(const string &condition) {
  TRACE("Field2D::applyBoundary(condition)");

  checkData(*this);

  /// Get the boundary factory (singleton)
  BoundaryFactory *bfact = BoundaryFactory::getInstance();

  /// Loop over the mesh boundary regions
  for(const auto& reg : fieldmesh->getBoundaries()) {
    BoundaryOp* op = static_cast<BoundaryOp*>(bfact->create(condition, reg));
    op->apply(*this);
    delete op;
  }

  // Set the corners to zero
  for(int jx=0;jx<fieldmesh->xstart;jx++) {
    for(int jy=0;jy<fieldmesh->ystart;jy++) {
      operator()(jx,jy) = 0.;
    }
    for(int jy=fieldmesh->yend+1;jy<fieldmesh->LocalNy;jy++) {
      operator()(jx,jy) = 0.;
    }
  }
  for(int jx=fieldmesh->xend+1;jx<fieldmesh->LocalNx;jx++) {
    for(int jy=0;jy<fieldmesh->ystart;jy++) {
      operator()(jx,jy) = 0.;
    }
    for(int jy=fieldmesh->yend+1;jy<fieldmesh->LocalNy;jy++) {
      operator()(jx,jy) = 0.;
    }
  }
}

void Field2D::applyBoundary(const string &region, const string &condition) {
  checkData(*this);

  /// Get the boundary factory (singleton)
  BoundaryFactory *bfact = BoundaryFactory::getInstance();

  /// Loop over the mesh boundary regions
  for(const auto& reg : fieldmesh->getBoundaries()) {
    if(reg->label.compare(region) == 0) {
      BoundaryOp* op = static_cast<BoundaryOp*>(bfact->create(condition, reg));
      op->apply(*this);
      delete op;
      break;
    }
  }

  // Set the corners to zero
  for(int jx=0;jx<fieldmesh->xstart;jx++) {
    for(int jy=0;jy<fieldmesh->ystart;jy++) {
      operator()(jx,jy) = 0.;
    }
    for(int jy=fieldmesh->yend+1;jy<fieldmesh->LocalNy;jy++) {
      operator()(jx,jy) = 0.;
    }
  }
  for(int jx=fieldmesh->xend+1;jx<fieldmesh->LocalNx;jx++) {
    for(int jy=0;jy<fieldmesh->ystart;jy++) {
      operator()(jx,jy) = 0.;
    }
    for(int jy=fieldmesh->yend+1;jy<fieldmesh->LocalNy;jy++) {
      operator()(jx,jy) = 0.;
    }
  }
}

void Field2D::applyTDerivBoundary() {
  TRACE("Field2D::applyTDerivBoundary()");

  checkData(*this);
  ASSERT1(deriv != nullptr);
  checkData(*deriv);

  for(const auto& bndry : bndry_op)
    bndry->apply_ddt(*this);
}

void Field2D::setBoundaryTo(const Field2D &f2d) {
  TRACE("Field2D::setBoundary(const Field2D&)");

  checkData(f2d);

  allocate(); // Make sure data allocated

  /// Loop over boundary regions
  for(const auto& reg : fieldmesh->getBoundaries()) {
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

// Unary minus
Field2D operator-(const Field2D &f) { return -1.0 * f; }

//////////////// NON-MEMBER FUNCTIONS //////////////////

BoutReal min(const Field2D &f, bool allpe, REGION rgn) {
  TRACE("Field2D::Min() %s",allpe? "over all PEs" : "");

  checkData(f);

  BoutReal result = f[f.region(rgn).begin()];

  for(const auto& i : f.region(rgn))
    if(f[i] < result)
      result = f[i];

  if(allpe) {
    // MPI reduce
    BoutReal localresult = result;
    MPI_Allreduce(&localresult, &result, 1, MPI_DOUBLE, MPI_MIN, BoutComm::get());
  }

  return result;
}

BoutReal max(const Field2D &f, bool allpe,REGION rgn) {
  TRACE("Field2D::Max() %s",allpe? "over all PEs" : "");

  checkData(f);

  BoutReal result = f[f.region(rgn).begin()];

  for(const auto& i : f.region(rgn))
    if(f[i] > result)
      result = f[i];

  if(allpe) {
    // MPI reduce
    BoutReal localresult = result;
    MPI_Allreduce(&localresult, &result, 1, MPI_DOUBLE, MPI_MAX, BoutComm::get());
  }

  return result;
}

bool finite(const Field2D &f, REGION rgn) {
  TRACE("finite(Field2D)");

  if (!f.isAllocated()) {
    return false;
  }

  for (const auto &i : f.region(rgn)) {
    if (!::finite(f[i])) {
      return false;
    }
  }

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
 * and uses checkData() to, if CHECK >= 3, check
 * result for non-finite numbers
 *
 */
#define F2D_FUNC(name, func)                                                             \
  const Field2D name(const Field2D &f, REGION rgn) {                                     \
    TRACE(#name "(Field2D)");                                                            \
    /* Check if the input is allocated */                                                \
    checkData(f);                                                                        \
    /* Define and allocate the output result */                                          \
    Field2D result(f.getMesh());                                                         \
    result.allocate();                                                                   \
    /* Loop over domain */                                                               \
    for (const auto &d : result.region(rgn)) {                                           \
      result[d] = func(f[d]);                                                            \
    }                                                                                    \
    result.setLocation(f.getLocation());                                                 \
    checkData(result);                                                                   \
    return result;                                                                       \
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

const Field2D floor(const Field2D &var, BoutReal f, REGION rgn) {
  checkData(var);

  Field2D result = copy(var);

  for(const auto& d : result.region(rgn))
    if(result[d] < f)
      result[d] = f;

  return result;
}

Field2D pow(const Field2D &lhs, const Field2D &rhs, REGION rgn) {
  TRACE("pow(Field2D, Field2D)");
  // Check if the inputs are allocated
  checkData(lhs);
  checkData(rhs);
  ASSERT1(lhs.getLocation() == rhs.getLocation());

  // Define and allocate the output result
  ASSERT1(lhs.getMesh() == rhs.getMesh());
  Field2D result(lhs.getMesh());
  result.allocate();

  // Loop over domain
  for (const auto &i : result.region(rgn)) {
    result[i] = ::pow(lhs[i], rhs[i]);
  }

  result.setLocation(lhs.getLocation());

  checkData(result);
  return result;
}

Field2D pow(const Field2D &lhs, BoutReal rhs, REGION rgn) {
  TRACE("pow(Field2D, BoutReal)");
  // Check if the inputs are allocated
  checkData(lhs);
  checkData(rhs);

  // Define and allocate the output result
  Field2D result(lhs.getMesh());
  result.allocate();

  // Loop over domain
  for(const auto& i: result.region(rgn)) {
    result[i] = ::pow(lhs[i], rhs);
  }

  result.setLocation(lhs.getLocation());

  checkData(result);
  return result;
}

Field2D pow(BoutReal lhs, const Field2D &rhs, REGION rgn) {
  TRACE("pow(lhs, Field2D)");
  // Check if the inputs are allocated
  checkData(lhs);
  checkData(rhs);

  // Define and allocate the output result
  Field2D result(rhs.getMesh());
  result.allocate();

  // Loop over domain
  for(const auto& i: result.region(rgn)) {
    result[i] = ::pow(lhs, rhs[i]);
  }

  result.setLocation(rhs.getLocation());

  checkData(result);
  return result;
}

#if CHECK > 0
/// Check if the data is valid
void checkData(const Field2D &f, REGION region) {
  if(!f.isAllocated()) {
    throw BoutException("Field2D: Operation on empty data\n");
  }

#if CHECK > 2
  // Do full checks
  for(const auto& i : f.region(region)){
    if(!::finite(f[i])) {
      throw BoutException("Field2D: Operation on non-finite data at [%d][%d]\n", i.x, i.y);
    }
  }
#endif
}
#endif

void invalidateGuards(Field2D &var){
#if CHECK > 2
  Mesh *localmesh = var.getMesh();

  // Inner x -- all y and all z
  for (int ix = 0; ix < localmesh->xstart; ix++) {
    for (int iy = 0; iy < localmesh->LocalNy; iy++) {
      var(ix, iy) = std::nan("");
    }
  }

  // Outer x -- all y and all z
  for (int ix = localmesh->xend + 1; ix < localmesh->LocalNx; ix++) {
    for (int iy = 0; iy < localmesh->LocalNy; iy++) {
      var(ix, iy) = std::nan("");
    }
  }

  // Remaining boundary point
  for (int ix = localmesh->xstart; ix <= localmesh->xend; ix++) {
    // Lower y -- non-boundary x and all z (could be all x but already set)
    for (int iy = 0; iy < localmesh->ystart; iy++) {
      var(ix, iy) = std::nan("");
    }

    // Lower y -- non-boundary x and all z (could be all x but already set)
    for (int iy = localmesh->yend + 1; iy < localmesh->LocalNy; iy++) {
      var(ix, iy) = std::nan("");
    }
  }
#endif  
  return;
}
