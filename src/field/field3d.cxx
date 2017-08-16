/*!*************************************************************************
 * \file field3d.cxx
 *
 * Class for 3D X-Y-Z scalar fields
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
#include <globals.hxx>

#include <cmath>

#include <field3d.hxx>
#include <utils.hxx>
#include <fft.hxx>
#include <dcomplex.hxx>
#include <interpolation.hxx>
#include <boundary_op.hxx>
#include <boundary_factory.hxx>
#include <boutexception.hxx>
#include <output.hxx>
#include <msg_stack.hxx>
#include <bout/constants.hxx>
#include <bout/assert.hxx>

/// Constructor
Field3D::Field3D(Mesh *msh)
    : Field(msh), background(nullptr), deriv(nullptr), yup_field(nullptr),
      ydown_field(nullptr) {
#ifdef TRACK
  name = "<F3D>";
#endif

  if (fieldmesh) {
    nx = fieldmesh->LocalNx;
    ny = fieldmesh->LocalNy;
    nz = fieldmesh->LocalNz;
  }
#if CHECK > 0
  else {
    nx = -1;
    ny = -1;
    nz = -1;
  }
#endif

  location = CELL_CENTRE; // Cell centred variable by default

  boundaryIsSet = false;
}

/// Doesn't copy any data, just create a new reference to the same data (copy on change
/// later)
Field3D::Field3D(const Field3D &f)
    : Field(f.fieldmesh),                // The mesh containing array sizes
      background(nullptr), data(f.data), // This handles references to the data array
      deriv(nullptr), yup_field(nullptr), ydown_field(nullptr) {

  TRACE("Field3D(Field3D&)");

#if CHECK > 2
  checkData(f);
#endif

  if (fieldmesh) {
    nx = fieldmesh->LocalNx;
    ny = fieldmesh->LocalNy;
    nz = fieldmesh->LocalNz;
  }
#if CHECK > 0
  else {
    nx = -1;
    ny = -1;
    nz = -1;
  }
#endif

  location = f.location;

  boundaryIsSet = false;
}

Field3D::Field3D(const Field2D &f)
    : Field(nullptr), background(nullptr), deriv(nullptr), yup_field(nullptr),
      ydown_field(nullptr) {

  TRACE("Field3D: Copy constructor from Field2D");

  location = CELL_CENTRE; // Cell centred variable by default

  boundaryIsSet = false;

  fieldmesh = mesh;
  nx = fieldmesh->LocalNx;
  ny = fieldmesh->LocalNy;
  nz = fieldmesh->LocalNz;

  *this = f;
}

Field3D::Field3D(const BoutReal val)
    : Field(nullptr), background(nullptr), deriv(nullptr), yup_field(nullptr),
      ydown_field(nullptr) {

  TRACE("Field3D: Copy constructor from value");

  location = CELL_CENTRE; // Cell centred variable by default

  boundaryIsSet = false;

  fieldmesh = mesh;
  nx = fieldmesh->LocalNx;
  ny = fieldmesh->LocalNy;
  nz = fieldmesh->LocalNz;

  *this = val;
}

Field3D::~Field3D() {
  /// Delete the time derivative variable if allocated
  if(deriv != NULL) {
    // The ddt of the yup/ydown_fields point to the same place as ddt.yup_field
    // only delete once
    // Also need to check that separate yup_field exists
    if ((yup_field != this) && (yup_field != nullptr))
      yup_field->deriv = nullptr;
    if ((ydown_field != this) && (ydown_field != nullptr))
      ydown_field->deriv = nullptr;

    // Now delete them as part of the deriv vector
    delete deriv;
  }
  
  if((yup_field != this) && (yup_field != nullptr))
    delete yup_field;
  
  if((ydown_field != this) && (ydown_field != nullptr))
    delete ydown_field;
}

void Field3D::allocate() {
  if(data.empty()) {
    if(!fieldmesh) {
      /// If no mesh, use the global
      fieldmesh = mesh;
      nx = fieldmesh->LocalNx;
      ny = fieldmesh->LocalNy;
      nz = fieldmesh->LocalNz;
    }
    data = Array<BoutReal>(nx*ny*nz);
  }else
    data.ensureUnique();
}

Field3D* Field3D::timeDeriv() {
  if(deriv == nullptr) {
    deriv = new Field3D(fieldmesh);
  }
  return deriv;
}

void Field3D::splitYupYdown() {
  TRACE("Field3D::splitYupYdown");
  
  if((yup_field != this) && (yup_field != nullptr))
    return;

  // yup_field and ydown_field null
  yup_field = new Field3D(fieldmesh);
  ydown_field = new Field3D(fieldmesh);
}

void Field3D::mergeYupYdown() {
  TRACE("Field3D::mergeYupYdown");
  
  if(yup_field == this && ydown_field == this)
    return;

  if(yup_field != nullptr){
    delete yup_field;
  }

  if(ydown_field != nullptr) {
    delete ydown_field;
  }

  yup_field = this;
  ydown_field = this;
}

Field3D& Field3D::ynext(int dir) {
  switch(dir) {
  case +1:
    return yup();
  case -1:
    return ydown();
  default:
    throw BoutException("Field3D: Call to ynext with strange direction %d. Only +/-1 currently supported", dir);
  }
}

const Field3D& Field3D::ynext(int dir) const {
  switch(dir) {
  case +1:
    return yup();
  case -1:
    return ydown();
  default:
    throw BoutException("Field3D: Call to ynext with strange direction %d. Only +/-1 currently supported", dir);
  }
}

void Field3D::setLocation(CELL_LOC loc) {
  if(loc == CELL_VSHIFT)
    throw BoutException("Field3D: CELL_VSHIFT cell location only makes sense for vectors");
  
  if(loc == CELL_DEFAULT)
    loc = CELL_CENTRE;
  
  location = loc;
}

CELL_LOC Field3D::getLocation() const {
  return location;
}

/***************************************************************
 *                         OPERATORS 
 ***************************************************************/

const DataIterator Field3D::iterator() const {
  return DataIterator(0, nx-1, 
                      0, ny-1,
                      0, nz-1);
}

const DataIterator Field3D::begin() const {
  return DataIterator(0, nx-1, 
                      0, ny-1,
                      0, nz-1);
}

const DataIterator Field3D::end() const {
  // end() iterator should be one past the last element
  return DataIterator(0, nx-1, 
                      0, ny-1,
                      0, nz-1,DI_GET_END);
}

const IndexRange Field3D::region(REGION rgn) const {
  switch(rgn) {
  case RGN_ALL: {
    return IndexRange{0, nx-1,
        0, ny-1,
        0, nz-1};
    break;
  }
  case RGN_NOBNDRY: {
    return IndexRange{fieldmesh->xstart, fieldmesh->xend,
        fieldmesh->ystart, fieldmesh->yend,
        0, nz-1};
    break;
  }
  case RGN_NOX: {
    return IndexRange{fieldmesh->xstart, fieldmesh->xend,
        0, ny-1,
        0, nz-1};
    break;
  }
  case RGN_NOY: {
    return IndexRange{0, nx-1,
        fieldmesh->ystart, fieldmesh->yend,
        0, nz-1};
    break;
  }
  default: {
    throw BoutException("Field3D::region() : Requested region not implemented");
  }
  };
}

/////////////////// ASSIGNMENT ////////////////////

Field3D & Field3D::operator=(const Field3D &rhs) {
  /// Check for self-assignment
  if(this == &rhs)
    return(*this); // skip this assignment

  TRACE("Field3D: Assignment from Field3D");
  
  /// Check that the data is valid
  checkData(rhs);
  
  // Copy the data and data sizes
  fieldmesh = rhs.fieldmesh;
  nx = rhs.nx; ny = rhs.ny; nz = rhs.nz; 
  
  data = rhs.data;
  
  location = rhs.location;
  
  return *this;
}

Field3D & Field3D::operator=(const Field2D &rhs) {
  TRACE("Field3D = Field2D");
  
  ASSERT1(rhs.isAllocated());
  
  /// Check that the data is valid
  checkData(rhs);
 
  /// Make sure there's a unique array to copy data into
  allocate();

  /// Copy data
  for(const auto& i : (*this))
    (*this)[i] = rhs[i];
  
  /// Only 3D fields have locations for now
  //location = CELL_CENTRE;
  
  return *this;
}

void Field3D::operator=(const FieldPerp &rhs) {
  ASSERT1(rhs.isAllocated());
  
  /// Make sure there's a unique array to copy data into
  allocate();

  /// Copy data
  for(const auto& i : rhs) {
    (*this)[i] = rhs[i];
  }
}

void Field3D::operator=(const bvalue &bv) {
  TRACE("Field3D = bvalue");
  
  allocate();

#if CHECK > 0
  if(!finite(bv.val))
    throw BoutException("Field3D: assignment from non-finite value at (%d,%d,%d)\n", 
			bv.jx, bv.jy,bv.jz);
#endif

  operator()(bv.jx, bv.jy,bv.jz) = bv.val;
}

Field3D & Field3D::operator=(const BoutReal val) {
  TRACE("Field3D = BoutReal");
  allocate();

#if CHECK > 0
  if(!finite(val))
    throw BoutException("Field3D: Assignment from non-finite BoutReal\n");
#endif
  for(const auto& i : (*this))
    (*this)[i] = val;

  // Only 3D fields have locations
  //location = CELL_CENTRE;
  // DON'T RE-SET LOCATION

  return *this;
}

/////////////////////////////////////////////////////////////////////

#define F3D_UPDATE_FIELD(op,bop,ftype)                       \
  Field3D & Field3D::operator op(const ftype &rhs) {         \
    TRACE("Field3D: %s %s", #op, #ftype);           \
    checkData(rhs) ;                                         \
    checkData(*this);                                        \
    if(data.unique()) {                                      \
      /* This is the only reference to this data */          \
      for(const auto& i : (*this))                                  \
        (*this)[i] op rhs[i];                                \
    }else {                                                  \
      /* Shared data */                                      \
      (*this) = (*this) bop rhs;                             \
    }                                                        \
    return *this;                                            \
  }

F3D_UPDATE_FIELD(+=, +, Field3D);    // operator+= Field3D
F3D_UPDATE_FIELD(-=, -, Field3D);    // operator-= Field3D
F3D_UPDATE_FIELD(*=, *, Field3D);    // operator*= Field3D
F3D_UPDATE_FIELD(/=, /, Field3D);    // operator/= Field3D

F3D_UPDATE_FIELD(+=, +, Field2D);    // operator+= Field2D
F3D_UPDATE_FIELD(-=, -, Field2D);    // operator-= Field2D
F3D_UPDATE_FIELD(*=, *, Field2D);    // operator*= Field2D
F3D_UPDATE_FIELD(/=, /, Field2D);    // operator/= Field2D

#define F3D_UPDATE_REAL(op,bop)                              \
  Field3D & Field3D::operator op(BoutReal rhs) {      \
    TRACE("Field3D: %s Field3D", #op);              \
    if(!finite(rhs))                                         \
      throw BoutException("Field3D: %s operator passed non-finite BoutReal number", #op); \
    checkData(*this);                                        \
                                                             \
    if(data.unique()) {                                      \
      /* This is the only reference to this data */          \
      for(const auto& i : (*this))                                  \
        (*this)[i] op rhs;                                   \
    }else {                                                  \
      /* Need to put result in a new block */                \
      (*this) = (*this) bop rhs;                             \
    }                                                        \
    return *this;                                            \
  }

F3D_UPDATE_REAL(+=,+);    // operator+= BoutReal
F3D_UPDATE_REAL(-=,-);    // operator-= BoutReal
F3D_UPDATE_REAL(*=,*);    // operator*= BoutReal
F3D_UPDATE_REAL(/=,/);    // operator/= BoutReal

/***************************************************************
 *                         STENCILS
 ***************************************************************/

void Field3D::setXStencil(stencil &fval, const bindex &bx, CELL_LOC loc) const {
  fval.jx = bx.jx;
  fval.jy = bx.jy;
  fval.jz = bx.jz;

  ASSERT1(isAllocated());
  
  fval.c  = operator()(bx.jx,  bx.jy, bx.jz);
  fval.p  = operator()(bx.jxp, bx.jy, bx.jz);
  fval.m  = operator()(bx.jxm, bx.jy, bx.jz);
  fval.pp = operator()(bx.jx2p, bx.jy, bx.jz);
  fval.mm = operator()(bx.jx2m, bx.jy, bx.jz);

  if(fieldmesh->StaggerGrids && (loc != CELL_DEFAULT) && (loc != location)) {
    // Non-centred stencil

    if((location == CELL_CENTRE) && (loc == CELL_XLOW)) {
      // Producing a stencil centred around a lower X value
      fval.pp = fval.p;
      fval.p  = fval.c;
      
    }else if(location == CELL_XLOW) {
      // Stencil centred around a cell centre
      
      fval.mm = fval.m;
      fval.m  = fval.c;
    }
    // Shifted in one direction -> shift in another
    // Could produce warning
  }
}

void Field3D::setYStencil(stencil &fval, const bindex &bx, CELL_LOC loc) const
{
  fval.jx = bx.jx;
  fval.jy = bx.jy;
  fval.jz = bx.jz;

  ASSERT1(isAllocated());
  
  fval.c = (*this)(bx.jx,bx.jy,bx.jz);
  fval.p = yup()(bx.jx,bx.jyp,bx.jz);
  fval.m = ydown()(bx.jx,bx.jym,bx.jz);
  if (yup_field == this && ydown_field == this){
    fval.pp = (*this)(bx.jx,bx.jy2p,bx.jz);
    fval.mm = (*this)(bx.jx,bx.jy2m,bx.jz);
  } else {
    fval.pp = nan("");
    fval.mm = nan("");
  }

  if(fieldmesh->StaggerGrids && (loc != CELL_DEFAULT) && (loc != location)) {
    // Non-centred stencil

    if((location == CELL_CENTRE) && (loc == CELL_YLOW)) {
      // Producing a stencil centred around a lower Y value
      fval.pp = fval.p;
      fval.p  = fval.c;
    }else if(location == CELL_YLOW) {
      // Stencil centred around a cell centre
      
      fval.mm = fval.m;
      fval.m  = fval.c;
    }
    // Shifted in one direction -> shift in another
    // Could produce warning
  }
}

void Field3D::setZStencil(stencil &fval, const bindex &bx, CELL_LOC loc) const {
  fval.jx = bx.jx;
  fval.jy = bx.jy;
  fval.jz = bx.jz;

  ASSERT1(isAllocated());

  fval.c = operator()(bx.jx,bx.jy,bx.jz);

  fval.p = operator()(bx.jx,bx.jy,bx.jzp);
  fval.m = operator()(bx.jx,bx.jy,bx.jzm);
  fval.pp = operator()(bx.jx,bx.jy,bx.jz2p);
  fval.mm = operator()(bx.jx,bx.jy,bx.jz2m);

  if(fieldmesh->StaggerGrids && (loc != CELL_DEFAULT) && (loc != location)) {
    // Non-centred stencil

    if((location == CELL_CENTRE) && (loc == CELL_ZLOW)) {
      // Producing a stencil centred around a lower Z value
      fval.pp = fval.p;
      fval.p  = fval.c;
      
    }else if(location == CELL_ZLOW) {
      // Stencil centred around a cell centre
      
      fval.mm = fval.m;
      fval.m  = fval.c;
    }
    // Shifted in one direction -> shift in another
    // Could produce warning
  }
}

///////////////////// BOUNDARY CONDITIONS //////////////////

void Field3D::setBackground(const Field2D &f2d) {
  background = &f2d;
}

void Field3D::applyBoundary(bool init) {
  TRACE("Field3D::applyBoundary()");

#if CHECK > 0
  if (init) {

    if(!boundaryIsSet)
      output << "WARNING: Call to Field3D::applyBoundary(), but no boundary set" << endl;
  }
#endif

  ASSERT1(isAllocated());
  
  if(background != NULL) {
    // Apply boundary to the total of this and background
    
    Field3D tot = *this + (*background);
    tot.copyBoundary(*this);
    tot.applyBoundary(init);
    *this = tot - (*background);
  } else {
    // Apply boundary to this field
    for(const auto& bndry : bndry_op)
      if ( !bndry->apply_to_ddt || init) // Always apply to the values when initialising fields, otherwise apply only if wanted
        bndry->apply(*this);
  }
}

void Field3D::applyBoundary(BoutReal t) {
  TRACE("Field3D::applyBoundary()");
  
#if CHECK > 0
  if(!boundaryIsSet)
    output << "WARNING: Call to Field3D::applyBoundary(t), but no boundary set." << endl;
#endif

  ASSERT1(isAllocated())

  if(background != NULL) {
    // Apply boundary to the total of this and background

    Field3D tot = *this + (*background);
    tot.copyBoundary(*this);
    tot.applyBoundary(t);
    *this = tot - (*background);
  }else {
    // Apply boundary to this field
    for(const auto& bndry : bndry_op)
      bndry->apply(*this,t);
  }
}

void Field3D::applyBoundary(const string &condition) {
  TRACE("Field3D::applyBoundary(condition)");
  
  ASSERT1(isAllocated());
  
  if(background != NULL) {
    // Apply boundary to the total of this and background
    
    Field3D tot = *this + (*background);
    tot.applyBoundary(condition);
    *this = tot - (*background);
    return;
  }

  /// Get the boundary factory (singleton)
  BoundaryFactory *bfact = BoundaryFactory::getInstance();
  
  /// Loop over the mesh boundary regions
  for(const auto& reg : fieldmesh->getBoundaries()) {
    BoundaryOp* op = static_cast<BoundaryOp*>(bfact->create(condition, reg));
    op->apply(*this);
    delete op;
  }

  //Field2D sets the corners to zero here, should we do the same here?
}

void Field3D::applyBoundary(const string &region, const string &condition) {
  ASSERT1(isAllocated());

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

  //Field2D sets the corners to zero here, should we do the same here?
}

void Field3D::applyTDerivBoundary() {
  TRACE("Field3D::applyTDerivBoundary()");
  
  ASSERT1(isAllocated());
  ASSERT1(deriv != NULL);
  ASSERT1(deriv->isAllocated());
  
  if(background != NULL)
    *this += *background;
    
  for(const auto& bndry : bndry_op)
    bndry->apply_ddt(*this);
  
  if(background != NULL)
    *this -= *background;
}

void Field3D::setBoundaryTo(const Field3D &f3d) {
  TRACE("Field3D::setBoundary(const Field3D&)");
  
  allocate(); // Make sure data allocated

  ASSERT1(f3d.isAllocated());

  /// Loop over boundary regions
  for(const auto& reg : fieldmesh->getBoundaries()) {
    /// Loop within each region
    for(reg->first(); !reg->isDone(); reg->next()) {
      for(int z=0;z<nz;z++) {
        // Get value half-way between cells
        BoutReal val = 0.5*(f3d(reg->x,reg->y,z) + f3d(reg->x-reg->bx, reg->y-reg->by, z));
        // Set to this value
        (*this)(reg->x,reg->y,z) = 2.*val - (*this)(reg->x-reg->bx, reg->y-reg->by, z);
      }
    }
  }
}

void Field3D::applyParallelBoundary() {

  TRACE("Field3D::applyParallelBoundary()");

  ASSERT1(isAllocated());

  if(background != NULL) {
    // Apply boundary to the total of this and background
    Field3D tot = *this + (*background);
    tot.applyParallelBoundary();
    *this = tot - (*background);
  } else {
    // Apply boundary to this field
    for(const auto& bndry : bndry_op_par) {
      bndry->apply(*this);
    }
  }
}

void Field3D::applyParallelBoundary(BoutReal t) {

  TRACE("Field3D::applyParallelBoundary(t)");

  ASSERT1(isAllocated());

  if(background != NULL) {
    // Apply boundary to the total of this and background
    Field3D tot = *this + (*background);
    tot.applyParallelBoundary(t);
    *this = tot - (*background);
  } else {
    // Apply boundary to this field
    for(const auto& bndry : bndry_op_par) {
      bndry->apply(*this, t);
    }
  }
}

void Field3D::applyParallelBoundary(const string &condition) {

  TRACE("Field3D::applyParallelBoundary(condition)");

  ASSERT1(isAllocated());

  if(background != NULL) {
    // Apply boundary to the total of this and background
    Field3D tot = *this + (*background);
    tot.applyParallelBoundary(condition);
    *this = tot - (*background);
  } else {
    /// Get the boundary factory (singleton)
    BoundaryFactory *bfact = BoundaryFactory::getInstance();

    /// Loop over the mesh boundary regions
    for(const auto& reg : fieldmesh->getBoundariesPar()) {
      BoundaryOpPar* op = static_cast<BoundaryOpPar*>(bfact->create(condition, reg));
      op->apply(*this);
      delete op;
    }
  }
}

void Field3D::applyParallelBoundary(const string &region, const string &condition) {

  TRACE("Field3D::applyParallelBoundary(region, condition)");

  ASSERT1(isAllocated());

  if(background != NULL) {
    // Apply boundary to the total of this and background
    Field3D tot = *this + (*background);
    tot.applyParallelBoundary(region, condition);
    *this = tot - (*background);
  } else {
    /// Get the boundary factory (singleton)
    BoundaryFactory *bfact = BoundaryFactory::getInstance();

    /// Loop over the mesh boundary regions
    for(const auto& reg : fieldmesh->getBoundariesPar()) {
      if(reg->label.compare(region) == 0) {
        BoundaryOpPar* op = static_cast<BoundaryOpPar*>(bfact->create(condition, reg));
        op->apply(*this);
        delete op;
        break;
      }
    }
  }
}

void Field3D::applyParallelBoundary(const string &region, const string &condition, Field3D *f) {

  TRACE("Field3D::applyParallelBoundary(region, condition, f)");

  ASSERT1(isAllocated());

  if(background != NULL) {
    // Apply boundary to the total of this and background
    Field3D tot = *this + (*background);
    tot.applyParallelBoundary(region, condition, f);
    *this = tot - (*background);
  } else {
    /// Get the boundary factory (singleton)
    BoundaryFactory *bfact = BoundaryFactory::getInstance();

    /// Loop over the mesh boundary regions
    for(const auto& reg : fieldmesh->getBoundariesPar()) {
      if(reg->label.compare(region) == 0) {
        // BoundaryFactory can't create boundaries using Field3Ds, so get temporary
        // boundary of the right type
        BoundaryOpPar* tmp = static_cast<BoundaryOpPar*>(bfact->create(condition, reg));
        // then clone that with the actual argument
        BoundaryOpPar* op = tmp->clone(reg, f);
        op->apply(*this);
        delete tmp;
        delete op;
        break;
      }
    }
  }
}


/***************************************************************
 *               NON-MEMBER OVERLOADED OPERATORS
 ***************************************************************/


const Field3D operator-(const Field3D &f) {
  return -1.0*f;
}

#define F3D_OP_FPERP(op)                     	                          \
  const FieldPerp operator op(const Field3D &lhs, const FieldPerp &rhs) { \
    FieldPerp result;                                                     \
    result.allocate();                                                    \
    result.setIndex(rhs.getIndex());                                      \
    for(const auto& i : rhs)                                                     \
      result[i] = lhs[i] op rhs[i];                                       \
    return result;                                                        \
  }

F3D_OP_FPERP(+);
F3D_OP_FPERP(-);
F3D_OP_FPERP(/);
F3D_OP_FPERP(*);

#define F3D_OP_FIELD(op, ftype)                                     \
  const Field3D operator op(const Field3D &lhs, const ftype &rhs) { \
    Field3D result;                                                 \
    result.allocate();                                              \
    for(const auto& i : lhs)                                               \
      result[i] = lhs[i] op rhs[i];                                 \
    result.setLocation( lhs.getLocation() );                        \
    return result;                                                  \
  }

F3D_OP_FIELD(+, Field3D);   // Field3D + Field3D
F3D_OP_FIELD(-, Field3D);   // Field3D - Field3D
F3D_OP_FIELD(*, Field3D);   // Field3D * Field3D
F3D_OP_FIELD(/, Field3D);   // Field3D / Field3D

F3D_OP_FIELD(+, Field2D);   // Field3D + Field2D
F3D_OP_FIELD(-, Field2D);   // Field3D - Field2D
F3D_OP_FIELD(*, Field2D);   // Field3D * Field2D
F3D_OP_FIELD(/, Field2D);   // Field3D / Field2D

#define F3D_OP_REAL(op)                                         \
  const Field3D operator op(const Field3D &lhs, BoutReal rhs) { \
    Field3D result;                                             \
    result.allocate();                                          \
    for(const auto& i : lhs)                                           \
      result[i] = lhs[i] op rhs;                                \
    result.setLocation( lhs.getLocation() );                    \
    return result;                                              \
  }

F3D_OP_REAL(+); // Field3D + BoutReal
F3D_OP_REAL(-); // Field3D - BoutReal
F3D_OP_REAL(*); // Field3D * BoutReal
F3D_OP_REAL(/); // Field3D / BoutReal

#define REAL_OP_F3D(op)                                         \
  const Field3D operator op(BoutReal lhs, const Field3D &rhs) { \
    Field3D result;                                             \
    result.allocate();                                          \
    for(const auto& i : rhs)                                           \
      result[i] = lhs op rhs[i];                                \
    result.setLocation( rhs.getLocation() );                    \
    return result;                                              \
  }

REAL_OP_F3D(+); // BoutReal + Field3D
REAL_OP_F3D(-); // BoutReal - Field3D
REAL_OP_F3D(*); // BoutReal * Field3D
REAL_OP_F3D(/); // BoutReal / Field3D

//////////////// NON-MEMBER FUNCTIONS //////////////////

Field3D pow(const Field3D &lhs, const Field3D &rhs) {
  TRACE("pow(Field3D, Field3D)");

  if(mesh->StaggerGrids && (lhs.getLocation() != rhs.getLocation())) {
    // Interpolate and call again
    return pow(lhs, interp_to(rhs, lhs.getLocation()));
  }
  
  Field3D result;
  result.allocate();

  // Iterate over indices
  for(const auto& i : result) {
    result[i] = ::pow(lhs[i], rhs[i]);
    ASSERT2( ::finite( result[i] ) );
  }
  
  result.setLocation( lhs.getLocation() );
  
  return result;
}

Field3D pow(const Field3D &lhs, const Field2D &rhs) {
  TRACE("pow(Field3D, Field2D)");
  // Check if the inputs are allocated
  ASSERT1(lhs.isAllocated());
  ASSERT1(rhs.isAllocated());

  // Define and allocate the output result  
  Field3D result;
  result.allocate();

  // Iterate over indices
  for(const auto& i : result) {
    result[i] = ::pow(lhs[i], rhs[i]);
    ASSERT3( ::finite( result[i] ) );
  }

  result.setLocation( lhs.getLocation() );
  
  return result;
}

Field3D pow(const Field3D &lhs, const FieldPerp &rhs) {
  TRACE("pow(Field3D, FieldPerp)");
  
  Field3D result;
  result.allocate();

  // Iterate over indices
  for(const auto& i : result) {
    result[i] = ::pow(lhs[i], rhs[i]);
    ASSERT2( ::finite( result[i] ) );
  }

  result.setLocation( lhs.getLocation() );
  return result;
}

Field3D pow(const Field3D &lhs, BoutReal rhs) {
  TRACE("pow(Field3D, BoutReal)");
  // Check if the inputs are allocated
  ASSERT1(lhs.isAllocated());

  Field3D result;
  result.allocate();
  for(const auto& i : result){
    result[i] = ::pow(lhs[i], rhs);
    ASSERT3(finite(result[i]));
  }
  
  result.setLocation( lhs.getLocation() );
  return result;
}

Field3D pow(BoutReal lhs, const Field3D &rhs) {
  TRACE("pow(lhs, Field3D)");
  // Check if the inputs are allocated
  ASSERT1(rhs.isAllocated());

  // Define and allocate the output result
  Field3D result;
  result.allocate();

  for(const auto& i : result){
    result[i] = ::pow(lhs, rhs[i]);
    ASSERT3(finite(result[i]));
  }
  
  result.setLocation( rhs.getLocation() );
  return result;
}

BoutReal min(const Field3D &f, bool allpe) {
  TRACE("Field3D::Min() %s",allpe? "over all PEs" : "");

  ASSERT2(f.isAllocated());

  BoutReal result = f[f.region(RGN_NOBNDRY).begin()];
  
  for(const auto& i: f.region(RGN_NOBNDRY))
    if(f[i] < result)
      result = f[i];
  
  if(allpe) {
    // MPI reduce
    BoutReal localresult = result;
    MPI_Allreduce(&localresult, &result, 1, MPI_DOUBLE, MPI_MIN, BoutComm::get());
  }

  return result;
}

BoutReal max(const Field3D &f, bool allpe) {
  TRACE("Field3D::Max() %s",allpe? "over all PEs" : "");

  ASSERT2(f.isAllocated());
  
  BoutReal result = f[f.region(RGN_NOBNDRY).begin()];
  
  for(const auto& i: f.region(RGN_NOBNDRY))
    if(f[i] > result)
      result = f[i];
  
  if(allpe) {
    // MPI reduce
    BoutReal localresult = result;
    MPI_Allreduce(&localresult, &result, 1, MPI_DOUBLE, MPI_MAX, BoutComm::get());
  }
  
  return result;
}

/////////////////////////////////////////////////////////////////////
// Friend functions

#define F3D_FUNC(name, func)                               \
  const Field3D name(const Field3D &f) {                   \
    TRACE(#name "(Field3D)");                     \
    /* Check if the input is allocated */                  \
    ASSERT1(f.isAllocated());                              \
    /* Define and allocate the output result */            \
    Field3D result;                                        \
    result.allocate();                                     \
    /* Loop over domain */                                 \
    for(const auto& d : result) {                                 \
      result[d] = func(f[d]);                              \
      /* If checking is set to 3 or higher, test result */ \
      ASSERT3(finite(result[d]));                          \
    }                                                      \
    result.setLocation(f.getLocation());                   \
    return result;                                         \
  }

F3D_FUNC(sqrt, ::sqrt);
F3D_FUNC(abs, ::fabs);

F3D_FUNC(exp, ::exp);
F3D_FUNC(log, ::log);

F3D_FUNC(sin, ::sin);
F3D_FUNC(cos, ::cos);
F3D_FUNC(tan, ::tan);

F3D_FUNC(sinh, ::sinh);
F3D_FUNC(cosh, ::cosh);
F3D_FUNC(tanh, ::tanh);

const Field3D filter(const Field3D &var, int N0) {
  TRACE("filter(Field3D, int)");
  
  ASSERT1(var.isAllocated());
  
  int ncz = mesh->LocalNz;
  Array<dcomplex> f(ncz/2 + 1);
  
  Field3D result;
  result.allocate();
  
  for(int jx=0;jx<mesh->LocalNx;jx++) {
    for(int jy=0;jy<mesh->LocalNy;jy++) {

      rfft(&(var(jx, jy, 0)), ncz, f.begin()); // Forward FFT

      for(int jz=0;jz<=ncz/2;jz++) {
	
	if(jz != N0) {
	  // Zero this component
	  f[jz] = 0.0;
	}
      }

      irfft(f.begin(), ncz, &(result(jx, jy, 0))); // Reverse FFT
    }
  }
  
#ifdef TRACK
  result.name = "filter("+var.name+")";
#endif
  
  result.setLocation(var.getLocation());

  return result;
}

// Fourier filter in z
const Field3D lowPass(const Field3D &var, int zmax) {
  TRACE("lowPass(Field3D, %d)", zmax);

  ASSERT1(var.isAllocated());
  
  int ncz = mesh->LocalNz;
  
  // Create an array 
  Array<dcomplex> f(ncz/2 + 1);
  
  if((zmax >= ncz/2) || (zmax < 0)) {
    // Removing nothing
    return var;
  }

  Field3D result;
  result.allocate();
  
  for(int jx=0;jx<mesh->LocalNx;jx++) {
    for(int jy=0;jy<mesh->LocalNy;jy++) {
      // Take FFT in the Z direction
      rfft(&(var(jx,jy,0)), ncz, f.begin());
      
      // Filter in z
      for(int jz=zmax+1;jz<=ncz/2;jz++)
	f[jz] = 0.0;

      irfft(f.begin(), ncz, &(result(jx,jy,0))); // Reverse FFT
    }
  }
  
  result.setLocation(var.getLocation());

  return result;
}

// Fourier filter in z with zmin
const Field3D lowPass(const Field3D &var, int zmax, int zmin) {
  TRACE("lowPass(Field3D, %d, %d)", zmax, zmin);

  ASSERT1(var.isAllocated());

  int ncz = mesh->LocalNz;
  Array<dcomplex> f(ncz/2 + 1);
 
  if(((zmax >= ncz/2) || (zmax < 0)) && (zmin < 0)) {
    // Removing nothing
    return var;
  }

  Field3D result;
  result.allocate();
  
  for(int jx=0;jx<mesh->LocalNx;jx++) {
    for(int jy=0;jy<mesh->LocalNy;jy++) {
      // Take FFT in the Z direction
      rfft(&(var(jx,jy,0)), ncz, f.begin());
      
      // Filter in z
      for(int jz=zmax+1;jz<=ncz/2;jz++)
	f[jz] = 0.0;

      // Filter zonal mode
      if(zmin==0) {
	f[0] = 0.0;
      }
      irfft(f.begin(), ncz, &(result(jx,jy,0))); // Reverse FFT
    }
  }
  
  result.setLocation(var.getLocation());
  
  return result;
}

/* 
 * Use FFT to shift by an angle in the Z direction
 */
void shiftZ(Field3D &var, int jx, int jy, double zangle) {
  TRACE("shiftZ");
  ASSERT1(var.isAllocated()); // Check that var has some data
  var.allocate(); // Ensure that var is unique
  
  int ncz = mesh->LocalNz;
  if(ncz == 1)
    return; // Shifting doesn't do anything
  
  Array<dcomplex> v(ncz/2 + 1);
  
  rfft(&(var(jx,jy,0)), ncz, v.begin()); // Forward FFT

  BoutReal zlength = mesh->coordinates()->zlength();
  // Apply phase shift
  for(int jz=1;jz<=ncz/2;jz++) {
    BoutReal kwave=jz*2.0*PI/zlength; // wave number is 1/[rad]
    v[jz] *= dcomplex(cos(kwave*zangle) , -sin(kwave*zangle));
  }

  irfft(v.begin(), ncz, &(var(jx,jy,0))); // Reverse FFT
}

void shiftZ(Field3D &var, double zangle) {
  for(int x=0;x<mesh->LocalNx;x++) 
    for(int y=0;y<mesh->LocalNy;y++)
      shiftZ(var, x, y, zangle);
}

bool finite(const Field3D &f) {
  TRACE("finite( Field3D )");

  if (!f.isAllocated()) {
    return false;
  }

  for (const auto &i : f) {
    if (!finite(f[i])) {
      return false;
    }
  }

  return true;
}

#if CHECK > 0
/// Check if the data is valid
void checkData(const Field3D &f)  {
  if(!f.isAllocated())
    throw BoutException("Field3D: Operation on empty data\n");
  
#if CHECK > 2
  //Do full checks
  for(const auto& d : f.region(RGN_NOBNDRY)) {
    if(!finite(f[d])) {
      throw BoutException("Field3D: Operation on non-finite data at [%d][%d][%d]\n", d.x, d.y, d.z);
    }
  }
#endif
}
#endif

const Field3D copy(const Field3D &f) {
  Field3D result = f;
  result.allocate();
  return result;
}

const Field3D floor(const Field3D &var, BoutReal f) {
  Field3D result = copy(var);
  
  for(const auto& d : result)
    if(result[d] < f)
      result[d] = f;
  
  return result;
}

Field2D DC(const Field3D &f) {
  TRACE("DC(Field3D)");
  
  Field2D result;
  result.allocate();

  for(int i=0;i<mesh->LocalNx;i++)
    for(int j=0;j<mesh->LocalNy;j++) {
      result(i,j) = 0.0;
      for(int k=0;k<mesh->LocalNz;k++)
	result(i,j) += f(i,j,k);
      result(i,j) /= (mesh->LocalNz);
    }
  
  return result;
}

