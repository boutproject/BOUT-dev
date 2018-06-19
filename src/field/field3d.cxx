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
Field3D::Field3D(Mesh *localmesh)
    : Field(localmesh), background(nullptr), deriv(nullptr), yup_field(nullptr),
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
    : Field(f.getMesh()), background(nullptr), deriv(nullptr), yup_field(nullptr),
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

Field3D::Field3D(const BoutReal val, Mesh *localmesh)
    : Field(localmesh), background(nullptr), deriv(nullptr), yup_field(nullptr),
      ydown_field(nullptr) {

  TRACE("Field3D: Copy constructor from value");

  location = CELL_CENTRE; // Cell centred variable by default

  boundaryIsSet = false;

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
#if CHECK > 2
    invalidateGuards(*this);
#endif
  } else
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

void Field3D::setLocation(CELL_LOC new_location) {
  if (getMesh()->StaggerGrids) {
    if (new_location == CELL_VSHIFT) {
      throw BoutException(
          "Field3D: CELL_VSHIFT cell location only makes sense for vectors");
    }
    if (new_location == CELL_DEFAULT) {
      new_location = CELL_CENTRE;
    }
    location = new_location;
  } else {
#if CHECK > 0
    if (new_location != CELL_CENTRE && new_location != CELL_DEFAULT) {
      throw BoutException("Field3D: Trying to set off-centre location on "
                          "non-staggered grid\n"
                          "         Did you mean to enable staggered grids?");
    }
#endif
    location = CELL_CENTRE;
  }
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
  }
  case RGN_NOBNDRY: {
    return IndexRange{fieldmesh->xstart, fieldmesh->xend,
        fieldmesh->ystart, fieldmesh->yend,
        0, nz-1};
  }
  case RGN_NOX: {
    return IndexRange{fieldmesh->xstart, fieldmesh->xend,
        0, ny-1,
        0, nz-1};
  }
  case RGN_NOY: {
    return IndexRange{0, nx-1,
        fieldmesh->ystart, fieldmesh->yend,
        0, nz-1};
  }
  default: {
    throw BoutException("Field3D::region() : Requested region not implemented");
  }
  };
}

const IndexRange Field3D::region2D(REGION rgn) const {
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

///////////////////// BOUNDARY CONDITIONS //////////////////

void Field3D::setBackground(const Field2D &f2d) {
  background = &f2d;
}

void Field3D::applyBoundary(bool init) {
  TRACE("Field3D::applyBoundary()");

#if CHECK > 0
  if (init) {

    if(!boundaryIsSet)
      output_warn << "WARNING: Call to Field3D::applyBoundary(), but no boundary set" << endl;
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
    output_warn << "WARNING: Call to Field3D::applyBoundary(t), but no boundary set." << endl;
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

Field3D operator-(const Field3D &f) { return -1.0 * f; }

#define F3D_OP_FPERP(op)                                                                 \
  FieldPerp operator op(const Field3D &lhs, const FieldPerp &rhs) {                      \
    FieldPerp result;                                                                    \
    result.allocate();                                                                   \
    result.setIndex(rhs.getIndex());                                                     \
    for (const auto &i : rhs)                                                            \
      result[i] = lhs[i] op rhs[i];                                                      \
    return result;                                                                       \
  }

//////////////// NON-MEMBER FUNCTIONS //////////////////

Field3D pow(const Field3D &lhs, const Field3D &rhs, REGION rgn) {
  TRACE("pow(Field3D, Field3D)");

  ASSERT1(lhs.getLocation() == rhs.getLocation());

  ASSERT1(lhs.getMesh() == rhs.getMesh());
  Field3D result(lhs.getMesh());
  result.allocate();

  // Iterate over indices
  for(const auto& i : result.region(rgn)) {
    result[i] = ::pow(lhs[i], rhs[i]);
  }
  
  result.setLocation( lhs.getLocation() );
  
  checkData(result);
  return result;
}

Field3D pow(const Field3D &lhs, const Field2D &rhs, REGION rgn) {
  TRACE("pow(Field3D, Field2D)");
  // Check if the inputs are allocated
  ASSERT1(lhs.isAllocated());
  ASSERT1(rhs.isAllocated());
  ASSERT1(lhs.getMesh() == rhs.getMesh());

  // Define and allocate the output result
  Field3D result(lhs.getMesh());
  result.allocate();

  // Iterate over indices
  for(const auto& i : result.region(rgn)) {
    result[i] = ::pow(lhs[i], rhs[i]);
  }

  result.setLocation( lhs.getLocation() );
  
  checkData(result);
  return result;
}

Field3D pow(const Field3D &lhs, const FieldPerp &rhs, REGION rgn) {
  TRACE("pow(Field3D, FieldPerp)");

  ASSERT1(lhs.getMesh() == rhs.getMesh());
  Field3D result(lhs.getMesh());
  result.allocate();

  // Iterate over indices
  for(const auto& i : result.region(rgn)) {
    result[i] = ::pow(lhs[i], rhs[i]);
  }

  result.setLocation( lhs.getLocation() );

  checkData(result);
  return result;
}

Field3D pow(const Field3D &lhs, BoutReal rhs, REGION rgn) {
  TRACE("pow(Field3D, BoutReal)");
  // Check if the inputs are allocated
  ASSERT1(lhs.isAllocated());

  Field3D result(lhs.getMesh());
  result.allocate();
  for(const auto& i : result.region(rgn)) {
    result[i] = ::pow(lhs[i], rhs);
  }
  
  result.setLocation( lhs.getLocation() );

  checkData(result);
  return result;
}

Field3D pow(BoutReal lhs, const Field3D &rhs, REGION rgn) {
  TRACE("pow(lhs, Field3D)");
  // Check if the inputs are allocated
  ASSERT1(rhs.isAllocated());

  // Define and allocate the output result
  Field3D result(rhs.getMesh());
  result.allocate();

  for(const auto& i : result.region(rgn)) {
    result[i] = ::pow(lhs, rhs[i]);
  }
  
  result.setLocation( rhs.getLocation() );

  checkData(result);
  return result;
}

BoutReal min(const Field3D &f, bool allpe, REGION rgn) {
  TRACE("Field3D::Min() %s",allpe? "over all PEs" : "");

  ASSERT2(f.isAllocated());

  BoutReal result = f[f.region(rgn).begin()];
  
  for(const auto& i: f.region(rgn))
    if(f[i] < result)
      result = f[i];
  
  if(allpe) {
    // MPI reduce
    BoutReal localresult = result;
    MPI_Allreduce(&localresult, &result, 1, MPI_DOUBLE, MPI_MIN, BoutComm::get());
  }

  return result;
}

BoutReal max(const Field3D &f, bool allpe, REGION rgn) {
  TRACE("Field3D::Max() %s",allpe? "over all PEs" : "");

  ASSERT2(f.isAllocated());
  
  BoutReal result = f[f.region(rgn).begin()];
  
  for(const auto& i: f.region(rgn))
    if(f[i] > result)
      result = f[i];
  
  if(allpe) {
    // MPI reduce
    BoutReal localresult = result;
    MPI_Allreduce(&localresult, &result, 1, MPI_DOUBLE, MPI_MAX, BoutComm::get());
  }
  
  return result;
}

BoutReal mean(const Field3D &f, bool allpe, REGION rgn) {
  TRACE("Field3D::mean() %s",allpe? "over all PEs" : "");

  ASSERT2(f.isAllocated());

  Mesh *localmesh = f.getMesh();

  // use first element for sum of values of f, second element for number of points
  BoutReal result[2] = {0., 0.};
  
  for(const auto& i: f.region(rgn)) {
    result[0] += f[i];
    result[1] += 1.;
  }

  if(allpe) {
    // MPI reduce
    BoutReal localresult[2] = {result[0], result[1]};
    MPI_Allreduce(&localresult, &result, 2, MPI_DOUBLE, MPI_SUM, BoutComm::get());
  }
  
  return result[0]/result[1];
}

/////////////////////////////////////////////////////////////////////
// Friend functions

/*!
 * This macro takes a function \p func, which is
 * assumed to operate on a single BoutReal and return
 * a single BoutReal, and wraps it up into a function
 * of a Field3D called \p name.
 *
 * @param name  The name of the function to define
 * @param func  The function to apply to each value
 *
 * If CHECK >= 1, checks if the Field3D is allocated
 *
 * Loops over the entire domain, applies function,
 * and uses checkData() to, if CHECK >= 3, check
 * result for non-finite numbers
 *
 */
#define F3D_FUNC(name, func)                                                             \
  const Field3D name(const Field3D &f, REGION rgn) {                                     \
    TRACE(#name "(Field3D)");                                                            \
    /* Check if the input is allocated */                                                \
    ASSERT1(f.isAllocated());                                                            \
    /* Define and allocate the output result */                                          \
    Field3D result(f.getMesh());                                                         \
    result.allocate();                                                                   \
    /* Loop over domain */                                                               \
    for (const auto &d : result.region(rgn)) {                                           \
      result[d] = func(f[d]);                                                            \
    }                                                                                    \
    result.setLocation(f.getLocation());                                                 \
    checkData(result);                                                                   \
    return result;                                                                       \
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

const Field3D filter(const Field3D &var, int N0, REGION rgn) {
  TRACE("filter(Field3D, int)");
  
  ASSERT1(var.isAllocated());

  Mesh *localmesh = var.getMesh();

  int ncz = localmesh->LocalNz;
  Array<dcomplex> f(ncz/2 + 1);

  Field3D result(localmesh);
  result.allocate();

  for (const auto &i : result.region2D(rgn)) {

    rfft(&(var(i.x, i.y, 0)), ncz, f.begin()); // Forward FFT

    for(int jz=0;jz<=ncz/2;jz++) {
      
      if(jz != N0) {
        // Zero this component
        f[jz] = 0.0;
      }
    }

    irfft(f.begin(), ncz, &(result(i.x, i.y, 0))); // Reverse FFT
  }
  
#ifdef TRACK
  result.name = "filter("+var.name+")";
#endif
  
  result.setLocation(var.getLocation());

  checkData(result);
  return result;
}

// Fourier filter in z
const Field3D lowPass(const Field3D &var, int zmax, REGION rgn) {
  TRACE("lowPass(Field3D, %d)", zmax);

  ASSERT1(var.isAllocated());

  Mesh *localmesh = var.getMesh();
  int ncz = localmesh->LocalNz;

  // Create an array 
  Array<dcomplex> f(ncz/2 + 1);
  
  if((zmax >= ncz/2) || (zmax < 0)) {
    // Removing nothing
    return var;
  }

  Field3D result(localmesh);
  result.allocate();

  for (const auto &i : result.region2D(rgn)) {
    // Take FFT in the Z direction
    rfft(&(var(i.x,i.y,0)), ncz, f.begin());
    
    // Filter in z
    for(int jz=zmax+1;jz<=ncz/2;jz++)
      f[jz] = 0.0;

    irfft(f.begin(), ncz, &(result(i.x,i.y,0))); // Reverse FFT
  }
  
  result.setLocation(var.getLocation());

  checkData(result);
  return result;
}

// Fourier filter in z with zmin
const Field3D lowPass(const Field3D &var, int zmax, int zmin, REGION rgn) {
  TRACE("lowPass(Field3D, %d, %d)", zmax, zmin);

  ASSERT1(var.isAllocated());
  Mesh *localmesh = var.getMesh();

  int ncz = localmesh->LocalNz;
  Array<dcomplex> f(ncz/2 + 1);
 
  if(((zmax >= ncz/2) || (zmax < 0)) && (zmin < 0)) {
    // Removing nothing
    return var;
  }

  Field3D result(localmesh);
  result.allocate();

  for (const auto &i : result.region2D(rgn)) {
    // Take FFT in the Z direction
    rfft(&(var(i.x,i.y,0)), ncz, f.begin());
    
    // Filter in z
    for(int jz=zmax+1;jz<=ncz/2;jz++)
      f[jz] = 0.0;

    // Filter zonal mode
    if(zmin==0) {
      f[0] = 0.0;
    }
    irfft(f.begin(), ncz, &(result(i.x,i.y,0))); // Reverse FFT
  }
  
  result.setLocation(var.getLocation());
  
  checkData(result);
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

void shiftZ(Field3D &var, double zangle, REGION rgn) {
  for (const auto &i : var.region2D(rgn))
    shiftZ(var, i.x, i.y, zangle);
}

bool finite(const Field3D &f, REGION rgn) {
  TRACE("finite( Field3D )");

  if (!f.isAllocated()) {
    return false;
  }

  for (const auto &i : f.region(rgn)) {
    if (!finite(f[i])) {
      return false;
    }
  }

  return true;
}

#if CHECK > 0
void checkData(const Field3D &f, REGION region) {
  if (!f.isAllocated())
    throw BoutException("Field3D: Operation on empty data\n");

#if CHECK > 2
  // Do full checks
  for (const auto &d : f.region(region)) {
    if (!finite(f[d])) {
      throw BoutException("Field3D: Operation on non-finite data at [%d][%d][%d]\n", d.x,
                          d.y, d.z);
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

const Field3D floor(const Field3D &var, BoutReal f, REGION rgn) {
  Field3D result = copy(var);
  
  for(const auto& d : result.region(rgn))
    if(result[d] < f)
      result[d] = f;
  
  return result;
}

Field2D DC(const Field3D &f, REGION rgn) {
  TRACE("DC(Field3D)");

  Mesh *localmesh = f.getMesh();
  Field2D result(localmesh);
  result.allocate();

  for (const auto &i : result.region(rgn)) {
    result(i.x,i.y) = 0.0;
    for (int k = 0; k < localmesh->LocalNz; k++) {
      result(i.x, i.y) += f(i.x, i.y, k);
    }
    result(i.x, i.y) /= (localmesh->LocalNz);
  }

  checkData(result);
  return result;
}

void invalidateGuards(Field3D &var){
#if CHECK > 2 // Strip out if not checking
  Mesh *localmesh = var.getMesh();

  // Inner x -- all y and all z
  for(int ix=0; ix<localmesh->xstart; ix++){
    for (int iy = 0; iy < localmesh->LocalNy; iy++) {
      for(int iz=0; iz<localmesh->LocalNz; iz++){
        var(ix, iy, iz) = std::nan("");
      }
    }
  }

  // Outer x -- all y and all z
  for (int ix = localmesh->xend + 1; ix < localmesh->LocalNx; ix++) {
    for (int iy = 0; iy < localmesh->LocalNy; iy++) {
      for(int iz=0; iz<localmesh->LocalNz; iz++){
        var(ix, iy, iz) = std::nan("");
      }
    }
  }

  // Remaining boundary point
  for (int ix = localmesh->xstart; ix <= localmesh->xend; ix++) {
    // Lower y -- non-boundary x and all z (could be all x but already set)
    for(int iy=0; iy<localmesh->ystart; iy++){
      for(int iz=0; iz<localmesh->LocalNz; iz++){
        var(ix, iy, iz) = std::nan("");
      }
    }
    // Lower y -- non-boundary x and all z (could be all x but already set)
    for(int iy=localmesh->yend+1; iy<localmesh->LocalNy; iy++){
      for(int iz=0; iz<localmesh->LocalNz; iz++){
        var(ix, iy, iz) = std::nan("");
      }
    }
  }
#endif  
  return;
}
