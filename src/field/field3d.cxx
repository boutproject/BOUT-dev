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
Field3D::Field3D(Mesh* localmesh, CELL_LOC location_in,
                 DirectionTypes directions_in)
    : Field(localmesh, location_in, directions_in) {
#ifdef TRACK
  name = "<F3D>";
#endif

  if (fieldmesh) {
    nx = fieldmesh->LocalNx;
    ny = fieldmesh->LocalNy;
    nz = fieldmesh->LocalNz;
  }
}

/// Doesn't copy any data, just create a new reference to the same data (copy on change
/// later)
Field3D::Field3D(const Field3D& f) : Field(f), data(f.data) {

  TRACE("Field3D(Field3D&)");

  if (fieldmesh) {
    nx = fieldmesh->LocalNx;
    ny = fieldmesh->LocalNy;
    nz = fieldmesh->LocalNz;
  }

  location = f.location;
  fieldCoordinates = f.fieldCoordinates;
}

Field3D::Field3D(const Field2D& f) : Field(f) {

  TRACE("Field3D: Copy constructor from Field2D");

  nx = fieldmesh->LocalNx;
  ny = fieldmesh->LocalNy;
  nz = fieldmesh->LocalNz;

  *this = f;
}

Field3D::Field3D(const BoutReal val, Mesh* localmesh) : Field3D(localmesh) {

  TRACE("Field3D: Copy constructor from value");

  *this = val;
}

Field3D::Field3D(Array<BoutReal> data, Mesh* localmesh, CELL_LOC datalocation,
                 DirectionTypes directions_in)
    : Field(localmesh, datalocation, directions_in), data(data) {
  TRACE("Field3D: Copy constructor from Array and Mesh");

  nx = fieldmesh->LocalNx;
  ny = fieldmesh->LocalNy;
  nz = fieldmesh->LocalNz;

  ASSERT1(data.size() == nx * ny * nz);

  setLocation(datalocation);
}

Field3D::~Field3D() {
  /// Delete the time derivative variable if allocated
  if (deriv != nullptr) {
    delete deriv;
  }
}

Field3D& Field3D::allocate() {
  if(data.empty()) {
    if(!fieldmesh) {
      // fieldmesh was not initialized when this field was initialized, so use
      // the global mesh and set some members to default values
      fieldmesh = bout::globals::mesh;
      nx = fieldmesh->LocalNx;
      ny = fieldmesh->LocalNy;
      nz = fieldmesh->LocalNz;
    }
    data.reallocate(nx * ny * nz);
#if CHECK > 2
    invalidateGuards(*this);
#endif
  } else
    data.ensureUnique();

  return *this;
}

Field3D* Field3D::timeDeriv() {
  if(deriv == nullptr) {
    deriv = new Field3D{emptyFrom(*this)};
  }
  return deriv;
}

void Field3D::splitParallelSlices() {
  TRACE("Field3D::splitParallelSlices");
  
#if CHECK > 2
  if (yup_fields.size() != ydown_fields.size()) {
    throw BoutException("Field3D::splitParallelSlices: forward/backward parallel slices not in sync.\n"
                        "    This is an internal library error");
  }
#endif

  if (!yup_fields.empty()) {
    return;
  }

  for (int i = 0; i < fieldmesh->ystart; ++i) {
    // Note the fields constructed here will be fully overwritten by the
    // ParallelTransform, so we don't need a full constructor
    yup_fields.emplace_back(fieldmesh);
    ydown_fields.emplace_back(fieldmesh);
  }
}

void Field3D::clearParallelSlices() {
  TRACE("Field3D::clearParallelSlices");

#if CHECK > 2
  if (yup_fields.size() != ydown_fields.size()) {
    throw BoutException("Field3D::mergeYupYdown: forward/backward parallel slices not in sync.\n"
                        "    This is an internal library error");
  }
#endif

  if (yup_fields.empty() && ydown_fields.empty()) {
    return;
  }

  yup_fields.clear();
  ydown_fields.clear();
}

const Field3D& Field3D::ynext(int dir) const {
#if CHECK > 0
  // Asked for more than yguards
  if (std::abs(dir) > fieldmesh->ystart) {
    throw BoutException(
        "Field3D: Call to ynext with %d which is more than number of yguards (%d)", dir,
        fieldmesh->ystart);
  }
#endif

  // ynext uses 1-indexing, but yup wants 0-indexing
  if (dir > 0) {
    return yup(dir - 1);
  } else if (dir < 0) {
    return ydown(std::abs(dir) - 1);
  } else {
    return *this;
  }
}

Field3D &Field3D::ynext(int dir) {
  // Call the `const` version: need to add `const` to `this` to call
  // it, then throw it away after. This is ok because `this` wasn't
  // `const` to begin with.
  // See Effective C++, Scott Meyers, p23, for a better explanation
  return const_cast<Field3D&>(static_cast<const Field3D&>(*this).ynext(dir));
}

// Not in header because we need to access fieldmesh
BoutReal &Field3D::operator()(const IndPerp &d, int jy) {
  return operator[](fieldmesh->indPerpto3D(d, jy));
}

const BoutReal &Field3D::operator()(const IndPerp &d, int jy) const {
  return operator[](fieldmesh->indPerpto3D(d, jy));
}

BoutReal &Field3D::operator()(const Ind2D &d, int jz) {
  return operator[](fieldmesh->ind2Dto3D(d, jz));
}

const BoutReal &Field3D::operator()(const Ind2D &d, int jz) const {
  return operator[](fieldmesh->ind2Dto3D(d, jz));
}

const Region<Ind3D> &Field3D::getRegion(REGION region) const {
  return fieldmesh->getRegion3D(toString(region));
};
const Region<Ind3D> &Field3D::getRegion(const std::string &region_name) const {
  return fieldmesh->getRegion3D(region_name);
};

const Region<Ind2D> &Field3D::getRegion2D(REGION region) const {
  return fieldmesh->getRegion2D(toString(region));
};
const Region<Ind2D> &Field3D::getRegion2D(const std::string &region_name) const {
  return fieldmesh->getRegion2D(region_name);
};

/***************************************************************
 *                         OPERATORS 
 ***************************************************************/

/////////////////// ASSIGNMENT ////////////////////

Field3D & Field3D::operator=(const Field3D &rhs) {
  /// Check for self-assignment
  if(this == &rhs)
    return(*this); // skip this assignment

  TRACE("Field3D: Assignment from Field3D");

  // Delete existing parallel slices. We don't copy parallel slices, so any
  // that currently exist will be incorrect.
  clearParallelSlices();

  copyFieldMembers(rhs);

  // Copy the data and data sizes
  nx = rhs.nx;
  ny = rhs.ny;
  nz = rhs.nz;

  data = rhs.data;

  return *this;
}

Field3D & Field3D::operator=(const Field2D &rhs) {
  TRACE("Field3D = Field2D");

  /// Check that the data is allocated
  ASSERT1(rhs.isAllocated());

  // Delete existing parallel slices. We don't copy parallel slices, so any
  // that currently exist will be incorrect.
  clearParallelSlices();

  setLocation(rhs.getLocation());

  /// Make sure there's a unique array to copy data into
  allocate();
  ASSERT1(areFieldsCompatible(*this, rhs));

  /// Copy data
  BOUT_FOR(i, getRegion("RGN_ALL")) { (*this)[i] = rhs[i]; }

  return *this;
}

void Field3D::operator=(const FieldPerp &rhs) {
  TRACE("Field3D = FieldPerp");

  ASSERT1(areFieldsCompatible(*this, rhs));
  /// Check that the data is allocated
  ASSERT1(rhs.isAllocated());

  // Delete existing parallel slices. We don't copy parallel slices, so any
  // that currently exist will be incorrect.
  clearParallelSlices();

  /// Make sure there's a unique array to copy data into
  allocate();

  /// Copy data
  BOUT_FOR(i, rhs.getRegion("RGN_ALL")) { (*this)(i, rhs.getIndex()) = rhs[i]; }
}

Field3D & Field3D::operator=(const BoutReal val) {
  TRACE("Field3D = BoutReal");

  // Delete existing parallel slices. We don't copy parallel slices, so any
  // that currently exist will be incorrect.
  clearParallelSlices();

  allocate();

  BOUT_FOR(i, getRegion("RGN_ALL")) { (*this)[i] = val; }

  return *this;
}

Field3D& Field3D::calcParallelSlices() {
  getCoordinates()->getParallelTransform().calcParallelSlices(*this);
  return *this;
}

///////////////////// BOUNDARY CONDITIONS //////////////////

void Field3D::applyBoundary(bool init) {
  TRACE("Field3D::applyBoundary()");

#if CHECK > 0
  if (init) {

    if(!boundaryIsSet)
      output_warn << "WARNING: Call to Field3D::applyBoundary(), but no boundary set" << endl;
  }
#endif

  checkData(*this);

  if (background != nullptr) {
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

  checkData(*this);

  if (background != nullptr) {
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

void Field3D::applyBoundary(const std::string &condition) {
  TRACE("Field3D::applyBoundary(condition)");
  
  checkData(*this);

  if (background != nullptr) {
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

void Field3D::applyBoundary(const std::string &region, const std::string &condition) {
  TRACE("Field3D::applyBoundary(string, string)");
  checkData(*this);

  /// Get the boundary factory (singleton)
  BoundaryFactory *bfact = BoundaryFactory::getInstance();

  bool region_found = false;
  /// Loop over the mesh boundary regions
  for (const auto &reg : fieldmesh->getBoundaries()) {
    if (reg->label.compare(region) == 0) {
      region_found = true;
      BoundaryOp *op = static_cast<BoundaryOp *>(bfact->create(condition, reg));
      op->apply(*this);
      delete op;
      break;
    }
  }

  if (!region_found) {
    throw BoutException("Region '%s' not found", region.c_str());
  }

  //Field2D sets the corners to zero here, should we do the same here?
}

void Field3D::applyTDerivBoundary() {
  TRACE("Field3D::applyTDerivBoundary()");
  
  checkData(*this);
  ASSERT1(deriv != nullptr);
  checkData(*deriv);

  if (background != nullptr)
    *this += *background;
    
  for(const auto& bndry : bndry_op)
    bndry->apply_ddt(*this);

  if (background != nullptr)
    *this -= *background;
}

void Field3D::setBoundaryTo(const Field3D &f3d) {
  TRACE("Field3D::setBoundary(const Field3D&)");
  
  checkData(f3d);

  allocate(); // Make sure data allocated

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

  checkData(*this);

  if (background != nullptr) {
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

  checkData(*this);

  if (background != nullptr) {
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

void Field3D::applyParallelBoundary(const std::string &condition) {

  TRACE("Field3D::applyParallelBoundary(condition)");

  checkData(*this);

  if (background != nullptr) {
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

void Field3D::applyParallelBoundary(const std::string &region, const std::string &condition) {

  TRACE("Field3D::applyParallelBoundary(region, condition)");

  checkData(*this);

  if (background != nullptr) {
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

void Field3D::applyParallelBoundary(const std::string &region, const std::string &condition, Field3D *f) {

  TRACE("Field3D::applyParallelBoundary(region, condition, f)");

  checkData(*this);

  if (background != nullptr) {
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

//////////////// NON-MEMBER FUNCTIONS //////////////////

Field3D toFieldAligned(const Field3D& f, const REGION region) {
  return f.getCoordinates()->getParallelTransform().toFieldAligned(f, region);
}

Field3D fromFieldAligned(const Field3D& f, const REGION region) {
  return f.getCoordinates()->getParallelTransform().fromFieldAligned(f, region);
}

Field3D pow(const Field3D &lhs, const Field3D &rhs, REGION rgn) {
  TRACE("pow(Field3D, Field3D)");

  ASSERT1(areFieldsCompatible(lhs, rhs));

  Field3D result{emptyFrom(lhs)};

  BOUT_FOR(i, result.getRegion(rgn)) { result[i] = ::pow(lhs[i], rhs[i]); }

  checkData(result);
  return result;
}

Field3D pow(const Field3D &lhs, const Field2D &rhs, REGION rgn) {
  TRACE("pow(Field3D, Field2D)");
  // Check if the inputs are allocated
  checkData(lhs);
  checkData(rhs);
  ASSERT1(areFieldsCompatible(lhs, rhs));

  // Define and allocate the output result
  Field3D result{emptyFrom(lhs)};

  BOUT_FOR(i, result.getRegion(rgn)) { result[i] = ::pow(lhs[i], rhs[i]); }

  checkData(result);
  return result;
}

FieldPerp pow(const Field3D &lhs, const FieldPerp &rhs, REGION rgn) {
  TRACE("pow(Field3D, FieldPerp)");

  checkData(lhs);
  checkData(rhs);
  ASSERT1(areFieldsCompatible(lhs, rhs));

  FieldPerp result{emptyFrom(rhs)};
  result.allocate();
  
  BOUT_FOR(i, result.getRegion(rgn)) {
    result[i] = ::pow(lhs(i, rhs.getIndex()), rhs[i]);
  }

  checkData(result);
  return result;
}

Field3D pow(const Field3D &lhs, BoutReal rhs, REGION rgn) {
  TRACE("pow(Field3D, BoutReal)");
  // Check if the inputs are allocated
  checkData(lhs);
  checkData(rhs);

  Field3D result{emptyFrom(lhs)};

  BOUT_FOR(i, result.getRegion(rgn)) { result[i] = ::pow(lhs[i], rhs); }

  checkData(result);
  return result;
}

Field3D pow(BoutReal lhs, const Field3D &rhs, REGION rgn) {
  TRACE("pow(lhs, Field3D)");
  // Check if the inputs are allocated
  checkData(lhs);
  checkData(rhs);

  // Define and allocate the output result
  Field3D result{emptyFrom(rhs)};

  BOUT_FOR(i, result.getRegion(rgn)) { result[i] = ::pow(lhs, rhs[i]); }

  checkData(result);
  return result;
}

BoutReal min(const Field3D &f, bool allpe, REGION rgn) {
  TRACE("Field3D::Min() %s",allpe? "over all PEs" : "");

  checkData(f);

  const auto region = f.getRegion(rgn);
  BoutReal result = f[*region.cbegin()];

  BOUT_FOR_OMP(i, region, parallel for reduction(min:result)) {
    if(f[i] < result) {
      result = f[i];
    }
  }
  
  if(allpe) {
    // MPI reduce
    BoutReal localresult = result;
    MPI_Allreduce(&localresult, &result, 1, MPI_DOUBLE, MPI_MIN, BoutComm::get());
  }

  return result;
}

BoutReal max(const Field3D &f, bool allpe, REGION rgn) {
  TRACE("Field3D::Max() %s",allpe? "over all PEs" : "");

  checkData(f);

  const auto region = f.getRegion(rgn);
  BoutReal result = f[*region.cbegin()];

  BOUT_FOR_OMP(i, region, parallel for reduction(max:result)) {
    if(f[i] > result) {
      result = f[i];
    }
  }
  
  if(allpe) {
    // MPI reduce
    BoutReal localresult = result;
    MPI_Allreduce(&localresult, &result, 1, MPI_DOUBLE, MPI_MAX, BoutComm::get());
  }
  
  return result;
}

BoutReal mean(const Field3D &f, bool allpe, REGION rgn) {
  TRACE("Field3D::mean() %s",allpe? "over all PEs" : "");

  checkData(f);

  // Intitialise the cummulative sum and counter
  BoutReal result = 0.;
  int count = 0;

  BOUT_FOR_OMP(i, f.getRegion(rgn), parallel for reduction(+:result,count)) {
    result += f[i];
    count += 1;
  }

  if(allpe) {
    // MPI reduce
    BoutReal localresult = result;
    MPI_Allreduce(&localresult, &result, 1, MPI_DOUBLE, MPI_SUM, BoutComm::get());
    int localcount = count;
    MPI_Allreduce(&localcount, &count, 1, MPI_INT, MPI_SUM, BoutComm::get());
  }

  return result / static_cast<BoutReal>(count);
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
  Field3D name(const Field3D &f, REGION rgn) {                                           \
    TRACE(#name "(Field3D)");                                                            \
    /* Check if the input is allocated */                                                \
    checkData(f);                                                                        \
    /* Define and allocate the output result */                                          \
    Field3D result{emptyFrom(f)};                                                        \
    BOUT_FOR(d, result.getRegion(rgn)) { result[d] = func(f[d]); }                       \
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

Field3D filter(const Field3D &var, int N0, REGION rgn) {
  TRACE("filter(Field3D, int)");
  
  checkData(var);

  Field3D result{emptyFrom(var)};

  const auto region_str = toString(rgn);

  // Only allow a whitelist of regions for now
  ASSERT2(region_str == "RGN_ALL" || region_str == "RGN_NOBNDRY" ||
          region_str == "RGN_NOX" || region_str == "RGN_NOY");

  const Region<Ind2D> &region = var.getRegion2D(region_str);

  auto localmesh = var.getMesh();
  int ncz = localmesh->zend + 1 - localmesh->zstart;

  BOUT_OMP(parallel)
  {
    Array<dcomplex> f(ncz / 2 + 1);

    BOUT_FOR_INNER(i, region) {
      // Forward FFT
      rfft(&var(i.x(), i.y(), localmesh->zstart), ncz, f.begin());

      for (int jz = 0; jz <= ncz / 2; jz++) {
        if (jz != N0) {
          // Zero this component
          f[jz] = 0.0;
        }
      }

      // Reverse FFT
      irfft(f.begin(), ncz, &result(i.x(), i.y(), localmesh->zstart));
    }
  }

#ifdef TRACK
  result.name = "filter(" + var.name + ")";
#endif

  checkData(result);
  return result;
}

// Fourier filter in z with zmin
Field3D lowPass(const Field3D &var, int zmax, bool keep_zonal, REGION rgn) {
  TRACE("lowPass(Field3D, %d, %d)", zmax, keep_zonal);

  checkData(var);
  Mesh *localmesh = var.getMesh();
  const int ncz = localmesh->zend + 1 - localmesh->zstart;

  if (((zmax >= ncz / 2) || (zmax < 0)) && keep_zonal) {
    // Removing nothing
    return var;
  }

  Field3D result{emptyFrom(var)};

  const auto region_str = toString(rgn);

  // Only allow a whitelist of regions for now
  ASSERT2(region_str == "RGN_ALL" || region_str == "RGN_NOBNDRY" ||
          region_str == "RGN_NOX" || region_str == "RGN_NOY");

  const Region<Ind2D> &region = var.getRegion2D(region_str);

  BOUT_OMP(parallel) {
    Array<dcomplex> f(ncz / 2 + 1);

    BOUT_FOR_INNER(i, region) {
      // Take FFT in the Z direction
      rfft(&var(i.x(), i.y(), localmesh->zstart), ncz, f.begin());

      // Filter in z
      for (int jz = zmax + 1; jz <= ncz / 2; jz++)
        f[jz] = 0.0;

      // Filter zonal mode
      if (!keep_zonal) {
        f[0] = 0.0;
      }
      // Reverse FFT
      irfft(f.begin(), ncz, &result(i.x(), i.y(), localmesh->zstart));
    }
  }

  checkData(result);
  return result;
}

/* 
 * Use FFT to shift by an angle in the Z direction
 */
void shiftZ(Field3D &var, int jx, int jy, double zangle) {
  TRACE("shiftZ");
  checkData(var);
  var.allocate(); // Ensure that var is unique
  Mesh *localmesh = var.getMesh();

  int ncz = localmesh->zend + 1 - localmesh->zstart;
  if(ncz == 1)
    return; // Shifting doesn't do anything
  
  Array<dcomplex> v(ncz/2 + 1);

  rfft(&(var(jx, jy, localmesh->zstart)), ncz, v.begin()); // Forward FFT

  BoutReal zlength = var.getCoordinates()->zlength();

  // Apply phase shift
  for(int jz=1;jz<=ncz/2;jz++) {
    BoutReal kwave=jz*2.0*PI/zlength; // wave number is 1/[rad]
    v[jz] *= dcomplex(cos(kwave*zangle) , -sin(kwave*zangle));
  }

  irfft(v.begin(), ncz, &(var(jx, jy, localmesh->zstart))); // Reverse FFT
}

void shiftZ(Field3D &var, double zangle, REGION rgn) {
  const auto region_str = toString(rgn);

  // Only allow a whitelist of regions for now
  ASSERT2(region_str == "RGN_ALL" || region_str == "RGN_NOBNDRY" ||
          region_str == "RGN_NOX" || region_str == "RGN_NOY");

  const Region<Ind2D> &region = var.getRegion2D(region_str);

  // Could be OpenMP if shiftZ(Field3D, int, int, double) didn't throw
  BOUT_FOR_SERIAL(i, region) {
    shiftZ(var, i.x(), i.y(), zangle);
  }
}

bool finite(const Field3D &f, REGION rgn) {
  TRACE("finite( Field3D )");

  if (!f.isAllocated()) {
    return false;
  }

  BOUT_FOR_SERIAL(i, f.getRegion(rgn)) {
    if (!finite(f[i])) {
      return false;
    }
  }

  return true;
}

namespace {
  // Internal routine to avoid ugliness with interactions between CHECK
  // levels and UNUSED parameters
#if CHECK > 2
void checkDataIsFiniteOnRegion(const Field3D& f, REGION region) {
  // Do full checks
  BOUT_FOR_SERIAL(i, f.getRegion(region)) {
    if (!finite(f[i])) {
      throw BoutException("Field3D: Operation on non-finite data at [%d][%d][%d]\n",
                          i.x(), i.y(), i.z());
    }
  }
}
#elif CHECK > 0
// No-op for no checking
void checkDataIsFiniteOnRegion(const Field3D &UNUSED(f), REGION UNUSED(region)) {}
#endif
}

#if CHECK > 0
void checkData(const Field3D &f, REGION region) {
  if (!f.isAllocated())
    throw BoutException("Field3D: Operation on empty data\n");

  checkDataIsFiniteOnRegion(f, region);
}
#endif

Field3D copy(const Field3D &f) {
  Field3D result = f;
  result.allocate();
  return result;
}

Field3D floor(const Field3D &var, BoutReal f, REGION rgn) {
  checkData(var);
  Field3D result = copy(var);

  BOUT_FOR(d, var.getRegion(rgn)) {
    if (result[d] < f) {
      result[d] = f;
    }
  }

  return result;
}

Field2D DC(const Field3D &f, REGION rgn) {
  TRACE("DC(Field3D)");

  checkData(f);

  Mesh *localmesh = f.getMesh();
  Field2D result(localmesh, f.getLocation());
  result.allocate();

  BOUT_FOR(i, result.getRegion(rgn)) {
    result[i] = 0.0;
    for (int k = localmesh->zstart; k <= localmesh->zend; k++) {
      result[i] += f[localmesh->ind2Dto3D(i, k)];
    }
    result[i] /= (localmesh->zend + 1 - localmesh->zstart);
  }

  checkData(result);
  return result;
}

#if CHECK > 2
void invalidateGuards(Field3D &var) {
  BOUT_FOR(i, var.getRegion("RGN_GUARDS")) { var[i] = BoutNaN; }
}
#endif

bool operator==(const Field3D &a, const Field3D &b) {
  if (!a.isAllocated() || !b.isAllocated()) {
    return false;
  }
  return min(abs(a - b)) < 1e-10;
}

std::ostream& operator<<(std::ostream &out, const Field3D &value) {
  out << toString(value);
  return out;
}
