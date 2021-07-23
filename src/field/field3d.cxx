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

#include "bout/build_config.hxx"

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
#if BOUT_USE_TRACK
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
Field3D::Field3D(const Field3D& f)
    : Field(f), data(f.data), yup_fields(f.yup_fields), ydown_fields(f.ydown_fields) {

  TRACE("Field3D(Field3D&)");

  if (fieldmesh) {
    nx = fieldmesh->LocalNx;
    ny = fieldmesh->LocalNy;
    nz = fieldmesh->LocalNz;
  }
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

Field3D::Field3D(Array<BoutReal> data_in, Mesh* localmesh, CELL_LOC datalocation,
                 DirectionTypes directions_in)
    : Field(localmesh, datalocation, directions_in), data(std::move(data_in)) {
  TRACE("Field3D: Copy constructor from Array and Mesh");

  nx = fieldmesh->LocalNx;
  ny = fieldmesh->LocalNy;
  nz = fieldmesh->LocalNz;

  ASSERT1(data.size() == nx * ny * nz);
}

Field3D::~Field3D() { delete deriv; }

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
        "Field3D: Call to ynext with {:d} which is more than number of yguards ({:d})",
        dir, fieldmesh->ystart);
  }
#endif

  // ynext uses 1-indexing, but yup wants 0-indexing
  if (dir > 0) {
    return yup(dir - 1);
  } else if (dir < 0) {
    return ydown(- dir - 1);
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

bool Field3D::requiresTwistShift(bool twist_shift_enabled) {
  // Workaround for 3D coordinates.
  // We need to communicate in the coordinates constructor in that
  // case a Field3D, but coordinates isn't valid yet. As such we
  // disable twist-shift in that case.
  if (getCoordinates() == nullptr) {
    return false;
  }
  return getCoordinates()->getParallelTransform().requiresTwistShift(twist_shift_enabled,
      getDirectionY());
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
}
const Region<Ind3D> &Field3D::getRegion(const std::string &region_name) const {
  return fieldmesh->getRegion3D(region_name);
}

const Region<Ind2D> &Field3D::getRegion2D(REGION region) const {
  return fieldmesh->getRegion2D(toString(region));
}
const Region<Ind2D> &Field3D::getRegion2D(const std::string &region_name) const {
  return fieldmesh->getRegion2D(region_name);
}

/***************************************************************
 *                         OPERATORS 
 ***************************************************************/

/////////////////// ASSIGNMENT ////////////////////

Field3D & Field3D::operator=(const Field3D &rhs) {
  /// Check for self-assignment
  if(this == &rhs)
    return(*this); // skip this assignment

  TRACE("Field3D: Assignment from Field3D");

  // Copy base slice
  Field::operator=(rhs);

  // Copy parallel slices or delete existing ones.
  yup_fields = rhs.yup_fields;
  ydown_fields = rhs.ydown_fields;

  // Copy the data and data sizes
  nx = rhs.nx;
  ny = rhs.ny;
  nz = rhs.nz;

  data = rhs.data;

  return *this;
}

Field3D& Field3D::operator=(Field3D&& rhs) {
  TRACE("Field3D: Assignment from Field3D");

  // Move parallel slices or delete existing ones.
  yup_fields = std::move(rhs.yup_fields);
  ydown_fields = std::move(rhs.ydown_fields);

  // Move the data and data sizes
  nx = rhs.nx;
  ny = rhs.ny;
  nz = rhs.nz;

  data = std::move(rhs.data);

  // Move base slice last
  Field::operator=(std::move(rhs));

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
  ASSERT1_FIELDS_COMPATIBLE(*this, rhs);

  /// Copy data
  BOUT_FOR(i, rhs.getRegion("RGN_ALL")) {
    for (int iz = 0; iz < nz; iz++) {
      (*this)(i, iz) = rhs[i];
    }
  }

  return *this;
}

void Field3D::operator=(const FieldPerp &rhs) {
  TRACE("Field3D = FieldPerp");

  ASSERT1_FIELDS_COMPATIBLE(*this, rhs);
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
    if (not isBoundarySet()) {
      output_warn << "WARNING: Call to Field3D::applyBoundary(), but no boundary set\n";
    }
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
    for (const auto& bndry : getBoundaryOps()) {
      // Always apply to the values when initialising
      // fields, otherwise apply only if wanted
      if (!bndry->apply_to_ddt || init) {
        bndry->apply(*this);
      }
    }
  }
}

void Field3D::applyBoundary(BoutReal t) {
  TRACE("Field3D::applyBoundary()");

#if CHECK > 0
  if (not isBoundarySet()) {
    output_warn << "WARNING: Call to Field3D::applyBoundary(t), but no boundary set.\n";
  }
#endif

  checkData(*this);

  if (background != nullptr) {
    // Apply boundary to the total of this and background

    Field3D tot = *this + (*background);
    tot.copyBoundary(*this);
    tot.applyBoundary(t);
    *this = tot - (*background);
  } else {
    // Apply boundary to this field
    for (const auto& bndry : getBoundaryOps()) {
      bndry->apply(*this, t);
    }
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
    auto op = std::unique_ptr<BoundaryOp>{
        dynamic_cast<BoundaryOp*>(bfact->create(condition, reg))};
    op->apply(*this);
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
    if (reg->label == region) {
      region_found = true;
      auto op = std::unique_ptr<BoundaryOp>{
          dynamic_cast<BoundaryOp*>(bfact->create(condition, reg))};
      op->apply(*this);
      break;
    }
  }

  if (!region_found) {
    throw BoutException("Region '{:s}' not found", region);
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

  for (const auto& bndry : getBoundaryOps()) {
    bndry->apply_ddt(*this);
  }

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
    for (const auto& bndry : getBoundaryOpPars()) {
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
    for (const auto& bndry : getBoundaryOpPars()) {
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
      auto op = std::unique_ptr<BoundaryOpPar>{
          dynamic_cast<BoundaryOpPar*>(bfact->create(condition, reg))};
      op->apply(*this);
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
      if (reg->label == region) {
        auto op = std::unique_ptr<BoundaryOpPar>{
            dynamic_cast<BoundaryOpPar*>(bfact->create(condition, reg))};
        op->apply(*this);
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
      if (reg->label == region) {
        // BoundaryFactory can't create boundaries using Field3Ds, so get temporary
        // boundary of the right type
        auto tmp = std::unique_ptr<BoundaryOpPar>{
            dynamic_cast<BoundaryOpPar*>(bfact->create(condition, reg))};
        // then clone that with the actual argument
        auto op = std::unique_ptr<BoundaryOpPar>{tmp->clone(reg, f)};
        op->apply(*this);
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

Field3D pow(const Field3D &lhs, const Field2D &rhs, const std::string& rgn) {
  TRACE("pow(Field3D, Field2D)");
  // Check if the inputs are allocated
  checkData(lhs);
  checkData(rhs);
  ASSERT1_FIELDS_COMPATIBLE(lhs, rhs);

  // Define and allocate the output result
  Field3D result{emptyFrom(lhs)};

  BOUT_FOR(i, result.getRegion(rgn)) { result[i] = ::pow(lhs[i], rhs[i]); }

  checkData(result);
  return result;
}

FieldPerp pow(const Field3D &lhs, const FieldPerp &rhs, const std::string& rgn) {
  TRACE("pow(Field3D, FieldPerp)");

  checkData(lhs);
  checkData(rhs);
  ASSERT1_FIELDS_COMPATIBLE(lhs, rhs);

  FieldPerp result{emptyFrom(rhs)};
  
  BOUT_FOR(i, result.getRegion(rgn)) {
    result[i] = ::pow(lhs(i, rhs.getIndex()), rhs[i]);
  }

  checkData(result);
  return result;
}

/////////////////////////////////////////////////////////////////////
// Friend functions

Field3D filter(const Field3D &var, int N0, const std::string& rgn) {
  TRACE("filter(Field3D, int)");
  
  checkData(var);

  int ncz = var.getNz();

  Field3D result{emptyFrom(var)};

  const auto region_str = toString(rgn);

  // Only allow a whitelist of regions for now
  ASSERT2(region_str == "RGN_ALL" || region_str == "RGN_NOBNDRY" ||
          region_str == "RGN_NOX" || region_str == "RGN_NOY");

  const Region<Ind2D> &region = var.getRegion2D(region_str);

  BOUT_OMP(parallel)
  {
    Array<dcomplex> f(ncz / 2 + 1);

    BOUT_FOR_INNER(i, region) {
      // Forward FFT
      rfft(var(i.x(), i.y()), ncz, f.begin());

      for (int jz = 0; jz <= ncz / 2; jz++) {
        if (jz != N0) {
          // Zero this component
          f[jz] = 0.0;
        }
      }

      // Reverse FFT
      irfft(f.begin(), ncz, result(i.x(), i.y()));
    }
  }

#if BOUT_USE_TRACK
  result.name = "filter(" + var.name + ")";
#endif

  checkData(result);
  return result;
}

// Fourier filter in z with zmin
Field3D lowPass(const Field3D &var, int zmax, bool keep_zonal, const std::string& rgn) {
  TRACE("lowPass(Field3D, {}, {})", zmax, keep_zonal);

  checkData(var);
  int ncz = var.getNz();

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
      rfft(var(i.x(), i.y()), ncz, f.begin());

      // Filter in z
      for (int jz = zmax + 1; jz <= ncz / 2; jz++)
        f[jz] = 0.0;

      // Filter zonal mode
      if (!keep_zonal) {
        f[0] = 0.0;
      }
      // Reverse FFT
      irfft(f.begin(), ncz, result(i.x(), i.y()));
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

  int ncz = localmesh->LocalNz;
  if(ncz == 1)
    return; // Shifting doesn't do anything
  
  Array<dcomplex> v(ncz/2 + 1);
  
  rfft(&(var(jx,jy,0)), ncz, v.begin()); // Forward FFT

  BoutReal zlength = var.getCoordinates()->zlength()(jx, jy);

  // Apply phase shift
  for(int jz=1;jz<=ncz/2;jz++) {
    BoutReal kwave=jz*2.0*PI/zlength; // wave number is 1/[rad]
    v[jz] *= dcomplex(cos(kwave*zangle) , -sin(kwave*zangle));
  }

  irfft(v.begin(), ncz, &(var(jx,jy,0))); // Reverse FFT
}

void shiftZ(Field3D &var, double zangle, const std::string& rgn) {
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

namespace {
  // Internal routine to avoid ugliness with interactions between CHECK
  // levels and UNUSED parameters
#if CHECK > 2
void checkDataIsFiniteOnRegion(const Field3D& f, const std::string& region) {
  // Do full checks
  BOUT_FOR_SERIAL(i, f.getDefaultRegion(region)) {
    if (!finite(f[i])) {
      throw BoutException("Field3D: Operation on non-finite data at [{:d}][{:d}][{:d}]\n",
                          i.x(), i.y(), i.z());
    }
  }
}
#elif CHECK > 0
// No-op for no checking
void checkDataIsFiniteOnRegion(const Field3D &UNUSED(f), const std::string& UNUSED(region)) {}
#endif
}

#if CHECK > 0
void checkData(const Field3D &f, const std::string& region) {
  if (!f.isAllocated())
    throw BoutException("Field3D: Operation on empty data\n");

  checkDataIsFiniteOnRegion(f, region);
}
#endif

Field2D DC(const Field3D &f, const std::string& rgn) {
  TRACE("DC(Field3D)");

  checkData(f);

  Mesh *localmesh = f.getMesh();
  Field2D result(localmesh, f.getLocation());
  result.allocate();

  BOUT_FOR(i, result.getRegion(rgn)) {
    result[i] = 0.0;
    for (int k = 0; k < localmesh->LocalNz; k++) {
      result[i] += f[localmesh->ind2Dto3D(i, k)];
    }
    result[i] /= (localmesh->LocalNz);
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

const Region<Ind3D>& Field3D::getDefaultRegion(const std::string& region_name) const {
  if (regionID != -1) {
    return fieldmesh->getRegion(regionID);
  }
  return fieldmesh->getRegion(region_name);
}

void Field3D::setRegion(const std::string& region_name) {
  regionID = fieldmesh->getRegionID(region_name);
}

void swap(Field3D& first, Field3D& second) noexcept {
  using std::swap;

  // Swap base class members
  swap(static_cast<Field&>(first), static_cast<Field&>(second));

  swap(first.data, second.data);
  swap(first.background, second.background);
  swap(first.nx, second.nx);
  swap(first.ny, second.ny);
  swap(first.nz, second.nz);
  swap(first.deriv, second.deriv);
  swap(first.yup_fields, second.yup_fields);
  swap(first.ydown_fields, second.ydown_fields);
}
