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

#include "bout/build_config.hxx"

#include <boutcomm.hxx>
#include <bout/rvec.hxx>

#include <globals.hxx> // for mesh

#include <field2d.hxx>

#include <utils.hxx>

#include <boundary_op.hxx>
#include <boundary_factory.hxx>

#include <boutexception.hxx>
#include <msg_stack.hxx>
#include <bout/mesh.hxx>

#include <cmath>
#include <output.hxx>

#include <bout/assert.hxx>

Field2D::Field2D(Mesh* localmesh, CELL_LOC location_in,
      DirectionTypes directions_in)
    : Field(localmesh, location_in, directions_in) {

  if (fieldmesh) {
    nx = fieldmesh->LocalNx;
    ny = fieldmesh->LocalNy;
  }

#if BOUT_USE_TRACK
  name = "<F2D>";
#endif
}

Field2D::Field2D(const Field2D& f) : Field(f), data(f.data) {
  TRACE("Field2D(Field2D&)");

#if BOUT_USE_TRACK
  name = f.name;
#endif

  if (fieldmesh) {
    nx = fieldmesh->LocalNx;
    ny = fieldmesh->LocalNy;
  }

  location = f.location;
  fieldCoordinates = f.fieldCoordinates;
}

Field2D::Field2D(BoutReal val, Mesh* localmesh) : Field2D(localmesh) {
  *this = val;
}

Field2D::Field2D(Array<BoutReal> data_in, Mesh* localmesh, CELL_LOC datalocation,
                 DirectionTypes directions_in)
    : Field(localmesh, datalocation, directions_in), data(std::move(data_in)) {

  ASSERT1(fieldmesh != nullptr);

  nx = fieldmesh->LocalNx;
  ny = fieldmesh->LocalNy;

  ASSERT1(data.size() == nx * ny);

  setLocation(datalocation);
}

Field2D::~Field2D() { delete deriv; }

Field2D& Field2D::allocate() {
  if(data.empty()) {
    if(!fieldmesh) {
      // fieldmesh was not initialized when this field was initialized, so use
      // the global mesh and set some members to default values
      fieldmesh = bout::globals::mesh;
      nx = fieldmesh->LocalNx;
      ny = fieldmesh->LocalNy;
    }
    data.reallocate(nx*ny);
#if CHECK > 2
    invalidateGuards(*this);
#endif
  }else
    data.ensureUnique();

  return *this;
}

__host__ __device__ Field2D* Field2D::timeDeriv() {
  if(deriv == nullptr)
    deriv = new Field2D{emptyFrom(*this)};
  return deriv;
}

////////////// Indexing ///////////////////

const Region<Ind2D> &Field2D::getRegion(REGION region) const {
  return fieldmesh->getRegion2D(toString(region));
}
const Region<Ind2D> &Field2D::getRegion(const std::string &region_name) const {
  return fieldmesh->getRegion2D(region_name);
}

// Not in header because we need to access fieldmesh
__host__ __device__ BoutReal& Field2D::operator[](const Ind3D &d) {
  return operator[](fieldmesh->map3Dto2D(d));
}

__host__ __device__ const BoutReal& Field2D::operator[](const Ind3D &d) const {
  return operator[](fieldmesh->map3Dto2D(d));
}

///////////// OPERATORS ////////////////

Field2D &Field2D::operator=(const Field2D &rhs) {
  // Check for self-assignment
  if (this == &rhs)
    return (*this); // skip this assignment

  TRACE("Field2D: Assignment from Field2D");

#if BOUT_USE_TRACK
  name = rhs.name;
#endif

  copyFieldMembers(rhs);

  // Copy the data and data sizes
  nx = rhs.nx;
  ny = rhs.ny;

  // Copy reference to data
  data = rhs.data;

  return *this;
}

Field2D &Field2D::operator=(const BoutReal rhs) {
#if BOUT_USE_TRACK
  name = "<r2D>";
#endif

  TRACE("Field2D = BoutReal");
  allocate();

  BOUT_FOR(i, getRegion("RGN_ALL")) { (*this)[i] = rhs; }

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

void Field2D::applyBoundary(BoutReal time) {
  TRACE("Field2D::applyBoundary(time)");

#if CHECK > 0
  if (!boundaryIsSet) {
    output_warn << "WARNING: Call to Field2D::applyBoundary(time), but no boundary set\n";
  }
#endif

  checkData(*this);

  for (const auto& bndry : bndry_op) {
    bndry->apply(*this, time);
  }
}

void Field2D::applyBoundary(const std::string &condition) {
  TRACE("Field2D::applyBoundary(condition)");

  checkData(*this);

  /// Get the boundary factory (singleton)
  BoundaryFactory *bfact = BoundaryFactory::getInstance();

  /// Loop over the mesh boundary regions
  for(const auto& reg : fieldmesh->getBoundaries()) {
    auto op = std::unique_ptr<BoundaryOp>{
        dynamic_cast<BoundaryOp*>(bfact->create(condition, reg))};
    op->apply(*this);
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

void Field2D::applyBoundary(const std::string &region, const std::string &condition) {
  TRACE("Field2D::applyBoundary(string, string)");
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

namespace {
  // Internal routine to avoid ugliness with interactions between CHECK
  // levels and UNUSED parameters
#if CHECK > 2
void checkDataIsFiniteOnRegion(const Field2D& f, const std::string& region) {
  // Do full checks
  BOUT_FOR_SERIAL(i, f.getRegion(region)) {
    if (!::finite(f[i])) {
      throw BoutException("Field2D: Operation on non-finite data at [{:d}][{:d}]\n",
                          i.x(), i.y());
    }
  }
}
#elif CHECK > 0
// No-op for no checking
void checkDataIsFiniteOnRegion(const Field2D &UNUSED(f), const std::string& UNUSED(region)) {}
#endif
}

#if CHECK > 0
/// Check if the data is valid
void checkData(const Field2D &f, const std::string& region) {
  if (!f.isAllocated()) {
    throw BoutException("Field2D: Operation on empty data\n");
  }

  checkDataIsFiniteOnRegion(f, region);
}
#endif

#if CHECK > 2
void invalidateGuards(Field2D &var) {
  BOUT_FOR(i, var.getRegion("RGN_GUARDS")) { var[i] = BoutNaN; }
}
#endif

bool operator==(const Field2D &a, const Field2D &b) {
  if (!a.isAllocated() || !b.isAllocated()) {
    return false;
  }
  return min(abs(a - b)) < 1e-10;
}

std::ostream& operator<<(std::ostream &out, const Field2D &value) {
  out << toString(value);
  return out;
}
