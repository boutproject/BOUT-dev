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

#ifdef TRACK
  name = "<F2D>";
#endif
}

Field2D::Field2D(const Field2D& f) : Field(f), data(f.data) {
  TRACE("Field2D(Field2D&)");

#ifdef TRACK
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

Field2D::Field2D(Array<BoutReal> data, Mesh* localmesh, CELL_LOC datalocation,
                 DirectionTypes directions_in)
    : Field(localmesh, datalocation, directions_in), data(data) {

  ASSERT1(fieldmesh != nullptr);

  nx = fieldmesh->LocalNx;
  ny = fieldmesh->LocalNy;

  ASSERT1(data.size() == nx * ny);

  setLocation(datalocation);
}

Field2D::~Field2D() {
  if(deriv)
    delete deriv;
}

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

Field2D* Field2D::timeDeriv() {
  if(deriv == nullptr)
    deriv = new Field2D{emptyFrom(*this)};
  return deriv;
}

////////////// Indexing ///////////////////

const Region<Ind2D> &Field2D::getRegion(REGION region) const {
  return fieldmesh->getRegion2D(toString(region));
};
const Region<Ind2D> &Field2D::getRegion(const std::string &region_name) const {
  return fieldmesh->getRegion2D(region_name);
};

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

#ifdef TRACK
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
#ifdef TRACK
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
    if (reg->label.compare(region) == 0) {
      region_found = true;
      auto op = std::unique_ptr<BoundaryOp>{
          dynamic_cast<BoundaryOp*>(bfact->create(condition, reg))};
      op->apply(*this);
      break;
    }
  }

  if (!region_found) {
    throw BoutException("Region '%s' not found", region.c_str());
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

  const auto region = f.getRegion(rgn);
  BoutReal result = f[*region.cbegin()];

  BOUT_FOR_OMP(i, region, parallel for reduction(min:result)) {
    if (f[i] < result) {
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

BoutReal max(const Field2D &f, bool allpe,REGION rgn) {
  TRACE("Field2D::Max() %s",allpe? "over all PEs" : "");

  checkData(f);

  const auto region = f.getRegion(rgn);
  BoutReal result = f[*region.cbegin()];

  BOUT_FOR_OMP(i, region, parallel for reduction(max:result)) {
    if (f[i] > result) {
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

bool finite(const Field2D &f, REGION rgn) {
  TRACE("finite(Field2D)");

  if (!f.isAllocated()) {
    return false;
  }

  BOUT_FOR_SERIAL(i, f.getRegion(rgn)) {
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
    Field2D result{emptyFrom(f)};                                                        \
    BOUT_FOR(d, result.getRegion(rgn)) { result[d] = func(f[d]); }                       \
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

  BOUT_FOR(d, result.getRegion(rgn)) {
    if (result[d] < f) {
      result[d] = f;
    }
  }

  return result;
}

Field2D pow(const Field2D &lhs, const Field2D &rhs, REGION rgn) {
  TRACE("pow(Field2D, Field2D)");
  // Check if the inputs are allocated
  checkData(lhs);
  checkData(rhs);
  ASSERT1(areFieldsCompatible(lhs, rhs));

  // Define and allocate the output result
  Field2D result{emptyFrom(lhs)};

  BOUT_FOR(i, result.getRegion(rgn)) { result[i] = ::pow(lhs[i], rhs[i]); }

  checkData(result);
  return result;
}

Field2D pow(const Field2D &lhs, BoutReal rhs, REGION rgn) {
  TRACE("pow(Field2D, BoutReal)");
  // Check if the inputs are allocated
  checkData(lhs);
  checkData(rhs);

  // Define and allocate the output result
  Field2D result{emptyFrom(lhs)};

  BOUT_FOR(i, result.getRegion(rgn)) { result[i] = ::pow(lhs[i], rhs); }

  checkData(result);
  return result;
}

Field2D pow(BoutReal lhs, const Field2D &rhs, REGION rgn) {
  TRACE("pow(lhs, Field2D)");
  // Check if the inputs are allocated
  checkData(lhs);
  checkData(rhs);

  // Define and allocate the output result
  Field2D result{emptyFrom(rhs)};

  BOUT_FOR(i, result.getRegion(rgn)) { result[i] = ::pow(lhs, rhs[i]); }

  checkData(result);
  return result;
}

namespace {
  // Internal routine to avoid ugliness with interactions between CHECK
  // levels and UNUSED parameters
#if CHECK > 2
void checkDataIsFiniteOnRegion(const Field2D& f, REGION region) {
  // Do full checks
  BOUT_FOR_SERIAL(i, f.getRegion(region)) {
    if (!::finite(f[i])) {
      throw BoutException("Field2D: Operation on non-finite data at [%d][%d]\n", i.x(),
                          i.y());
    }
  }
}
#elif CHECK > 0
// No-op for no checking
void checkDataIsFiniteOnRegion(const Field2D &UNUSED(f), REGION UNUSED(region)) {}
#endif
}

#if CHECK > 0
/// Check if the data is valid
void checkData(const Field2D &f, REGION region) {
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
