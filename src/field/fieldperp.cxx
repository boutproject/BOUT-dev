/**************************************************************************
 * Class for 2D X-Z slices
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

#include <bout/mesh.hxx>
#include <fieldperp.hxx>
#include <utils.hxx>
#include <boutexception.hxx>
#include <msg_stack.hxx>

FieldPerp::FieldPerp(Mesh *localmesh, CELL_LOC location_in, int yindex_in,
      DirectionTypes directions)
    : Field(localmesh, location_in, directions),
      yindex(yindex_in) {
  if (fieldmesh) {
    nx = fieldmesh->LocalNx;
    nz = fieldmesh->LocalNz;
  }
}

FieldPerp::FieldPerp(BoutReal val, Mesh *localmesh) : FieldPerp(localmesh) {
  *this = val;
}

FieldPerp::FieldPerp(Array<BoutReal> data_in, Mesh* localmesh, CELL_LOC location_in,
                     int yindex_in, DirectionTypes directions)
    : Field(localmesh, location_in, directions), yindex(yindex_in),
      nx(fieldmesh->LocalNx), nz(fieldmesh->LocalNz), data(std::move(data_in)) {

  ASSERT1(data.size() == nx * nz);

  setLocation(location_in);
}

FieldPerp& FieldPerp::allocate() {
  if (data.empty()) {
    if (!fieldmesh) {
      // fieldmesh was not initialized when this field was initialized, so use
      // the global mesh and set some members to default values
      fieldmesh = bout::globals::mesh;
      nx = fieldmesh->LocalNx;
      nz = fieldmesh->LocalNz;
    }
    data.reallocate(nx * nz);
#if CHECK > 2
    invalidateGuards(*this);
#endif
  } else
    data.ensureUnique();

  return *this;
}

/***************************************************************
 *                         ASSIGNMENT 
 ***************************************************************/

FieldPerp &FieldPerp::operator=(const FieldPerp &rhs) {
  /// Check for self-assignment
  if (this == &rhs) {
    return (*this); // skip this assignment
  }

  copyFieldMembers(rhs);

  nx = rhs.nx;
  nz = rhs.nz;
  yindex = rhs.yindex;
  data = rhs.data;

  return *this;
}

FieldPerp & FieldPerp::operator=(const BoutReal rhs) {
  TRACE("FieldPerp = BoutReal");

  allocate();

  BOUT_FOR(i, getRegion("RGN_ALL")) { (*this)[i] = rhs; }

  return *this;
}

const Region<IndPerp> &FieldPerp::getRegion(REGION region) const {
  return fieldmesh->getRegionPerp(toString(region));
};
const Region<IndPerp> &FieldPerp::getRegion(const std::string &region_name) const {
  return fieldmesh->getRegionPerp(region_name);
};

int FieldPerp::getGlobalIndex() const {
  auto& fieldmesh = *getMesh();
  const int start = fieldmesh.hasBndryLowerY() ? 0 : fieldmesh.ystart;
  const int end = fieldmesh.hasBndryUpperY() ? fieldmesh.LocalNy : fieldmesh.yend + 1;

  // Only use the global y index if it's either an interior (grid)
  // point, or a boundary point. Otherwise, use -1 to indicate a guard
  // cell or an invalid value. The actual FieldPerp value is still
  // written to file
  return (yindex >= start and yindex < end) ? fieldmesh.getGlobalYIndex(yindex) : -1;
}

FieldPerp& FieldPerp::setIndexFromGlobal(int y_global) {
  auto& fieldmesh = *getMesh();
  const int start = fieldmesh.hasBndryLowerY() ? 0 : fieldmesh.ystart;
  const int end = fieldmesh.hasBndryUpperY() ? fieldmesh.LocalNy : fieldmesh.yend + 1;

  // Only use the global y index if it's either an interior (grid)
  // point, or a boundary point. Otherwise, use -1 to indicate a
  // guard cell or an invalid value
  const int yindex_local = fieldmesh.getLocalYIndex(y_global);
  yindex = (yindex_local >= start and yindex_local < end) ? yindex_local : -1;
  return *this;
}

//////////////// NON-MEMBER FUNCTIONS //////////////////

FieldPerp toFieldAligned(const FieldPerp& f, const std::string& region) {
  return f.getCoordinates()->getParallelTransform().toFieldAligned(f, region);
}

FieldPerp fromFieldAligned(const FieldPerp& f, const std::string& region) {
  return f.getCoordinates()->getParallelTransform().fromFieldAligned(f, region);
}

////////////// NON-MEMBER OVERLOADED OPERATORS //////////////

// Unary minus
FieldPerp operator-(const FieldPerp &f) { return -1.0 * f; }

/////////////////////////////////////////////////
// functions

const FieldPerp sliceXZ(const Field3D& f, int y) {
  // Source field should be valid
  checkData(f);

  FieldPerp result(f.getMesh(), f.getLocation(), y,
                   {f.getDirectionY(), f.getDirectionZ()});

  // Allocate memory
  result.allocate();
  BOUT_FOR(i, result.getRegion("RGN_ALL")) { result[i] = f(i, y); }

  checkData(result);
  return result;
}

#if CHECK > 2
void checkDataIsFiniteOnRegion(const FieldPerp &f, const std::string& region) {
  // Do full checks
  BOUT_FOR_SERIAL(i, f.getRegion(region)) {
    if (!::finite(f[i])) {
      throw BoutException("FieldPerp: Operation on non-finite data at [%d][%d]\n", i.x(),
                          i.z());
    }
  }
}
#else
void checkDataIsFiniteOnRegion(const FieldPerp &UNUSED(f), const std::string& UNUSED(region)) {}
#endif


#if CHECK > 0
/// Check if the data is valid
void checkData(const FieldPerp &f, const std::string& region) {
  if (!f.isAllocated()) {
    throw BoutException("FieldPerp: Operation on empty data\n");
  }

  ASSERT3(f.getIndex() >= 0 && f.getIndex() < f.getMesh()->LocalNy);

  checkDataIsFiniteOnRegion(f, region);
}
#endif

#if CHECK > 2
void invalidateGuards(FieldPerp &var) {
  BOUT_FOR(i, var.getRegion("RGN_GUARDS")) { var[i] = BoutNaN; }
}
#endif

bool operator==(const FieldPerp &a, const FieldPerp &b) {
  if (!a.isAllocated() || !b.isAllocated()) {
    return false;
  }
  return min(abs(a - b)) < 1e-10;
}

std::ostream& operator<<(std::ostream &out, const FieldPerp &value) {
  out << toString(value);
  return out;
}
