/**************************************************************************
 * Base class for fields
 *
 **************************************************************************
 * Copyright 2010 - 2026 BOUT++ contributors
 *
 * Contact: Ben Dudson, dudson2@llnl.gov
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

#include <bout/bout_types.hxx>
#include <bout/boutexception.hxx>
#include <bout/coordinates.hxx>
#include <bout/field.hxx>
#include <bout/field2d.hxx>
#include <bout/field3d.hxx>
#include <bout/field_data.hxx>
#include <bout/output_bout_types.hxx>
#include <string>

Field::Field(Mesh* localmesh, CELL_LOC location_in, DirectionTypes directions_in)
    : FieldData(localmesh, location_in), directions(directions_in) {}

int Field::getNx() const { return getMesh()->LocalNx; }

int Field::getNy() const { return getMesh()->LocalNy; }

int Field::getNz() const { return getMesh()->LocalNz; }

bool Field::isFci() const {
  const auto* coords = this->getCoordinates();
  if (coords == nullptr) {
    return false;
  }
  if (not coords->hasParallelTransform()) {
    return false;
  }
  return not coords->getParallelTransform().canToFromFieldAligned();
}

namespace bout {

template <typename T>
void checkFinite(const T& f, const std::string& name, const std::string& rgn) {

  if (!f.isAllocated()) {
    throw BoutException("{:s} is not allocated", name);
  }

  BOUT_FOR_SERIAL(i, f.getRegion(rgn)) {
    if (!std::isfinite(f[i])) {
      throw BoutException("{:s} is not finite at {}", name, i);
    }
  }
}

template <typename T>
void checkPositive(const T& f, const std::string& name, const std::string& rgn) {

  if (!f.isAllocated()) {
    throw BoutException("{:s} is not allocated", name);
  }

  BOUT_FOR_SERIAL(i, f.getRegion(rgn)) {
    if (f[i] <= 0.) {
      throw BoutException("{:s} ({:s} {:s}) is {:e} (not positive) at {}", name,
                          f.getLocation(), f.getDirections(), f[i], i);
    }
  }
}

template void checkFinite<Field2D>(const Field2D&, const std::string&,
                                   const std::string&);

template void checkFinite<Field3D>(const Field3D&, const std::string&,
                                   const std::string&);

template void checkFinite<FieldPerp>(const FieldPerp&, const std::string&,
                                     const std::string&);

template void checkPositive<Field2D>(const Field2D&, const std::string&,
                                     const std::string&);

template void checkPositive<Field3D>(const Field3D&, const std::string&,
                                     const std::string&);

template void checkPositive<FieldPerp>(const FieldPerp&, const std::string&,
                                       const std::string&);

} // namespace bout
