/*!************************************************************************
 * \file index_derivs_interface.hxx
 *
 * Definition of main derivative kernels
 *
 **************************************************************************
 * Copyright 2018
 *    D.Dickinson
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

#ifndef __INDEX_DERIVS_INTERFACE_HXX__
#define __INDEX_DERIVS_INTERFACE_HXX__

#include <bout_types.hxx>
#include <msg_stack.hxx>

class Field3D;
class Field2D;

/// The main kernel used for all upwind and flux derivatives
template <typename T, DIRECTION direction, DERIV derivType>
T indexFlowDerivative(const T& vel, const T& f, CELL_LOC outloc,
                      const std::string& method, REGION region) {
  AUTO_TRACE();

  // Checks
  static_assert(std::is_base_of<Field2D, T>::value || std::is_base_of<Field3D, T>::value,
                "indexFlowDerivative only works on Field2D or Field3D input");

  static_assert(derivType == DERIV::Upwind || derivType == DERIV::Flux,
                "indexFlowDerivative only works for derivType in {Upwind, Flux}.");

  auto* localmesh = f.getMesh();

  // Check that the mesh is correct
  ASSERT1(vel.getMesh() == localmesh);
  // Check that the input variable has data
  ASSERT1(f.isAllocated());
  ASSERT1(vel.isAllocated());

  // Check the input data is valid
  {
    TRACE("Checking inputs");
    checkData(f);
    checkData(vel);
  }

  // Define properties of this approach
  const CELL_LOC allowedStaggerLoc = localmesh->getAllowedStaggerLoc(direction);

  // Handle the staggering
  const CELL_LOC inloc = f.getLocation(); // Input locations
  const CELL_LOC vloc = vel.getLocation();
  if (outloc == CELL_DEFAULT) {
    outloc = inloc;
  }
  const STAGGER stagger = localmesh->getStagger(vloc, inloc, outloc, allowedStaggerLoc);

  // Check for early exit
  const int nPoint = localmesh->getNpoints(direction);

  if (nPoint == 1) {
    auto tmp = T(0., localmesh);
    tmp.setLocation(outloc);
    return tmp;
  }

  // Lookup the method
  auto derivativeMethod = DerivativeStore<T>::getInstance().getFlowDerivative(
      method, direction, stagger, derivType);

  // Create the result field
  T result(localmesh);
  result.allocate(); // Make sure data allocated
  result.setLocation(outloc);

  // Apply method
  derivativeMethod(vel, f, result, region);

  // Check the result is valid
  {
    TRACE("Checking result");
    checkData(result);
  }

  return result;
}

/// The main kernel used for all standard derivatives
template <typename T, DIRECTION direction, DERIV derivType>
T indexStandardDerivative(const T& f, CELL_LOC outloc, const std::string& method,
                          REGION region) {
  AUTO_TRACE();

  // Checks
  static_assert(std::is_base_of<Field2D, T>::value || std::is_base_of<Field3D, T>::value,
                "indexStandardDerivative only works on Field2D or Field3D input");

  static_assert(derivType == DERIV::Standard || derivType == DERIV::StandardSecond
                    || derivType == DERIV::StandardFourth,
                "indexStandardDerivative only works for derivType in {Standard, "
                "StandardSecond, StandardFourth}");

  auto* localmesh = f.getMesh();

  // Check that the input variable has data
  ASSERT1(f.isAllocated());

  // Check the input data is valid
  {
    TRACE("Checking input");
    checkData(f);
  }

  // Define properties of this approach
  const CELL_LOC allowedStaggerLoc = localmesh->getAllowedStaggerLoc(direction);

  // Handle the staggering
  const CELL_LOC inloc = f.getLocation(); // Input location
  if (outloc == CELL_DEFAULT) {
    outloc = inloc;
  }
  const STAGGER stagger = localmesh->getStagger(inloc, outloc, allowedStaggerLoc);

  // Check for early exit
  const int nPoint = localmesh->getNpoints(direction);

  if (nPoint == 1) {
    auto tmp = T(0., localmesh);
    tmp.setLocation(outloc);
    return tmp;
  }

  // Lookup the method
  auto derivativeMethod = DerivativeStore<T>::getInstance().getStandardDerivative(
      method, direction, stagger, derivType);

  // Create the result field
  T result(localmesh);
  result.allocate(); // Make sure data allocated
  result.setLocation(outloc);

  // Apply method
  derivativeMethod(f, result, region);

  // Check the result is valid
  {
    TRACE("Checking result");
    checkData(result);
  }

  return result;
}

#endif
