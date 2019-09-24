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

#include <bout/deriv_store.hxx>
#include <bout_types.hxx>
#include <msg_stack.hxx>
#include "bout/traits.hxx"

class Field3D;
class Field2D;

namespace bout {
namespace derivatives {
namespace index {

/// The main kernel used for all upwind and flux derivatives
template <typename T, DIRECTION direction, DERIV derivType>
T flowDerivative(const T& vel, const T& f, CELL_LOC outloc, const std::string& method,
                 const std::string& region) {
  AUTO_TRACE();

  // Checks
  static_assert(bout::utils::is_Field2D<T>::value || bout::utils::is_Field3D<T>::value,
                "flowDerivative only works on Field2D or Field3D input");

  static_assert(derivType == DERIV::Upwind || derivType == DERIV::Flux,
                "flowDerivative only works for derivType in {Upwind, Flux}.");

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
    return zeroFrom(f).setLocation(outloc);
  }

  // Lookup the method
  auto derivativeMethod = DerivativeStore<T>::getInstance().getFlowDerivative(
      method, direction, stagger, derivType);

  // Create the result field
  T result{emptyFrom(f).setLocation(outloc)};

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
T standardDerivative(const T& f, CELL_LOC outloc, const std::string& method,
                     const std::string& region) {
  AUTO_TRACE();

  // Checks
  static_assert(bout::utils::is_Field2D<T>::value || bout::utils::is_Field3D<T>::value,
                "standardDerivative only works on Field2D or Field3D input");

  static_assert(derivType == DERIV::Standard || derivType == DERIV::StandardSecond
                    || derivType == DERIV::StandardFourth,
                "standardDerivative only works for derivType in {Standard, "
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
    return zeroFrom(f).setLocation(outloc);
  }

  // Lookup the method
  auto derivativeMethod = DerivativeStore<T>::getInstance().getStandardDerivative(
      method, direction, stagger, derivType);

  // Create the result field
  T result{emptyFrom(f).setLocation(outloc)};

  // Apply method
  derivativeMethod(f, result, region);

  // Check the result is valid
  {
    TRACE("Checking result");
    checkData(result);
  }

  return result;
}

////// STANDARD OPERATORS

////////////// X DERIVATIVE /////////////////
template <typename T>
T DDX(const T& f, CELL_LOC outloc = CELL_DEFAULT, const std::string& method = "DEFAULT",
      const std::string& region = "RGN_NOBNDRY") {
  AUTO_TRACE();
  return standardDerivative<T, DIRECTION::X, DERIV::Standard>(f, outloc, method, region);
}

template <typename T>
T D2DX2(const T& f, CELL_LOC outloc = CELL_DEFAULT, const std::string& method = "DEFAULT",
        const std::string& region = "RGN_NOBNDRY") {
  AUTO_TRACE();
  return standardDerivative<T, DIRECTION::X, DERIV::StandardSecond>(f, outloc, method,
                                                                    region);
}

template <typename T>
T D4DX4(const T& f, CELL_LOC outloc = CELL_DEFAULT, const std::string& method = "DEFAULT",
        const std::string& region = "RGN_NOBNDRY") {
  AUTO_TRACE();
  return standardDerivative<T, DIRECTION::X, DERIV::StandardFourth>(f, outloc, method,
                                                                    region);
}

////////////// Y DERIVATIVE /////////////////

template <typename T>
T DDY(const T& f, CELL_LOC outloc = CELL_DEFAULT, const std::string& method = "DEFAULT",
      const std::string& region = "RGN_NOBNDRY") {
  AUTO_TRACE();
  if (f.hasParallelSlices()) {
    ASSERT1(f.getDirectionY() == YDirectionType::Standard);
    return standardDerivative<T, DIRECTION::YOrthogonal, DERIV::Standard>(f, outloc,
                                                                          method, region);
  } else {
    const bool is_unaligned = (f.getDirectionY() == YDirectionType::Standard);
    const T f_aligned = is_unaligned ? toFieldAligned(f, "RGN_NOX") : f;
    T result = standardDerivative<T, DIRECTION::Y, DERIV::Standard>(f_aligned, outloc,
                                                                    method, region);
    return is_unaligned ? fromFieldAligned(result, region) : result;
  }
}

template <typename T>
T D2DY2(const T& f, CELL_LOC outloc = CELL_DEFAULT, const std::string& method = "DEFAULT",
        const std::string& region = "RGN_NOBNDRY") {
  AUTO_TRACE();
  if (f.hasParallelSlices()) {
    ASSERT1(f.getDirectionY() == YDirectionType::Standard);
    return standardDerivative<T, DIRECTION::YOrthogonal, DERIV::StandardSecond>(
        f, outloc, method, region);
  } else {
    const bool is_unaligned = (f.getDirectionY() == YDirectionType::Standard);
    const T f_aligned = is_unaligned ? toFieldAligned(f, "RGN_NOX") : f;
    T result = standardDerivative<T, DIRECTION::Y, DERIV::StandardSecond>(
        f_aligned, outloc, method, region);
    return is_unaligned ? fromFieldAligned(result, region) : result;
  }
}

template <typename T>
T D4DY4(const T& f, CELL_LOC outloc = CELL_DEFAULT, const std::string& method = "DEFAULT",
        const std::string& region = "RGN_NOBNDRY") {
  AUTO_TRACE();
  if (f.hasParallelSlices()) {
    ASSERT1(f.getDirectionY() == YDirectionType::Standard);
    return standardDerivative<T, DIRECTION::YOrthogonal, DERIV::StandardFourth>(
        f, outloc, method, region);
  } else {
    const bool is_unaligned = (f.getDirectionY() == YDirectionType::Aligned);
    const T f_aligned = is_unaligned ? toFieldAligned(f, "RGN_NOX") : f;
    T result = standardDerivative<T, DIRECTION::Y, DERIV::StandardFourth>(
        f_aligned, outloc, method, region);
    return is_unaligned ? fromFieldAligned(result, region) : result;
  }
}

////////////// Z DERIVATIVE /////////////////
template <typename T>
T DDZ(const T& f, CELL_LOC outloc = CELL_DEFAULT, const std::string& method = "DEFAULT",
      const std::string& region = "RGN_NOBNDRY") {
  AUTO_TRACE();
  return standardDerivative<T, DIRECTION::Z, DERIV::Standard>(f, outloc, method, region);
}

template <typename T>
T D2DZ2(const T& f, CELL_LOC outloc = CELL_DEFAULT, const std::string& method = "DEFAULT",
        const std::string& region = "RGN_NOBNDRY") {
  AUTO_TRACE();
  return standardDerivative<T, DIRECTION::Z, DERIV::StandardSecond>(f, outloc, method,
                                                                    region);
}

template <typename T>
T D4DZ4(const T& f, CELL_LOC outloc = CELL_DEFAULT, const std::string& method = "DEFAULT",
        const std::string& region = "RGN_NOBNDRY") {
  AUTO_TRACE();
  return standardDerivative<T, DIRECTION::Z, DERIV::StandardFourth>(f, outloc, method,
                                                                    region);
}

////// ADVECTION AND FLUX OPERATORS

/// Advection operator in index space in [] direction
///
/// \f[
///   v \frac{d}{di} f
/// \f]
///
/// @param[in] v  The velocity in the Y direction
/// @param[in] f  The field being advected
/// @param[in] outloc The cell location where the result is desired. The default is the
/// same as \p f
/// @param[in] method  The differencing method to use
/// @param[in] region  The region of the grid for which the result is calculated.

////////////// X DERIVATIVE /////////////////

template <typename T>
T VDDX(const T& vel, const T& f, CELL_LOC outloc = CELL_DEFAULT,
       const std::string& method = "DEFAULT", const std::string& region = "RGN_NOBNDRY") {
  AUTO_TRACE();
  return flowDerivative<T, DIRECTION::X, DERIV::Upwind>(vel, f, outloc, method, region);
}

template <typename T>
T FDDX(const T& vel, const T& f, CELL_LOC outloc = CELL_DEFAULT,
       const std::string& method = "DEFAULT", const std::string& region = "RGN_NOBNDRY") {
  AUTO_TRACE();
  return flowDerivative<T, DIRECTION::X, DERIV::Flux>(vel, f, outloc, method, region);
}

////////////// Y DERIVATIVE /////////////////

template <typename T>
T VDDY(const T& vel, const T& f, CELL_LOC outloc = CELL_DEFAULT,
       const std::string& method = "DEFAULT", const std::string& region = "RGN_NOBNDRY") {
  AUTO_TRACE();
  const bool fHasParallelSlices = (f.hasParallelSlices());
  const bool velHasParallelSlices = (vel.hasParallelSlices());
  if (fHasParallelSlices && velHasParallelSlices) {
    ASSERT1(vel.getDirectionY() == YDirectionType::Standard);
    ASSERT1(f.getDirectionY() == YDirectionType::Standard);
    return flowDerivative<T, DIRECTION::YOrthogonal, DERIV::Upwind>(vel, f, outloc,
                                                                    method, region);
  } else {
    ASSERT2(f.getDirectionY() == vel.getDirectionY());
    const bool are_unaligned = ((f.getDirectionY() == YDirectionType::Standard)
                                and (vel.getDirectionY() == YDirectionType::Standard));

    const T f_aligned = are_unaligned ? toFieldAligned(f, "RGN_NOX") : f;
    const T vel_aligned = are_unaligned ? toFieldAligned(vel, "RGN_NOX") : vel;
    T result = flowDerivative<T, DIRECTION::Y, DERIV::Upwind>(vel_aligned, f_aligned,
                                                              outloc, method, region);
    return are_unaligned ? fromFieldAligned(result, region) : result;
  }
}

template <typename T>
T FDDY(const T& vel, const T& f, CELL_LOC outloc = CELL_DEFAULT,
       const std::string& method = "DEFAULT", const std::string& region = "RGN_NOBNDRY") {
  AUTO_TRACE();
  const bool fHasParallelSlices = (f.hasParallelSlices());
  const bool velHasParallelSlices = (vel.hasParallelSlices());
  if (fHasParallelSlices && velHasParallelSlices) {
    ASSERT1(vel.getDirectionY() == YDirectionType::Standard);
    ASSERT1(f.getDirectionY() == YDirectionType::Standard);
    return flowDerivative<T, DIRECTION::YOrthogonal, DERIV::Flux>(vel, f, outloc, method,
                                                                  region);
  } else {
    ASSERT2(f.getDirectionY() == vel.getDirectionY());
    const bool are_unaligned = ((f.getDirectionY() == YDirectionType::Standard)
                                and (vel.getDirectionY() == YDirectionType::Standard));

    const T f_aligned = are_unaligned ? toFieldAligned(f, "RGN_NOX") : f;
    const T vel_aligned = are_unaligned ? toFieldAligned(vel, "RGN_NOX") : vel;
    T result = flowDerivative<T, DIRECTION::Y, DERIV::Flux>(vel_aligned, f_aligned,
                                                            outloc, method, region);
    return are_unaligned ? fromFieldAligned(result, region) : result;
  }
}

////////////// Z DERIVATIVE /////////////////

template <typename T>
T VDDZ(const T& vel, const T& f, CELL_LOC outloc = CELL_DEFAULT,
       const std::string& method = "DEFAULT", const std::string& region = "RGN_NOBNDRY") {
  AUTO_TRACE();
  return flowDerivative<T, DIRECTION::Z, DERIV::Upwind>(vel, f, outloc, method, region);
}

template <typename T>
T FDDZ(const T& vel, const T& f, CELL_LOC outloc = CELL_DEFAULT,
       const std::string& method = "DEFAULT", const std::string& region = "RGN_NOBNDRY") {
  AUTO_TRACE();
  return flowDerivative<T, DIRECTION::Z, DERIV::Flux>(vel, f, outloc, method, region);
}

} // Namespace index
} // Namespace derivatives
} // Namespace bout
#endif
