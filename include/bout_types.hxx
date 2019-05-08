/**************************************************************************
 * Copyright 2010 B.D.Dudson, S.Farley, M.V.Umansky, X.Q.Xu
 *
 * Contact Ben Dudson, bd512@york.ac.uk
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
 **************************************************************************/

#ifndef __BOUT_TYPES_H__
#define __BOUT_TYPES_H__

#include "bout/deprecated.hxx"

#include <limits>
#include <string>

/// Size of real numbers
typedef double BoutReal;

/// Quiet NaN
const BoutReal BoutNaN = std::numeric_limits<BoutReal>::quiet_NaN();

#define ENUMSTR(val) {val, #val}
#define STRENUM(val) {#val, val}

/// 4 possible variable locations. Default is for passing to functions
enum CELL_LOC {CELL_DEFAULT=0, CELL_CENTRE=1, CELL_CENTER=1, CELL_XLOW=2, CELL_YLOW=3, CELL_ZLOW=4, CELL_VSHIFT=5};

std::string toString(CELL_LOC location);
CELL_LOC CELL_LOCFromString(std::string location_string);
DEPRECATED(inline std::string CELL_LOC_STRING(CELL_LOC location)) {
  return toString(location);
}

/// Differential methods. Both central and upwind
enum DIFF_METHOD {DIFF_DEFAULT, DIFF_U1, DIFF_U2, DIFF_C2, DIFF_W2, DIFF_W3, DIFF_C4, DIFF_U3, DIFF_FFT, DIFF_SPLIT, DIFF_S2};

std::string toString(DIFF_METHOD location);
DEPRECATED(inline std::string DIFF_METHOD_STRING(DIFF_METHOD location)) {
  return toString(location);
}

/// Specify grid region for looping
enum REGION { RGN_ALL, RGN_NOBNDRY, RGN_NOX, RGN_NOY, RGN_NOZ };

std::string toString(REGION region);
DEPRECATED(inline std::string REGION_STRING(REGION region)) { return toString(region); }

/// To identify particular directions (in index space):
///   - X, Y, Z are the coordinate directions
///   - YAligned is a special case of Y, indicating a field-aligned grid, where
///     the x- and z- axes are not necessarily orthogonal
///   - YOrthogonal is a special case of Y, indicating a grid where the x and z
///     axes are orthogonal but the y-direction is not necessarily
///     field-aligned
enum class DIRECTION { X, Y, Z, YAligned, YOrthogonal };

std::string toString(DIRECTION direction);
DEPRECATED(inline std::string DIRECTION_STRING(DIRECTION direction)) {
  return toString(direction);
}

/// Identify kind of a field's y-direction
/// - Standard is the default for the Mesh/Coordinates/ParallelTransform
/// - Aligned indicates that the field has been transformed to field-aligned
///   coordinates
enum class YDirectionType { Standard, Aligned };

std::string toString(YDirectionType d);
YDirectionType YDirectionTypeFromString(std::string y_direction_string);

/// Identify kind of a field's z-direction
/// - Standard is the default
/// - Average indicates that the field represents an average over the
///   z-direction, rather than having a particular z-position (i.e. is a
///   Field2D)
enum class ZDirectionType { Standard, Average };

std::string toString(ZDirectionType d);
ZDirectionType ZDirectionTypeFromString(std::string z_direction_string);

/// Container for direction types
struct DirectionTypes {
  YDirectionType y;
  ZDirectionType z;
};

/// Check whether direction types are compatible, so two fields with attributes
/// d1 and d2 respectively can be added, subtracted, etc.
bool areDirectionsCompatible(const DirectionTypes& d1, const DirectionTypes& d2);

void swap(const DirectionTypes& first, const DirectionTypes& second);

/// To identify valid staggering combinations
enum class STAGGER { None, C2L, L2C };

std::string toString(STAGGER stagger);
DEPRECATED(inline std::string STAGGER_STRING(STAGGER stagger)) {
  return toString(stagger);
}

/// To identify types of derivative method combinations
enum class DERIV { Standard, StandardSecond, StandardFourth, Upwind, Flux };

std::string toString(DERIV deriv);
DEPRECATED(inline std::string DERIV_STRING(DERIV deriv)) { return toString(deriv); }

// A small struct that can be used to wrap a specific enum value, giving
// it a unique type that can be passed as a valid type to templates and
// which can be inspected to provide the actual value of the enum
template<typename T, T val>
struct enumWrapper {
  using type = T;
  static const type value = val;
  T lookup(){return val;};
};

/// Boundary condition function
using FuncPtr = BoutReal(*)(BoutReal t, BoutReal x, BoutReal y, BoutReal z);

#endif // __BOUT_TYPES_H__
