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

#include <string>
#include <limits>

/// Size of real numbers
typedef double BoutReal;

/// Quiet NaN
const BoutReal BoutNaN = std::numeric_limits<BoutReal>::quiet_NaN();

#define ENUMSTR(val) {val, #val}

/// 4 possible variable locations. Default is for passing to functions
enum CELL_LOC {CELL_DEFAULT=0, CELL_CENTRE=1, CELL_CENTER=1, CELL_XLOW=2, CELL_YLOW=3, CELL_ZLOW=4, CELL_VSHIFT=5};

const std::string& CELL_LOC_STRING(CELL_LOC location);

/// Differential methods. Both central and upwind
enum DIFF_METHOD {DIFF_DEFAULT, DIFF_U1, DIFF_U2, DIFF_C2, DIFF_W2, DIFF_W3, DIFF_C4, DIFF_U3, DIFF_FFT, DIFF_SPLIT, DIFF_S2};

const std::string& DIFF_METHOD_STRING(DIFF_METHOD location);

/// Specify grid region for looping
enum REGION {RGN_ALL, RGN_NOBNDRY, RGN_NOX, RGN_NOY, RGN_NOZ};

const std::string& REGION_STRING(REGION region);

/// To identify particular directions (in index space)
enum class DIRECTION { X = 0, Y = 1, Z = 3, YAligned = 4, YOrthogonal = 5, Special = 6, Null = 7 };

const std::string& DIRECTION_STRING(DIRECTION direction);

/// To identify valid staggering combinations
enum class STAGGER { None = 0, C2L = 1, L2C = 2};

const std::string& STAGGER_STRING(STAGGER stagger);

/// To identify types of derivative method combinations
enum class DERIV { Standard = 0, StandardSecond = 1, StandardFourth = 2,
		   Upwind = 3, Flux = 4 };

const std::string& DERIV_STRING(DERIV deriv);

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

bool compatibleDirections(DIRECTION d1, DIRECTION d2);
bool isXDirectionType(DIRECTION x);
bool isYDirectionType(DIRECTION y);
bool isZDirectionType(DIRECTION z);

#endif // __BOUT_TYPES_H__
