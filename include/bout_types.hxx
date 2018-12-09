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

#include <vector>
#include <string>
#include <map>
#include <limits>

/// Size of real numbers
typedef double BoutReal;

/// Quiet NaN
const BoutReal BoutNaN = std::numeric_limits<BoutReal>::quiet_NaN();

#define ENUMSTR(val) {val, #val}

/// 4 possible variable locations. Default is for passing to functions
enum CELL_LOC {CELL_DEFAULT=0, CELL_CENTRE=1, CELL_CENTER=1, CELL_XLOW=2, CELL_YLOW=3, CELL_ZLOW=4, CELL_VSHIFT=5};

const std::map<CELL_LOC, std::string> CELL_LOCtoString = {
  ENUMSTR(CELL_DEFAULT),
  ENUMSTR(CELL_CENTRE),
  ENUMSTR(CELL_XLOW),
  ENUMSTR(CELL_YLOW),
  ENUMSTR(CELL_ZLOW),
  ENUMSTR(CELL_VSHIFT)
};

inline const std::string& CELL_LOC_STRING(CELL_LOC location) {
  return CELL_LOCtoString.at(location);
}

/// Differential methods. Both central and upwind
enum DIFF_METHOD {DIFF_DEFAULT, DIFF_U1, DIFF_U2, DIFF_C2, DIFF_W2, DIFF_W3, DIFF_C4, DIFF_U3, DIFF_FFT, DIFF_SPLIT, DIFF_S2};

const std::map<DIFF_METHOD, std::string> DIFF_METHODtoString = {
  {DIFF_DEFAULT, "DEFAULT"},
  {DIFF_U1, "U1"}, {DIFF_U2, "U2"}, {DIFF_U3, "U3"},
  {DIFF_C2, "C2"}, {DIFF_C4, "C4"}, {DIFF_S2, "S2"},
  {DIFF_W2, "W2"}, {DIFF_W3, "W3"}, {DIFF_FFT, "FFT"},
  {DIFF_SPLIT, "SPLIT"}
};

inline const std::string& DIFF_METHOD_STRING(DIFF_METHOD location) {
  return DIFF_METHODtoString.at(location);
}

/// Specify grid region for looping
enum REGION {RGN_ALL, RGN_NOBNDRY, RGN_NOX, RGN_NOY, RGN_NOZ};

const std::map<REGION, std::string> REGIONtoString = {
  ENUMSTR(RGN_ALL),
  ENUMSTR(RGN_NOBNDRY),
  ENUMSTR(RGN_NOX),
  ENUMSTR(RGN_NOY),
  ENUMSTR(RGN_NOZ)
};

inline const std::string& REGION_STRING(REGION region) {
  return REGIONtoString.at(region);
}

/// To identify particular directions (in index space)
enum class DIRECTION { X = 0, Y = 1, Z = 3, YAligned = 4, YOrthogonal = 5 };

const std::map<DIRECTION, std::string> DIRECTIONtoString = {
  {DIRECTION::X, "X"},
  {DIRECTION::Y, "Y"},
  {DIRECTION::Z, "Z"},
  {DIRECTION::YAligned, "Y - field aligned"},
  {DIRECTION::YOrthogonal, "Y - orthogonal"}

};

inline const std::string& DIRECTION_STRING(DIRECTION direction) {
  return DIRECTIONtoString.at(direction);
}

/// To identify valid staggering combinations
enum class STAGGER { None = 0, C2L = 1, L2C = 2};

const std::map<STAGGER, std::string> STAGGERtoString = {
  {STAGGER::None, "No staggering"},
  {STAGGER::C2L, "Centre to Low"},
  {STAGGER::L2C, "Low to Centre"}

};

inline const std::string& STAGGER_STRING(STAGGER stagger) {
  return STAGGERtoString.at(stagger);
}

/// To identify types of derivative method combinations
enum class DERIV { Standard = 0, StandardSecond = 1, StandardFourth = 2,
		   Upwind = 3, Flux = 4 };

static std::map<DERIV, std::string> DERIVtoString = {
  {DERIV::Standard, "Standard"},
  {DERIV::StandardSecond, "Standard -- second order"},
  {DERIV::StandardFourth, "Standard -- fourth order"},
  {DERIV::Upwind, "Upwind"},
  {DERIV::Flux, "Flux"}  
};

inline const std::string& DERIV_STRING(DERIV deriv) {
  return DERIVtoString.at(deriv);
}

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
typedef BoutReal (*FuncPtr)(BoutReal t, BoutReal x, BoutReal y, BoutReal z);

#endif // __BOUT_TYPES_H__
