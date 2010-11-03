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

using std::string; // String class 
using std::vector;

typedef double BoutReal;

typedef vector<BoutReal> rvec;  // Vector of BoutReals

/// 4 possible variable locations. Default is for passing to functions
enum CELL_LOC {CELL_DEFAULT, CELL_CENTRE, CELL_XLOW, CELL_YLOW, CELL_ZLOW, CELL_VSHIFT};

/// Differential methods. Both central and upwind
enum DIFF_METHOD {DIFF_DEFAULT, DIFF_U1, DIFF_C2, DIFF_W2, DIFF_W3, DIFF_C4, DIFF_U4, DIFF_FFT, DIFF_SPLIT, DIFF_NND};

/// Specify grid region for looping
enum REGION {RGN_ALL, RGN_NOBNDRY, RGN_NOX, RGN_NOY, RGN_NOZ};

#endif // __BOUT_TYPES_H__
