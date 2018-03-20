/*!************************************************************************
 * \file derivs.hxx
 *
 * Basic differential functions
 *
 **************************************************************************
 * Copyright 2010,2017
 *    B.D.Dudson, S.Farley, M.V.Umansky, X.Q.Xu, D. Schwörer
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

#ifndef __DERIVS_H__
#define __DERIVS_H__

#include "field2d.hxx"
#include "field3d.hxx"
#include "vector2d.hxx"
#include "vector3d.hxx"

#include "bout_types.hxx" // See this for code

////////// FIRST DERIVATIVES //////////

/// Calculate first partial derivative in X
///
///   $\partial / \partial x$
///
/// @param[in] f       The field to be differentiated
/// @param[in] outloc  The cell location where the result is desired. If
///                    staggered grids is not enabled then this has no effect
///                    If not given, defaults to CELL_DEFAULT
/// @param[in] method  Differencing method to use. This overrides the default
///                    If not given, defaults to DIFF_DEFAULT
/// @param[in] region  What region is expected to be calculated
///                    If not given, defaults to RGN_NOBNDRY
///
///
///
const Field3D DDX(const Field3D &f, CELL_LOC outloc = CELL_DEFAULT,
                  REGION region = RGN_NOBNDRY, DIFF_METHOD method = DIFF_DEFAULT);

DEPRECATED(const Field3D DDX(const Field3D &f, DIFF_METHOD method, CELL_LOC outloc,
                             REGION region));
inline const Field3D DDX(const Field3D &f, DIFF_METHOD method,
                         CELL_LOC outloc = CELL_DEFAULT, REGION region = RGN_NOBNDRY) {
  return DDX(f, outloc, region, method);
}

DEPRECATED(const Field3D DDX(const Field3D &f, CELL_LOC outloc, DIFF_METHOD method,
                             REGION region));
inline const Field3D DDX(const Field3D &f, CELL_LOC outloc, DIFF_METHOD method,
                         REGION region = RGN_NOBNDRY) {
  return DDX(f, outloc, region, method);
}

/// Calculate first partial derivative in X
///
///   $\partial / \partial x$
///
/// @param[in] f       The field to be differentiated
/// @param[in] outloc  The cell location where the result is desired. If
///                    staggered grids is not enabled then this has no effect
///                    If not given, defaults to CELL_DEFAULT
/// @param[in] method  Differencing method to use. This overrides the default
///                    If not given, defaults to DIFF_DEFAULT
/// @param[in] region  What region is expected to be calculated
///                    If not given, defaults to RGN_NOBNDRY
///
///
///
const Field2D DDX(const Field2D &f, CELL_LOC outloc = CELL_DEFAULT,
                  REGION region = RGN_NOBNDRY, DIFF_METHOD method = DIFF_DEFAULT);

DEPRECATED(const Field2D DDX(const Field2D &f, DIFF_METHOD method, CELL_LOC outloc,
                             REGION region));
inline const Field2D DDX(const Field2D &f, DIFF_METHOD method,
                         CELL_LOC outloc = CELL_DEFAULT, REGION region = RGN_NOBNDRY) {
  return DDX(f, outloc, region, method);
}

DEPRECATED(const Field2D DDX(const Field2D &f, CELL_LOC outloc, DIFF_METHOD method,
                             REGION region));
inline const Field2D DDX(const Field2D &f, CELL_LOC outloc, DIFF_METHOD method,
                         REGION region = RGN_NOBNDRY) {
  return DDX(f, outloc, region, method);
}

/// Calculate first partial derivative in Y
///
///   $\partial / \partial y$
///
/// @param[in] f       The field to be differentiated
/// @param[in] outloc  The cell location where the result is desired. If
///                    staggered grids is not enabled then this has no effect
///                    If not given, defaults to CELL_DEFAULT
/// @param[in] method  Differencing method to use. This overrides the default
///                    If not given, defaults to DIFF_DEFAULT
/// @param[in] region  What region is expected to be calculated
///                    If not given, defaults to RGN_NOBNDRY
///
///
///
const Field3D DDY(const Field3D &f, CELL_LOC outloc = CELL_DEFAULT,
                  REGION region = RGN_NOBNDRY, DIFF_METHOD method = DIFF_DEFAULT);

DEPRECATED(const Field3D DDY(const Field3D &f, DIFF_METHOD method, CELL_LOC outloc,
                             REGION region));
inline const Field3D DDY(const Field3D &f, DIFF_METHOD method,
                         CELL_LOC outloc = CELL_DEFAULT, REGION region = RGN_NOBNDRY) {
  return DDY(f, outloc, region, method);
}

DEPRECATED(const Field3D DDY(const Field3D &f, CELL_LOC outloc, DIFF_METHOD method,
                             REGION region));
inline const Field3D DDY(const Field3D &f, CELL_LOC outloc, DIFF_METHOD method,
                         REGION region = RGN_NOBNDRY) {
  return DDY(f, outloc, region, method);
}

/// Calculate first partial derivative in Y
///
///   $\partial / \partial y$
///
/// @param[in] f       The field to be differentiated
/// @param[in] outloc  The cell location where the result is desired. If
///                    staggered grids is not enabled then this has no effect
///                    If not given, defaults to CELL_DEFAULT
/// @param[in] method  Differencing method to use. This overrides the default
///                    If not given, defaults to DIFF_DEFAULT
/// @param[in] region  What region is expected to be calculated
///                    If not given, defaults to RGN_NOBNDRY
///
///
///
const Field2D DDY(const Field2D &f, CELL_LOC outloc = CELL_DEFAULT,
                  REGION region = RGN_NOBNDRY, DIFF_METHOD method = DIFF_DEFAULT);

DEPRECATED(const Field2D DDY(const Field2D &f, DIFF_METHOD method, CELL_LOC outloc,
                             REGION region));
inline const Field2D DDY(const Field2D &f, DIFF_METHOD method,
                         CELL_LOC outloc = CELL_DEFAULT, REGION region = RGN_NOBNDRY) {
  return DDY(f, outloc, region, method);
}

DEPRECATED(const Field2D DDY(const Field2D &f, CELL_LOC outloc, DIFF_METHOD method,
                             REGION region));
inline const Field2D DDY(const Field2D &f, CELL_LOC outloc, DIFF_METHOD method,
                         REGION region = RGN_NOBNDRY) {
  return DDY(f, outloc, region, method);
}

/// Calculate first partial derivative in Z
///
///   $\partial / \partial z$
///
/// @param[in] f       The field to be differentiated
/// @param[in] outloc  The cell location where the result is desired. If
///                    staggered grids is not enabled then this has no effect
///                    If not given, defaults to CELL_DEFAULT
/// @param[in] method  Differencing method to use. This overrides the default
///                    If not given, defaults to DIFF_DEFAULT
/// @param[in] region  What region is expected to be calculated
///                    If not given, defaults to RGN_NOBNDRY
///
///
///
const Field3D DDZ(const Field3D &f, CELL_LOC outloc = CELL_DEFAULT,
                  REGION region = RGN_NOBNDRY, DIFF_METHOD method = DIFF_DEFAULT);

DEPRECATED(const Field3D DDZ(const Field3D &f, DIFF_METHOD method, CELL_LOC outloc,
                             REGION region));
inline const Field3D DDZ(const Field3D &f, DIFF_METHOD method,
                         CELL_LOC outloc = CELL_DEFAULT, REGION region = RGN_NOBNDRY) {
  return DDZ(f, outloc, region, method);
}

DEPRECATED(const Field3D DDZ(const Field3D &f, CELL_LOC outloc, DIFF_METHOD method,
                             REGION region));
inline const Field3D DDZ(const Field3D &f, CELL_LOC outloc, DIFF_METHOD method,
                         REGION region = RGN_NOBNDRY) {
  return DDZ(f, outloc, region, method);
}

/// Calculate first partial derivative in Z
///
///   $\partial / \partial z$
///
/// @param[in] f       The field to be differentiated
/// @param[in] outloc  The cell location where the result is desired. If
///                    staggered grids is not enabled then this has no effect
///                    If not given, defaults to CELL_DEFAULT
/// @param[in] method  Differencing method to use. This overrides the default
///                    If not given, defaults to DIFF_DEFAULT
/// @param[in] region  What region is expected to be calculated
///                    If not given, defaults to RGN_NOBNDRY
///
///
///
const Field2D DDZ(const Field2D &f, CELL_LOC outloc = CELL_DEFAULT,
                  REGION region = RGN_NOBNDRY, DIFF_METHOD method = DIFF_DEFAULT);

DEPRECATED(const Field2D DDZ(const Field2D &f, DIFF_METHOD method, CELL_LOC outloc,
                             REGION region));
inline const Field2D DDZ(const Field2D &f, DIFF_METHOD method,
                         CELL_LOC outloc = CELL_DEFAULT, REGION region = RGN_NOBNDRY) {
  return DDZ(f, outloc, region, method);
}

DEPRECATED(const Field2D DDZ(const Field2D &f, CELL_LOC outloc, DIFF_METHOD method,
                             REGION region));
inline const Field2D DDZ(const Field2D &f, CELL_LOC outloc, DIFF_METHOD method,
                         REGION region = RGN_NOBNDRY) {
  return DDZ(f, outloc, region, method);
}
////////// SECOND DERIVATIVES //////////

/// Calculate second partial derivative in X
///
///   $\partial^2 / \partial x^2$
///
/// @param[in] f       The field to be differentiated
/// @param[in] outloc  The cell location where the result is desired. If
///                    staggered grids is not enabled then this has no effect
///                    If not given, defaults to CELL_DEFAULT
/// @param[in] method  Differencing method to use. This overrides the default
///                    If not given, defaults to DIFF_DEFAULT
/// @param[in] region  What region is expected to be calculated
///                    If not given, defaults to RGN_NOBNDRY
///
///
///
const Field3D D2DX2(const Field3D &f, CELL_LOC outloc = CELL_DEFAULT,
                    REGION region = RGN_NOBNDRY, DIFF_METHOD method = DIFF_DEFAULT);

DEPRECATED(const Field3D D2DX2(const Field3D &f, DIFF_METHOD method, CELL_LOC outloc,
                               REGION region));
inline const Field3D D2DX2(const Field3D &f, DIFF_METHOD method,
                           CELL_LOC outloc = CELL_DEFAULT, REGION region = RGN_NOBNDRY) {
  return D2DX2(f, outloc, region, method);
}

DEPRECATED(const Field3D D2DX2(const Field3D &f, CELL_LOC outloc, DIFF_METHOD method,
                               REGION region));
inline const Field3D D2DX2(const Field3D &f, CELL_LOC outloc, DIFF_METHOD method,
                           REGION region = RGN_NOBNDRY) {
  return D2DX2(f, outloc, region, method);
}

/// Calculate second partial derivative in X
///
///   $\partial^2 / \partial x^2$
///
/// @param[in] f       The field to be differentiated
/// @param[in] outloc  The cell location where the result is desired. If
///                    staggered grids is not enabled then this has no effect
///                    If not given, defaults to CELL_DEFAULT
/// @param[in] method  Differencing method to use. This overrides the default
///                    If not given, defaults to DIFF_DEFAULT
/// @param[in] region  What region is expected to be calculated
///                    If not given, defaults to RGN_NOBNDRY
///
///
///
const Field2D D2DX2(const Field2D &f, CELL_LOC outloc = CELL_DEFAULT,
                    REGION region = RGN_NOBNDRY, DIFF_METHOD method = DIFF_DEFAULT);

DEPRECATED(const Field2D D2DX2(const Field2D &f, DIFF_METHOD method, CELL_LOC outloc,
                               REGION region));
inline const Field2D D2DX2(const Field2D &f, DIFF_METHOD method,
                           CELL_LOC outloc = CELL_DEFAULT, REGION region = RGN_NOBNDRY) {
  return D2DX2(f, outloc, region, method);
}

DEPRECATED(const Field2D D2DX2(const Field2D &f, CELL_LOC outloc, DIFF_METHOD method,
                               REGION region));
inline const Field2D D2DX2(const Field2D &f, CELL_LOC outloc, DIFF_METHOD method,
                           REGION region = RGN_NOBNDRY) {
  return D2DX2(f, outloc, region, method);
}

/// Calculate second partial derivative in Y
///
///   $\partial^2 / \partial y^2$
///
/// @param[in] f       The field to be differentiated
/// @param[in] outloc  The cell location where the result is desired. If
///                    staggered grids is not enabled then this has no effect
///                    If not given, defaults to CELL_DEFAULT
/// @param[in] method  Differencing method to use. This overrides the default
///                    If not given, defaults to DIFF_DEFAULT
/// @param[in] region  What region is expected to be calculated
///                    If not given, defaults to RGN_NOBNDRY
///
///
///
const Field3D D2DY2(const Field3D &f, CELL_LOC outloc = CELL_DEFAULT,
                    REGION region = RGN_NOBNDRY, DIFF_METHOD method = DIFF_DEFAULT);

DEPRECATED(const Field3D D2DY2(const Field3D &f, DIFF_METHOD method, CELL_LOC outloc,
                               REGION region));
inline const Field3D D2DY2(const Field3D &f, DIFF_METHOD method,
                           CELL_LOC outloc = CELL_DEFAULT, REGION region = RGN_NOBNDRY) {
  return D2DY2(f, outloc, region, method);
}

DEPRECATED(const Field3D D2DY2(const Field3D &f, CELL_LOC outloc, DIFF_METHOD method,
                               REGION region));
inline const Field3D D2DY2(const Field3D &f, CELL_LOC outloc, DIFF_METHOD method,
                           REGION region = RGN_NOBNDRY) {
  return D2DY2(f, outloc, region, method);
}

/// Calculate second partial derivative in Y
///
///   $\partial^2 / \partial y^2$
///
/// @param[in] f       The field to be differentiated
/// @param[in] outloc  The cell location where the result is desired. If
///                    staggered grids is not enabled then this has no effect
///                    If not given, defaults to CELL_DEFAULT
/// @param[in] method  Differencing method to use. This overrides the default
///                    If not given, defaults to DIFF_DEFAULT
/// @param[in] region  What region is expected to be calculated
///                    If not given, defaults to RGN_NOBNDRY
///
///
///
const Field2D D2DY2(const Field2D &f, CELL_LOC outloc = CELL_DEFAULT,
                    REGION region = RGN_NOBNDRY, DIFF_METHOD method = DIFF_DEFAULT);

DEPRECATED(const Field2D D2DY2(const Field2D &f, DIFF_METHOD method, CELL_LOC outloc,
                               REGION region));
inline const Field2D D2DY2(const Field2D &f, DIFF_METHOD method,
                           CELL_LOC outloc = CELL_DEFAULT, REGION region = RGN_NOBNDRY) {
  return D2DY2(f, outloc, region, method);
}

DEPRECATED(const Field2D D2DY2(const Field2D &f, CELL_LOC outloc, DIFF_METHOD method,
                               REGION region));
inline const Field2D D2DY2(const Field2D &f, CELL_LOC outloc, DIFF_METHOD method,
                           REGION region = RGN_NOBNDRY) {
  return D2DY2(f, outloc, region, method);
}

/// Calculate second partial derivative in Z
///
///   $\partial^2 / \partial z^2$
///
/// @param[in] f       The field to be differentiated
/// @param[in] outloc  The cell location where the result is desired. If
///                    staggered grids is not enabled then this has no effect
///                    If not given, defaults to CELL_DEFAULT
/// @param[in] method  Differencing method to use. This overrides the default
///                    If not given, defaults to DIFF_DEFAULT
/// @param[in] region  What region is expected to be calculated
///                    If not given, defaults to RGN_NOBNDRY
///
///
///
const Field3D D2DZ2(const Field3D &f, CELL_LOC outloc = CELL_DEFAULT,
                    REGION region = RGN_NOBNDRY, DIFF_METHOD method = DIFF_DEFAULT);

DEPRECATED(const Field3D D2DZ2(const Field3D &f, DIFF_METHOD method, CELL_LOC outloc,
                               REGION region));
inline const Field3D D2DZ2(const Field3D &f, DIFF_METHOD method,
                           CELL_LOC outloc = CELL_DEFAULT, REGION region = RGN_NOBNDRY) {
  return D2DZ2(f, outloc, region, method);
}

DEPRECATED(const Field3D D2DZ2(const Field3D &f, CELL_LOC outloc, DIFF_METHOD method,
                               REGION region));
inline const Field3D D2DZ2(const Field3D &f, CELL_LOC outloc, DIFF_METHOD method,
                           REGION region = RGN_NOBNDRY) {
  return D2DZ2(f, outloc, region, method);
}

/// Calculate second partial derivative in Z
///
///   $\partial^2 / \partial z^2$
///
/// @param[in] f       The field to be differentiated
/// @param[in] outloc  The cell location where the result is desired. If
///                    staggered grids is not enabled then this has no effect
///                    If not given, defaults to CELL_DEFAULT
/// @param[in] method  Differencing method to use. This overrides the default
///                    If not given, defaults to DIFF_DEFAULT
/// @param[in] region  What region is expected to be calculated
///                    If not given, defaults to RGN_NOBNDRY
///
///
///
const Field2D D2DZ2(const Field2D &f, CELL_LOC outloc = CELL_DEFAULT,
                    REGION region = RGN_NOBNDRY, DIFF_METHOD method = DIFF_DEFAULT);

DEPRECATED(const Field2D D2DZ2(const Field2D &f, DIFF_METHOD method, CELL_LOC outloc,
                               REGION region));
inline const Field2D D2DZ2(const Field2D &f, DIFF_METHOD method,
                           CELL_LOC outloc = CELL_DEFAULT, REGION region = RGN_NOBNDRY) {
  return D2DZ2(f, outloc, region, method);
}

DEPRECATED(const Field2D D2DZ2(const Field2D &f, CELL_LOC outloc, DIFF_METHOD method,
                               REGION region));
inline const Field2D D2DZ2(const Field2D &f, CELL_LOC outloc, DIFF_METHOD method,
                           REGION region = RGN_NOBNDRY) {
  return D2DZ2(f, outloc, region, method);
}
////////// FORTH DERIVATIVES //////////

/// Calculate forth partial derivative in X
///
///   $\partial^4 / \partial x^4$
///
/// @param[in] f       The field to be differentiated
/// @param[in] outloc  The cell location where the result is desired. If
///                    staggered grids is not enabled then this has no effect
///                    If not given, defaults to CELL_DEFAULT
/// @param[in] method  Differencing method to use. This overrides the default
///                    If not given, defaults to DIFF_DEFAULT
/// @param[in] region  What region is expected to be calculated
///                    If not given, defaults to RGN_NOBNDRY
///
///
///
const Field3D D4DX4(const Field3D &f, CELL_LOC outloc = CELL_DEFAULT,
                    REGION region = RGN_NOBNDRY, DIFF_METHOD method = DIFF_DEFAULT);

DEPRECATED(const Field3D D4DX4(const Field3D &f, DIFF_METHOD method, CELL_LOC outloc,
                               REGION region));
inline const Field3D D4DX4(const Field3D &f, DIFF_METHOD method,
                           CELL_LOC outloc = CELL_DEFAULT, REGION region = RGN_NOBNDRY) {
  return D4DX4(f, outloc, region, method);
}

DEPRECATED(const Field3D D4DX4(const Field3D &f, CELL_LOC outloc, DIFF_METHOD method,
                               REGION region));
inline const Field3D D4DX4(const Field3D &f, CELL_LOC outloc, DIFF_METHOD method,
                           REGION region = RGN_NOBNDRY) {
  return D4DX4(f, outloc, region, method);
}

/// Calculate forth partial derivative in X
///
///   $\partial^4 / \partial x^4$
///
/// @param[in] f       The field to be differentiated
/// @param[in] outloc  The cell location where the result is desired. If
///                    staggered grids is not enabled then this has no effect
///                    If not given, defaults to CELL_DEFAULT
/// @param[in] method  Differencing method to use. This overrides the default
///                    If not given, defaults to DIFF_DEFAULT
/// @param[in] region  What region is expected to be calculated
///                    If not given, defaults to RGN_NOBNDRY
///
///
///
const Field2D D4DX4(const Field2D &f, CELL_LOC outloc = CELL_DEFAULT,
                    REGION region = RGN_NOBNDRY, DIFF_METHOD method = DIFF_DEFAULT);

DEPRECATED(const Field2D D4DX4(const Field2D &f, DIFF_METHOD method, CELL_LOC outloc,
                               REGION region));
inline const Field2D D4DX4(const Field2D &f, DIFF_METHOD method,
                           CELL_LOC outloc = CELL_DEFAULT, REGION region = RGN_NOBNDRY) {
  return D4DX4(f, outloc, region, method);
}

DEPRECATED(const Field2D D4DX4(const Field2D &f, CELL_LOC outloc, DIFF_METHOD method,
                               REGION region));
inline const Field2D D4DX4(const Field2D &f, CELL_LOC outloc, DIFF_METHOD method,
                           REGION region = RGN_NOBNDRY) {
  return D4DX4(f, outloc, region, method);
}

/// Calculate forth partial derivative in Y
///
///   $\partial^4 / \partial y^4$
///
/// @param[in] f       The field to be differentiated
/// @param[in] outloc  The cell location where the result is desired. If
///                    staggered grids is not enabled then this has no effect
///                    If not given, defaults to CELL_DEFAULT
/// @param[in] method  Differencing method to use. This overrides the default
///                    If not given, defaults to DIFF_DEFAULT
/// @param[in] region  What region is expected to be calculated
///                    If not given, defaults to RGN_NOBNDRY
///
///
///
const Field3D D4DY4(const Field3D &f, CELL_LOC outloc = CELL_DEFAULT,
                    REGION region = RGN_NOBNDRY, DIFF_METHOD method = DIFF_DEFAULT);

DEPRECATED(const Field3D D4DY4(const Field3D &f, DIFF_METHOD method, CELL_LOC outloc,
                               REGION region));
inline const Field3D D4DY4(const Field3D &f, DIFF_METHOD method,
                           CELL_LOC outloc = CELL_DEFAULT, REGION region = RGN_NOBNDRY) {
  return D4DY4(f, outloc, region, method);
}

DEPRECATED(const Field3D D4DY4(const Field3D &f, CELL_LOC outloc, DIFF_METHOD method,
                               REGION region));
inline const Field3D D4DY4(const Field3D &f, CELL_LOC outloc, DIFF_METHOD method,
                           REGION region = RGN_NOBNDRY) {
  return D4DY4(f, outloc, region, method);
}

/// Calculate forth partial derivative in Y
///
///   $\partial^4 / \partial y^4$
///
/// @param[in] f       The field to be differentiated
/// @param[in] outloc  The cell location where the result is desired. If
///                    staggered grids is not enabled then this has no effect
///                    If not given, defaults to CELL_DEFAULT
/// @param[in] method  Differencing method to use. This overrides the default
///                    If not given, defaults to DIFF_DEFAULT
/// @param[in] region  What region is expected to be calculated
///                    If not given, defaults to RGN_NOBNDRY
///
///
///
const Field2D D4DY4(const Field2D &f, CELL_LOC outloc = CELL_DEFAULT,
                    REGION region = RGN_NOBNDRY, DIFF_METHOD method = DIFF_DEFAULT);

DEPRECATED(const Field2D D4DY4(const Field2D &f, DIFF_METHOD method, CELL_LOC outloc,
                               REGION region));
inline const Field2D D4DY4(const Field2D &f, DIFF_METHOD method,
                           CELL_LOC outloc = CELL_DEFAULT, REGION region = RGN_NOBNDRY) {
  return D4DY4(f, outloc, region, method);
}

DEPRECATED(const Field2D D4DY4(const Field2D &f, CELL_LOC outloc, DIFF_METHOD method,
                               REGION region));
inline const Field2D D4DY4(const Field2D &f, CELL_LOC outloc, DIFF_METHOD method,
                           REGION region = RGN_NOBNDRY) {
  return D4DY4(f, outloc, region, method);
}

/// Calculate forth partial derivative in Z
///
///   $\partial^4 / \partial z^4$
///
/// @param[in] f       The field to be differentiated
/// @param[in] outloc  The cell location where the result is desired. If
///                    staggered grids is not enabled then this has no effect
///                    If not given, defaults to CELL_DEFAULT
/// @param[in] method  Differencing method to use. This overrides the default
///                    If not given, defaults to DIFF_DEFAULT
/// @param[in] region  What region is expected to be calculated
///                    If not given, defaults to RGN_NOBNDRY
///
///
///
const Field3D D4DZ4(const Field3D &f, CELL_LOC outloc = CELL_DEFAULT,
                    REGION region = RGN_NOBNDRY, DIFF_METHOD method = DIFF_DEFAULT);

DEPRECATED(const Field3D D4DZ4(const Field3D &f, DIFF_METHOD method, CELL_LOC outloc,
                               REGION region));
inline const Field3D D4DZ4(const Field3D &f, DIFF_METHOD method,
                           CELL_LOC outloc = CELL_DEFAULT, REGION region = RGN_NOBNDRY) {
  return D4DZ4(f, outloc, region, method);
}

DEPRECATED(const Field3D D4DZ4(const Field3D &f, CELL_LOC outloc, DIFF_METHOD method,
                               REGION region));
inline const Field3D D4DZ4(const Field3D &f, CELL_LOC outloc, DIFF_METHOD method,
                           REGION region = RGN_NOBNDRY) {
  return D4DZ4(f, outloc, region, method);
}

/// Calculate forth partial derivative in Z
///
///   $\partial^4 / \partial z^4$
///
/// @param[in] f       The field to be differentiated
/// @param[in] outloc  The cell location where the result is desired. If
///                    staggered grids is not enabled then this has no effect
///                    If not given, defaults to CELL_DEFAULT
/// @param[in] method  Differencing method to use. This overrides the default
///                    If not given, defaults to DIFF_DEFAULT
/// @param[in] region  What region is expected to be calculated
///                    If not given, defaults to RGN_NOBNDRY
///
///
///
const Field2D D4DZ4(const Field2D &f, CELL_LOC outloc = CELL_DEFAULT,
                    REGION region = RGN_NOBNDRY, DIFF_METHOD method = DIFF_DEFAULT);

DEPRECATED(const Field2D D4DZ4(const Field2D &f, DIFF_METHOD method, CELL_LOC outloc,
                               REGION region));
inline const Field2D D4DZ4(const Field2D &f, DIFF_METHOD method,
                           CELL_LOC outloc = CELL_DEFAULT, REGION region = RGN_NOBNDRY) {
  return D4DZ4(f, outloc, region, method);
}

DEPRECATED(const Field2D D4DZ4(const Field2D &f, CELL_LOC outloc, DIFF_METHOD method,
                               REGION region));
inline const Field2D D4DZ4(const Field2D &f, CELL_LOC outloc, DIFF_METHOD method,
                           REGION region = RGN_NOBNDRY) {
  return D4DZ4(f, outloc, region, method);
}
///////// UPWINDING METHODS /////////////

/// For terms of form v * grad(f)
///
///   $v \cdot \partial f / \partial x$
///
/// @param[in] v       The velocity field
/// @param[in] f       The field of the advected quantity
/// @param[in] outloc  The cell location where the result is desired. If
///                    staggered grids is not enabled then this has no effect
///                    If not given, defaults to CELL_DEFAULT
/// @param[in] method  Differencing method to use. This overrides the default
///                    If not given, defaults to DIFF_DEFAULT
/// @param[in] region  What region is expected to be calculated
///                    If not given, defaults to RGN_NOBNDRY
///
///
///
const Field3D VDDX(const Field3D &v, const Field3D &f, CELL_LOC outloc = CELL_DEFAULT,
                   REGION region = RGN_NOBNDRY, DIFF_METHOD method = DIFF_DEFAULT);

DEPRECATED(const Field3D VDDX(const Field3D &v, const Field3D &f, DIFF_METHOD method,
                              CELL_LOC outloc, REGION region));
inline const Field3D VDDX(const Field3D &v, const Field3D &f, DIFF_METHOD method,
                          CELL_LOC outloc = CELL_DEFAULT, REGION region = RGN_NOBNDRY) {
  return VDDX(v, f, outloc, region, method);
}

DEPRECATED(const Field3D VDDX(const Field3D &v, const Field3D &f, CELL_LOC outloc,
                              DIFF_METHOD method, REGION region));
inline const Field3D VDDX(const Field3D &v, const Field3D &f, CELL_LOC outloc,
                          DIFF_METHOD method, REGION region = RGN_NOBNDRY) {
  return VDDX(v, f, outloc, region, method);
}

/// For terms of form v * grad(f)
///
///   $v \cdot \partial f / \partial x$
///
/// @param[in] v       The velocity field
/// @param[in] f       The field of the advected quantity
/// @param[in] outloc  The cell location where the result is desired. If
///                    staggered grids is not enabled then this has no effect
///                    If not given, defaults to CELL_DEFAULT
/// @param[in] method  Differencing method to use. This overrides the default
///                    If not given, defaults to DIFF_DEFAULT
/// @param[in] region  What region is expected to be calculated
///                    If not given, defaults to RGN_NOBNDRY
///
///
///
const Field2D VDDX(const Field2D &v, const Field2D &f, CELL_LOC outloc = CELL_DEFAULT,
                   REGION region = RGN_NOBNDRY, DIFF_METHOD method = DIFF_DEFAULT);

DEPRECATED(const Field2D VDDX(const Field2D &v, const Field2D &f, DIFF_METHOD method,
                              CELL_LOC outloc, REGION region));
inline const Field2D VDDX(const Field2D &v, const Field2D &f, DIFF_METHOD method,
                          CELL_LOC outloc = CELL_DEFAULT, REGION region = RGN_NOBNDRY) {
  return VDDX(v, f, outloc, region, method);
}

DEPRECATED(const Field2D VDDX(const Field2D &v, const Field2D &f, CELL_LOC outloc,
                              DIFF_METHOD method, REGION region));
inline const Field2D VDDX(const Field2D &v, const Field2D &f, CELL_LOC outloc,
                          DIFF_METHOD method, REGION region = RGN_NOBNDRY) {
  return VDDX(v, f, outloc, region, method);
}

/// For terms of form v * grad(f)
///
///   $v \cdot \partial f / \partial y$
///
/// @param[in] v       The velocity field
/// @param[in] f       The field of the advected quantity
/// @param[in] outloc  The cell location where the result is desired. If
///                    staggered grids is not enabled then this has no effect
///                    If not given, defaults to CELL_DEFAULT
/// @param[in] method  Differencing method to use. This overrides the default
///                    If not given, defaults to DIFF_DEFAULT
/// @param[in] region  What region is expected to be calculated
///                    If not given, defaults to RGN_NOBNDRY
///
///
///
const Field3D VDDY(const Field3D &v, const Field3D &f, CELL_LOC outloc = CELL_DEFAULT,
                   REGION region = RGN_NOBNDRY, DIFF_METHOD method = DIFF_DEFAULT);

DEPRECATED(const Field3D VDDY(const Field3D &v, const Field3D &f, DIFF_METHOD method,
                              CELL_LOC outloc, REGION region));
inline const Field3D VDDY(const Field3D &v, const Field3D &f, DIFF_METHOD method,
                          CELL_LOC outloc = CELL_DEFAULT, REGION region = RGN_NOBNDRY) {
  return VDDY(v, f, outloc, region, method);
}

DEPRECATED(const Field3D VDDY(const Field3D &v, const Field3D &f, CELL_LOC outloc,
                              DIFF_METHOD method, REGION region));
inline const Field3D VDDY(const Field3D &v, const Field3D &f, CELL_LOC outloc,
                          DIFF_METHOD method, REGION region = RGN_NOBNDRY) {
  return VDDY(v, f, outloc, region, method);
}

/// For terms of form v * grad(f)
///
///   $v \cdot \partial f / \partial y$
///
/// @param[in] v       The velocity field
/// @param[in] f       The field of the advected quantity
/// @param[in] outloc  The cell location where the result is desired. If
///                    staggered grids is not enabled then this has no effect
///                    If not given, defaults to CELL_DEFAULT
/// @param[in] method  Differencing method to use. This overrides the default
///                    If not given, defaults to DIFF_DEFAULT
/// @param[in] region  What region is expected to be calculated
///                    If not given, defaults to RGN_NOBNDRY
///
///
///
const Field2D VDDY(const Field2D &v, const Field2D &f, CELL_LOC outloc = CELL_DEFAULT,
                   REGION region = RGN_NOBNDRY, DIFF_METHOD method = DIFF_DEFAULT);

DEPRECATED(const Field2D VDDY(const Field2D &v, const Field2D &f, DIFF_METHOD method,
                              CELL_LOC outloc, REGION region));
inline const Field2D VDDY(const Field2D &v, const Field2D &f, DIFF_METHOD method,
                          CELL_LOC outloc = CELL_DEFAULT, REGION region = RGN_NOBNDRY) {
  return VDDY(v, f, outloc, region, method);
}

DEPRECATED(const Field2D VDDY(const Field2D &v, const Field2D &f, CELL_LOC outloc,
                              DIFF_METHOD method, REGION region));
inline const Field2D VDDY(const Field2D &v, const Field2D &f, CELL_LOC outloc,
                          DIFF_METHOD method, REGION region = RGN_NOBNDRY) {
  return VDDY(v, f, outloc, region, method);
}

/// For terms of form v * grad(f)
///
///   $v \cdot \partial f / \partial z$
///
/// @param[in] v       The velocity field
/// @param[in] f       The field of the advected quantity
/// @param[in] outloc  The cell location where the result is desired. If
///                    staggered grids is not enabled then this has no effect
///                    If not given, defaults to CELL_DEFAULT
/// @param[in] method  Differencing method to use. This overrides the default
///                    If not given, defaults to DIFF_DEFAULT
/// @param[in] region  What region is expected to be calculated
///                    If not given, defaults to RGN_NOBNDRY
///
///
///
const Field3D VDDZ(const Field3D &v, const Field3D &f, CELL_LOC outloc = CELL_DEFAULT,
                   REGION region = RGN_NOBNDRY, DIFF_METHOD method = DIFF_DEFAULT);

DEPRECATED(const Field3D VDDZ(const Field3D &v, const Field3D &f, DIFF_METHOD method,
                              CELL_LOC outloc, REGION region));
inline const Field3D VDDZ(const Field3D &v, const Field3D &f, DIFF_METHOD method,
                          CELL_LOC outloc = CELL_DEFAULT, REGION region = RGN_NOBNDRY) {
  return VDDZ(v, f, outloc, region, method);
}

DEPRECATED(const Field3D VDDZ(const Field3D &v, const Field3D &f, CELL_LOC outloc,
                              DIFF_METHOD method, REGION region));
inline const Field3D VDDZ(const Field3D &v, const Field3D &f, CELL_LOC outloc,
                          DIFF_METHOD method, REGION region = RGN_NOBNDRY) {
  return VDDZ(v, f, outloc, region, method);
}

/// For terms of form v * grad(f)
///
///   $v \cdot \partial f / \partial z$
///
/// @param[in] v       The velocity field
/// @param[in] f       The field of the advected quantity
/// @param[in] outloc  The cell location where the result is desired. If
///                    staggered grids is not enabled then this has no effect
///                    If not given, defaults to CELL_DEFAULT
/// @param[in] method  Differencing method to use. This overrides the default
///                    If not given, defaults to DIFF_DEFAULT
/// @param[in] region  What region is expected to be calculated
///                    If not given, defaults to RGN_NOBNDRY
///
///
///
const Field2D VDDZ(const Field2D &v, const Field2D &f, CELL_LOC outloc = CELL_DEFAULT,
                   REGION region = RGN_NOBNDRY, DIFF_METHOD method = DIFF_DEFAULT);

DEPRECATED(const Field2D VDDZ(const Field2D &v, const Field2D &f, DIFF_METHOD method,
                              CELL_LOC outloc, REGION region));
inline const Field2D VDDZ(const Field2D &v, const Field2D &f, DIFF_METHOD method,
                          CELL_LOC outloc = CELL_DEFAULT, REGION region = RGN_NOBNDRY) {
  return VDDZ(v, f, outloc, region, method);
}

DEPRECATED(const Field2D VDDZ(const Field2D &v, const Field2D &f, CELL_LOC outloc,
                              DIFF_METHOD method, REGION region));
inline const Field2D VDDZ(const Field2D &v, const Field2D &f, CELL_LOC outloc,
                          DIFF_METHOD method, REGION region = RGN_NOBNDRY) {
  return VDDZ(v, f, outloc, region, method);
}
///////// FLUX METHODS /////////////

/// for terms of form div(v * f)
///
///   $\partial (v f) / \partial x$
///
/// @param[in] v       The velocity field
/// @param[in] f       The field of the advected quantity
/// @param[in] outloc  The cell location where the result is desired. If
///                    staggered grids is not enabled then this has no effect
///                    If not given, defaults to CELL_DEFAULT
/// @param[in] method  Differencing method to use. This overrides the default
///                    If not given, defaults to DIFF_DEFAULT
/// @param[in] region  What region is expected to be calculated
///                    If not given, defaults to RGN_NOBNDRY
///
///
///
const Field3D FDDX(const Field3D &v, const Field3D &f, CELL_LOC outloc = CELL_DEFAULT,
                   REGION region = RGN_NOBNDRY, DIFF_METHOD method = DIFF_DEFAULT);

DEPRECATED(const Field3D FDDX(const Field3D &v, const Field3D &f, DIFF_METHOD method,
                              CELL_LOC outloc, REGION region));
inline const Field3D FDDX(const Field3D &v, const Field3D &f, DIFF_METHOD method,
                          CELL_LOC outloc = CELL_DEFAULT, REGION region = RGN_NOBNDRY) {
  return FDDX(v, f, outloc, region, method);
}

DEPRECATED(const Field3D FDDX(const Field3D &v, const Field3D &f, CELL_LOC outloc,
                              DIFF_METHOD method, REGION region));
inline const Field3D FDDX(const Field3D &v, const Field3D &f, CELL_LOC outloc,
                          DIFF_METHOD method, REGION region = RGN_NOBNDRY) {
  return FDDX(v, f, outloc, region, method);
}

/// for terms of form div(v * f)
///
///   $\partial (v f) / \partial x$
///
/// @param[in] v       The velocity field
/// @param[in] f       The field of the advected quantity
/// @param[in] outloc  The cell location where the result is desired. If
///                    staggered grids is not enabled then this has no effect
///                    If not given, defaults to CELL_DEFAULT
/// @param[in] method  Differencing method to use. This overrides the default
///                    If not given, defaults to DIFF_DEFAULT
/// @param[in] region  What region is expected to be calculated
///                    If not given, defaults to RGN_NOBNDRY
///
///
///
const Field2D FDDX(const Field2D &v, const Field2D &f, CELL_LOC outloc = CELL_DEFAULT,
                   REGION region = RGN_NOBNDRY, DIFF_METHOD method = DIFF_DEFAULT);

DEPRECATED(const Field2D FDDX(const Field2D &v, const Field2D &f, DIFF_METHOD method,
                              CELL_LOC outloc, REGION region));
inline const Field2D FDDX(const Field2D &v, const Field2D &f, DIFF_METHOD method,
                          CELL_LOC outloc = CELL_DEFAULT, REGION region = RGN_NOBNDRY) {
  return FDDX(v, f, outloc, region, method);
}

DEPRECATED(const Field2D FDDX(const Field2D &v, const Field2D &f, CELL_LOC outloc,
                              DIFF_METHOD method, REGION region));
inline const Field2D FDDX(const Field2D &v, const Field2D &f, CELL_LOC outloc,
                          DIFF_METHOD method, REGION region = RGN_NOBNDRY) {
  return FDDX(v, f, outloc, region, method);
}

/// for terms of form div(v * f)
///
///   $\partial (v f) / \partial y$
///
/// @param[in] v       The velocity field
/// @param[in] f       The field of the advected quantity
/// @param[in] outloc  The cell location where the result is desired. If
///                    staggered grids is not enabled then this has no effect
///                    If not given, defaults to CELL_DEFAULT
/// @param[in] method  Differencing method to use. This overrides the default
///                    If not given, defaults to DIFF_DEFAULT
/// @param[in] region  What region is expected to be calculated
///                    If not given, defaults to RGN_NOBNDRY
///
///
///
const Field3D FDDY(const Field3D &v, const Field3D &f, CELL_LOC outloc = CELL_DEFAULT,
                   REGION region = RGN_NOBNDRY, DIFF_METHOD method = DIFF_DEFAULT);

DEPRECATED(const Field3D FDDY(const Field3D &v, const Field3D &f, DIFF_METHOD method,
                              CELL_LOC outloc, REGION region));
inline const Field3D FDDY(const Field3D &v, const Field3D &f, DIFF_METHOD method,
                          CELL_LOC outloc = CELL_DEFAULT, REGION region = RGN_NOBNDRY) {
  return FDDY(v, f, outloc, region, method);
}

DEPRECATED(const Field3D FDDY(const Field3D &v, const Field3D &f, CELL_LOC outloc,
                              DIFF_METHOD method, REGION region));
inline const Field3D FDDY(const Field3D &v, const Field3D &f, CELL_LOC outloc,
                          DIFF_METHOD method, REGION region = RGN_NOBNDRY) {
  return FDDY(v, f, outloc, region, method);
}

/// for terms of form div(v * f)
///
///   $\partial (v f) / \partial y$
///
/// @param[in] v       The velocity field
/// @param[in] f       The field of the advected quantity
/// @param[in] outloc  The cell location where the result is desired. If
///                    staggered grids is not enabled then this has no effect
///                    If not given, defaults to CELL_DEFAULT
/// @param[in] method  Differencing method to use. This overrides the default
///                    If not given, defaults to DIFF_DEFAULT
/// @param[in] region  What region is expected to be calculated
///                    If not given, defaults to RGN_NOBNDRY
///
///
///
const Field2D FDDY(const Field2D &v, const Field2D &f, CELL_LOC outloc = CELL_DEFAULT,
                   REGION region = RGN_NOBNDRY, DIFF_METHOD method = DIFF_DEFAULT);

DEPRECATED(const Field2D FDDY(const Field2D &v, const Field2D &f, DIFF_METHOD method,
                              CELL_LOC outloc, REGION region));
inline const Field2D FDDY(const Field2D &v, const Field2D &f, DIFF_METHOD method,
                          CELL_LOC outloc = CELL_DEFAULT, REGION region = RGN_NOBNDRY) {
  return FDDY(v, f, outloc, region, method);
}

DEPRECATED(const Field2D FDDY(const Field2D &v, const Field2D &f, CELL_LOC outloc,
                              DIFF_METHOD method, REGION region));
inline const Field2D FDDY(const Field2D &v, const Field2D &f, CELL_LOC outloc,
                          DIFF_METHOD method, REGION region = RGN_NOBNDRY) {
  return FDDY(v, f, outloc, region, method);
}

/// for terms of form div(v * f)
///
///   $\partial (v f) / \partial z$
///
/// @param[in] v       The velocity field
/// @param[in] f       The field of the advected quantity
/// @param[in] outloc  The cell location where the result is desired. If
///                    staggered grids is not enabled then this has no effect
///                    If not given, defaults to CELL_DEFAULT
/// @param[in] method  Differencing method to use. This overrides the default
///                    If not given, defaults to DIFF_DEFAULT
/// @param[in] region  What region is expected to be calculated
///                    If not given, defaults to RGN_NOBNDRY
///
///
///
const Field3D FDDZ(const Field3D &v, const Field3D &f, CELL_LOC outloc = CELL_DEFAULT,
                   REGION region = RGN_NOBNDRY, DIFF_METHOD method = DIFF_DEFAULT);

DEPRECATED(const Field3D FDDZ(const Field3D &v, const Field3D &f, DIFF_METHOD method,
                              CELL_LOC outloc, REGION region));
inline const Field3D FDDZ(const Field3D &v, const Field3D &f, DIFF_METHOD method,
                          CELL_LOC outloc = CELL_DEFAULT, REGION region = RGN_NOBNDRY) {
  return FDDZ(v, f, outloc, region, method);
}

DEPRECATED(const Field3D FDDZ(const Field3D &v, const Field3D &f, CELL_LOC outloc,
                              DIFF_METHOD method, REGION region));
inline const Field3D FDDZ(const Field3D &v, const Field3D &f, CELL_LOC outloc,
                          DIFF_METHOD method, REGION region = RGN_NOBNDRY) {
  return FDDZ(v, f, outloc, region, method);
}

/// for terms of form div(v * f)
///
///   $\partial (v f) / \partial z$
///
/// @param[in] v       The velocity field
/// @param[in] f       The field of the advected quantity
/// @param[in] outloc  The cell location where the result is desired. If
///                    staggered grids is not enabled then this has no effect
///                    If not given, defaults to CELL_DEFAULT
/// @param[in] method  Differencing method to use. This overrides the default
///                    If not given, defaults to DIFF_DEFAULT
/// @param[in] region  What region is expected to be calculated
///                    If not given, defaults to RGN_NOBNDRY
///
///
///
const Field2D FDDZ(const Field2D &v, const Field2D &f, CELL_LOC outloc = CELL_DEFAULT,
                   REGION region = RGN_NOBNDRY, DIFF_METHOD method = DIFF_DEFAULT);

DEPRECATED(const Field2D FDDZ(const Field2D &v, const Field2D &f, DIFF_METHOD method,
                              CELL_LOC outloc, REGION region));
inline const Field2D FDDZ(const Field2D &v, const Field2D &f, DIFF_METHOD method,
                          CELL_LOC outloc = CELL_DEFAULT, REGION region = RGN_NOBNDRY) {
  return FDDZ(v, f, outloc, region, method);
}

DEPRECATED(const Field2D FDDZ(const Field2D &v, const Field2D &f, CELL_LOC outloc,
                              DIFF_METHOD method, REGION region));
inline const Field2D FDDZ(const Field2D &v, const Field2D &f, CELL_LOC outloc,
                          DIFF_METHOD method, REGION region = RGN_NOBNDRY) {
  return FDDZ(v, f, outloc, region, method);
}

/// Calculate first partial derivative in Z
///
///   $\partial / \partial z$
///
/// @param[in] f       The field to be differentiated
/// @param[in] outloc  The cell location where the result is desired. If
///                    staggered grids is not enabled then this has no effect
///                    If not given, defaults to CELL_DEFAULT
/// @param[in] method  Differencing method to use. This overrides the default
///                    If not given, defaults to DIFF_DEFAULT
/// @param[in] region  What region is expected to be calculated
///                    If not given, defaults to RGN_NOBNDRY
///
///
///
const Vector3D DDZ(const Vector3D &f, CELL_LOC outloc = CELL_DEFAULT,
                   REGION region = RGN_NOBNDRY, DIFF_METHOD method = DIFF_DEFAULT);

DEPRECATED(const Vector3D DDZ(const Vector3D &f, DIFF_METHOD method, CELL_LOC outloc,
                              REGION region));
inline const Vector3D DDZ(const Vector3D &f, DIFF_METHOD method,
                          CELL_LOC outloc = CELL_DEFAULT, REGION region = RGN_NOBNDRY) {
  return DDZ(f, outloc, region, method);
}

DEPRECATED(const Vector3D DDZ(const Vector3D &f, CELL_LOC outloc, DIFF_METHOD method,
                              REGION region));
inline const Vector3D DDZ(const Vector3D &f, CELL_LOC outloc, DIFF_METHOD method,
                          REGION region = RGN_NOBNDRY) {
  return DDZ(f, outloc, region, method);
}

/// Calculate mixed partial derivative in x and y
///
///   $\partial^2 / \partial x \partial y$
///
/// @param[in] f       The field to be differentiated
/// @param[in] outloc  The cell location where the result is desired. If
///                    staggered grids is not enabled then this has no effect
///                    If not given, defaults to CELL_DEFAULT
/// @param[in] method  Differencing method to use. This overrides the default
///                    If not given, defaults to DIFF_DEFAULT
/// @param[in] region  What region is expected to be calculated
///                    If not given, defaults to RGN_NOBNDRY
///
///
///
const Field2D D2DXDY(const Field2D &f, CELL_LOC outloc = CELL_DEFAULT,
                     REGION region = RGN_NOBNDRY, DIFF_METHOD method = DIFF_DEFAULT);

DEPRECATED(const Field2D D2DXDY(const Field2D &f, DIFF_METHOD method, CELL_LOC outloc,
                                REGION region));
inline const Field2D D2DXDY(const Field2D &f, DIFF_METHOD method,
                            CELL_LOC outloc = CELL_DEFAULT, REGION region = RGN_NOBNDRY) {
  return D2DXDY(f, outloc, region, method);
}

DEPRECATED(const Field2D D2DXDY(const Field2D &f, CELL_LOC outloc, DIFF_METHOD method,
                                REGION region));
inline const Field2D D2DXDY(const Field2D &f, CELL_LOC outloc, DIFF_METHOD method,
                            REGION region = RGN_NOBNDRY) {
  return D2DXDY(f, outloc, region, method);
}

/// Calculate mixed partial derivative in x and y
///
///   $\partial^2 / \partial x \partial y$
///
/// @param[in] f       The field to be differentiated
/// @param[in] outloc  The cell location where the result is desired. If
///                    staggered grids is not enabled then this has no effect
///                    If not given, defaults to CELL_DEFAULT
/// @param[in] method  Differencing method to use. This overrides the default
///                    If not given, defaults to DIFF_DEFAULT
/// @param[in] region  What region is expected to be calculated
///                    If not given, defaults to RGN_NOBNDRY
///
///
///
const Field3D D2DXDY(const Field3D &f, CELL_LOC outloc = CELL_DEFAULT,
                     REGION region = RGN_NOBNDRY, DIFF_METHOD method = DIFF_DEFAULT);

DEPRECATED(const Field3D D2DXDY(const Field3D &f, DIFF_METHOD method, CELL_LOC outloc,
                                REGION region));
inline const Field3D D2DXDY(const Field3D &f, DIFF_METHOD method,
                            CELL_LOC outloc = CELL_DEFAULT, REGION region = RGN_NOBNDRY) {
  return D2DXDY(f, outloc, region, method);
}

DEPRECATED(const Field3D D2DXDY(const Field3D &f, CELL_LOC outloc, DIFF_METHOD method,
                                REGION region));
inline const Field3D D2DXDY(const Field3D &f, CELL_LOC outloc, DIFF_METHOD method,
                            REGION region = RGN_NOBNDRY) {
  return D2DXDY(f, outloc, region, method);
}

/// Calculate mixed partial derivative in x and z
///
///   $\partial^2 / \partial x \partial z$
///
/// @param[in] f       The field to be differentiated
/// @param[in] outloc  The cell location where the result is desired. If
///                    staggered grids is not enabled then this has no effect
///                    If not given, defaults to CELL_DEFAULT
/// @param[in] method  Differencing method to use. This overrides the default
///                    If not given, defaults to DIFF_DEFAULT
/// @param[in] region  What region is expected to be calculated
///                    If not given, defaults to RGN_NOBNDRY
///
///
///
const Field2D D2DXDZ(const Field2D &f, CELL_LOC outloc = CELL_DEFAULT,
                     REGION region = RGN_NOBNDRY, DIFF_METHOD method = DIFF_DEFAULT);

DEPRECATED(const Field2D D2DXDZ(const Field2D &f, DIFF_METHOD method, CELL_LOC outloc,
                                REGION region));
inline const Field2D D2DXDZ(const Field2D &f, DIFF_METHOD method,
                            CELL_LOC outloc = CELL_DEFAULT, REGION region = RGN_NOBNDRY) {
  return D2DXDZ(f, outloc, region, method);
}

DEPRECATED(const Field2D D2DXDZ(const Field2D &f, CELL_LOC outloc, DIFF_METHOD method,
                                REGION region));
inline const Field2D D2DXDZ(const Field2D &f, CELL_LOC outloc, DIFF_METHOD method,
                            REGION region = RGN_NOBNDRY) {
  return D2DXDZ(f, outloc, region, method);
}

/// Calculate mixed partial derivative in x and z
///
///   $\partial^2 / \partial x \partial z$
///
/// @param[in] f       The field to be differentiated
/// @param[in] outloc  The cell location where the result is desired. If
///                    staggered grids is not enabled then this has no effect
///                    If not given, defaults to CELL_DEFAULT
/// @param[in] method  Differencing method to use. This overrides the default
///                    If not given, defaults to DIFF_DEFAULT
/// @param[in] region  What region is expected to be calculated
///                    If not given, defaults to RGN_NOBNDRY
///
///
///
const Field3D D2DXDZ(const Field3D &f, CELL_LOC outloc = CELL_DEFAULT,
                     REGION region = RGN_NOBNDRY, DIFF_METHOD method = DIFF_DEFAULT);

DEPRECATED(const Field3D D2DXDZ(const Field3D &f, DIFF_METHOD method, CELL_LOC outloc,
                                REGION region));
inline const Field3D D2DXDZ(const Field3D &f, DIFF_METHOD method,
                            CELL_LOC outloc = CELL_DEFAULT, REGION region = RGN_NOBNDRY) {
  return D2DXDZ(f, outloc, region, method);
}

DEPRECATED(const Field3D D2DXDZ(const Field3D &f, CELL_LOC outloc, DIFF_METHOD method,
                                REGION region));
inline const Field3D D2DXDZ(const Field3D &f, CELL_LOC outloc, DIFF_METHOD method,
                            REGION region = RGN_NOBNDRY) {
  return D2DXDZ(f, outloc, region, method);
}

/// Calculate mixed partial derivative in y and z
///
///   $\partial^2 / \partial y \partial z$
///
/// @param[in] f       The field to be differentiated
/// @param[in] outloc  The cell location where the result is desired. If
///                    staggered grids is not enabled then this has no effect
///                    If not given, defaults to CELL_DEFAULT
/// @param[in] method  Differencing method to use. This overrides the default
///                    If not given, defaults to DIFF_DEFAULT
/// @param[in] region  What region is expected to be calculated
///                    If not given, defaults to RGN_NOBNDRY
///
///
///
const Field2D D2DYDZ(const Field2D &f, CELL_LOC outloc = CELL_DEFAULT,
                     REGION region = RGN_NOBNDRY, DIFF_METHOD method = DIFF_DEFAULT);

DEPRECATED(const Field2D D2DYDZ(const Field2D &f, DIFF_METHOD method, CELL_LOC outloc,
                                REGION region));
inline const Field2D D2DYDZ(const Field2D &f, DIFF_METHOD method,
                            CELL_LOC outloc = CELL_DEFAULT, REGION region = RGN_NOBNDRY) {
  return D2DYDZ(f, outloc, region, method);
}

DEPRECATED(const Field2D D2DYDZ(const Field2D &f, CELL_LOC outloc, DIFF_METHOD method,
                                REGION region));
inline const Field2D D2DYDZ(const Field2D &f, CELL_LOC outloc, DIFF_METHOD method,
                            REGION region = RGN_NOBNDRY) {
  return D2DYDZ(f, outloc, region, method);
}

/// Calculate mixed partial derivative in y and z
///
///   $\partial^2 / \partial y \partial z$
///
/// @param[in] f       The field to be differentiated
/// @param[in] outloc  The cell location where the result is desired. If
///                    staggered grids is not enabled then this has no effect
///                    If not given, defaults to CELL_DEFAULT
/// @param[in] method  Differencing method to use. This overrides the default
///                    If not given, defaults to DIFF_DEFAULT
/// @param[in] region  What region is expected to be calculated
///                    If not given, defaults to RGN_NOBNDRY
///
///
///
const Field3D D2DYDZ(const Field3D &f, CELL_LOC outloc = CELL_DEFAULT,
                     REGION region = RGN_NOBNDRY, DIFF_METHOD method = DIFF_DEFAULT);

DEPRECATED(const Field3D D2DYDZ(const Field3D &f, DIFF_METHOD method, CELL_LOC outloc,
                                REGION region));
inline const Field3D D2DYDZ(const Field3D &f, DIFF_METHOD method,
                            CELL_LOC outloc = CELL_DEFAULT, REGION region = RGN_NOBNDRY) {
  return D2DYDZ(f, outloc, region, method);
}

DEPRECATED(const Field3D D2DYDZ(const Field3D &f, CELL_LOC outloc, DIFF_METHOD method,
                                REGION region));
inline const Field3D D2DYDZ(const Field3D &f, CELL_LOC outloc, DIFF_METHOD method,
                            REGION region = RGN_NOBNDRY) {
  return D2DYDZ(f, outloc, region, method);
}

/// For terms of form v * grad(f)
///
///   $v \cdot \partial f / \partial z$
///
/// @param[in] v       The velocity field
/// @param[in] f       The field of the advected quantity
/// @param[in] outloc  The cell location where the result is desired. If
///                    staggered grids is not enabled then this has no effect
///                    If not given, defaults to CELL_DEFAULT
/// @param[in] method  Differencing method to use. This overrides the default
///                    If not given, defaults to DIFF_DEFAULT
/// @param[in] region  What region is expected to be calculated
///                    If not given, defaults to RGN_NOBNDRY
///
///
///
const Field2D VDDZ(const Field3D &v, const Field2D &f, CELL_LOC outloc = CELL_DEFAULT,
                   REGION region = RGN_NOBNDRY, DIFF_METHOD method = DIFF_DEFAULT);

DEPRECATED(const Field2D VDDZ(const Field3D &v, const Field2D &f, DIFF_METHOD method,
                              CELL_LOC outloc, REGION region));
inline const Field2D VDDZ(const Field3D &v, const Field2D &f, DIFF_METHOD method,
                          CELL_LOC outloc = CELL_DEFAULT, REGION region = RGN_NOBNDRY) {
  return VDDZ(v, f, outloc, region, method);
}

DEPRECATED(const Field2D VDDZ(const Field3D &v, const Field2D &f, CELL_LOC outloc,
                              DIFF_METHOD method, REGION region));
inline const Field2D VDDZ(const Field3D &v, const Field2D &f, CELL_LOC outloc,
                          DIFF_METHOD method, REGION region = RGN_NOBNDRY) {
  return VDDZ(v, f, outloc, region, method);
}

// Deprecated methods
//
// Calculate first partial derivative in Z
//
//   $\partial / \partial z$
//
// @param[in] f       The field to be differentiated
// @param[in] outloc  The cell location where the result is desired.
//                    If staggered grids is not enabled then this has no effect
// @param[in] method  Differencing method to use. This overrides the default
// @param[in] inc_xbndry  DEPRECATED: use REGION flags
//                    Determines whether the derivative should be calculated in
//                    the X boundaries. This allows mixed operators (e.g.
//                    D2DXDZ) without additional communication

inline const Field3D DDZ(const Field3D &f, CELL_LOC outloc, DIFF_METHOD method,
                         bool inc_xbndry) {
  return DDZ(f, outloc, inc_xbndry ? RGN_NOY : RGN_NOBNDRY, method);
}

inline const Field3D DDZ(const Field3D &f, DIFF_METHOD method, CELL_LOC outloc,
                         bool inc_xbndry) {
  return DDZ(f, outloc, inc_xbndry ? RGN_NOY : RGN_NOBNDRY, method);
}

inline const Field3D DDZ(const Field3D &f, DIFF_METHOD method, bool inc_xbndry) {
  return DDZ(f, CELL_DEFAULT, inc_xbndry ? RGN_NOY : RGN_NOBNDRY, method);
}

inline const Field3D DDZ(const Field3D &f, bool inc_xbndry) {
  return DDZ(f, CELL_DEFAULT, inc_xbndry ? RGN_NOY : RGN_NOBNDRY, DIFF_DEFAULT);
}

#endif // __DERIVS_H__
