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

#include "bout_types.hxx"

#ifdef DERIV_FUNC_REGION_ENUM_TO_STRING
#error This utility macro should not clash with another one
#else
#define DERIV_FUNC_REGION_ENUM_TO_STRING(func, T) \
[[deprecated("Please use #func(const #T& f, CELL_LOC outloc = CELL_DEFAULT, " \
    "const std::string& method = \"DEFAULT\", const std::string& region = \"RGN_ALL\") " \
    "instead")]] \
inline T func(const T& f, CELL_LOC outloc, const std::string& method, \
    REGION region) { \
  return func(f, outloc, method, toString(region)); \
} \
[[deprecated("Please use #func(const #T& f, CELL_LOC outloc = CELL_DEFAULT, " \
    "const std::string& method = \"DEFAULT\", const std::string& region = \"RGN_ALL\") " \
    "instead")]] \
inline T func(const T& f, CELL_LOC outloc, DIFF_METHOD method, \
    REGION region = RGN_NOBNDRY) { \
  return func(f, outloc, toString(method), toString(region)); \
}
#endif

#ifdef VDERIV_FUNC_REGION_ENUM_TO_STRING
#error This utility macro should not clash with another one
#else
#define VDERIV_FUNC_REGION_ENUM_TO_STRING(func, T, T1, T2) \
[[deprecated("Please use #func(const #T1 v, const #T2& f, " \
    "CELL_LOC outloc = CELL_DEFAULT, const std::string& method = \"DEFAULT\", const " \
    "std::string& region = \"RGN_ALL\") instead")]] \
inline T func(const T1& v, const T2& f, CELL_LOC outloc, const std::string& method, \
    REGION region) { \
  return func(v, f, outloc, method, toString(region)); \
} \
[[deprecated("Please use #func(const #T1& v, const #T2& f, " \
    "CELL_LOC outloc = CELL_DEFAULT, const std::string& method = \"DEFAULT\", " \
    "const std::string& region = \"RGN_ALL\") instead")]] \
inline T func(const T1& v, const T2& f, CELL_LOC outloc, DIFF_METHOD method, \
    REGION region = RGN_NOBNDRY) { \
  return func(v, f, outloc, toString(method), toString(region)); \
}
#endif

////////// FIRST DERIVATIVES //////////

/// Calculate first partial derivative in X
///
///   \f$\partial / \partial x\f$
///
/// @param[in] f       The field to be differentiated
/// @param[in] outloc  The cell location where the result is desired. If
///                    staggered grids is not enabled then this has no effect
///                    If not given, defaults to CELL_DEFAULT
/// @param[in] method  Differencing method to use. This overrides the default
///                    If not given, defaults to DIFF_DEFAULT
/// @param[in] region  What region is expected to be calculated
///                    If not given, defaults to RGN_NOBNDRY
Field3D DDX(const Field3D& f, CELL_LOC outloc = CELL_DEFAULT, const std::string&
    method = "DEFAULT", const std::string& region = "RGN_NOBNDRY");
DERIV_FUNC_REGION_ENUM_TO_STRING(DDX, Field3D)

/// Calculate first partial derivative in X
///
///   \f$\partial / \partial x\f$
///
/// @param[in] f       The field to be differentiated
/// @param[in] outloc  The cell location where the result is desired. If
///                    staggered grids is not enabled then this has no effect
///                    If not given, defaults to CELL_DEFAULT
/// @param[in] method  Differencing method to use. This overrides the default
///                    If not given, defaults to DIFF_DEFAULT
/// @param[in] region  What region is expected to be calculated
///                    If not given, defaults to RGN_NOBNDRY
Field2D DDX(const Field2D& f, CELL_LOC outloc = CELL_DEFAULT, const std::string&
    method = "DEFAULT", const std::string& region = "RGN_NOBNDRY");
DERIV_FUNC_REGION_ENUM_TO_STRING(DDX, Field2D)

/// Calculate first partial derivative in Y
///
///   \f$\partial / \partial y\f$
///
/// @param[in] f       The field to be differentiated
/// @param[in] outloc  The cell location where the result is desired. If
///                    staggered grids is not enabled then this has no effect
///                    If not given, defaults to CELL_DEFAULT
/// @param[in] method  Differencing method to use. This overrides the default
///                    If not given, defaults to DIFF_DEFAULT
/// @param[in] region  What region is expected to be calculated
///                    If not given, defaults to RGN_NOBNDRY
Field3D DDY(const Field3D& f, CELL_LOC outloc = CELL_DEFAULT, const std::string&
    method = "DEFAULT", const std::string& region = "RGN_NOBNDRY");
DERIV_FUNC_REGION_ENUM_TO_STRING(DDY, Field3D)

/// Calculate first partial derivative in Y
///
///   \f$\partial / \partial y\f$
///
/// @param[in] f       The field to be differentiated
/// @param[in] outloc  The cell location where the result is desired. If
///                    staggered grids is not enabled then this has no effect
///                    If not given, defaults to CELL_DEFAULT
/// @param[in] method  Differencing method to use. This overrides the default
///                    If not given, defaults to DIFF_DEFAULT
/// @param[in] region  What region is expected to be calculated
///                    If not given, defaults to RGN_NOBNDRY
Field2D DDY(const Field2D& f, CELL_LOC outloc = CELL_DEFAULT, const std::string&
    method = "DEFAULT", const std::string& region = "RGN_NOBNDRY");
DERIV_FUNC_REGION_ENUM_TO_STRING(DDY, Field2D)

/// Calculate first partial derivative in Z
///
///   \f$\partial / \partial z\f$
///
/// @param[in] f       The field to be differentiated
/// @param[in] outloc  The cell location where the result is desired. If
///                    staggered grids is not enabled then this has no effect
///                    If not given, defaults to CELL_DEFAULT
/// @param[in] method  Differencing method to use. This overrides the default
///                    If not given, defaults to DIFF_DEFAULT
/// @param[in] region  What region is expected to be calculated
///                    If not given, defaults to RGN_NOBNDRY
Field3D DDZ(const Field3D& f, CELL_LOC outloc = CELL_DEFAULT, const std::string&
    method = "DEFAULT", const std::string& region = "RGN_NOBNDRY");
DERIV_FUNC_REGION_ENUM_TO_STRING(DDZ, Field3D)

/// Calculate first partial derivative in Z
///
///   \f$\partial / \partial z\f$
///
/// @param[in] f       The field to be differentiated
/// @param[in] outloc  The cell location where the result is desired. If
///                    staggered grids is not enabled then this has no effect
///                    If not given, defaults to CELL_DEFAULT
/// @param[in] method  Differencing method to use. This overrides the default
///                    If not given, defaults to DIFF_DEFAULT
/// @param[in] region  What region is expected to be calculated
///                    If not given, defaults to RGN_NOBNDRY
Field2D DDZ(const Field2D& f, CELL_LOC outloc = CELL_DEFAULT, const std::string&
    method = "DEFAULT", const std::string& region = "RGN_NOBNDRY");
DERIV_FUNC_REGION_ENUM_TO_STRING(DDZ, Field2D)

/// Calculate first partial derivative in Z
///
///   \f$\partial / \partial z\f$
///
/// @param[in] f       The field to be differentiated
/// @param[in] outloc  The cell location where the result is desired. If
///                    staggered grids is not enabled then this has no effect
///                    If not given, defaults to CELL_DEFAULT
/// @param[in] method  Differencing method to use. This overrides the default
///                    If not given, defaults to DIFF_DEFAULT
/// @param[in] region  What region is expected to be calculated
///                    If not given, defaults to RGN_NOBNDRY
Vector3D DDZ(const Vector3D& f, CELL_LOC outloc = CELL_DEFAULT, const std::string&
    method = "DEFAULT", const std::string& region = "RGN_NOBNDRY");
DERIV_FUNC_REGION_ENUM_TO_STRING(DDZ, Vector3D)

/// Calculate first partial derivative in Z
///
///   \f$\partial / \partial z\f$
///
/// @param[in] f       The field to be differentiated
/// @param[in] outloc  The cell location where the result is desired. If
///                    staggered grids is not enabled then this has no effect
///                    If not given, defaults to CELL_DEFAULT
/// @param[in] method  Differencing method to use. This overrides the default
///                    If not given, defaults to DIFF_DEFAULT
/// @param[in] region  What region is expected to be calculated
///                    If not given, defaults to RGN_NOBNDRY
Vector2D DDZ(const Vector2D& f, CELL_LOC outloc = CELL_DEFAULT, const std::string&
    method = "DEFAULT", const std::string& region = "RGN_NOBNDRY");
DERIV_FUNC_REGION_ENUM_TO_STRING(DDZ, Vector2D)

////////// SECOND DERIVATIVES //////////

/// Calculate second partial derivative in X
///
///   \f$\partial^2 / \partial x^2\f$
///
/// @param[in] f       The field to be differentiated
/// @param[in] outloc  The cell location where the result is desired. If
///                    staggered grids is not enabled then this has no effect
///                    If not given, defaults to CELL_DEFAULT
/// @param[in] method  Differencing method to use. This overrides the default
///                    If not given, defaults to DIFF_DEFAULT
/// @param[in] region  What region is expected to be calculated
///                    If not given, defaults to RGN_NOBNDRY
Field3D D2DX2(const Field3D& f, CELL_LOC outloc = CELL_DEFAULT, const std::string&
    method = "DEFAULT", const std::string& region = "RGN_NOBNDRY");
DERIV_FUNC_REGION_ENUM_TO_STRING(D2DX2, Field3D)

/// Calculate second partial derivative in X
///
///   \f$\partial^2 / \partial x^2\f$
///
/// @param[in] f       The field to be differentiated
/// @param[in] outloc  The cell location where the result is desired. If
///                    staggered grids is not enabled then this has no effect
///                    If not given, defaults to CELL_DEFAULT
/// @param[in] method  Differencing method to use. This overrides the default
///                    If not given, defaults to DIFF_DEFAULT
/// @param[in] region  What region is expected to be calculated
///                    If not given, defaults to RGN_NOBNDRY
Field2D D2DX2(const Field2D& f, CELL_LOC outloc = CELL_DEFAULT, const std::string&
    method = "DEFAULT", const std::string& region = "RGN_NOBNDRY");
DERIV_FUNC_REGION_ENUM_TO_STRING(D2DX2, Field2D)

/// Calculate second partial derivative in Y
///
///   \f$\partial^2 / \partial y^2\f$
///
/// @param[in] f       The field to be differentiated
/// @param[in] outloc  The cell location where the result is desired. If
///                    staggered grids is not enabled then this has no effect
///                    If not given, defaults to CELL_DEFAULT
/// @param[in] method  Differencing method to use. This overrides the default
///                    If not given, defaults to DIFF_DEFAULT
/// @param[in] region  What region is expected to be calculated
///                    If not given, defaults to RGN_NOBNDRY
Field3D D2DY2(const Field3D& f, CELL_LOC outloc = CELL_DEFAULT, const std::string&
    method = "DEFAULT", const std::string& region = "RGN_NOBNDRY");
DERIV_FUNC_REGION_ENUM_TO_STRING(D2DY2, Field3D)

/// Calculate second partial derivative in Y
///
///   \f$\partial^2 / \partial y^2\f$
///
/// @param[in] f       The field to be differentiated
/// @param[in] outloc  The cell location where the result is desired. If
///                    staggered grids is not enabled then this has no effect
///                    If not given, defaults to CELL_DEFAULT
/// @param[in] method  Differencing method to use. This overrides the default
///                    If not given, defaults to DIFF_DEFAULT
/// @param[in] region  What region is expected to be calculated
///                    If not given, defaults to RGN_NOBNDRY
Field2D D2DY2(const Field2D& f, CELL_LOC outloc = CELL_DEFAULT, const std::string&
    method = "DEFAULT", const std::string& region = "RGN_NOBNDRY");
DERIV_FUNC_REGION_ENUM_TO_STRING(D2DY2, Field2D)

/// Calculate second partial derivative in Z
///
///   \f$\partial^2 / \partial z^2\f$
///
/// @param[in] f       The field to be differentiated
/// @param[in] outloc  The cell location where the result is desired. If
///                    staggered grids is not enabled then this has no effect
///                    If not given, defaults to CELL_DEFAULT
/// @param[in] method  Differencing method to use. This overrides the default
///                    If not given, defaults to DIFF_DEFAULT
/// @param[in] region  What region is expected to be calculated
///                    If not given, defaults to RGN_NOBNDRY
Field3D D2DZ2(const Field3D& f, CELL_LOC outloc = CELL_DEFAULT, const std::string&
    method = "DEFAULT", const std::string& region = "RGN_NOBNDRY");
DERIV_FUNC_REGION_ENUM_TO_STRING(D2DZ2, Field3D)

/// Calculate second partial derivative in Z
///
///   \f$\partial^2 / \partial z^2\f$
///
/// @param[in] f       The field to be differentiated
/// @param[in] outloc  The cell location where the result is desired. If
///                    staggered grids is not enabled then this has no effect
///                    If not given, defaults to CELL_DEFAULT
/// @param[in] method  Differencing method to use. This overrides the default
///                    If not given, defaults to DIFF_DEFAULT
/// @param[in] region  What region is expected to be calculated
///                    If not given, defaults to RGN_NOBNDRY
Field2D D2DZ2(const Field2D& f, CELL_LOC outloc = CELL_DEFAULT, const std::string&
    method = "DEFAULT", const std::string& region = "RGN_NOBNDRY");
DERIV_FUNC_REGION_ENUM_TO_STRING(D2DZ2, Field2D)

////////// FOURTH DERIVATIVES //////////

/// Calculate forth partial derivative in X
///
///   \f$\partial^4 / \partial x^4\f$
///
/// @param[in] f       The field to be differentiated
/// @param[in] outloc  The cell location where the result is desired. If
///                    staggered grids is not enabled then this has no effect
///                    If not given, defaults to CELL_DEFAULT
/// @param[in] method  Differencing method to use. This overrides the default
///                    If not given, defaults to DIFF_DEFAULT
/// @param[in] region  What region is expected to be calculated
///                    If not given, defaults to RGN_NOBNDRY
Field3D D4DX4(const Field3D& f, CELL_LOC outloc = CELL_DEFAULT, const std::string&
    method = "DEFAULT", const std::string& region = "RGN_NOBNDRY");
DERIV_FUNC_REGION_ENUM_TO_STRING(D4DX4, Field3D)

/// Calculate forth partial derivative in X
///
///   \f$\partial^4 / \partial x^4\f$
///
/// @param[in] f       The field to be differentiated
/// @param[in] outloc  The cell location where the result is desired. If
///                    staggered grids is not enabled then this has no effect
///                    If not given, defaults to CELL_DEFAULT
/// @param[in] method  Differencing method to use. This overrides the default
///                    If not given, defaults to DIFF_DEFAULT
/// @param[in] region  What region is expected to be calculated
///                    If not given, defaults to RGN_NOBNDRY
Field2D D4DX4(const Field2D& f, CELL_LOC outloc = CELL_DEFAULT, const std::string&
    method = "DEFAULT", const std::string& region = "RGN_NOBNDRY");
DERIV_FUNC_REGION_ENUM_TO_STRING(D4DX4, Field2D)

/// Calculate forth partial derivative in Y
///
///   \f$\partial^4 / \partial y^4\f$
///
/// @param[in] f       The field to be differentiated
/// @param[in] outloc  The cell location where the result is desired. If
///                    staggered grids is not enabled then this has no effect
///                    If not given, defaults to CELL_DEFAULT
/// @param[in] method  Differencing method to use. This overrides the default
///                    If not given, defaults to DIFF_DEFAULT
/// @param[in] region  What region is expected to be calculated
///                    If not given, defaults to RGN_NOBNDRY
Field3D D4DY4(const Field3D& f, CELL_LOC outloc = CELL_DEFAULT, const std::string&
    method = "DEFAULT", const std::string& region = "RGN_NOBNDRY");
DERIV_FUNC_REGION_ENUM_TO_STRING(D4DY4, Field3D)

/// Calculate forth partial derivative in Y
///
///   \f$\partial^4 / \partial y^4\f$
///
/// @param[in] f       The field to be differentiated
/// @param[in] outloc  The cell location where the result is desired. If
///                    staggered grids is not enabled then this has no effect
///                    If not given, defaults to CELL_DEFAULT
/// @param[in] method  Differencing method to use. This overrides the default
///                    If not given, defaults to DIFF_DEFAULT
/// @param[in] region  What region is expected to be calculated
///                    If not given, defaults to RGN_NOBNDRY
Field2D D4DY4(const Field2D& f, CELL_LOC outloc = CELL_DEFAULT, const std::string&
    method = "DEFAULT", const std::string& region = "RGN_NOBNDRY");
DERIV_FUNC_REGION_ENUM_TO_STRING(D4DY4, Field2D)

/// Calculate forth partial derivative in Z
///
///   \f$\partial^4 / \partial z^4\f$
///
/// @param[in] f       The field to be differentiated
/// @param[in] outloc  The cell location where the result is desired. If
///                    staggered grids is not enabled then this has no effect
///                    If not given, defaults to CELL_DEFAULT
/// @param[in] method  Differencing method to use. This overrides the default
///                    If not given, defaults to DIFF_DEFAULT
/// @param[in] region  What region is expected to be calculated
///                    If not given, defaults to RGN_NOBNDRY
Field3D D4DZ4(const Field3D& f, CELL_LOC outloc = CELL_DEFAULT, const std::string&
    method = "DEFAULT", const std::string& region = "RGN_NOBNDRY");
DERIV_FUNC_REGION_ENUM_TO_STRING(D4DZ4, Field3D)

/// Calculate forth partial derivative in Z
///
///   \f$\partial^4 / \partial z^4\f$
///
/// @param[in] f       The field to be differentiated
/// @param[in] outloc  The cell location where the result is desired. If
///                    staggered grids is not enabled then this has no effect
///                    If not given, defaults to CELL_DEFAULT
/// @param[in] method  Differencing method to use. This overrides the default
///                    If not given, defaults to DIFF_DEFAULT
/// @param[in] region  What region is expected to be calculated
///                    If not given, defaults to RGN_NOBNDRY
Field2D D4DZ4(const Field2D& f, CELL_LOC outloc = CELL_DEFAULT, const std::string&
    method = "DEFAULT", const std::string& region = "RGN_NOBNDRY");
DERIV_FUNC_REGION_ENUM_TO_STRING(D4DZ4, Field2D)

/// For terms of form v * grad(f)
///
///   \f$v \cdot \partial f / \partial x\f$
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
Field3D VDDX(const Field3D& v, const Field3D& f, CELL_LOC outloc = CELL_DEFAULT,
    const std::string& method = "DEFAULT", const std::string& region = "RGN_NOBNDRY");
VDERIV_FUNC_REGION_ENUM_TO_STRING(VDDX, Field3D, Field3D, Field3D)

/// For terms of form v * grad(f)
///
///   \f$v \cdot \partial f / \partial x\f$
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
Field2D VDDX(const Field2D& v, const Field2D& f, CELL_LOC outloc = CELL_DEFAULT,
    const std::string& method = "DEFAULT", const std::string& region = "RGN_NOBNDRY");
VDERIV_FUNC_REGION_ENUM_TO_STRING(VDDX, Field2D, Field2D, Field2D)

/// For terms of form v * grad(f)
///
///   \f$v \cdot \partial f / \partial y\f$
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
Field3D VDDY(const Field3D& v, const Field3D& f, CELL_LOC outloc = CELL_DEFAULT,
    const std::string& method = "DEFAULT", const std::string& region = "RGN_NOBNDRY");
VDERIV_FUNC_REGION_ENUM_TO_STRING(VDDY, Field3D, Field3D, Field3D)

/// For terms of form v * grad(f)
///
///   \f$v \cdot \partial f / \partial y\f$
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
Field2D VDDY(const Field2D& v, const Field2D& f, CELL_LOC outloc = CELL_DEFAULT,
    const std::string& method = "DEFAULT", const std::string& region = "RGN_NOBNDRY");
VDERIV_FUNC_REGION_ENUM_TO_STRING(VDDY, Field2D, Field2D, Field2D)

/// For terms of form v * grad(f)
///
///   \f$v \cdot \partial f / \partial z\f$
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
Field3D VDDZ(const Field3D& v, const Field3D& f, CELL_LOC outloc = CELL_DEFAULT,
    const std::string& method = "DEFAULT", const std::string& region = "RGN_NOBNDRY");
VDERIV_FUNC_REGION_ENUM_TO_STRING(VDDZ, Field3D, Field3D, Field3D)

/// For terms of form v * grad(f)
///
///   \f$v \cdot \partial f / \partial z\f$
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
Field2D VDDZ(const Field2D& v, const Field2D& f, CELL_LOC outloc = CELL_DEFAULT,
    const std::string& method = "DEFAULT", const std::string& region = "RGN_NOBNDRY");
VDERIV_FUNC_REGION_ENUM_TO_STRING(VDDZ, Field2D, Field2D, Field2D)

/// For terms of form v * grad(f)
///
///   \f$v \cdot \partial f / \partial z\f$
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
Field2D VDDZ(const Field3D& v, const Field2D& f, CELL_LOC outloc = CELL_DEFAULT,
    const std::string& method = "DEFAULT", const std::string& region = "RGN_NOBNDRY");
VDERIV_FUNC_REGION_ENUM_TO_STRING(VDDZ, Field2D, Field3D, Field2D)

/// for terms of form div(v * f)
///
///   \f$\partial (v f) / \partial x\f$
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
Field3D FDDX(const Field3D& v, const Field3D& f, CELL_LOC outloc = CELL_DEFAULT,
    const std::string& method = "DEFAULT", const std::string& region = "RGN_NOBNDRY");
VDERIV_FUNC_REGION_ENUM_TO_STRING(FDDX, Field3D, Field3D, Field3D)

/// for terms of form div(v * f)
///
///   \f$\partial (v f) / \partial x\f$
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
Field2D FDDX(const Field2D& v, const Field2D& f, CELL_LOC outloc = CELL_DEFAULT,
    const std::string& method = "DEFAULT", const std::string& region = "RGN_NOBNDRY");
VDERIV_FUNC_REGION_ENUM_TO_STRING(FDDX, Field2D, Field2D, Field2D)

/// for terms of form div(v * f)
///
///   \f$\partial (v f) / \partial y\f$
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
Field3D FDDY(const Field3D& v, const Field3D& f, CELL_LOC outloc = CELL_DEFAULT,
    const std::string& method = "DEFAULT", const std::string& region = "RGN_NOBNDRY");
VDERIV_FUNC_REGION_ENUM_TO_STRING(FDDY, Field3D, Field3D, Field3D)

/// for terms of form div(v * f)
///
///   \f$\partial (v f) / \partial y\f$
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
Field2D FDDY(const Field2D& v, const Field2D& f, CELL_LOC outloc = CELL_DEFAULT,
    const std::string& method = "DEFAULT", const std::string& region = "RGN_NOBNDRY");
VDERIV_FUNC_REGION_ENUM_TO_STRING(FDDY, Field2D, Field2D, Field2D)

/// for terms of form div(v * f)
///
///   \f$\partial (v f) / \partial z\f$
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
Field3D FDDZ(const Field3D& v, const Field3D& f, CELL_LOC outloc = CELL_DEFAULT,
    const std::string& method = "DEFAULT", const std::string& region = "RGN_NOBNDRY");
VDERIV_FUNC_REGION_ENUM_TO_STRING(FDDZ, Field3D, Field3D, Field3D)

/// for terms of form div(v * f)
///
///   \f$\partial (v f) / \partial z\f$
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
Field2D FDDZ(const Field2D& v, const Field2D& f, CELL_LOC outloc = CELL_DEFAULT,
    const std::string& method = "DEFAULT", const std::string& region = "RGN_NOBNDRY");
VDERIV_FUNC_REGION_ENUM_TO_STRING(FDDZ, Field2D, Field2D, Field2D)

/// Calculate mixed partial derivative in x and y
///
///   \f$\partial^2 / \partial x \partial y\f$
///
/// @param[in] f       The field to be differentiated
/// @param[in] outloc  The cell location where the result is desired. If
///                    staggered grids is not enabled then this has no effect
///                    If not given, defaults to CELL_DEFAULT
/// @param[in] method  Differencing method to use. This overrides the default
///                    If not given, defaults to DIFF_DEFAULT
/// @param[in] region  What region is expected to be calculated
///                    If not given, defaults to RGN_NOBNDRY
/// @param[in] dfdy_boundary_condition Boundary condition to use to set the guard cells of
///                                    df/dy, before calculating the x-derivative.
Field3D D2DXDY(const Field3D& f, CELL_LOC outloc = CELL_DEFAULT,
                     const std::string& method = "DEFAULT",
                     const std::string& region = "RGN_NOBNDRY",
                     const std::string& dfdy_boundary_condition = "free_o3");
[[deprecated("Please use D2DXDY(const Field3D& f, CELL_LOC outloc = CELL_DEFAULT, " 
    "const std::string& method = \"DEFAULT\", const std::string& region = \"RGN_ALL\", "
    "const std::string& dfdy_boundary_condition) instead")]]
inline Field3D D2DXDY(const Field3D& f, CELL_LOC outloc, const std::string& method,
                            REGION region,
                            const std::string& dfdy_boundary_condition = "free_o3") {
  return D2DXDY(f, outloc, method, toString(region), dfdy_boundary_condition);
}
[[deprecated("Please use D2DXDY(const Field3D& f, CELL_LOC outloc = CELL_DEFAULT, " 
    "const std::string& method = \"DEFAULT\", const std::string& region = \"RGN_ALL\", "
    "const std::string& dfdy_boundary_condition) instead")]]
inline Field3D D2DXDY(const Field3D& f, CELL_LOC outloc, DIFF_METHOD method,
                            REGION region = RGN_NOBNDRY,
                            const std::string& dfdy_boundary_condition = "free_o3") {
  return D2DXDY(f, outloc, toString(method), toString(region), dfdy_boundary_condition);
};

/// Calculate mixed partial derivative in x and y
///
///   \f$\partial^2 / \partial x \partial y\f$
///
/// @param[in] f       The field to be differentiated
/// @param[in] outloc  The cell location where the result is desired. If
///                    staggered grids is not enabled then this has no effect
///                    If not given, defaults to CELL_DEFAULT
/// @param[in] method  Differencing method to use. This overrides the default
///                    If not given, defaults to DIFF_DEFAULT
/// @param[in] region  What region is expected to be calculated
///                    If not given, defaults to RGN_NOBNDRY
/// @param[in] dfdy_boundary_condition Boundary condition to use to set the guard cells of
///                                    df/dy, before calculating the x-derivative.
Field2D D2DXDY(const Field2D& f, CELL_LOC outloc = CELL_DEFAULT,
                     const std::string& method = "DEFAULT",
                     const std::string& region = "RGN_NOBNDRY",
                     const std::string& dfdy_boundary_condition = "free_o3");
[[deprecated("Please use D2DXDY(const Field2D& f, CELL_LOC outloc = CELL_DEFAULT, " 
    "const std::string& method = \"DEFAULT\", const std::string& region = \"RGN_ALL\", "
    "const std::string& dfdy_boundary_condition) instead")]]
inline Field2D D2DXDY(const Field2D& f, CELL_LOC outloc, const std::string& method,
                            REGION region,
                            const std::string& dfdy_boundary_condition = "free_o3") {
  return D2DXDY(f, outloc, method, toString(region), dfdy_boundary_condition);
}
[[deprecated("Please use D2DXDY(const Field2D& f, CELL_LOC outloc = CELL_DEFAULT, " 
    "const std::string& method = \"DEFAULT\", const std::string& region = \"RGN_ALL\", "
    "const std::string& dfdy_boundary_condition) instead")]]
inline Field2D D2DXDY(const Field2D& f, CELL_LOC outloc, DIFF_METHOD method,
                            REGION region = RGN_NOBNDRY,
                            const std::string& dfdy_boundary_condition = "free_o3") {
  return D2DXDY(f, outloc, toString(method), toString(region), dfdy_boundary_condition);
};

/// Calculate mixed partial derivative in x and z
///
///   \f$\partial^2 / \partial x \partial z\f$
///
/// @param[in] f       The field to be differentiated
/// @param[in] outloc  The cell location where the result is desired. If
///                    staggered grids is not enabled then this has no effect
///                    If not given, defaults to CELL_DEFAULT
/// @param[in] method  Differencing method to use. This overrides the default
///                    If not given, defaults to DIFF_DEFAULT
/// @param[in] region  What region is expected to be calculated
///                    If not given, defaults to RGN_NOBNDRY
Field3D D2DXDZ(const Field3D& f, CELL_LOC outloc = CELL_DEFAULT, const std::string&
    method = "DEFAULT", const std::string& region = "RGN_NOBNDRY");
DERIV_FUNC_REGION_ENUM_TO_STRING(D2DXDZ, Field3D)

/// Calculate mixed partial derivative in x and z
///
///   \f$\partial^2 / \partial x \partial z\f$
///
/// @param[in] f       The field to be differentiated
/// @param[in] outloc  The cell location where the result is desired. If
///                    staggered grids is not enabled then this has no effect
///                    If not given, defaults to CELL_DEFAULT
/// @param[in] method  Differencing method to use. This overrides the default
///                    If not given, defaults to DIFF_DEFAULT
/// @param[in] region  What region is expected to be calculated
///                    If not given, defaults to RGN_NOBNDRY
Field2D D2DXDZ(const Field2D& f, CELL_LOC outloc = CELL_DEFAULT, const std::string&
    method = "DEFAULT", const std::string& region = "RGN_NOBNDRY");
DERIV_FUNC_REGION_ENUM_TO_STRING(D2DXDZ, Field2D)

/// Calculate mixed partial derivative in y and z
///
///   \f$\partial^2 / \partial y \partial z\f$
///
/// @param[in] f       The field to be differentiated
/// @param[in] outloc  The cell location where the result is desired. If
///                    staggered grids is not enabled then this has no effect
///                    If not given, defaults to CELL_DEFAULT
/// @param[in] method  Differencing method to use. This overrides the default
///                    If not given, defaults to DIFF_DEFAULT
/// @param[in] region  What region is expected to be calculated
///                    If not given, defaults to RGN_NOBNDRY
Field3D D2DYDZ(const Field3D& f, CELL_LOC outloc = CELL_DEFAULT, const std::string&
    method = "DEFAULT", const std::string& region = "RGN_NOBNDRY");
DERIV_FUNC_REGION_ENUM_TO_STRING(D2DYDZ, Field3D)

/// Calculate mixed partial derivative in y and z
///
///   \f$\partial^2 / \partial y \partial z\f$
///
/// @param[in] f       The field to be differentiated
/// @param[in] outloc  The cell location where the result is desired. If
///                    staggered grids is not enabled then this has no effect
///                    If not given, defaults to CELL_DEFAULT
/// @param[in] method  Differencing method to use. This overrides the default
///                    If not given, defaults to DIFF_DEFAULT
/// @param[in] region  What region is expected to be calculated
///                    If not given, defaults to RGN_NOBNDRY
Field2D D2DYDZ(const Field2D& f, CELL_LOC outloc = CELL_DEFAULT, const std::string&
    method = "DEFAULT", const std::string& region = "RGN_NOBNDRY");
DERIV_FUNC_REGION_ENUM_TO_STRING(D2DYDZ, Field2D)

#undef DERIV_FUNC_REGION_ENUM_TO_STRING
#undef VDERIV_FUNC_REGION_ENUM_TO_STRING

#endif // __DERIVS_H__
