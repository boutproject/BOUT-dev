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
/// @param[in] dfdy_region Region in which to calculate df/dy. If an empty string (default)
///                        then the same as the region for the calculation as a whole.
///                        If dfdy_region < region in size then this will cause errors.
Field3D D2DXDY(const Field3D& f, CELL_LOC outloc = CELL_DEFAULT,
               const std::string& method = "DEFAULT",
               const std::string& region = "RGN_NOBNDRY",
               const std::string& dfdy_boundary_condition = "free_o3",
               const std::string& dfdy_region = "");

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
/// @param[in] dfdy_region Region in which to calculate df/dy. If an empty string
///                        (default) then the same as the region for the calculation as a
///                        whole. If dfdy_region < region in size then this will cause
///                        errors.
Field2D D2DXDY(const Field2D& f, CELL_LOC outloc = CELL_DEFAULT,
               const std::string& method = "DEFAULT",
               const std::string& region = "RGN_NOBNDRY",
               const std::string& dfdy_boundary_condition = "free_o3",
               const std::string& dfdy_region = "");

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


#endif // __DERIVS_H__
