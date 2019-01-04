/*!*************************************************************************
 * \file where.hxx
 * 
 * A set of functions which choose between two values
 *
 **************************************************************************
 * Copyright 2010 B.D.Dudson, S.Farley, M.V.Umansky, X.Q.Xu
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

#ifndef __WHERE_H__
#define __WHERE_H__

#include "bout_types.hxx"

class Field2D;
class Field3D;

/// For each point, choose between two inputs based on a third input
///
/// @param[in] test   The value which determines which input to use
/// @param[in] gt0    Uses this value if test > 0.0
/// @param[in] le0    Uses this value if test <= 0.0
const Field3D where(const Field2D &test, const Field3D &gt0, const Field3D &le0);
const Field3D where(const Field2D &test, const Field3D &gt0, BoutReal le0);
const Field3D where(const Field2D &test, BoutReal gt0, const Field3D &le0);
const Field3D where(const Field2D &test, const Field3D &gt0, const Field2D &le0);
const Field3D where(const Field2D &test, const Field2D &gt0, const Field3D &le0);

/// For each point, choose between two inputs based on a third input
///
/// @param[in] test   The value which determines which input to use
/// @param[in] gt0    Uses this value if test > 0.0
/// @param[in] le0    Uses this value if test <= 0.0
const Field2D where(const Field2D &test, const Field2D &gt0, const Field2D &le0);
const Field2D where(const Field2D &test, const Field2D &gt0, BoutReal le0);
const Field2D where(const Field2D &test, BoutReal gt0, const Field2D &le0);
const Field2D where(const Field2D &test, BoutReal gt0, BoutReal le0);

/// For each point, choose between two inputs based on a third input
///
/// @param[in] test   The value which determines which input to use
/// @param[in] gt0    Uses this value if test > 0.0
/// @param[in] le0    Uses this value if test <= 0.0
const Field3D where(const Field3D &test, BoutReal gt0, const Field3D &le0);

#endif // __WHERE_H__

