/**************************************************************************
 * Sets initial profiles
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

#ifndef __INITIALPROF_H__
#define __INITIALPROF_H__

#include <string>

class Field3D;
class Field2D;
class Vector2D;
class Vector3D;

/*!
 * Set a field from options
 *
 * This is called by Solver for each evolving field at the beginning
 * of a simulation.
 *
 *
 * @param[in] name   The name of the field. This will be used
 *                   as the section name in the options
 * @param[out] var   The field, which will be set to a value
 *                   depending on the options
 *
 *
 * To create the value, it looks for a setting "function"
 * in a section called name. If that is not found, then it looks
 * for "function" in a section called "All". If that is also not
 * found, then the value defaults to zero.
 *
 * A second variable, "scale", can be used to multiply the function,
 * and defaults to 1.0
 *
 * Example
 * -------
 * Given the input file:
 *
 * [All]
 * function = sin(y)
 *
 * [pressure]
 *
 * [density]
 * scale = 0.2
 *
 * [vorticity]
 * function = cos(y)
 *
 * initial_profile would generate:
 *
 * o pressure -> sin(y)
 * o density  -> 0.2*sin(y)
 * o vorticity -> cos(y)
 *
 */
void initial_profile(const std::string& name, Field3D& var);

/*!
 * Set a Field2D from options
 *
 * @param[in] name   The name of the field, used as a section name
 *
 * @param[out] var   The field which will be set to a value
 */
void initial_profile(const std::string& name, Field2D& var);

/*!
 * Set a vector to a value. The options used depend
 * on whether the vector is covariant or contravariant.
 *
 * If covariant, then each component will be initialised
 * by adding "_x", "_y", "_z" to the name.
 *
 * If contravariant, then each component will be initialised
 * by adding "x", "y" and "z" to the name.
 */
void initial_profile(const std::string& name, Vector2D& var);

/*!
 * Set a vector to a value. The options used depend
 * on whether the vector is covariant or contravariant.
 *
 * If covariant, then each component will be initialised
 * by adding "_x", "_y", "_z" to the name.
 *
 * If contravariant, then each component will be initialised
 * by adding "x", "y" and "z" to the name.
 *
 *
 */
void initial_profile(const std::string& name, Vector3D& var);

#endif // __INITIALPROF_H__
