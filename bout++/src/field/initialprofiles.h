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

#include "field3d.h"
#include "field2d.h"
#include "vector2d.h"
#include "vector3d.h"

int initial_profile(const char *name, Field3D &var);
int initial_profile(const char *name, Field2D &var);

int initial_profile(const char *name, Vector2D &var);
int initial_profile(const char *name, Vector3D &var);

// Generate a 3D field with a given Z oscillation
const Field3D genZMode(int n, real phase = 0.0);

#endif // __INITIALPROF_H__
