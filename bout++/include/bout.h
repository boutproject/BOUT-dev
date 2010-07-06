/*!************************************************************************
 *
 * \file bout.h
 * \brief File included into the physics code
 *
 * Just includes commonly needed definitions from other include files
 *
 **************************************************************************
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
 *
 **************************************************************************/

#ifndef __BOUT_H__
#define __BOUT_H__

#include "globals.h"

#include "field2d.h"
#include "field3d.h"
#include "vector2d.h"
#include "vector3d.h"

#include "grid.h" // For reading uedge.grd files

#include "difops.h" // Differential operators

#include "vecops.h" // Vector differential operations

#include "smoothing.h" // Smoothing functions

#include "sourcex.h"     // source and mask functions

#include "solver.h"

#include "datafile.h"

#include "where.h"

const BoutReal BOUT_VERSION = 0.85;  ///< Version number

// BOUT++ functions (bout++.cpp). Call to add a variable to evolve

void bout_solve(Field2D &var, const char *name);
void bout_solve(Field3D &var, const char *name);
void bout_solve(Vector2D &var, const char *name);
void bout_solve(Vector3D &var, const char *name);

bool bout_constrain(Field3D &var, Field3D &F_var, const char *name);

// Physics functions

int physics_init(bool restarting);
int physics_run(BoutReal t);

#ifndef GLOBALORIGIN
#define GLOBAL extern
#else
#define GLOBAL
#endif

// Solver object
GLOBAL Solver *solver;    // Interface to PVODE

#undef GLOBAL

#endif // __BOUT_H__
