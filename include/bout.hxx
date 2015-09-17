/*!************************************************************************
 *
 * @mainpage BOUT++
 * 
 * @version 3.0
 * 
 * @par Description 
 * Framework for the solution of partial differential
 * equations, in particular fluid models in plasma physics.
 *
 * @par Include files
 * - bout++.hxx includes commonly used routines and classes
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

#include "boutcomm.hxx"

#include "globals.hxx"

#include "field2d.hxx"
#include "field3d.hxx"
#include "vector2d.hxx"
#include "vector3d.hxx"

#include "difops.hxx" // Differential operators

#include "vecops.hxx" // Vector differential operations

#include "smoothing.hxx" // Smoothing functions

#include "sourcex.hxx"     // source and mask functions

#include "bout/solver.hxx"

#include "datafile.hxx"

#include "where.hxx"

#include "output.hxx"

#include "utils.hxx"

const BoutReal BOUT_VERSION = 3;  ///< Version number

// BOUT++ main functions

/*!
 * BOUT++ initialisation. This function must be
 * called first, passing command-line arguments.
 * 
 * This will call MPI_Initialize, and if BOUT++
 * has been configured with external libraries such as
 * PETSc then these will be initialised as well.
 * 
 */
int BoutInitialise(int &argc, char **&argv);

int bout_run(Solver *solver, rhsfunc physics_run);

int bout_monitor(Solver *solver, BoutReal t, int iter, int NOUT); 

/*!
 * BOUT++ finalisation. This should be called at the
 * end of the program.
 *
 * Frees memory, flushes buffers, and closes files.
 * If BOUT++ initialised MPI or external libraries,
 * then these are also finalised.
 */
int BoutFinalise();

#endif // __BOUT_H__
