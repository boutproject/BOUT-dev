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

const BoutReal BOUT_VERSION = BOUT_VERSION_DOUBLE;  ///< Version number

// BOUT++ main functions

/*!
 * BOUT++ initialisation. This function must be
 * called first, passing command-line arguments.
 * 
 * This will call MPI_Initialize, and if BOUT++
 * has been configured with external libraries such as
 * PETSc then these will be initialised as well.
 * 
 * Example
 * -------
 *
 * A minimal BOUT++ program consists of:
 *
 *     int main(int argc, char** argv) {
 *       BoutInitialise(argc, argv);
 *       
 *       BoutFinalise();
 *     }
 *
 * Usually this function is called in a standard main() function,
 * either by including boutmain.hxx or by including bout/physicsmodel.hxx
 * and using the BOUTMAIN macro.
 *     
 */
int BoutInitialise(int &argc, char **&argv);

/*!
 * Run the given solver. This function is only used
 * for old-style physics models with standalone C functions
 * The main() function in boutmain.hxx calls this function
 * to set up the RHS function and add bout_monitor.
 * 
 */
int bout_run(Solver *solver, rhsfunc physics_run);

/*!
 * Monitor class for output. Called by the solver every output timestep.
 * 
 * This is added to the solver in bout_run (for C-style models)
 * or in bout/physicsmodel.hxx
 */
class BoutMonitor: public Monitor{
  int call(Solver *solver, BoutReal t, int iter, int NOUT) override;
};

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
