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

#include <bout/mesh.hxx>
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

#ifndef BOUT_NO_USING_NAMESPACE_BOUTGLOBALS
// Include using statement by default in user code.
// Macro allows us to include bout.hxx or physicsmodel.hxx without the using
// statement in library code.
using namespace bout::globals;
#endif // BOUT_NO_USING_NAMESPACE_BOUTGLOBALS

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

namespace bout {
namespace experimental {
/// Function type for handling signals
using SignalHandler = void(*)(int);

/// Set a signal handler for user-requested clean exit, and
/// (optionally) segmentation faults and floating point errors
///
/// - For segmentation faults, compile with `--enable-signal`.
/// - For floating point errors, compile with `--enable-sigfpe`
void setupSignalHandler(SignalHandler signal_handler);

/// The default BOUT++ signal handler: throw an exception with an
/// appropriate message
void defaultSignalHandler(int sig);

/// Set up the i18n environment
void setupGetText();

/// Results of parsing the command line arguments
struct CommandLineArgs {
  int verbosity{4};
  bool color_output{false};
  std::string data_dir{"data"};          ///< Directory for data input/output
  std::string opt_file{"BOUT.inp"};      ///< Filename for the options file
  std::string set_file{"BOUT.settings"}; ///< Filename for the options file
  std::string log_file{"BOUT.log"};      ///< File name for the log file
  /// The original set of command line arguments
  std::vector<std::string> original_argv;
};

/// Parse the "fixed" command line arguments, like --help and -d
CommandLineArgs parseCommandLineArgs(int argc, char** argv);

/// Throw an exception if \p data_dir is either not a directory or not
/// accessible. We do not check whether we can write, as it is
/// sufficient that the files we need are writeable
void checkDataDirectoryIsAccessible(const std::string& data_dir);

/// Print the initial header
void printStartupHeader(int MYPE, int NPES);

/// Print the compile-time options
void printCompileTimeOptions();

/// Print the arguments given on the command line
void printCommandLineArguments(const std::vector<std::string>& original_argv);

/// Setup the pipe etc and run stdout through bout-log-color. Return
/// true if it was successful
bool setupBoutLogColor(bool color_output, int MYPE);
} // namespace experimental
} // namespace bout

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
class BoutMonitor: public Monitor {
public:
  BoutMonitor(BoutReal timestep = -1) : Monitor(timestep) {
    // Add wall clock time etc to dump file
    run_data.outputVars(bout::globals::dump);
  }
private:
  int call(Solver* solver, BoutReal t, int iter, int NOUT) override;
  RunMetrics run_data;
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
