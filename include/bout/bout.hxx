/*!************************************************************************
 *
 * @mainpage BOUT++
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

#ifndef BOUT_H
#define BOUT_H

#include "bout/build_config.hxx"

#include "bout/boutcomm.hxx"
#include "bout/difops.hxx" // Differential operators
#include "bout/field2d.hxx"
#include "bout/field3d.hxx"
#include "bout/globals.hxx"
#include "bout/mesh.hxx"
#include "bout/options_io.hxx"
#include "bout/output.hxx"
#include "bout/smoothing.hxx" // Smoothing functions
#include "bout/solver.hxx"
#include "bout/sourcex.hxx" // source and mask functions
#include "bout/utils.hxx"
#include "bout/vecops.hxx" // Vector differential operations
#include "bout/vector2d.hxx"
#include "bout/vector3d.hxx"
#include "bout/version.hxx"
#include "bout/where.hxx"

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
int BoutInitialise(int& argc, char**& argv);

namespace bout {
namespace experimental {
/// Function type for handling signals
using SignalHandler = void (*)(int);

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
  /// The "canonicalised" command line arguments, with single-letter
  /// arguments expanded
  std::vector<std::string> argv;
};

/// Parse the "fixed" command line arguments, like --help and -d
CommandLineArgs parseCommandLineArgs(int argc, char** argv);

/// Throw an exception if \p data_dir is either not a directory or not
/// accessible. We do not check whether we can write, as it is
/// sufficient that the files we need are writeable
void checkDataDirectoryIsAccessible(const std::string& data_dir);

/// Set up the output: open the log file for each processor, enable or
/// disable the default outputs based on \p verbosity, disable writing
/// to stdout for \p MYPE != 0
void setupOutput(const std::string& data_dir, const std::string& log_file, int verbosity,
                 int MYPE = 0);

/// Save the process ID for processor N = \p MYPE to file ".BOUT.pid.N"
/// in \p data_dir, so it can be shut down by user signal
///
/// Throws if it was not possible to create the file
void savePIDtoFile(const std::string& data_dir, int MYPE);

/// Print the initial header
void printStartupHeader(int MYPE, int NPES);

/// Print the compile-time options
void printCompileTimeOptions();

/// Print the arguments given on the command line
void printCommandLineArguments(const std::vector<std::string>& original_argv);

/// Setup the pipe etc and run stdout through bout-log-color. Return
/// true if it was successful
bool setupBoutLogColor(bool color_output, int MYPE);

/// Set BOUT++ version information, along with current time (as
/// `started`), into `run` section of \p options
void setRunStartInfo(Options& options);

/// Set the current time (as `finished`) into `run` section of \p
/// options
void setRunFinishInfo(Options& options);

/// Write \p options to \p settings_file in directory \p data_dir
void writeSettingsFile(Options& options, const std::string& data_dir,
                       const std::string& settings_file);

/// Add the configure-time build options to \p options
void addBuildFlagsToOptions(Options& options);
} // namespace experimental
} // namespace bout

/*!
 * Monitor class for output. Called by the solver every output timestep.
 *
 * This is added to the solver in bout/physicsmodel.hxx
 */
class BoutMonitor : public Monitor {
public:
  BoutMonitor(BoutReal timestep = -1);
  BoutMonitor(BoutReal timestep, Options& options);

private:
  int call(Solver* solver, BoutReal t, int iter, int NOUT) override;
  RunMetrics run_data;
  /// Wall time limit in seconds
  BoutReal wall_limit;
  /// Starting time
  BoutReal mpi_start_time;
  /// Stop if file `stop_check_name` exists
  bool stop_check;
  /// Filename for `stop_check`
  std::string stop_check_name;
};

/*!
 * BOUT++ finalisation. This should be called at the
 * end of the program.
 *
 * Frees memory, flushes buffers, and closes files.
 * If BOUT++ initialised MPI or external libraries,
 * then these are also finalised.
 *
 * If \p write_settings is true, output the settings, showing which
 * options were used. This overwrites the file written during
 * initialisation (BOUT.settings by default)
 *
 */
int BoutFinalise(bool write_settings = true);

#endif // BOUT_H
