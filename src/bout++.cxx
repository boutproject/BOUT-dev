/**************************************************************************
 *
 * Main BOUT++ functions
 * Adapted from the BOUT code by B.Dudson, University of York, Oct 2007
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

const char DEFAULT_DIR[] = "data";
const char DEFAULT_OPT[] = "BOUT.inp";

// MD5 Checksum passed at compile-time
#define CHECKSUM1_(x) #x
#define CHECKSUM_(x) CHECKSUM1_(x)
#define CHECKSUM CHECKSUM_(MD5SUM)

// Revision passed at compile time
#define REV1_(x) #x
#define REV_(x) REV1_(x)
#define REV REV_(REVISION)

#define GLOBALORIGIN

#include "mpi.h"

#include <boutcomm.hxx>
#include <bout.hxx>
#include <datafile.hxx>
#include <bout/solver.hxx>
#include <boutexception.hxx>
#include <optionsreader.hxx>
#include <derivs.hxx>
#include <msg_stack.hxx>

#include <bout/sys/timer.hxx>

#include <boundary_factory.hxx>

#include <bout/petsclib.hxx>

#include <time.h>

#include <strings.h>
#include <string>
#include <list>
using std::string;
using std::list;
#include <sys/types.h>
#include <unistd.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef SIGHANDLE
#include <signal.h>
void bout_signal_handler(int sig);  // Handles segmentation faults
#endif

#include <output.hxx>

BoutReal simtime;
int iteration;

const string time_to_hms(BoutReal t);   // Converts to h:mm:ss.s format
char get_spin();                    // Produces a spinning bar

/*!
  Initialise BOUT++
  
  Inputs
  ------
  
  The command-line arguments argc and argv are passed by
  reference, and pointers to these will be stored in various
  places in BOUT++.
  
 */
void BoutInitialise(int &argc, char **&argv) {

  string dump_ext; ///< Extensions for restart and dump files

  const char *data_dir; ///< Directory for data input/output
  const char *opt_file; ///< Filename for the options file
  
#ifdef SIGHANDLE
  /// Set a signal handler for segmentation faults
  signal(SIGSEGV, bout_signal_handler);
#endif

  // Set default data directory
  data_dir = DEFAULT_DIR;
  opt_file = DEFAULT_OPT;

  /// Check command-line arguments
  /// NB: "restart" and "append" are now caught by options
  for (int i=1;i<argc;i++) {
    if (strncasecmp(argv[i], "-d", 2) == 0) {
      // Set data directory
      if (i+1 >= argc) {
        fprintf(stderr, "Useage is %s -d <data directory>\n", argv[0]);
        return;
      }
      i++;
      data_dir = argv[i];
    }
    if (strncasecmp(argv[i], "-f", 2) == 0) {
      // Set options file
      if (i+1 >= argc) {
        fprintf(stderr, "Useage is %s -f <options filename>\n", argv[0]);
        return;
      }
      i++;
      opt_file = argv[i];
    }
  }
  
  // Set options
  Options::getRoot()->set("datadir", string(data_dir));
  Options::getRoot()->set("optionfile", string(opt_file));

  // Set the command-line arguments
  PetscLib::setArgs(argc, argv); // PETSc initialisation
  Solver::setArgs(argc, argv);   // Solver initialisation
  BoutComm::setArgs(argc, argv); // MPI initialisation

  int NPES = BoutComm::size();
  int MYPE = BoutComm::rank();

  /// Set up the output
  if (MYPE == 0) output.enable(); // Enable writing to stdout
  else output.disable(); // No writing to stdout

  /// Open an output file to echo everything to
  /// On processor 0 anything written to output will go to stdout and the file
  output.open("%s/BOUT.log.%d", data_dir, MYPE);

  /// Print intro
  output.write("\nBOUT++ version %.2f\n", BOUT_VERSION);
#ifdef REVISION
  output.write("Revision: %s\n", REV);
#endif
#ifdef MD5SUM
  output.write("MD5 checksum: %s\n", CHECKSUM);
#endif
  output.write("Code compiled on %s at %s\n\n", __DATE__, __TIME__);
  output.write("B.Dudson (University of York), M.Umansky (LLNL) 2007\n");
  output.write("Based on BOUT by Xueqiao Xu, 1999\n\n");

  output.write("Processor number: %d of %d\n\n", MYPE, NPES);

  output.write("pid: %d\n\n",getpid());

  /// Print compile-time options

  output.write("Compile-time options:\n");

#ifdef CHECK
  output.write("\tChecking enabled, level %d\n", CHECK);
#else
  output.write("\tChecking disabled\n");
#endif

#ifdef SIGHANDLE
  output.write("\tSignal handling enabled\n");
#else
  output.write("\tSignal handling disabled\n");
#endif

#ifdef PDBF
  output.write("\tPDB support enabled\n");
#else
  output.write("\tPDB support disabled\n");
#endif

#ifdef NCDF
  output.write("\tnetCDF support enabled\n");
#else
  output.write("\tnetCDF support disabled\n");
#endif

#ifdef PNCDF
  output.write("\tParallel NetCDF support enabled\n");
#else
  output.write("\tParallel NetCDF support disabled\n");
#endif

#ifdef _OPENMP
  output.write("\tOpenMP parallelisation enabled\n");
#else
  output.write("\tOpenMP parallelisation disabled\n");
#endif

#ifdef METRIC3D
  output.write("\tRUNNING IN 3D-METRIC MODE\n");
#endif

  /// Get the options tree
  Options *options = Options::getRoot();

  try {
    /// Load settings file
    OptionsReader *reader = OptionsReader::getInstance();
    reader->read(options, "%s/%s", data_dir, opt_file);

    // Get options override from command-line
    reader->parseCommandLine(options, argc, argv);
  }catch(BoutException *e) {
    output << "Error encountered during initialisation\n";
    output << e->what() << endl;
    return;
  }

  try {
    /////////////////////////////////////////////
    /// Get some settings

    // Check if restarting
    bool append;
    OPTION(options, append, false);

    /// Get file extensions
    options->get("dump_format", dump_ext, "nc");

    /// Setup derivative methods
    if (derivs_init()) {
      output.write("Failed to initialise derivative methods. Aborting\n");
      return;
    }

    ////////////////////////////////////////////

    // Set up the "dump" data output file
    output << "Setting up output (dump) file\n";

    if(!options->getSection("output")->isSet("floats"))
      options->getSection("output")->set("floats", true, "default"); // by default output floats

    dump = Datafile(options->getSection("output"));
    
    /// Open a file for the output
    if(append) {
      dump.opena("%s/BOUT.dmp.%s", data_dir, dump_ext.c_str());
    }else {
      dump.openw("%s/BOUT.dmp.%s", data_dir, dump_ext.c_str());
    }

    /// Add book-keeping variables to the output files
    dump.writeVar(BOUT_VERSION, "BOUT_VERSION");
    dump.add(simtime, "t_array", 1); // Appends the time of dumps into an array
    dump.add(iteration, "iteration", 0);

    ///////////////////////////////////////////////
    
    mesh = Mesh::create();  ///< Create the mesh
    mesh->load();           ///< Load from sources. Required for Field initialisation
    mesh->outputVars(dump); ///< Save mesh configuration into output file
    
  }catch(BoutException &e) {
    output.write("Error encountered during initialisation: %s\n", e.what());
    BoutComm::cleanup();
    throw;
  }
}

int bout_run(Solver *solver, rhsfunc physics_run) {
  
  /// Set the RHS function
  solver->setRHS(physics_run);
  
  /// Add the monitor function
  solver->addMonitor(bout_monitor);

  /// Run the simulation
  return solver->solve();
}

int BoutFinalise() {
  // Delete the mesh
  delete mesh;

  // Close the output file
  dump.close();

  // Delete field memory
  Field2D::cleanup();
  Field3D::cleanup();

  // Cleanup boundary factory
  BoundaryFactory::cleanup();
  
  // Cleanup timer
  Timer::cleanup();

  // Options tree
  Options::cleanup();
  OptionsReader::cleanup();

  // Debugging message stack
  msg_stack.clear();

  // Call PetscFinalize if not already called
  PetscLib::cleanup();
    
  // Logging output
  Output::cleanup();

  // MPI communicator, including MPI_Finalize()
  BoutComm::cleanup();
  
  return 0;
}

/*!*************************************************************************
 * SOLUTION MONITOR FUNCTION
 *
 * Called each timestep by the solver
 **************************************************************************/

int bout_monitor(Solver *solver, BoutReal t, int iter, int NOUT) {
  // Data used for timing
  static bool first_time = true;
  static BoutReal wall_limit, mpi_start_time; // Keep track of remaining wall time

#ifdef CHECK
  int msg_point = msg_stack.push("bout_monitor(%e, %d, %d)", t, iter, NOUT);
#endif

  // Set the global variables. This is done because they need to be
  // written to the output file before the first step (initial condition)
  simtime = t;
  iteration = iter;

  /// Write dump file
  dump.write();

  /// Collect timing information
  BoutReal wtime        = Timer::resetTime("run");
  int ncalls            = solver->rhs_ncalls;
  BoutReal wtime_rhs    = Timer::resetTime("rhs");
  BoutReal wtime_invert = Timer::resetTime("invert");
  BoutReal wtime_comms  = Timer::resetTime("comms");  // Time spent communicating (part of RHS)
  BoutReal wtime_io     = Timer::resetTime("io");      // Time spend on I/O

  output.print("\r"); // Only goes to screen

  if (first_time) {
    /// First time the monitor has been called

    /// Get some options
    Options *options = Options::getRoot();
    OPTION(options, wall_limit, -1.0); // Wall time limit. By default, no limit
    wall_limit *= 60.0*60.0;  // Convert from hours to seconds

    /// Record the starting time
    mpi_start_time = MPI_Wtime(); // NB: Miss time for first step (can be big!)

    first_time = false;

    /// Print the column header for timing info
    output.write("Sim Time  |  RHS evals  | Wall Time |  Calc    Inv   Comm    I/O   SOLVER\n\n");

  }
  
  output.write("%.3e      %5d       %.2e   %5.1f  %5.1f  %5.1f  %5.1f  %5.1f\n", 
               simtime, ncalls, wtime,
               100.0*(wtime_rhs - wtime_comms - wtime_invert)/wtime,
               100.*wtime_invert/wtime,  // Inversions
               100.0*wtime_comms/wtime,  // Communications
               100.* wtime_io / wtime,      // I/O
               100.*(wtime - wtime_io - wtime_rhs)/wtime); // Everything else
  
  
  // This bit only to screen, not log file

  BoutReal t_elapsed = MPI_Wtime() - mpi_start_time;
  output.print("%c  Step %d of %d. Elapsed %s", get_spin(), iteration+1, NOUT, (time_to_hms(t_elapsed)).c_str());
  output.print(" ETA %s", (time_to_hms(wtime * ((BoutReal) (NOUT - iteration - 1)))).c_str());

  if (wall_limit > 0.0) {
    // Check if enough time left

    BoutReal t_remain = mpi_start_time + wall_limit - MPI_Wtime();
    if (t_remain < wtime) {
      // Less than 1 time-step left
      output.write("Only %e seconds left. Quitting\n", t_remain);

#ifdef CHECK
      msg_stack.pop(msg_point);
#endif
      return 1; // Return an error code to quit
    } else {
      output.print(" Wall %s", (time_to_hms(t_remain)).c_str());
    }
  }
  
#ifdef CHECK
  msg_stack.pop(msg_point);
#endif

  return 0;
}

/**************************************************************************
 * Global error handling
 **************************************************************************/

/// Print an error message and exit
void bout_error() {
  bout_error(NULL);
}

void bout_error(const char *str) {
  output.write("****** ERROR CAUGHT ******\n");

  if (str != NULL) output.write(str);

  output.write("\n");

#ifdef CHECK
  msg_stack.dump();
#endif

  MPI_Abort(BoutComm::get(), 1);

  exit(1);
}

#ifdef SIGHANDLE
/// Signal handler - catch segfaults
void bout_signal_handler(int sig) {
  /// Set signal handler back to default to prevent possible infinite loop
  signal(SIGSEGV, SIG_DFL);

  output.write("\n****** SEGMENTATION FAULT CAUGHT ******\n\n");
#ifdef CHECK
  /// Print out the message stack to help debugging
  msg_stack.dump();
#else
  output.write("Enable checking (-DCHECK flag) to get a trace\n");
#endif

  exit(sig);
}
#endif

/**************************************************************************
 * Utilities
 **************************************************************************/

/// Write a time in h:mm:ss.s format
const string time_to_hms(BoutReal t) {
  int h, m;

  h = (int) (t / 3600); t -= 3600.*((BoutReal) h);
  m = (int) (t / 60);   t -= 60 * ((BoutReal) m);

  char buffer[256];
  sprintf(buffer,"%d:%02d:%04.1f", h, m, t);

  return string(buffer);
}

/// Produce a spinning bar character
char get_spin() {
  static int i = 0;
  char c = '|'; // Doesn't need to be assigned; squash warning

  switch(i) {
  case 0: 
    c = '|'; break;
  case 1:
    c = '/'; break;
  case 2:
    c = '-'; break;
  case 3:
    c = '\\'; break;
  }
  i = (i+1) % 4;
  return c;
}
