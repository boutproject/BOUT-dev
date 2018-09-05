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
const char DEFAULT_SET[] = "BOUT.settings";
const char DEFAULT_LOG[] = "BOUT.log";

// MD5 Checksum passed at compile-time
#define CHECKSUM1_(x) #x
#define CHECKSUM_(x) CHECKSUM1_(x)
#define CHECKSUM CHECKSUM_(MD5SUM)

// Revision passed at compile time
#define REV1_(x) #x
#define REV_(x) REV1_(x)
#define REV REV_(REVISION)

#define GLOBALORIGIN


#define INDIRECT1_BOUTMAIN(a) #a
#define INDIRECT0_BOUTMAIN(...) INDIRECT1_BOUTMAIN(#__VA_ARGS__)
#define STRINGIFY(a) INDIRECT0_BOUTMAIN(a)

#include "mpi.h"

#include <boutcomm.hxx>
#include <bout.hxx>
#include <datafile.hxx>
#include <bout/solver.hxx>
#include <boutexception.hxx>
#include <optionsreader.hxx>
#include <msg_stack.hxx>

#include <bout/sys/timer.hxx>

#include <boundary_factory.hxx>

#include <invert_laplace.hxx>

#include <bout/slepclib.hxx>
#include <bout/petsclib.hxx>

#include <time.h>

#include <strings.h>
#include <string>
#include <list>
using std::string;
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <signal.h>
void bout_signal_handler(int sig);  // Handles signals
#ifdef BOUT_FPE
#include <fenv.h>
#endif


#include <output.hxx>

BoutReal simtime;
int iteration;
bool user_requested_exit=false;

const string time_to_hms(BoutReal t);   // Converts to h:mm:ss.s format
char get_spin();                    // Produces a spinning bar

/*!
  Initialise BOUT++

  Inputs
  ------

  The command-line arguments argc and argv are passed by
  reference, and pointers to these will be stored in various
  places in BOUT++.

  Outputs
  -------

  Any non-zero return value should halt the simulation. If the return value is
  less than zero, the exit status from BOUT++ is 0, otherwise it is the return
  value of BoutInitialise.

 */
int BoutInitialise(int &argc, char **&argv) {

  string dump_ext; ///< Extensions for restart and dump files

  const char *data_dir; ///< Directory for data input/output
  const char *opt_file; ///< Filename for the options file
  const char *set_file; ///< Filename for the options file
  const char *log_file; ///< File name for the log file

#ifdef SIGHANDLE
  /// Set a signal handler for segmentation faults
  signal(SIGSEGV, bout_signal_handler);
#endif
#ifdef BOUT_FPE
  signal(SIGFPE,  bout_signal_handler);
  feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
#endif

  /// Trap SIGUSR1 to allow a clean exit after next write
  signal(SIGUSR1, bout_signal_handler);

  // Set default data directory
  data_dir = DEFAULT_DIR;
  opt_file = DEFAULT_OPT;
  set_file = DEFAULT_SET;
  log_file = DEFAULT_LOG;

  int verbosity=4;
  /// Check command-line arguments
  /// NB: "restart" and "append" are now caught by options
  /// Check for help flag separately
  for (int i=1;i<argc;i++) {
    if (string(argv[i]) == "-h" ||
    	string(argv[i]) == "--help") {
      // Print help message -- note this will be displayed once per processor as we've not started MPI yet.
      fprintf(stdout, "Usage: %s [-d <data directory>] [-f <options filename>] [restart [append]] [VAR=VALUE]\n", argv[0]);
      fprintf(stdout,
              "\n"
              "  -d <data directory>\tLook in <data directory> for input/output files\n"
              "  -f <options filename>\tUse OPTIONS given in <options filename>\n"
              "  -o <settings filename>\tSave used OPTIONS given to <options filename>\n"
              "  -l, --log <log filename>\tPrint log to <log filename>\n"
              "  -v, --verbose\t\tIncrease verbosity\n"
              "  -q, --quiet\t\tDecrease verbosity\n"
#ifdef LOGCOLOR
              "  -c, --color\t\tColor output using bout-log-color\n"
#endif
              "  -h, --help\t\tThis message\n"
              "  restart [append]\tRestart the simulation. If append is specified, "
              "append to the existing output files, otherwise overwrite them\n"
              "  VAR=VALUE\t\tSpecify a VALUE for input parameter VAR\n"
              "\nFor all possible input parameters, see the user manual and/or the "
              "physics model source (e.g. %s.cxx)\n",
              argv[0]);

      return -1;
    }
  }
  bool color_output = false; // Will be set true if -c is in the options
  for (int i=1;i<argc;i++) {
    if (string(argv[i]) == "-d") {
      // Set data directory
      if (i+1 >= argc) {
        fprintf(stderr, "Usage is %s -d <data directory>\n", argv[0]);
        return 1;
      }
      i++;
      data_dir = argv[i];
      
    } else if (string(argv[i]) == "-f") {
      // Set options file
      if (i+1 >= argc) {
        fprintf(stderr, "Usage is %s -f <options filename>\n", argv[0]);
        return 1;
      }
      i++;
      opt_file = argv[i];
      
    } else if (string(argv[i]) == "-o") {
      // Set options file
      if (i+1 >= argc) {
        fprintf(stderr, "Usage is %s -o <settings filename>\n", argv[0]);
        return 1;
      }
      i++;
      set_file = argv[i];

    } else if ((string(argv[i]) == "-l") || (string(argv[i]) == "--log")) {
      if (i + 1 >= argc) {
        fprintf(stderr, "Usage is %s -l <log filename>\n", argv[0]);
        return 1;
      }
      i++;
      log_file = argv[i];

    } else if ( (string(argv[i]) == "-v") ||
                (string(argv[i]) == "--verbose") ){
      verbosity++;
      
    } else if ( (string(argv[i]) == "-q") ||
                (string(argv[i]) == "--quiet")) {
      verbosity--;
      
    } else if ( (string(argv[i]) == "-c") ||
                (string(argv[i]) == "--color") ) {
      // Add color to the output by piping through bout-log-color
      // This is done after checking all command-line inputs
      // in case -c is set multiple times
      color_output = true;
    }
  }
  
  if (std::string(set_file) == std::string(opt_file)){
    throw BoutException("Input and output file for settings must be different.\nProvide -o <settings file> to avoid this issue.\n");
  }

  // Check that data_dir exists. We do not check whether we can write, as it is
  // sufficient that the files we need are writeable ...
  struct stat test;
  if (stat(data_dir, &test) == 0){
    if (!S_ISDIR(test.st_mode)){
      throw BoutException("DataDir \"%s\" is not a directory\n",data_dir);
    }
  } else {
    throw BoutException("DataDir \"%s\" does not exist or is not accessible\n",data_dir);
  }
  
  // Set options
  Options::getRoot()->set("datadir", string(data_dir));
  Options::getRoot()->set("optionfile", string(opt_file));
  Options::getRoot()->set("settingsfile", string(set_file));

  // Set the command-line arguments
  SlepcLib::setArgs(argc, argv); // SLEPc initialisation
  PetscLib::setArgs(argc, argv); // PETSc initialisation
  Solver::setArgs(argc, argv);   // Solver initialisation
  BoutComm::setArgs(argc, argv); // MPI initialisation

  int NPES = BoutComm::size();
  int MYPE = BoutComm::rank();
  
#ifdef LOGCOLOR
  if (color_output && (MYPE == 0)) {
    // Color stdout by piping through bout-log-color script
    // Only done on processor 0, since this is the only processor which writes to stdout
    // This uses popen, fileno and dup2 functions, which are POSIX
    bool success = false;

    // Run bout-log-color through the shell. This should share stdout with BOUT++,
    // and read stdin from the pipe
    FILE *outpipe = popen("bout-log-color", "w");

    if (outpipe != nullptr) {
      // Valid pipe
      // Get the integer file descriptor
      int fno = fileno(outpipe);
      if (fno != -1) {
        // Valid file descriptor

        // Note: We can get to here if bout-log-color failed to run
        // This seems to cause code to fail later
        
        // Replace stdout with the pipe.
        int status = dup2(fno, STDOUT_FILENO);
        if (status != -1) {
          success = true;
        }
      }
    }
    if (!success) {
      // Failed . Probably not important enough to stop the simulation
      fprintf(stderr, "Could not run bout-log-color. Make sure it is in your PATH\n");
    }
  }
#endif // LOGCOLOR
  
  /// Set up the output, accessing underlying Output object
  {
    Output &output = *Output::getInstance();
    if (MYPE == 0) output.enable(); // Enable writing to stdout
    else output.disable(); // No writing to stdout

    /// Open an output file to echo everything to
    /// On processor 0 anything written to output will go to stdout and the file
    if (output.open("%s/%s.%d", data_dir, log_file, MYPE)) {
      return 1;
    }
  }
  
  output_error.enable(verbosity>0);
  output_warn.enable(verbosity>1);
  output_progress.enable(verbosity>2);
  output_info.enable(verbosity>3);
  output_debug.enable(verbosity>4); //Only actually enabled if also compiled with DEBUG
  
  // The backward-compatible output object same as output_progress
  output.enable(verbosity>2);

  // Save the PID of this process to file, so it can be shut down by user signal
  {
    std::string filename;
    std::stringstream(filename) << data_dir << "/.BOUT.pid." << MYPE;
    std::ofstream pid_file;
    pid_file.open(filename, std::ios::out);
    if (pid_file.is_open()) {
      pid_file << getpid() << "\n";
      pid_file.close();
    }
  }
  
  /// Print intro
  output_progress.write("BOUT++ version %s\n", BOUT_VERSION_STRING);
#ifdef REVISION
  output_progress.write("Revision: %s\n", REV);
#endif
#ifdef MD5SUM
  output_progress.write("MD5 checksum: %s\n", CHECKSUM);
#endif
  output_progress.write("Code compiled on %s at %s\n\n", __DATE__, __TIME__);
  output_info.write("B.Dudson (University of York), M.Umansky (LLNL) 2007\n");
  output_info.write("Based on BOUT by Xueqiao Xu, 1999\n\n");

  output_info.write("Processor number: %d of %d\n\n", MYPE, NPES);

  output_info.write("pid: %d\n\n",getpid());

  /// Print compile-time options

  output_info.write("Compile-time options:\n");

#if CHECK > 0
  output_info.write("\tChecking enabled, level %d\n", CHECK);
#else
  output_info.write("\tChecking disabled\n");
#endif

#ifdef SIGHANDLE
  output_info.write("\tSignal handling enabled\n");
#else
  output_info.write("\tSignal handling disabled\n");
#endif

#ifdef NCDF
  output_info.write("\tnetCDF support enabled\n");
#else
#ifdef NCDF4
  output_info.write("\tnetCDF4 support enabled\n");
#else
  output_info.write("\tnetCDF support disabled\n");
#endif
#endif

#ifdef PNCDF
  output_info.write("\tParallel NetCDF support enabled\n");
#else
  output_info.write("\tParallel NetCDF support disabled\n");
#endif

#ifdef _OPENMP
  output_info.write("\tOpenMP parallelisation enabled, using %d threads\n",omp_get_max_threads());
#else
  output_info.write("\tOpenMP parallelisation disabled\n");
#endif

#ifdef METRIC3D
  output_info.write("\tRUNNING IN 3D-METRIC MODE\n");
#endif

#ifdef BOUT_FPE
  output_info.write("\tFloatingPointExceptions enabled\n");
#endif

  //The stringify is needed here as BOUT_FLAGS_STRING may already contain quoted strings
  //which could cause problems (e.g. terminate strings).
  output_info.write("\tCompiled with flags : %s\n",STRINGIFY(BOUT_FLAGS_STRING));
  
  /// Get the options tree
  Options *options = Options::getRoot();

  try {
    /// Load settings file
    OptionsReader *reader = OptionsReader::getInstance();
    reader->read(options, "%s/%s", data_dir, opt_file);

    // Get options override from command-line
    reader->parseCommandLine(options, argc, argv);

    // Save settings
    if (BoutComm::rank() == 0) {
      reader->write(options, "%s/%s", data_dir, set_file);
    }
  } catch (BoutException &e) {
    output << "Error encountered during initialisation\n";
    output << e.what() << endl;
    return 1;
  }

  try {
    /////////////////////////////////////////////
    
    mesh = Mesh::create();  ///< Create the mesh
    mesh->load();           ///< Load from sources. Required for Field initialisation
    mesh->setParallelTransform(); ///< Set the parallel transform from options
    /////////////////////////////////////////////
    /// Get some settings

    // Check if restarting
    bool append;
    OPTION(options, append, false);

    /// Get file extensions
    options->get("dump_format", dump_ext, "nc");
    
    ////////////////////////////////////////////

    // Set up the "dump" data output file
    output << "Setting up output (dump) file\n";

    dump = Datafile(options->getSection("output"));
    
    /// Open a file for the output
    if(append) {
      dump.opena("%s/BOUT.dmp.%s", data_dir, dump_ext.c_str());
    }else {
      dump.openw("%s/BOUT.dmp.%s", data_dir, dump_ext.c_str());
    }

    /// Add book-keeping variables to the output files
    dump.add(const_cast<BoutReal&>(BOUT_VERSION), "BOUT_VERSION", false);
    dump.add(simtime, "t_array", true); // Appends the time of dumps into an array
    dump.add(iteration, "iteration", false);

    ////////////////////////////////////////////

    mesh->outputVars(dump); ///< Save mesh configuration into output file
    
  }catch(BoutException &e) {
    output_error.write("Error encountered during initialisation: %s\n", e.what());
    throw;
  }
  return 0;
}

int bout_run(Solver *solver, rhsfunc physics_run) {
  
  /// Set the RHS function
  solver->setRHS(physics_run);
  
  /// Add the monitor function
  Monitor * bout_monitor = new BoutMonitor();
  solver->addMonitor(bout_monitor, Solver::BACK);

  /// Run the simulation
  return solver->solve();
}

int BoutFinalise() {

  // Output the settings, showing which options were used
  // This overwrites the file written during initialisation
  try {
    if (BoutComm::rank() == 0) {
      string data_dir;
      Options::getRoot()->get("datadir", data_dir, "data");

      OptionsReader *reader = OptionsReader::getInstance();
      std::string settingsfile;
      OPTION(Options::getRoot(), settingsfile, "");
      reader->write(Options::getRoot(), "%s/%s", data_dir.c_str(), settingsfile.c_str());
    }
  } catch (BoutException &e) {
    output_error << "Error whilst writing settings" << endl;
    output_error << e.what() << endl;
  }

  // Delete the mesh
  delete mesh;

  // Close the output file
  dump.close();

  // Make sure all processes have finished writing before exit
  MPI_Barrier(BoutComm::get());

  // Laplacian inversion
  Laplacian::cleanup();

  // Delete field memory
  Array<BoutReal>::cleanup();
  Array<dcomplex>::cleanup();
  Array<fcmplx>::cleanup();
  Array<int>::cleanup();
  Array<unsigned long>::cleanup();
  
  // Cleanup boundary factory
  BoundaryFactory::cleanup();
  
  // Cleanup timer
  Timer::cleanup();

  // Options tree
  Options::cleanup();
  OptionsReader::cleanup();

  // Debugging message stack
  msg_stack.clear();

  // Call SlepcFinalize if not already called
  SlepcLib::cleanup();

  // Call PetscFinalize if not already called
  PetscLib::cleanup();

  // MPI communicator, including MPI_Finalize()
  BoutComm::cleanup();
  
  return 0;
}

/*!*************************************************************************
 * SOLUTION MONITOR FUNCTION
 *
 * Called each timestep by the solver
 **************************************************************************/

int BoutMonitor::call(Solver *solver, BoutReal t, int iter, int NOUT) {
  TRACE("BoutMonitor::call(%e, %d, %d)", t, iter, NOUT);

  // Data used for timing
  static bool first_time = true;
  static BoutReal wall_limit, mpi_start_time; // Keep track of remaining wall time
  
  static bool stopCheck;       // Check for file, exit if exists?
  static std::string stopCheckName; // File checked, whose existence triggers a stop
  
  // Set the global variables. This is done because they need to be
  // written to the output file before the first step (initial condition)
  simtime = t;
  iteration = iter;

  /// Write dump file
  dump.write();

  /// Collect timing information
  BoutReal wtime        = Timer::resetTime("run");
  int ncalls            = solver->resetRHSCounter();
  int ncalls_e		= solver->resetRHSCounter_e();
  int ncalls_i		= solver->resetRHSCounter_i();

  bool output_split     = solver->splitOperator();
  BoutReal wtime_rhs    = Timer::resetTime("rhs");
  BoutReal wtime_invert = Timer::resetTime("invert");
  BoutReal wtime_comms  = Timer::resetTime("comms");  // Time spent communicating (part of RHS)
  BoutReal wtime_io     = Timer::resetTime("io");      // Time spend on I/O

  output_progress.print("\r"); // Only goes to screen

  if (first_time) {
    /// First time the monitor has been called

    /// Get some options
    Options *options = Options::getRoot();
    OPTION(options, wall_limit, -1.0); // Wall time limit. By default, no limit
    wall_limit *= 60.0 * 60.0;         // Convert from hours to seconds

    OPTION(options, stopCheck, false);
    if (stopCheck) {
      // Get name of file whose existence triggers a stop
      OPTION(options, stopCheckName, "BOUT.stop");
      // Now add data directory to start of name to ensure we look in a run specific location
      std::string data_dir;
      Options::getRoot()->get("datadir", data_dir, string(DEFAULT_DIR));
      stopCheckName = data_dir + "/" + stopCheckName;
    }

    /// Record the starting time
    mpi_start_time = MPI_Wtime() - wtime;

    first_time = false;

    /// Print the column header for timing info
    if (!output_split) {
      output_progress.write("Sim Time  |  RHS evals  | Wall Time |  Calc    Inv   Comm    I/O   "
                            "SOLVER\n\n");
    } else {
      output_progress.write("Sim Time  |  RHS_e evals  | RHS_I evals  | Wall Time |  Calc    Inv  "
                            " Comm    I/O   SOLVER\n\n");
    }
  }

  if (!output_split) {
    output_progress.write("%.3e      %5d       %.2e   %5.1f  %5.1f  %5.1f  %5.1f  %5.1f\n", 
               simtime, ncalls, wtime,
               100.0*(wtime_rhs - wtime_comms - wtime_invert)/wtime,
               100.*wtime_invert/wtime,  // Inversions
               100.0*wtime_comms/wtime,  // Communications
               100.* wtime_io / wtime,      // I/O
               100.*(wtime - wtime_io - wtime_rhs)/wtime); // Everything else

  } else {
    output_progress.write("%.3e      %5d            %5d       %.2e   %5.1f  %5.1f  %5.1f  %5.1f  %5.1f\n",
               simtime, ncalls_e, ncalls_i, wtime,
               100.0*(wtime_rhs - wtime_comms - wtime_invert)/wtime,
               100.*wtime_invert/wtime,  // Inversions
               100.0*wtime_comms/wtime,  // Communications
               100.* wtime_io / wtime,      // I/O
               100.*(wtime - wtime_io - wtime_rhs)/wtime); // Everything else
  }

  // This bit only to screen, not log file

  BoutReal t_elapsed = MPI_Wtime() - mpi_start_time;
  output_progress.print("%c  Step %d of %d. Elapsed %s", get_spin(), iteration+1, NOUT, (time_to_hms(t_elapsed)).c_str());
  output_progress.print(" ETA %s", (time_to_hms(wtime * static_cast<BoutReal>(NOUT - iteration - 1))).c_str());

  if (wall_limit > 0.0) {
    // Check if enough time left

    BoutReal t_remain = mpi_start_time + wall_limit - MPI_Wtime();
    if (t_remain < wtime) {
      // Less than 1 time-step left
      output_warn.write("Only %e seconds left. Quitting\n", t_remain);

      return 1; // Return an error code to quit
    } else {
      output_progress.print(" Wall %s", (time_to_hms(t_remain)).c_str());
    }
  }

  // Check if the user has created the stop file and if so trigger an exit
  if (stopCheck) {
    std::ifstream f(stopCheckName);
    if (f.good()) {
      output << "\n" << "File " << stopCheckName
             << " exists -- triggering exit." << endl;
      return 1;
    }
  }

  return 0;
}

/**************************************************************************
 * Global error handling
 **************************************************************************/

/// Signal handler - handles all signals
void bout_signal_handler(int sig) {
  /// Set signal handler back to default to prevent possible infinite loop
  signal(SIGSEGV, SIG_DFL);
  // print number of process to stderr, so the user knows which log to check
  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  fprintf(stderr,"\nSighandler called on process %d with sig %d\n"
          ,world_rank,sig);

  switch (sig){
  case SIGSEGV:
    throw BoutException("\n****** SEGMENTATION FAULT CAUGHT ******\n\n");
    break;
  case SIGFPE:
    throw BoutException("\n****** Floating Point Exception "
                        "(FPE) caught ******\n\n");
    break;
  case SIGINT:
    throw BoutException("\n****** SigInt caught ******\n\n");
    break;
  case SIGKILL:
    throw BoutException("\n****** SigKill caught ******\n\n");
    break;
  case SIGUSR1:
    user_requested_exit=true;
    break;
  default:
    throw BoutException("\n****** Signal %d  caught ******\n\n",sig);
    break;
  }
}

/**************************************************************************
 * Utilities
 **************************************************************************/

/// Write a time in h:mm:ss.s format
const string time_to_hms(BoutReal t) {
  int h, m;

  h = static_cast<int>(t / 3600);
  t -= 3600. * static_cast<BoutReal>(h);
  m = static_cast<int>(t / 60);
  t -= 60 * static_cast<BoutReal>(m);

  char buffer[256];
  sprintf(buffer, "%d:%02d:%04.1f", h, m, t);

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
