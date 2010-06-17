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
static char DEFAULT_GRID[] = "data/bout.grd.pdb";
#ifdef PETSC
static char help[] = "BOUT++: Uses finite difference methods to solve plasma fluid problems in curvilinear coordinates";
#endif

// MD5 Checksum passed at compile-time
#define CHECKSUM1_(x) #x
#define CHECKSUM_(x) CHECKSUM1_(x)
#define CHECKSUM CHECKSUM_(MD5SUM)

#define GLOBALORIGIN

#include "bout.h"
#include "datafile.h"
#include "grid.h"
#include "solver.h"
#include "field2d.h"
#include "field3d.h"
#include "vector2d.h"
#include "vector3d.h"
#include "initialprofiles.h"
#include "derivs.h"
#include "utils.h"
#include "invert_laplace.h"
#include "interpolation.h"
#include "boutmesh.h"

#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include <string>
using std::string;

#ifdef SIGHANDLE
#include <signal.h>
void bout_signal_handler(int sig);  // Handles segmentation faults
#endif

bool appending = false;
char dumpname[512];

real simtime;
int iteration;

void setup_files();                 // Adds basic quantities to the output files 
const string time_to_hms(real t);   // Converts to h:mm:ss.s format
char get_spin();                    // Produces a spinning bar

int bout_monitor(real t, int iter, int NOUT); // Function called by the solver each timestep

/*!************************************************************************
 * Main function
 **************************************************************************/

int main(int argc, char **argv)
{
  int i, NOUT;
  real TIMESTEP;
  char *grid_name;
  time_t start_time, end_time;
  bool dump_float; // Output dump files as floats

  char *grid_ext, *dump_ext; ///< Extensions for restart and dump files
  
  const char *data_dir; ///< Directory for data input/output

#ifdef CHECK
  int msg_point; ///< Used to return the message stack to a fixed point
#endif

#ifdef SIGHANDLE
  /// Set a signal handler for segmentation faults
  signal(SIGSEGV, bout_signal_handler);
#endif

  // Set default data directory
  data_dir = DEFAULT_DIR;

  /// Check command-line arguments
  for(i=1;i<argc;i++) {
    if(strncasecmp(argv[i], "re", 2) == 0) {
      restarting = true;
    }
    if(strncasecmp(argv[i], "ap", 2) == 0) {
      appending = true;
    }
    if(strncasecmp(argv[i], "no", 2) == 0) {
      // No output. Used for scaling studies
      Datafile::enabled = false;
    }
    if(strncasecmp(argv[i], "-d", 2) == 0) {
      // Set data directory
      if(i+1 >= argc) {
		output.write("Useage is %s -d <data directory>\n");
		return 1;
      }
      i++;
      data_dir = argv[i];
    }
  }
  
  /// Start MPI
#ifdef PETSC
  PetscInitialize(&argc,&argv,"../petscopt",help);
#else
  MPI_Init(&argc,&argv);
#endif
  MPI_Comm_size(MPI_COMM_WORLD, &NPES);
  MPI_Comm_rank(MPI_COMM_WORLD, &MYPE);

  /// Set up the output
  if(MYPE == 0) {
    output.enable(); // Enable writing to stdout
  }else {
    output.disable(); // No writing to stdout
  }
  /// Open an output file to echo everything to
  /// On processor 0 anything written to output will go to stdout and the file
  output.open("%s/BOUT.log.%d", data_dir, MYPE);

  /// Print intro
  output.write("\nBOUT++ version %.2f\n", BOUT_VERSION);
#ifdef MD5SUM
  output.write("MD5 checksum: %s\n", CHECKSUM);
#endif
  output.write("Code compiled on %s at %s\n", __DATE__, __TIME__);
  output.write("B.Dudson (University of York), M.Umansky (LLNL) 2007\n");
  output.write("Based on BOUT by Xueqiao Xu, 1999\n\n");

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

#ifdef METRIC3D
  output.write("\tRUNNING IN 3D-METRIC MODE\n");
#endif

  start_time = time((time_t*) NULL);
  output.write("\nRun started at  : %s\n", ctime(&start_time));
  output.write("Processor number: %d of %d\n\n", MYPE, NPES);

  /// Load settings file
  if(options.read("%s/BOUT.inp", data_dir)) {
    output.write("\tFailed to read settings file. Aborting\n");
    return(1);
  }

  /////////////////////////////////////////////
  /// Get some settings
  
  /// GET GLOBAL OPTIONS
  options.setSection(NULL);

  options.get("NOUT", NOUT, 1);
  options.get("TIMESTEP", TIMESTEP, 1.0);

  if((grid_name = options.getString("grid")) == (char*) NULL)
    grid_name = DEFAULT_GRID;
  
  OPTION(dump_float,   true);
  OPTION(IncIntShear,  false);
  OPTION(ShiftOrder,   0);
  OPTION(non_uniform,  false);

  options.get("BoundaryOnCell", BoundaryOnCell, false); // Determine location of boundary
  options.get("StaggerGrids",   StaggerGrids,   false); // Stagger grids
  
  /// Get file extensions
  if((dump_ext = options.getString("dump_format")) == NULL) {
    // Set default extension
    dump_ext = DEFAULT_FILE_EXT;
  }
  
  /// Setup derivative methods
  if(derivs_init()) {
    output.write("Failed to initialise derivative methods. Aborting\n");
    return(1);
  }
  
  ////////////////////////////////////////////

  /// Create the mesh
  mesh = new BoutMesh();
  
  output.write("Setting grid format\n");
  /// Load the grid
  if((grid_ext = options.getString("grid_format")) == NULL) {
    // Guess format based on grid filename
    mesh->add_source(new GridFile(data_format(grid_name), grid_name));
  }else {
    // User-specified format
    mesh->add_source(new GridFile(data_format(grid_ext), grid_name));
  }
  if(mesh->load()) {
    output << "Failed to read grid. Aborting\n";
    return 1;
  }

  /// Set the file names
  sprintf(dumpname, "%s/BOUT.dmp.%d.%s", data_dir, MYPE, dump_ext);

  // Set file formats
  output.write("Setting file formats\n");
  dump.setFormat(data_format(dumpname));

  if(dump_float)
    dump.setLowPrecision(); // Down-convert to floats

  /// Add book-keeping variables to the output files
  setup_files();

  /// initialise Laplacian inversion code
  invert_init();

  output.write("Initialising physics module\n");
  /// Initialise physics module
#ifdef CHECK
  msg_point = msg_stack.push("Initialising physics module");
#endif

  if(physics_init()) {
    output.write("Failed to initialise physics. Aborting\n");
    return(1);
  }

#ifdef CHECK
  // Can't trust that the user won't leave messages on the stack
  msg_stack.pop(msg_point);
#endif

  /// Initialise the solver
  solver.setRestartDir(data_dir);
  if(solver.init(physics_run, argc, argv, restarting, NOUT, TIMESTEP)) {
    output.write("Failed to initialise solver. Aborting\n");
    return(1);
  }
  
  if(!restarting) {
    /// Write initial state as time-point 0
    
    // Run RHS once to ensure all variables set
    if(physics_run(0.0)) {
      output.write("Physics RHS call failed\n");
      return(1);
    }

    if(appending) {
      dump.append(dumpname);
    }else {
      dump.write(dumpname);
      appending = true;
    }
  }

  //////// MAIN LOOP /////////

  output.write("Running simulation\n\n");

  /// Run the solver
  solver.run(bout_monitor);

  // close MPI
#ifdef PETSC
  PetscFinalize();
#else
  MPI_Finalize();
#endif

  end_time = time((time_t*) NULL);
  output.write("\nRun finished at  : %s\n", ctime(&end_time));
  output.write("Run time : ");

  int dt;
  dt = end_time - start_time;
  i = (int) (dt / (60.*60.));
  if(i > 0) {
    output.write("%d h ", i);
    dt -= i*60*60;
  }
  i = (int) (dt / 60.);
  if(i > 0) {
    output.write("%d m ", i);
    dt -= i*60;
  }
  output.write("%d s\n", dt);

  return(0);
}

/*!*************************************************************************
 * SOLUTION MONITOR FUNCTION
 *
 * Called each timestep by the solver
 **************************************************************************/

int bout_monitor(real t, int iter, int NOUT)
{
  // Data used for timing
  static bool first_time = true;
  static real wtime = 0.0;       ///< Wall-time since last output
  static real wall_limit, mpi_start_time; // Keep track of remaining wall time

#ifdef CHECK
  int msg_point = msg_stack.push("bout_monitor(%e, %d, %d)", t, iter, NOUT);
#endif

  // Set the global variables. This is done because they need to be
  // written to the output file before the first step (initial condition)
  simtime = t;
  iteration = iter;

  /// Write (append) dump file
  
  if(appending) {
    dump.append(dumpname);
  }else {
    dump.write(dumpname);
    appending = true;
  }
  
  /// Collect timing information
  int ncalls = solver.rhs_ncalls;
  real wtime_rhs   = solver.rhs_wtime;
  //real wtime_invert = 0.0; // wtime_invert is a global
  real wtime_comms = Communicator::wtime;  // Time spent communicating (part of RHS)
  real wtime_io    = Datafile::wtime;      // Time spend on I/O

  output.print("\r"); // Only goes to screen
  
  if(first_time) {
    /// First time the monitor has been called
    
    /// Get some options
    options.setSection(NULL);
    OPTION(wall_limit, -1.0); // Wall time limit. By default, no limit
    wall_limit *= 60.0*60.0;  // Convert from hours to seconds
    
    /// Record the starting time
    mpi_start_time = MPI_Wtime(); // NB: Miss time for first step (can be big!)

    first_time = false;

    /// Print the column header for timing info
    output.write("Sim Time  |  RHS evals  | Wall Time |  Calc    Inv   Comm    I/O   SOLVER\n\n");
    
    /// Don't know wall time, so don't print timings
    output.write("%.3e      %5d        -     -    -    -    -    - \n", 
	       simtime, ncalls); // Everything else
  }else {
    wtime = MPI_Wtime() - wtime;
    
    output.write("%.3e      %5d       %.2e   %5.1f  %5.1f  %5.1f  %5.1f  %5.1f\n", 
		 simtime, ncalls, wtime,
		 100.0*(wtime_rhs - wtime_comms - wtime_invert)/wtime,
		 100.*wtime_invert/wtime,  // Inversions
		 100.0*wtime_comms/wtime,  // Communications
		 100.*wtime_io/wtime,      // I/O
		 100.*(wtime - wtime_io - wtime_rhs)/wtime); // Everything else
  }
  
  // This bit only to screen, not log file
  
  real t_elapsed = MPI_Wtime() - mpi_start_time;
  output.print("%c  Step %d of %d. Elapsed %s", get_spin(), iteration+1, NOUT, (time_to_hms(t_elapsed)).c_str());
  output.print(" ETA %s", (time_to_hms(wtime * ((real) (NOUT - iteration - 1)))).c_str());
  
  if(wall_limit > 0.0) {
    // Check if enough time left
    
    real t_remain = mpi_start_time + wall_limit - MPI_Wtime();
    if(t_remain < wtime) {
      // Less than 1 time-step left
      output.write("Only %e seconds left. Quitting\n", t_remain);
      
#ifdef CHECK
      msg_stack.pop(msg_point);
#endif
      return 1; // Return an error code to quit
    }else {
      output.print(" Wall %s", (time_to_hms(t_remain)).c_str());
    } 
  }

  /// Reset clocks for next timestep
  
  Communicator::wtime = 0.0; // Reset communicator clock
  Datafile::wtime = 0.0;
  wtime_invert = 0.0;
  wtime = MPI_Wtime();
  
#ifdef CHECK
  msg_stack.pop(msg_point);
#endif
  
  return 0;
}

/*!*************************************************************************
 * SETUP RESTART & DUMP FILES
 **************************************************************************/

void setup_files()
{
  // Setup dump file

  // This is a temporary hack to get around datafile's limitations (fix soon)
  static real version = BOUT_VERSION;
  dump.add(version, "BOUT_VERSION", 0);
  dump.add(MXSUB, "MXSUB", 0);
  dump.add(MYSUB, "MYSUB", 0);
  dump.add(MXG,   "MXG",   0);
  dump.add(MYG,   "MYG",   0);
  dump.add(MZ,    "MZ",    0);
  dump.add(mesh->NXPE,  "NXPE",  0);
  dump.add(NYPE,  "NYPE",  0);
  dump.add(ZMAX,  "ZMAX",  0);
  dump.add(ZMIN,  "ZMIN",  0);

  dump.add(simtime, "t_array", 1); // Appends the time of dumps into an array
  dump.add(iteration, "iteration", 0);
}

/*!************************************************************************
 * Add variables to be solved
 **************************************************************************/

// NOTE: Here bout_solve is for backwards-compatibility. Eventually will be removed

void bout_solve(Field2D &var, Field2D &F_var, const char *name)
{
  // Add to solver
  solver.add(var, F_var, name);
}

void bout_solve(Field3D &var, Field3D &F_var, const char *name)
{
  solver.add(var, F_var, name);
}

void bout_solve(Vector2D &var, Vector2D &F_var, const char *name)
{
  solver.add(var, F_var, name);
}

void bout_solve(Vector3D &var, Vector3D &F_var, const char *name)
{
  solver.add(var, F_var, name);
}

/*!************************************************************************
 * Add constraints
 **************************************************************************/

bool bout_constrain(Field3D &var, Field3D &F_var, const char *name)
{
  // Add to solver
  solver.constraint(var, F_var, name);

  return true;
}

/**************************************************************************
 * Global error handling
 **************************************************************************/

/// Print an error message and exit
void bout_error()
{
  bout_error(NULL);
}

void bout_error(const char *str)
{
  output.write("****** ERROR CAUGHT ******\n");

  if(str != NULL)
    output.write(str);

#ifdef CHECK
  msg_stack.dump();
#endif

  MPI_Abort(MPI_COMM_WORLD, 1);

  exit(1);
}

#ifdef SIGHANDLE
/// Signal handler - catch segfaults
void bout_signal_handler(int sig)
{
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
const string time_to_hms(real t)
{
  int h, m;
  
  h = (int) (t / 3600); t -= 3600.*((real) h);
  m = (int) (t / 60);   t -= 60 * ((real) m);
  
  char buffer[256];
  sprintf(buffer,"%d:%02d:%04.1f", h, m, t);
  
  return string(buffer);
}

/// Produce a spinning bar character
char get_spin()
{
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
