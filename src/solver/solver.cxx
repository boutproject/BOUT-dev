/**************************************************************************
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

#include <boutcomm.hxx>
#include <bout/solver.hxx>
#include <string.h>
#include <time.h>

#include <initialprofiles.hxx>
#include <interpolation.hxx>
#include <boutexception.hxx>

#include <field_factory.hxx>

#include "solverfactory.hxx"

#include <bout/sys/timer.hxx>
#include <msg_stack.hxx>
#include <output.hxx>
#include <bout/assert.hxx>

#include <bout/array.hxx>

// Static member variables

int* Solver::pargc = 0;
char*** Solver::pargv = 0;

/**************************************************************************
 * Constructor
 **************************************************************************/

Solver::Solver(Options *opts) : options(opts), model(0), prefunc(0) {
  if(options == NULL)
    options = Options::getRoot()->getSection("solver");

  // Set flags to defaults
  has_constraints = false;
  initialised = false;
  canReset = false;

  // Zero timing
  rhs_ncalls = 0;
  rhs_ncalls_e = 0;
  rhs_ncalls_i = 0;
  // Restart directory
  if(options->isSet("restartdir")) {
    // Solver-specific restart directory
    options->get("restartdir", restartdir, "data");
  }else {
    // Use the root data directory
    Options::getRoot()->get("datadir", restartdir, "data");
  }
  
  // Restart option
  options->get("enablerestart", enablerestart, true);
  if(enablerestart) {
    Options::getRoot()->get("restart", restarting, false);
  }else
    restarting = false;

  // Set up restart options
  restart = Datafile(Options::getRoot()->getSection("restart"));
  
  // Split operator
  split_operator = false;
  max_dt = -1.0;

  // Output monitor
  options->get("monitor_timestep", monitor_timestep, false);
  
  // Method of Manufactured Solutions (MMS)
  options->get("mms", mms, false);
  options->get("mms_initialise", mms_initialise, mms);
}

/**************************************************************************
 * Add physics models
 **************************************************************************/

void Solver::setModel(PhysicsModel *m) {
  if(model)
    throw BoutException("Solver can only evolve one model");
  
  if(initialised)
    throw BoutException("Solver already initialised");
  
  if(m->initialise(this, restarting))
    throw BoutException("Couldn't initialise physics model");
  
  // Check if the model is split operator
  split_operator = m->splitOperator();
  
  model = m;
}

/**************************************************************************
 * Add fields
 **************************************************************************/

void Solver::add(Field2D &v, const char* name) {
  int msg_point = msg_stack.push("Adding 2D field: Solver::add(%s)", name);
  
  if(varAdded(string(name)))
    throw BoutException("Variable '%s' already added to Solver", name);

  if(initialised)
    throw BoutException("Error: Cannot add to solver after initialisation\n");
  
  // Set boundary conditions
  v.setBoundary(name);
  ddt(v).copyBoundary(v); // Set boundary to be the same as v

  VarStr<Field2D> d;
  
  d.constraint = false;
  d.var = &v;
  d.F_var = &ddt(v);
  d.name = string(name);

#ifdef TRACK
  v.name = name;
#endif

  /// Generate initial perturbations.
  /// NOTE: This could be done in init, but this would prevent the user
  ///       from modifying the initial perturbation (e.g. to prevent unphysical situations)
  ///       before it's loaded into the solver. If restarting, this perturbation
  ///       will be over-written anyway
  if(mms_initialise) {
    // Load solution at t = 0
    
    FieldFactory *fact = FieldFactory::get();
    
    v = fact->create2D("solution", Options::getRoot()->getSection(name), mesh);
  }else {
    initial_profile(name, v);
  }
  
  if(mms) {
    // Allocate storage for error variable
    d.MMS_err = new Field2D(0.0);
  }
  
  // Check if the boundary regions should be evolved
  // First get option from section "All"
  // then use that as default for specific section
  Options::getRoot()->getSection("all")->get("evolve_bndry", d.evolve_bndry, false);
  Options::getRoot()->getSection(name)->get("evolve_bndry", d.evolve_bndry, d.evolve_bndry);

  v.applyBoundary(true);

  f2d.push_back(d);

  msg_stack.pop(msg_point);
}

void Solver::add(Field3D &v, const char* name) {

#ifdef CHECK
  int msg_point = msg_stack.push("Adding 3D field: Solver::add(%s)", name);
  
  if(varAdded(string(name)))
    throw BoutException("Variable '%s' already added to Solver", name);
#endif

  if(initialised)
    throw BoutException("Error: Cannot add to solver after initialisation\n");

  // Set boundary conditions
  v.setBoundary(name);
  ddt(v).copyBoundary(v); // Set boundary to be the same as v

  if(mesh->StaggerGrids && (v.getLocation() != CELL_CENTRE)) {
    output.write("\tVariable %s shifted to %s\n", name, strLocation(v.getLocation()));
    ddt(v).setLocation(v.getLocation()); // Make sure both at the same location
  }

  VarStr<Field3D> d;
  
  d.constraint = false;
  d.var = &v;
  d.F_var = &ddt(v);
  d.location = v.getLocation();
  d.name = string(name);
  
#ifdef TRACK
  v.name = name;
#endif

  if(mms_initialise) {
    // Load solution at t = 0
    FieldFactory *fact = FieldFactory::get();
    
    v = fact->create3D("solution", Options::getRoot()->getSection(name), mesh, v.getLocation());
    
  }else {
    initial_profile(name, v);
  }
  
  if(mms) {
    d.MMS_err = new Field3D();
    (*d.MMS_err) = 0.0;
  }
  
  // Check if the boundary regions should be evolved
  // First get option from section "All"
  // then use that as default for specific section
  Options::getRoot()->getSection("all")->get("evolve_bndry", d.evolve_bndry, false);
  Options::getRoot()->getSection(name)->get("evolve_bndry", d.evolve_bndry, d.evolve_bndry);

  v.applyBoundary(true); // Make sure initial profile obeys boundary conditions
  v.setLocation(d.location); // Restore location if changed
  
  f3d.push_back(d);
              
#ifdef CHECK
  msg_stack.pop(msg_point);
#endif
}

void Solver::add(Vector2D &v, const char* name) {

  int msg_point = msg_stack.push("Adding 2D vector: Solver::add(%s)", name);
  
  if(varAdded(string(name)))
    throw BoutException("Variable '%s' already added to Solver", name);

  if(initialised)
    throw BoutException("Error: Cannot add to solver after initialisation\n");

  // Set boundary conditions
  v.setBoundary(name);
  ddt(v).copyBoundary(v); // Set boundary to be the same as v
  
  VarStr<Vector2D> d;
  
  d.constraint = false;
  d.var = &v;
  d.F_var = &ddt(v);
  d.covariant = v.covariant;
  d.name = string(name);

  v2d.push_back(d);

  /// NOTE: No initial_profile call, because this will be done for each
  ///       component individually.
  
  /// Add suffix, depending on co- /contravariance
  if(v.covariant) {
    add(v.x, (d.name+"_x").c_str());
    add(v.y, (d.name+"_y").c_str());
    add(v.z, (d.name+"_z").c_str());
  }else {
    add(v.x, (d.name+"x").c_str());
    add(v.y, (d.name+"y").c_str());
    add(v.z, (d.name+"z").c_str());
  }
  
  /// Make sure initial profile obeys boundary conditions
  v.applyBoundary(true);

  msg_stack.pop(msg_point);
}

void Solver::add(Vector3D &v, const char* name) {

  int msg_point = msg_stack.push("Adding 3D vector: Solver::add(%s)", name);
  
  if(varAdded(string(name)))
    throw BoutException("Variable '%s' already added to Solver", name);

  if(initialised)
    throw BoutException("Error: Cannot add to solver after initialisation\n");
  
  // Set boundary conditions
  v.setBoundary(name);
  ddt(v).copyBoundary(v); // Set boundary to be the same as v

  VarStr<Vector3D> d;
  
  d.constraint = false;
  d.var = &v;
  d.F_var = &ddt(v);
  d.covariant = v.covariant;
  d.name = string(name);
  
  v3d.push_back(d);

  // Add suffix, depending on co- /contravariance
  if(v.covariant) {
    add(v.x, (d.name+"_x").c_str());
    add(v.y, (d.name+"_y").c_str());
    add(v.z, (d.name+"_z").c_str());
  }else {
    add(v.x, (d.name+"x").c_str());
    add(v.y, (d.name+"y").c_str());
    add(v.z, (d.name+"z").c_str());
  }

  v.applyBoundary(true);

  msg_stack.pop(msg_point);
}

/**************************************************************************
 * Constraints
 **************************************************************************/

void Solver::constraint(Field2D &v, Field2D &C_v, const char* name) {

#ifdef CHECK
  int msg_point = msg_stack.push("Constrain 2D scalar: Solver::constraint(%s)", name);
  
  if(varAdded(string(name)))
    throw BoutException("Variable '%s' already added to Solver", name);
#endif

  if(!has_constraints)
    throw BoutException("ERROR: This solver doesn't support constraints\n");

  if(initialised)
    throw BoutException("Error: Cannot add constraints to solver after initialisation\n");

  if(name == NULL)
    throw BoutException("WARNING: Constraint requested for variable with NULL name\n");
  
  VarStr<Field2D> d;
  
  d.constraint = true;
  d.var = &v;
  d.F_var = &C_v;
  d.name = string(name);

  f2d.push_back(d);

#ifdef CHECK
  msg_stack.pop(msg_point);
#endif
}

void Solver::constraint(Field3D &v, Field3D &C_v, const char* name) {

#ifdef CHECK
  int msg_point = msg_stack.push("Constrain 3D scalar: Solver::constraint(%s)", name);

  if(varAdded(string(name)))
    throw BoutException("Variable '%s' already added to Solver", name);
#endif

  if(!has_constraints)
    throw BoutException("ERROR: This solver doesn't support constraints\n");

  if(initialised)
    throw BoutException("Error: Cannot add constraints to solver after initialisation\n");

  if(name == NULL)
    throw BoutException("WARNING: Constraint requested for variable with NULL name\n");

  VarStr<Field3D> d;
  
  d.constraint = true;
  d.var = &v;
  d.F_var = &C_v;
  d.location = v.getLocation();
  d.name = string(name);
  
  f3d.push_back(d);

#ifdef CHECK
  msg_stack.pop(msg_point);
#endif
}

void Solver::constraint(Vector2D &v, Vector2D &C_v, const char* name) {

#ifdef CHECK
  int msg_point = msg_stack.push("Constrain 2D vector: Solver::constraint(%s)", name);
  
  if(varAdded(string(name)))
    throw BoutException("Variable '%s' already added to Solver", name);
#endif

  if(!has_constraints)
    throw BoutException("ERROR: This solver doesn't support constraints\n");

  if(initialised)
    throw BoutException("Error: Cannot add constraints to solver after initialisation\n");

  if(name == NULL)
    throw BoutException("WARNING: Constraint requested for variable with NULL name\n");
    
  VarStr<Vector2D> d;
  
  d.constraint = true;
  d.var = &v;
  d.F_var = &C_v;
  d.covariant = v.covariant;
  d.name = string(name);
  
  v2d.push_back(d);

  // Add suffix, depending on co- /contravariance
  if(v.covariant) {
    constraint(v.x, C_v.x, (d.name+"_x").c_str());
    constraint(v.y, C_v.y, (d.name+"_x").c_str());
    constraint(v.z, C_v.z, (d.name+"_x").c_str());
  }else {
    constraint(v.x, C_v.x, (d.name+"x").c_str());
    constraint(v.y, C_v.y, (d.name+"x").c_str());
    constraint(v.z, C_v.z, (d.name+"x").c_str());
  }

#ifdef CHECK
  msg_stack.pop(msg_point);
#endif
}

void Solver::constraint(Vector3D &v, Vector3D &C_v, const char* name) {

#ifdef CHECK
  int msg_point = msg_stack.push("Constrain 3D vector: Solver::constraint(%s)", name);
  
  if(varAdded(string(name)))
    throw BoutException("Variable '%s' already added to Solver", name);
#endif

  if(!has_constraints)
    throw BoutException("ERROR: This solver doesn't support constraints\n");

  if(initialised)
    throw BoutException("Error: Cannot add constraints to solver after initialisation\n");

  if(name == NULL)
    throw BoutException("WARNING: Constraint requested for variable with NULL name\n");

  VarStr<Vector3D> d;
  
  d.constraint = true;
  d.var = &v;
  d.F_var = &C_v;
  d.covariant = v.covariant;
  d.name = string(name);
  
  v3d.push_back(d);

  // Add suffix, depending on co- /contravariance
  if(v.covariant) {
    constraint(v.x, C_v.x, (d.name+"_x").c_str());
    constraint(v.y, C_v.y, (d.name+"_x").c_str());
    constraint(v.z, C_v.z, (d.name+"_x").c_str());
  }else {
    constraint(v.x, C_v.x, (d.name+"x").c_str());
    constraint(v.y, C_v.y, (d.name+"x").c_str());
    constraint(v.z, C_v.z, (d.name+"x").c_str());
  }

#ifdef CHECK
  msg_stack.pop(msg_point);
#endif
}

/**************************************************************************
 * Solver main loop: Initialise, run, and finish
 **************************************************************************/

int Solver::solve(int NOUT, BoutReal TIMESTEP) {
  
  dump_on_restart = false;
  bool append = false;
  if(NOUT < 0) {
    /// Get options
    
    Options *globaloptions = Options::getRoot(); // Default from global options
    OPTION(globaloptions, NOUT, 1);
    OPTION(globaloptions, TIMESTEP, 1.0);
    OPTION(globaloptions, append, false);
    OPTION(globaloptions, dump_on_restart, !restarting || !append);
    
    // Check specific solver options, which override global options
    OPTION(options, NOUT, NOUT);
    options->get("output_step", TIMESTEP, TIMESTEP);
  }
  
  output.write("Solver running for %d outputs with output timestep of %e\n", NOUT, TIMESTEP);
  
  // Initialise
  if(init(restarting, NOUT, TIMESTEP)) {
    throw BoutException("Failed to initialise solver-> Aborting\n");
  }
  
  /// Run the solver
  output.write("Running simulation\n\n");
  
  time_t start_time = time((time_t*) NULL);
  output.write("\nRun started at  : %s\n", ctime(&start_time));
  
  Timer timer("run"); // Start timer
  
  if ( dump_on_restart ) {
    /// Write initial state as time-point 0
    
    // Run RHS once to ensure all variables set
    if (run_rhs(simtime)) {
      throw BoutException("Physics RHS call failed\n");
    }
    
    // Call monitors so initial values are written to output dump files
    if (call_monitors(simtime, 0, NOUT)){
      throw BoutException("Initial monitor call failed!");
    }
  }
  
  int status;
  try {
    status = run();

    time_t end_time = time((time_t*) NULL);
    output.write("\nRun finished at  : %s\n", ctime(&end_time));
    output.write("Run time : ");

    int dt = end_time - start_time;
    int i = (int) (dt / (60.*60.));
    if (i > 0) {
      output.write("%d h ", i);
      dt -= i*60*60;
    }
    i = (int) (dt / 60.);
    if (i > 0) {
      output.write("%d m ", i);
      dt -= i*60;
    }
    output.write("%d s\n", dt);
  }catch(BoutException &e) {
    output << "Error encountered in solver run\n";
    output << e.what() << endl;
    
    if(enablerestart) {
      // Write restart to a different file
      restart.write("%s/BOUT.failed.%s", restartdir.c_str(), restartext.c_str());
    }
    
    throw e;
  }

  return status;
}


/**************************************************************************
 * Initialisation
 **************************************************************************/

int Solver::init(bool restarting, int nout, BoutReal tstep) {
  
  TRACE("Solver::init()");

  if(initialised)
    throw BoutException("ERROR: Solver is already initialised\n");

  output.write("Initialising solver\n");

  MPI_Comm_size(BoutComm::get(), &NPES);
  MPI_Comm_rank(BoutComm::get(), &MYPE);
  
  if(enablerestart) {
    // Set up restart file
    
    options->get("archive", archive_restart, -1);

    if(archive_restart > 0) {
      output.write("Archiving restart files every %d iterations\n",
                   archive_restart);
    }
    
    /// Get restart file extension
    string dump_ext, restart_ext;
    
    options->get("dump_format", dump_ext, "nc");
    options->get("restart_format", restart_ext, dump_ext);
    restartext = string(restart_ext);
  
    /// Add basic variables to the restart file
    restart.add(simtime,  "tt",    0);
    restart.add(iteration, "hist_hi", 0);
    
    restart.add(NPES, "NPES", 0);
    restart.add(mesh->NXPE, "NXPE", 0);

    /// Add variables to the restart and dump files.
    /// NOTE: Since vector components are already in the field arrays,
    ///       only loop over scalars, not vectors
    for(const auto& f : f2d) {
      // Add to restart file (not appending)
      restart.add(*(f.var), f.name.c_str(), 0);
      
      /// NOTE: Initial perturbations have already been set in add()
      
      /// Make sure boundary condition is satisfied
      //f.var->applyBoundary();
      /// NOTE: boundary conditions on the initial profiles have also been set in add()
    }  
    for(const auto& f : f3d) {
      // Add to restart file (not appending)
      restart.add(*(f.var), f.name.c_str(), 0);
      
      /// Make sure boundary condition is satisfied
      //f.var->applyBoundary();
      /// NOTE: boundary conditions on the initial profiles have also been set in add()
    }
  }
  
  if(restarting) {
    /// Load state from the restart file
    
    // Copy processor numbers for comparison after. Very useful for checking
    // that the restart file is for the correct number of processors etc.
    int tmp_NP = NPES;
    int tmp_NX = mesh->NXPE;
    
    TRACE("Loading restart file");
    
    /// Load restart file
    if(!restart.openr("%s/BOUT.restart.%s", restartdir.c_str(), restartext.c_str()))
      throw BoutException("Error: Could not open restart file\n");
    if(!restart.read())
      throw BoutException("Error: Could not read restart file\n");
    restart.close();

    if(NPES == 0) {
      // Old restart file
      output.write("WARNING: Cannot verify processor numbers\n");
      NPES = tmp_NP;
      mesh->NXPE = tmp_NX;
    }else {
      // Check the processor numbers match
      if(NPES != tmp_NP) {
	output.write("ERROR: Number of processors (%d) doesn't match restart file number (%d)\n",
		     tmp_NP, NPES);
	return(1);
      }
      if(mesh->NXPE != tmp_NX) {
	output.write("ERROR: Number of X processors (%d) doesn't match restart file number (%d)\n",
		     tmp_NX, mesh->NXPE);
	return(1);
      }

      output.write("Restarting at iteration %d, simulation time %e\n", iteration, simtime);
    }
    
  }else {
    // Not restarting
    simtime = 0.0; iteration = 0;
  }
  
  if(enablerestart) {
    /// Open the restart file for writing
    if(!restart.openw("%s/BOUT.restart.%s", restartdir.c_str(), restartext.c_str()))
      throw BoutException("Error: Could not open restart file for writing\n");
  }
  
  /// Mark as initialised. No more variables can be added
  initialised = true;

  return 0;
}

void Solver::outputVars(Datafile &outputfile) {
  // Add 2D and 3D evolving fields to output file
  for(const auto& f : f2d) {
    // Add to dump file (appending)
    outputfile.add(*(f.var), f.name.c_str(), 1);
  }  
  for(const auto& f : f3d) {
    // Add to dump file (appending)
    outputfile.add(*(f.var), f.name.c_str(), 1);
    
    if(mms) {
      // Add an error variable
      dump.add(*(f.MMS_err), (string("E_")+f.name).c_str(), 1);
    }
  }
}

/////////////////////////////////////////////////////

void Solver::addMonitor(MonitorFunc f, MonitorPosition pos) {
  if(pos == Solver::FRONT) {
    monitors.push_front(f);
  }else
    monitors.push_back(f);
}

void Solver::removeMonitor(MonitorFunc f) {
  monitors.remove(f);
}

int Solver::call_monitors(BoutReal simtime, int iter, int NOUT) {
  if(mms) {
    // Calculate MMS errors
    calculate_mms_error(simtime);
  }
  
  if( enablerestart ) {
    /// Write the restart file
    restart.write();
    
    if((archive_restart > 0) && (iteration % archive_restart == 0)) {
      restart.write("%s/BOUT.restart_%04d.%s", restartdir.c_str(), iteration, restartext.c_str());
    }
  }
  
  try {
    // Call physics model monitor
    if(model) {
      if(model->runOutputMonitor(simtime, iter, NOUT))
        throw BoutException("Monitor signalled to quit");
    }
    
    // Call C function monitors
    for(const auto& monitor : monitors) {
      // Call each monitor one by one
      int ret = monitor(this, simtime,iter, NOUT);
      if(ret)
        throw BoutException("Monitor signalled to quit");
    }
  } catch (BoutException &e) {
    // User signalled to quit
    if( enablerestart ) {
      // Write restart to a different file
      restart.write("%s/BOUT.final.%s", restartdir.c_str(), restartext.c_str());
    }
    
    output.write("Monitor signalled to quit. Returning\n");
    return 1;
  }
  
  // Reset iteration and wall-time count
  rhs_ncalls = 0;
  rhs_ncalls_i = 0;
  rhs_ncalls_e = 0;
  
  return 0;
}

/////////////////////////////////////////////////////

void Solver::addTimestepMonitor(TimestepMonitorFunc f) {
  timestep_monitors.push_front(f);
}

void Solver::removeTimestepMonitor(TimestepMonitorFunc f) {
  timestep_monitors.remove(f);
}

int Solver::call_timestep_monitors(BoutReal simtime, BoutReal lastdt) {
  if(!monitor_timestep)
    return 0;
  
  for(const auto& monitor : timestep_monitors) {
    // Call each monitor one by one
    int ret = monitor(this, simtime, lastdt);
    if(ret)
      return ret; // Return first time an error is encountered
  }
  
  // Call physics model monitor
  if(model) {
    int ret = model->runTimestepMonitor(simtime, lastdt);
    if(ret)
      return ret; // Return first time an error is encountered
  }
  return 0;
}

void Solver::setRestartDir(const string &dir) {
  restartdir = dir;
}

/**************************************************************************
 * Useful routines (protected)
 **************************************************************************/

int Solver::getLocalN() {

  /// Cache the value, so this is not repeatedly called.
  /// This value should not change after initialisation
  static int cacheLocalN = -1;
  if(cacheLocalN != -1) {
    return cacheLocalN;
  }
  
  ASSERT0(initialised); // Must be initialised
  
  int n2d = n2Dvars();
  int n3d = n3Dvars();
  
  int ncz = mesh->LocalNz;
  int MYSUB = mesh->yend - mesh->ystart + 1;

  int local_N = (mesh->xend - mesh->xstart + 1) *
    (mesh->yend - mesh->ystart + 1)*(n2d + ncz*n3d); // NOTE: Not including extra toroidal point

  //////////// How many variables have evolving boundaries?
  
  int n2dbndry = 0;
  for(const auto& f : f2d) {
    if(f.evolve_bndry)
      n2dbndry++;
  }
  
  int n3dbndry = 0;
  for(const auto& f : f3d) {
    if(f.evolve_bndry)
      n3dbndry++;
  }

  //////////// Find boundary regions ////////////
  
  // Y up
  for(RangeIterator xi = mesh->iterateBndryUpperY(); !xi.isDone(); xi++) {
    local_N +=  (mesh->LocalNy - mesh->yend - 1) * (n2dbndry + ncz * n3dbndry);
  }
  
  // Y down
  for(RangeIterator xi = mesh->iterateBndryLowerY(); !xi.isDone(); xi++) {
    local_N +=  mesh->ystart * (n2dbndry + ncz * n3dbndry);
  }
  
  // X inner
  if(mesh->firstX() && !mesh->periodicX) {
    local_N += mesh->xstart * MYSUB * (n2dbndry + ncz * n3dbndry);
    output.write("\tBoundary region inner X\n");
  }

  // X outer
  if(mesh->lastX() && !mesh->periodicX) {
    local_N += (mesh->LocalNx - mesh->xend - 1) * MYSUB * (n2dbndry + ncz * n3dbndry);
    output.write("\tBoundary region outer X\n");
  }
  
  cacheLocalN = local_N;

  return local_N;
}

Solver* Solver::create(Options *opts) {  
  return SolverFactory::getInstance()->createSolver(opts);
}

Solver* Solver::create(SolverType &type, Options *opts) {  
  return SolverFactory::getInstance()->createSolver(type, opts);
}

/**************************************************************************
 * Looping over variables
 *
 * NOTE: This part is very inefficient, and should be replaced ASAP
 * Is the interleaving of variables needed or helpful to the solver?
 **************************************************************************/

/// Perform an operation at a given (jx,jy) location, moving data between BOUT++ and CVODE
void Solver::loop_vars_op(int jx, int jy, BoutReal *udata, int &p, SOLVER_VAR_OP op, bool bndry) {
  int jz;
 
  switch(op) {
  case LOAD_VARS: {
    /// Load variables from IDA into BOUT++
    
    // Loop over 2D variables
    for(const auto& f : f2d) {
      if(bndry && !f.evolve_bndry)
        continue;
      (*f.var)(jx, jy) = udata[p];
      p++;
    }
    
    for (jz=0; jz < mesh->LocalNz; jz++) {
      
      // Loop over 3D variables
      for(const auto& f : f3d) {
        if(bndry && !f.evolve_bndry)
          continue;
        (*f.var)(jx, jy, jz) = udata[p];
        p++;
      }  
    }
    break;
  }
  case LOAD_DERIVS: {
    /// Load derivatives from IDA into BOUT++
    /// Used for preconditioner
    
    // Loop over 2D variables
    for(const auto& f : f2d) {
      if(bndry && !f.evolve_bndry)
        continue;
      (*f.F_var)(jx, jy) = udata[p];
      p++;
    }
    
    for (jz=0; jz < mesh->LocalNz; jz++) {
      
      // Loop over 3D variables
      for(const auto& f : f3d) {
        if(bndry && !f.evolve_bndry)
          continue;
        (*f.F_var)(jx, jy, jz) = udata[p];
        p++;
      }  
    }
    
    break;
  }
  case SET_ID: {
    /// Set the type of equation (Differential or Algebraic)
    
    // Loop over 2D variables
    for(const auto& f : f2d) {
      if(bndry && !f.evolve_bndry)
	continue;
      if(f.constraint) {
	udata[p] = 0;
      }else {
	udata[p] = 1;
      }
      p++;
    }
    
    for (jz=0; jz < mesh->LocalNz; jz++) {
      
      // Loop over 3D variables
      for(const auto& f : f3d) {
        if(bndry && !f.evolve_bndry)
	  continue;
	if(f.constraint) {
	  udata[p] = 0;
	}else {
	  udata[p] = 1;
	}
	p++;
      }
    }
    
    break;
  }
  case SAVE_VARS: {
    /// Save variables from BOUT++ into IDA (only used at start of simulation)
    
    // Loop over 2D variables
    for(const auto& f : f2d) {
      if(bndry && !f.evolve_bndry)
        continue;
      udata[p] = (*f.var)(jx, jy);
      p++;
    }
    
    for (jz=0; jz < mesh->LocalNz; jz++) {
      
      // Loop over 3D variables
      for(const auto& f : f3d) {
        if(bndry && !f.evolve_bndry)
          continue;
        udata[p] = (*f.var)(jx, jy, jz);
        p++;
      }  
    }
    break;
  }
    /// Save time-derivatives from BOUT++ into CVODE (returning RHS result)
  case SAVE_DERIVS: {
    
    // Loop over 2D variables
    for(const auto& f : f2d) {
      if(bndry && !f.evolve_bndry)
        continue;
      udata[p] = (*f.F_var)(jx, jy);
      p++;
    }
    
    for (jz=0; jz < mesh->LocalNz; jz++) {
      
      // Loop over 3D variables
      for(const auto& f : f3d) {
        if(bndry && !f.evolve_bndry)
          continue;
        udata[p] = (*f.F_var)(jx, jy, jz);
        p++;
      }
    }
    break;
  }
  }
}

/// Loop over variables and domain. Used for all data operations for consistency
void Solver::loop_vars(BoutReal *udata, SOLVER_VAR_OP op) {
  int jx, jy;
  int p = 0; // Counter for location in udata array

  int MYSUB = mesh->yend - mesh->ystart + 1;

  // Inner X boundary
  if(mesh->firstX() && !mesh->periodicX) {
    for(jx=0;jx<mesh->xstart;jx++)
      for(jy=0;jy<MYSUB;jy++)
	loop_vars_op(jx, jy+mesh->ystart, udata, p, op, true);
  }

  // Lower Y boundary region
  for(RangeIterator xi = mesh->iterateBndryLowerY(); !xi.isDone(); xi++) {
    for(jy=0;jy<mesh->ystart;jy++)
      loop_vars_op(*xi, jy, udata, p, op, true);
  }

  // Bulk of points
  for (jx=mesh->xstart; jx <= mesh->xend; jx++)
    for (jy=mesh->ystart; jy <= mesh->yend; jy++)
      loop_vars_op(jx, jy, udata, p, op, false);
  
  // Upper Y boundary condition
  for(RangeIterator xi = mesh->iterateBndryUpperY(); !xi.isDone(); xi++) {
    for(jy=mesh->yend+1;jy<mesh->LocalNy;jy++)
      loop_vars_op(*xi, jy, udata, p, op, true);
  }

  // Outer X boundary
  if(mesh->lastX() && !mesh->periodicX) {
    for(jx=mesh->xend+1;jx<mesh->LocalNx;jx++)
      for(jy=mesh->ystart;jy<=mesh->yend;jy++)
	loop_vars_op(jx, jy, udata, p, op, true);
  }
}

void Solver::load_vars(BoutReal *udata) {
  // Make sure data is allocated
  for(const auto& f : f2d) 
    f.var->allocate();
  for(const auto& f : f3d) {
    f.var->allocate();
    f.var->setLocation(f.location);
  }

  loop_vars(udata, LOAD_VARS);

  // Mark each vector as either co- or contra-variant

  for(const auto& v : v2d) 
    v.var->covariant = v.covariant;
  for(const auto& v : v3d) 
    v.var->covariant = v.covariant;
}

void Solver::load_derivs(BoutReal *udata) {
  // Make sure data is allocated
  for(const auto& f : f2d) 
    f.F_var->allocate();
  for(const auto& f : f3d) {
    f.F_var->allocate();
    f.F_var->setLocation(f.location);
  }

  loop_vars(udata, LOAD_DERIVS);

  // Mark each vector as either co- or contra-variant

  for(const auto& v : v2d) 
    v.F_var->covariant = v.covariant;
  for(const auto& v : v3d) 
    v.F_var->covariant = v.covariant;
}

// This function only called during initialisation
void Solver::save_vars(BoutReal *udata) {
  for(const auto& f : f2d) 
    if(!f.var->isAllocated())
      throw BoutException("Variable '%s' not initialised", f.name.c_str());

  for(const auto& f : f3d) 
    if(!f.var->isAllocated())
      throw BoutException("Variable '%s' not initialised", f.name.c_str());
  
  // Make sure vectors in correct basis
  for(const auto& v : v2d) {
    if(v.covariant) {
      v.var->toCovariant();
    }else
      v.var->toContravariant();
  }
  for(const auto& v : v3d) {
    if(v.covariant) {
      v.var->toCovariant();
    }else
      v.var->toContravariant();
  }

  loop_vars(udata, SAVE_VARS);
}

void Solver::save_derivs(BoutReal *dudata) {
  // Make sure vectors in correct basis
  for(const auto& v : v2d) {
    if(v.covariant) {
      v.F_var->toCovariant();
    }else
      v.F_var->toContravariant();
  }
  for(const auto& v : v3d) {
    if(v.covariant) {
      v.F_var->toCovariant();
    }else
      v.F_var->toContravariant();
  }

  // Make sure 3D fields are at the correct cell location
  for(const auto& f : f3d) {
    if(f.location != (f.F_var)->getLocation()) {
      //output.write("SOLVER: Interpolating\n");
      *(f.F_var) = interp_to(*(f.F_var), f.location);
    }
  }

  loop_vars(dudata, SAVE_DERIVS);
}

void Solver::set_id(BoutReal *udata) {
  loop_vars(udata, SET_ID);
}


/*!
 * Returns a Field3D containing the global indices
 *
 */
const Field3D Solver::globalIndex(int localStart) {
  Field3D index = -1; // Set to -1, indicating out of domain

  int n2d = f2d.size();
  int n3d = f3d.size();

  int ind = localStart;

  // Find how many boundary cells are evolving
  int n2dbndry = 0;
  for(const auto& f : f2d) {
    if(f.evolve_bndry)
      ++n2dbndry;
  }
  int n3dbndry = 0;
  for(const auto& f : f3d) {
    if(f.evolve_bndry)
      n3dbndry++;
  }

  if(n2dbndry + n3dbndry > 0) {
    // Some boundary points evolving
    
    // Inner X boundary
    if(mesh->firstX() && !mesh->periodicX) {
      for(int jx=0;jx<mesh->xstart;jx++)
        for(int jy=mesh->ystart;jy<=mesh->yend;jy++) {
          // Zero index contains 2D and 3D variables
          index(jx, jy, 0) = ind;
          ind += n2dbndry + n3dbndry;
          for(int jz=1;jz<mesh->LocalNz; jz++) {
            index(jx, jy, jz) = ind;
            ind += n3dbndry;
          }
        }
    }
    
    // Lower Y boundary region
    for(RangeIterator xi = mesh->iterateBndryLowerY(); !xi.isDone(); xi++) {
      for(int jy=0;jy<mesh->ystart;jy++) {
        index(*xi, jy, 0) = ind;
        ind += n2dbndry + n3dbndry;
        for(int jz=1;jz<mesh->LocalNz; jz++) {
          index(*xi, jy, jz) = ind;
          ind += n3dbndry;
        }
      }
    }
  }
  
  // Bulk of points
  for (int jx=mesh->xstart; jx <= mesh->xend; jx++)
    for (int jy=mesh->ystart; jy <= mesh->yend; jy++) {
      index(jx, jy, 0) = ind;
      ind += n2d + n3d;
      for(int jz=1;jz<mesh->LocalNz; jz++) {
        index(jx, jy, jz) = ind;
        ind += n3d;
      }
    }

  if(n2dbndry + n3dbndry > 0) {
    // Some boundary points evolving
    
    // Upper Y boundary condition
    for(RangeIterator xi = mesh->iterateBndryUpperY(); !xi.isDone(); xi++) {
      for(int jy=mesh->yend+1;jy<mesh->LocalNy;jy++) {
        index(*xi, jy, 0) = ind;
        ind += n2dbndry + n3dbndry;
        for(int jz=1;jz<mesh->LocalNz; jz++) {
          index(*xi, jy, jz) = ind;
          ind += n3dbndry;
        }
      }
    }
    
    // Outer X boundary
    if(mesh->lastX() && !mesh->periodicX) {
      for(int jx=mesh->xend+1;jx<mesh->LocalNx;jx++)
        for(int jy=mesh->ystart;jy<=mesh->yend;jy++) {
          index(jx, jy, 0) = ind;
          ind += n2dbndry + n3dbndry;
          for(int jz=1;jz<mesh->LocalNz; jz++) {
            index(jx, jy, jz) = ind;
            ind += n3dbndry;
          }
        }
    }
  }
  
  // Should have included all evolving variables
  ASSERT1(ind == localStart + getLocalN());
  
  // Now swap guard cells
  mesh->communicate(index);
  
  return index;
}

/**************************************************************************
 * Running user-supplied functions
 **************************************************************************/

void Solver::setSplitOperator(rhsfunc fC, rhsfunc fD) {
  split_operator = true;
  phys_conv = fC;
  phys_diff = fD;
}

int Solver::run_rhs(BoutReal t) {
  int status;
  
  Timer timer("rhs");
  
  if(split_operator) {
    // Run both parts
    
    int nv = getLocalN();
    // Create two temporary arrays for system state
    Array<BoutReal> tmp(nv);
    Array<BoutReal> tmp2(nv);
    
    save_vars(tmp.begin()); // Copy variables into tmp
    pre_rhs(t);
    if(model) {
      status = model->runConvective(t);
    }else 
      status = (*phys_conv)(t);
    post_rhs(t); // Check variables, apply boundary conditions
    
    load_vars(tmp.begin()); // Reset variables
    save_derivs(tmp.begin()); // Save time derivatives
    pre_rhs(t);
    if(model) {
      status = model->runDiffusive(t, false);
    }else
      status = (*phys_diff)(t);
    post_rhs(t);
    save_derivs(tmp2.begin()); // Save time derivatives
    for(BoutReal *t = tmp.begin(), *t2 = tmp2.begin(); t != tmp.end(); ++t, ++t2)
      *t += *t2;
    load_derivs(tmp.begin()); // Put back time-derivatives
  }else {
    pre_rhs(t);
    if(model) {
      status = model->runRHS(t);
    }else
      status = (*phys_run)(t);
    post_rhs(t);
  }

  // If using Method of Manufactured Solutions
  add_mms_sources(t);

  rhs_ncalls++;
  rhs_ncalls_e++;
  rhs_ncalls_i++;
  return status;
}

/// NOTE: This calls add_mms_sources
int Solver::run_convective(BoutReal t) {
  int status;
  
  Timer timer("rhs");
  pre_rhs(t);
  if(split_operator) {
    if(model) {
      status = model->runConvective(t);
    }else
      status = (*phys_conv)(t);
  }else {
    // Zero if not split
    for(const auto& f : f3d)
      *(f.F_var) = 0.0;
    for(const auto& f : f2d)
      *(f.F_var) = 0.0;
    status = 0;
  }
  post_rhs(t);
  
  // If using Method of Manufactured Solutions
  add_mms_sources(t);
  
  rhs_ncalls++;
  rhs_ncalls_e++;
  return status;
}

int Solver::run_diffusive(BoutReal t, bool linear) {
  int status = 0;
  
  Timer timer("rhs");
  pre_rhs(t);
  if(split_operator) {

    if(model) {
      status = model->runDiffusive(t, linear);
    }else 
      status = (*phys_diff)(t);
    post_rhs(t);
  }else {
    // Return total
    if(model) {
      status = model->runRHS(t);
    }else
      status = (*phys_run)(t);
  }
  rhs_ncalls_i++;
  return status;
}

void Solver::pre_rhs(BoutReal t) {

  // Apply boundary conditions to the values
  for(const auto& f : f2d) {
    if(!f.constraint) // If it's not a constraint
      f.var->applyBoundary(t);
  }
  
  for(const auto& f : f3d) {
    if(!f.constraint)
      f.var->applyBoundary(t);
  }
  
}

void Solver::post_rhs(BoutReal t) {
#ifdef CHECK
  for(const auto& f : f3d) {
    if(!f.F_var->isAllocated())
      throw BoutException("Time derivative for '%s' not set", f.name.c_str());
  }
#endif
  // Make sure vectors in correct basis
  for(const auto& v : v2d) {
    if(v.covariant) {
      v.F_var->toCovariant();
    }else
      v.F_var->toContravariant();
  }
  for(const auto& v : v3d) {
    if(v.covariant) {
      v.F_var->toCovariant();
    }else
      v.F_var->toContravariant();
  }

  // Make sure 3D fields are at the correct cell location
  for(const auto& f : f3d) {
    if(f.location != (f.F_var)->getLocation()) {
      //output.write("SOLVER: Interpolating\n");
      *(f.F_var) = interp_to(*(f.F_var), f.location);
    }
  }

  // Apply boundary conditions to the time-derivatives
  for(const auto& f : f2d) {
    if(!f.constraint && f.evolve_bndry) // If it's not a constraint and if the boundary is evolving
      f.var->applyTDerivBoundary();
  }
  
  for(const auto& f : f3d) {
    if(!f.constraint && f.evolve_bndry)
      f.var->applyTDerivBoundary();
  }
#if CHECK > 2
  msg_stack.push("Solver checking time derivatives");
  for(const auto& f : f3d) {
    msg_stack.push("Variable: %s", f.name.c_str());
    checkData(*f.F_var);
    msg_stack.pop();
  }
  msg_stack.pop();
#endif
}

bool Solver::varAdded(const string &name) {
  for(const auto& f : f2d) {
    if(f.name == name)
      return true;
  }
  
  for(const auto& f : f3d) {
    if(f.name == name)
      return true;
  }
  
  for(const auto& f : v2d) {
    if(f.name == name)
      return true;
  }
  
  for(const auto& f : v3d) {
    if(f.name == name)
      return true;
  }
  
  return false;
}

bool Solver::have_user_precon() {
  if(model)
    return model->hasPrecon();
  
  return prefunc != 0;
}

int Solver::run_precon(BoutReal t, BoutReal gamma, BoutReal delta) {
  if(!have_user_precon())
    return 1;

  if(model)
    return model->runPrecon(t, gamma, delta);
  
  return (*prefunc)(t, gamma, delta);
}

// Add source terms to time derivatives
void Solver::add_mms_sources(BoutReal t) {
  if(!mms)
    return;

  FieldFactory *fact = FieldFactory::get();
    
  // Iterate over 2D variables
  for(const auto& f : f2d) {
    *f.F_var += fact->create2D("source", Options::getRoot()->getSection(f.name), mesh, (f.var)->getLocation(), t);
  }
  
  for(const auto& f : f3d) {
    *f.F_var += fact->create3D("source", Options::getRoot()->getSection(f.name), mesh, (f.var)->getLocation(), t);
  }
}

// Calculate 
void Solver::calculate_mms_error(BoutReal t) {
  FieldFactory *fact = FieldFactory::get();
  
  for(const auto& f : f3d) {
    Field3D solution = fact->create3D("solution", Options::getRoot()->getSection(f.name), mesh, (f.var)->getLocation(), t);
    
    *(f.MMS_err) = *(f.var) - solution;
  }
}
