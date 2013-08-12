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

#include "solverfactory.hxx"

#include <bout/sys/timer.hxx>
#include <msg_stack.hxx>
#include <output.hxx>

// Static member variables

int* Solver::pargc = 0;
char*** Solver::pargv = 0;

/**************************************************************************
 * Constructor
 **************************************************************************/

Solver::Solver(Options *opts) : options(opts), model(0) {
  if(options == NULL)
    options = Options::getRoot()->getSection("solver");

  // Set flags to defaults
  has_constraints = false;
  initialised = false;

  // Zero timing
  rhs_ncalls = 0;

  // Restart directory
  if(options->isSet("restartdir")) {
    // Solver-specific restart directory
    options->get("restartdir", restartdir, "data");
  }else {
    // Use the root data directory
    Options::getRoot()->get("datadir", restartdir, "data");
  }
  
  // Restart option
  Options::getRoot()->get("restart", restarting, false);

  // Split operator
  split_operator = false;
  max_dt = -1.0;
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

  VarStr<Field2D> d;
  
  d.constraint = false;
  d.var = &v;
  d.F_var = &ddt(v);
  d.name = string(name);

  f2d.push_back(d);

#ifdef TRACK
  var.name = name;
#endif

  /// Generate initial perturbations.
  /// NOTE: This could be done in init, but this would prevent the user
  ///       from modifying the initial perturbation (e.g. to prevent unphysical situations)
  ///       before it's loaded into the solver. If restarting, this perturbation
  ///       will be over-written anyway
  initial_profile(name, v);
  v.applyBoundary();

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
  
  f3d.push_back(d);

#ifdef TRACK
  var.name = name;
#endif

  initial_profile(name, v);
  v.applyBoundary(); // Make sure initial profile obeys boundary conditions
  v.setLocation(d.location); // Restore location if changed
                
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
  v.applyBoundary();

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

  v.applyBoundary();

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
    bout_error("ERROR: This solver doesn't support constraints\n");

  if(initialised)
    bout_error("Error: Cannot add constraints to solver after initialisation\n");

  if(name == NULL)
    bout_error("WARNING: Constraint requested for variable with NULL name\n");
  
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
    bout_error("ERROR: This solver doesn't support constraints\n");

  if(initialised)
    bout_error("Error: Cannot add constraints to solver after initialisation\n");

  if(name == NULL)
    bout_error("WARNING: Constraint requested for variable with NULL name\n");

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
    bout_error("ERROR: This solver doesn't support constraints\n");

  if(initialised)
    bout_error("Error: Cannot add constraints to solver after initialisation\n");

  if(name == NULL)
    bout_error("WARNING: Constraint requested for variable with NULL name\n");
    
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
    bout_error("ERROR: This solver doesn't support constraints\n");

  if(initialised)
    bout_error("Error: Cannot add constraints to solver after initialisation\n");

  if(name == NULL)
    bout_error("WARNING: Constraint requested for variable with NULL name\n");

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
 * Initialisation
 **************************************************************************/

int Solver::solve() {
  /// Get options
  Options *options = Options::getRoot();
  int NOUT;
  OPTION(options, NOUT, 1);
  BoutReal TIMESTEP;
  OPTION(options, TIMESTEP, 1.0);

  // Initialise
  if(init(restarting, NOUT, TIMESTEP)) {
    output.write("Failed to initialise solver-> Aborting\n");
    return 1;
  }
  
  if (!restarting) {
    /// Write initial state as time-point 0
    
    // Run RHS once to ensure all variables set
    if (run_rhs(0.0)) {
      output.write("Physics RHS call failed\n");
      return 1;
    }
    
    dump.write();
  }
  
  /// Run the solver
  output.write("Running simulation\n\n");
  int status;
  try {
    time_t start_time = time((time_t*) NULL);
    output.write("\nRun started at  : %s\n", ctime(&start_time));
   
    Timer timer("run"); // Start timer
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
  }catch(BoutException *e) {
    output << "Error encountered during initialisation\n";
    output << e->what() << endl;
    return 1;
  }

  return 0;
}


/**************************************************************************
 * Initialisation
 **************************************************************************/

int Solver::init(bool restarting, int nout, BoutReal tstep) {
  
#ifdef CHECK
  int msg_point = msg_stack.push("Solver::init()");
#endif

  if(initialised)
    throw BoutException("ERROR: Solver is already initialised\n");

  output.write("Initialising solver\n");

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
  
  // Set up restart options
  restart = Datafile(options->getSection("restart"));
  
  /// Add basic variables to the restart file
  restart.add(simtime,  "tt",    0);
  restart.add(iteration, "hist_hi", 0);
  
  MPI_Comm_size(BoutComm::get(), &NPES);
  MPI_Comm_rank(BoutComm::get(), &MYPE);
  
  restart.add(NPES, "NPES", 0);
  restart.add(mesh->NXPE, "NXPE", 0);

  /// Add variables to the restart and dump files.
  /// NOTE: Since vector components are already in the field arrays,
  ///       only loop over scalars, not vectors
  for(vector< VarStr<Field2D> >::iterator it = f2d.begin(); it != f2d.end(); it++) {
    // Add to restart file (not appending)
    restart.add(*(it->var), it->name.c_str(), 0);
    
    // Add to dump file (appending)
    dump.add(*(it->var), it->name.c_str(), 1);
    
    /// NOTE: Initial perturbations have already been set in add()
    
    /// Make sure boundary condition is satisfied
    it->var->applyBoundary();
  }  
  for(vector< VarStr<Field3D> >::iterator it = f3d.begin(); it != f3d.end(); it++) {
    // Add to restart file (not appending)
    restart.add(*(it->var), it->name.c_str(), 0);
    
    // Add to dump file (appending)
    dump.add(*(it->var), it->name.c_str(), 1);
    
    /// Make sure boundary condition is satisfied
    it->var->applyBoundary();
  }

  if(restarting) {
    /// Load state from the restart file
    
    // Copy processor numbers for comparison after. Very useful for checking
    // that the restart file is for the correct number of processors etc.
    int tmp_NP = NPES;
    int tmp_NX = mesh->NXPE;
    
#ifdef CHECK
    int msg_pt2 = msg_stack.push("Loading restart file");
#endif
    
    /// Load restart file
    if(!restart.openr("%s/BOUT.restart.%s", restartdir.c_str(), restartext.c_str()))
      throw new BoutException("Error: Could not open restart file\n");
    if(!restart.read())
      throw new BoutException("Error: Could not read restart file\n");
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

#ifdef CHECK
    msg_stack.pop(msg_pt2);
#endif
    
  }else {
    // Not restarting
    simtime = 0.0; iteration = 0;
  }
  
  /// Open the restart file for writing
  if(!restart.openw("%s/BOUT.restart.%s", restartdir.c_str(), restartext.c_str()))
    throw new BoutException("Error: Could not open restart file for writing\n");
  
  /// Mark as initialised. No more variables can be added
  initialised = true;

#ifdef CHECK
  msg_stack.pop(msg_point);
#endif

  return 0;
}

void Solver::addMonitor(MonitorFunc f) {
  monitors.push_front(f);
}

void Solver::removeMonitor(MonitorFunc f) {
  monitors.remove(f);
}

int Solver::call_monitors(BoutReal simtime, int iter, int NOUT) {
  for(std::list<MonitorFunc>::iterator it = monitors.begin(); it != monitors.end(); it++) {
    // Call each monitor one by one
    int ret = (*it)(this, simtime,iter, NOUT);
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
  int n2d = n2Dvars();
  int n3d = n3Dvars();
  
  int ncz = mesh->ngz-1;
  int MYSUB = mesh->yend - mesh->ystart + 1;

  int local_N = (mesh->xend - mesh->xstart + 1) *
    (mesh->yend - mesh->ystart + 1)*(n2d + ncz*n3d); // NOTE: Not including extra toroidal point

  //////////// Find boundary regions ////////////
  
  // Y up
  for(RangeIterator xi = mesh->iterateBndryUpperY(); !xi.isDone(); xi++) {
    local_N +=  (mesh->ngy - mesh->yend - 1) * (n2d + ncz * n3d);
  }
  
  // Y down
  for(RangeIterator xi = mesh->iterateBndryLowerY(); !xi.isDone(); xi++) {
    local_N +=  mesh->ystart * (n2d + ncz * n3d);
  }
  
  // X inner
  if(mesh->firstX() && !mesh->periodicX) {
    local_N += mesh->xstart * MYSUB * (n2d + ncz * n3d);
    output.write("\tBoundary region inner X\n");
  }

  // X outer
  if(mesh->lastX() && !mesh->periodicX) {
    local_N += (mesh->ngx - mesh->xend - 1) * MYSUB * (n2d + ncz * n3d);
    output.write("\tBoundary region outer X\n");
  }
  
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
 **************************************************************************/

/// Perform an operation at a given (jx,jy) location, moving data between BOUT++ and CVODE
void Solver::loop_vars_op(int jx, int jy, BoutReal *udata, int &p, SOLVER_VAR_OP op) {
  int i;
  int jz;
 
  int n2d = f2d.size();
  int n3d = f3d.size();

  switch(op) {
  case LOAD_VARS: {
    /// Load variables from IDA into BOUT++
    
    // Loop over 2D variables
    for(i=0;i<n2d;i++) {
      (*f2d[i].var)(jx, jy) = udata[p];
      p++;
    }
    
    for (jz=0; jz < mesh->ngz-1; jz++) {
      
      // Loop over 3D variables
      for(i=0;i<n3d;i++) {
	(*f3d[i].var)(jx, jy, jz) = udata[p];
	p++;
      }  
    }
    break;
  }
  case LOAD_DERIVS: {
    /// Load derivatives from IDA into BOUT++
    /// Used for preconditioner
    
    // Loop over 2D variables
    for(i=0;i<n2d;i++) {
      (*f2d[i].F_var)(jx, jy) = udata[p];
      p++;
    }
    
    for (jz=0; jz < mesh->ngz-1; jz++) {
      
      // Loop over 3D variables
      for(i=0;i<n3d;i++) {
	(*f3d[i].F_var)(jx, jy, jz) = udata[p];
	p++;
      }  
    }
    
    break;
  }
  case SAVE_VARS: {
    /// Save variables from BOUT++ into IDA (only used at start of simulation)
    
    // Loop over 2D variables
    for(i=0;i<n2d;i++) {
      udata[p] = (*f2d[i].var)(jx, jy);
      p++;
    }
    
    for (jz=0; jz < mesh->ngz-1; jz++) {
      
      // Loop over 3D variables
      for(i=0;i<n3d;i++) {
	udata[p] = (*f3d[i].var)(jx, jy, jz);
	p++;
      }  
    }
    break;
  }
    /// Save time-derivatives from BOUT++ into CVODE (returning RHS result)
  case SAVE_DERIVS: {
    
    // Loop over 2D variables
    for(i=0;i<n2d;i++) {
      udata[p] = (*f2d[i].F_var)(jx, jy);
      p++;
    }
    
    for (jz=0; jz < mesh->ngz-1; jz++) {
      
      // Loop over 3D variables
      for(i=0;i<n3d;i++) {
	udata[p] = (*f3d[i].F_var)(jx, jy, jz);
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
	loop_vars_op(jx, jy+mesh->ystart, udata, p, op);
  }

  // Lower Y boundary region
  for(RangeIterator xi = mesh->iterateBndryLowerY(); !xi.isDone(); xi++) {
    for(jy=0;jy<mesh->ystart;jy++)
      loop_vars_op(*xi, jy, udata, p, op);
  }

  // Bulk of points
  for (jx=mesh->xstart; jx <= mesh->xend; jx++)
    for (jy=mesh->ystart; jy <= mesh->yend; jy++)
      loop_vars_op(jx, jy, udata, p, op);
  
  // Upper Y boundary condition
  for(RangeIterator xi = mesh->iterateBndryUpperY(); !xi.isDone(); xi++) {
    for(jy=mesh->yend+1;jy<mesh->ngy;jy++)
      loop_vars_op(*xi, jy, udata, p, op);
  }

  // Outer X boundary
  if(mesh->lastX() && !mesh->periodicX) {
    for(jx=mesh->xend+1;jx<mesh->ngx;jx++)
      for(jy=mesh->ystart;jy<=mesh->yend;jy++)
	loop_vars_op(jx, jy, udata, p, op);
  }
}

void Solver::load_vars(BoutReal *udata) {
  unsigned int i;
  
  // Make sure data is allocated
  for(i=0;i<f2d.size();i++)
    f2d[i].var->allocate();
  for(i=0;i<f3d.size();i++) {
    f3d[i].var->allocate();
    f3d[i].var->setLocation(f3d[i].location);
  }

  loop_vars(udata, LOAD_VARS);

  // Mark each vector as either co- or contra-variant

  for(i=0;i<v2d.size();i++)
    v2d[i].var->covariant = v2d[i].covariant;
  for(i=0;i<v3d.size();i++)
    v3d[i].var->covariant = v3d[i].covariant;
}

void Solver::load_derivs(BoutReal *udata) {
  unsigned int i;
  
  // Make sure data is allocated
  for(i=0;i<f2d.size();i++)
    f2d[i].F_var->allocate();
  for(i=0;i<f3d.size();i++) {
    f3d[i].F_var->allocate();
    f3d[i].F_var->setLocation(f3d[i].location);
  }

  loop_vars(udata, LOAD_DERIVS);

  // Mark each vector as either co- or contra-variant

  for(i=0;i<v2d.size();i++)
    v2d[i].F_var->covariant = v2d[i].covariant;
  for(i=0;i<v3d.size();i++)
    v3d[i].F_var->covariant = v3d[i].covariant;
}

// This function only called during initialisation
int Solver::save_vars(BoutReal *udata) {
  unsigned int i;

  for(i=0;i<f2d.size();i++)
    if(f2d[i].var->getData() == (BoutReal**) NULL)
      return(1);

  for(i=0;i<f3d.size();i++)
    if(f3d[i].var->getData() == (BoutReal***) NULL)
      return(1);
  
  // Make sure vectors in correct basis
  for(i=0;i<v2d.size();i++) {
    if(v2d[i].covariant) {
      v2d[i].var->toCovariant();
    }else
      v2d[i].var->toContravariant();
  }
  for(i=0;i<v3d.size();i++) {
    if(v3d[i].covariant) {
      v3d[i].var->toCovariant();
    }else
      v3d[i].var->toContravariant();
  }

  loop_vars(udata, SAVE_VARS);

  return(0);
}

void Solver::save_derivs(BoutReal *dudata) {
  unsigned int i;

  // Make sure vectors in correct basis
  for(i=0;i<v2d.size();i++) {
    if(v2d[i].covariant) {
      v2d[i].F_var->toCovariant();
    }else
      v2d[i].F_var->toContravariant();
  }
  for(i=0;i<v3d.size();i++) {
    if(v3d[i].covariant) {
      v3d[i].F_var->toCovariant();
    }else
      v3d[i].F_var->toContravariant();
  }

  // Make sure 3D fields are at the correct cell location
  for(vector< VarStr<Field3D> >::iterator it = f3d.begin(); it != f3d.end(); it++) {
    if((*it).location != ((*it).F_var)->getLocation()) {
      //output.write("SOLVER: Interpolating\n");
      *((*it).F_var) = interp_to(*((*it).F_var), (*it).location);
    }
  }

  loop_vars(dudata, SAVE_DERIVS);
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
    
    static int nv;
    static BoutReal *tmp = NULL, *tmp2;
    if(tmp == NULL) {
      nv = getLocalN();
      tmp = new BoutReal[nv];
      tmp2 = new BoutReal[nv];
    }
    save_vars(tmp); // Copy variables into tmp
    if(model) {
      status = model->runConvective(t);
    }else 
      status = (*phys_conv)(t);
    post_rhs(); // Check variables, apply boundary conditions
    
    load_vars(tmp); // Reset variables
    save_derivs(tmp); // Save time derivatives
    if(model) {
      status = model->runDiffusive(t);
    }else
      status = (*phys_diff)(t);
    post_rhs();
    save_derivs(tmp2); // Save time derivatives
    for(int i=0;i<nv;i++)
      tmp[i] += tmp2[i];
    load_derivs(tmp); // Put back time-derivatives
  }else {
    if(model) {
      status = model->runRHS(t);
    }else
      status = (*phys_run)(t);
    post_rhs();
  }

  rhs_ncalls++;

  return status;
}

int Solver::run_convective(BoutReal t) {
  int status;
  
  Timer timer("rhs");
  
  if(split_operator) {
    if(model) {
      status = model->runConvective(t);
    }else
      status = (*phys_conv)(t);
  }else {
    // Return total
    if(model) {
      status = model->runRHS(t);
    }else
      status = (*phys_run)(t);
  }
  post_rhs();
  
  rhs_ncalls++;
  
  return status;
}

int Solver::run_diffusive(BoutReal t) {
  int status = 0;
  
  Timer timer("rhs");

  if(split_operator) {
    if(model) {
      status = model->runDiffusive(t);
    }else 
      status = (*phys_diff)(t);
    post_rhs();
  }else {
    // Zero if not split
    for(vector< VarStr<Field3D> >::iterator it = f3d.begin(); it != f3d.end(); it++)
      *((*it).F_var) = 0.0;
    for(vector< VarStr<Field2D> >::iterator it = f2d.begin(); it != f2d.end(); it++)
      *((*it).F_var) = 0.0;
  }
  
  return status;
}

void Solver::post_rhs() {

  // Make sure vectors in correct basis
  for(int i=0;i<v2d.size();i++) {
    if(v2d[i].covariant) {
      v2d[i].F_var->toCovariant();
    }else
      v2d[i].F_var->toContravariant();
  }
  for(int i=0;i<v3d.size();i++) {
    if(v3d[i].covariant) {
      v3d[i].F_var->toCovariant();
    }else
      v3d[i].F_var->toContravariant();
  }

  // Make sure 3D fields are at the correct cell location
  for(vector< VarStr<Field3D> >::iterator it = f3d.begin(); it != f3d.end(); it++) {
    if((*it).location != ((*it).F_var)->getLocation()) {
      //output.write("SOLVER: Interpolating\n");
      *((*it).F_var) = interp_to(*((*it).F_var), (*it).location);
    }
  }

  // Apply boundary conditions to the time-derivatives
  for(vector< VarStr<Field2D> >::iterator it = f2d.begin(); it != f2d.end(); it++) {
    if(!it->constraint) // If it's not a constraint
      it->var->applyTDerivBoundary();
  }
  
  for(vector< VarStr<Field3D> >::iterator it = f3d.begin(); it != f3d.end(); it++) {
    if(!it->constraint)
      it->var->applyTDerivBoundary();
  }

#ifdef CHECK
  msg_stack.push("Solver checking time derivatives");
  for(vector< VarStr<Field3D> >::iterator it = f3d.begin(); it != f3d.end(); it++)
    it->F_var->checkData();
  msg_stack.pop();
#endif
}

bool Solver::varAdded(const string &name) {
  for(vector< VarStr<Field2D> >::iterator it = f2d.begin(); it != f2d.end(); it++) {
    if(it->name == name)
      return true;
  }
  
  for(vector< VarStr<Field3D> >::iterator it = f3d.begin(); it != f3d.end(); it++) {
    if(it->name == name)
      return true;
  }
  
  for(vector< VarStr<Vector2D> >::iterator it = v2d.begin(); it != v2d.end(); it++) {
    if(it->name == name)
      return true;
  }
  
  for(vector< VarStr<Vector3D> >::iterator it = v3d.begin(); it != v3d.end(); it++) {
    if(it->name == name)
      return true;
  }
  
  return false;
}
