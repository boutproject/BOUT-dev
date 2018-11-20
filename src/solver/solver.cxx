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

#include <bout/boutcomm.hxx>
#include <bout/solver.hxx>
#include <string.h>
#include <time.h>

#include <bout/initialprofiles.hxx>
#include <bout/interpolation.hxx>
#include <bout/boutexception.hxx>

#include <bout/field_factory.hxx>

#include "bout/solverfactory.hxx"

#include <bout/sys/timer.hxx>
#include <bout/msg_stack.hxx>
#include <bout/output.hxx>
#include <bout/assert.hxx>

#include <bout/array.hxx>
#include "bout/region.hxx"

// Static member variables

int *Solver::pargc = nullptr;
char ***Solver::pargv = nullptr;

/**************************************************************************
 * Constructor
 **************************************************************************/

Solver::Solver(Options *opts) : options(opts), model(nullptr), prefunc(nullptr) {
  if(options == nullptr)
    options = Options::getRoot()->getSection("solver");

  // Set flags to defaults
  has_constraints = false;
  initialised = false;
  canReset = false;

  // Zero timing
  rhs_ncalls = 0;
  rhs_ncalls_e = 0;
  rhs_ncalls_i = 0;
  
  // Split operator
  split_operator = false;
  max_dt = -1.0;
  
  // Set simulation time and iteration count
  // This may be modified by restart
  simtime = 0.0; iteration = 0;
  
  // Output monitor
  options->get("monitor_timestep", monitor_timestep, false);
  
  // Method of Manufactured Solutions (MMS)
  options->get("mms", mms, false);
  options->get("mms_initialise", mms_initialise, mms);
}

/**************************************************************************
 * Destructor
 **************************************************************************/
Solver::~Solver(){
  //Ensure all MMS_err fields allocated here are destroyed etc.
  for(const auto& f : f3d) {
    if(f.MMS_err) {
      delete f.MMS_err;
    }
  }

  for(const auto& f : f2d) {
    if(f.MMS_err) {
      delete f.MMS_err;
    }
  }
}

/**************************************************************************
 * Add physics models
 **************************************************************************/

void Solver::setModel(PhysicsModel *m) {
  if(model)
    throw BoutException("Solver can only evolve one model");
  
  if(initialised)
    throw BoutException("Solver already initialised");
  
  // Initialise them model, which specifies which variables to evolve
  m->initialise(this);
  
  // Check if the model is split operator
  split_operator = m->splitOperator();
  
  model = m;
}

/**************************************************************************
 * Add fields
 **************************************************************************/

void Solver::add(Field2D &v, const std::string name) {
  TRACE("Adding 2D field: Solver::add(%s)", name.c_str());

  if (varAdded(name))
    throw BoutException("Variable '%s' already added to Solver", name.c_str());

  if (initialised)
    throw BoutException("Error: Cannot add to solver after initialisation\n");

  // Set boundary conditions
  v.setBoundary(name);
  ddt(v).copyBoundary(v); // Set boundary to be the same as v

  VarStr<Field2D> d;
  
  d.constraint = false;
  d.var = &v;
  d.F_var = &ddt(v);
  d.location = v.getLocation();
  d.covariant = false;
  d.name = name;
  
#ifdef TRACK
  v.name = name;
#endif

  /// Generate initial perturbations.
  /// NOTE: This could be done in init, but this would prevent the user
  ///       from modifying the initial perturbation (e.g. to prevent unphysical situations)
  ///       before it's loaded into the solver. If restarting, this perturbation
  ///       will be over-written anyway
  if (mms_initialise) {
    // Load solution at t = 0
    
    FieldFactory *fact = FieldFactory::get();
    
    v = fact->create2D("solution", Options::getRoot()->getSection(name), mesh);
  } else {
    initial_profile(name, v);
  }
  
  if (mms) {
    // Allocate storage for error variable
    d.MMS_err = new Field2D(0.0);
  } else {
    d.MMS_err = nullptr;
  }
  
  // Check if the boundary regions should be evolved
  // First get option from section "All"
  // then use that as default for specific section
  Options::getRoot()->getSection("all")->get("evolve_bndry", d.evolve_bndry, false);
  Options::getRoot()->getSection(name)->get("evolve_bndry", d.evolve_bndry, d.evolve_bndry);

  v.applyBoundary(true);

  f2d.push_back(d);
}

void Solver::add(Field3D &v, const std::string name) {
  TRACE("Adding 3D field: Solver::add(%s)", name.c_str());

#if CHECK > 0  
  if (varAdded(name))
    throw BoutException("Variable '%s' already added to Solver", name.c_str());
#endif

  if (initialised)
    throw BoutException("Error: Cannot add to solver after initialisation\n");

  // Set boundary conditions
  v.setBoundary(name);
  ddt(v).copyBoundary(v); // Set boundary to be the same as v

  if (mesh->StaggerGrids && (v.getLocation() != CELL_CENTRE)) {
    output_info.write("\tVariable %s shifted to %s\n", name.c_str(), strLocation(v.getLocation()));
    ddt(v).setLocation(v.getLocation()); // Make sure both at the same location
  }

  VarStr<Field3D> d;
  
  d.constraint = false;
  d.var = &v;
  d.F_var = &ddt(v);
  d.location = v.getLocation();
  d.covariant = false;
  d.name = name;
  
#ifdef TRACK
  v.name = name;
#endif

  if (mms_initialise) {
    // Load solution at t = 0
    FieldFactory *fact = FieldFactory::get();
    
    v = fact->create3D("solution", Options::getRoot()->getSection(name), mesh, v.getLocation());
    
  } else {
    initial_profile(name, v);
  }
  
  if (mms) {
    d.MMS_err = new Field3D(v.getMesh());
    (*d.MMS_err) = 0.0;
  } else {
    d.MMS_err = nullptr;
  }
  
  // Check if the boundary regions should be evolved
  // First get option from section "All"
  // then use that as default for specific section
  Options::getRoot()->getSection("all")->get("evolve_bndry", d.evolve_bndry, false);
  Options::getRoot()->getSection(name)->get("evolve_bndry", d.evolve_bndry, d.evolve_bndry);

  v.applyBoundary(true); // Make sure initial profile obeys boundary conditions
  v.setLocation(d.location); // Restore location if changed
  
  f3d.push_back(d);
}

void Solver::add(Vector2D &v, const std::string name) {
  TRACE("Adding 2D vector: Solver::add(%s)", name.c_str());
  
  if (varAdded(name))
    throw BoutException("Variable '%s' already added to Solver", name.c_str());

  if (initialised)
    throw BoutException("Error: Cannot add to solver after initialisation\n");

  // Set boundary conditions
  v.setBoundary(name);
  ddt(v).copyBoundary(v); // Set boundary to be the same as v
  
  VarStr<Vector2D> d;
  
  d.constraint = false;
  d.var = &v;
  d.F_var = &ddt(v);
  d.location = CELL_DEFAULT;
  d.covariant = v.covariant;
  d.name = name;
  // MMS errors set on individual components
  d.MMS_err = nullptr;

  v2d.push_back(d);

  /// NOTE: No initial_profile call, because this will be done for each
  ///       component individually.
  
  /// Add suffix, depending on co- /contravariance
  if (v.covariant) {
    add(v.x, d.name+"_x");
    add(v.y, d.name+"_y");
    add(v.z, d.name+"_z");
  } else {
    add(v.x, d.name+"x");
    add(v.y, d.name+"y");
    add(v.z, d.name+"z");
  }
  
  /// Make sure initial profile obeys boundary conditions
  v.applyBoundary(true);
}

void Solver::add(Vector3D &v, const std::string name) {
  TRACE("Adding 3D vector: Solver::add(%s)", name.c_str());
  
  if (varAdded(name))
    throw BoutException("Variable '%s' already added to Solver", name.c_str());

  if (initialised)
    throw BoutException("Error: Cannot add to solver after initialisation\n");
  
  // Set boundary conditions
  v.setBoundary(name);
  ddt(v).copyBoundary(v); // Set boundary to be the same as v

  VarStr<Vector3D> d;
  
  d.constraint = false;
  d.var = &v;
  d.F_var = &ddt(v);
  d.location = CELL_DEFAULT;
  d.covariant = v.covariant;
  d.name = name;
  // MMS errors set on individual components
  d.MMS_err = nullptr;
  
  v3d.push_back(d);

  // Add suffix, depending on co- /contravariance
  if (v.covariant) {
    add(v.x, d.name+"_x");
    add(v.y, d.name+"_y");
    add(v.z, d.name+"_z");
  } else {
    add(v.x, d.name+"x");
    add(v.y, d.name+"y");
    add(v.z, d.name+"z");
  }

  v.applyBoundary(true);
}

/**************************************************************************
 * Constraints
 **************************************************************************/

void Solver::constraint(Field2D &v, Field2D &C_v, const std::string name) {

  if (name.empty()) {
    throw BoutException("ERROR: Constraint requested for variable with empty name\n");
  }

  TRACE("Constrain 2D scalar: Solver::constraint(%s)", name.c_str());

#if CHECK > 0  
  if (varAdded(name))
    throw BoutException("Variable '%s' already added to Solver", name.c_str());
#endif

  if (!has_constraints)
    throw BoutException("ERROR: This solver doesn't support constraints\n");

  if (initialised)
    throw BoutException("Error: Cannot add constraints to solver after initialisation\n");

  VarStr<Field2D> d;
  
  d.constraint = true;
  d.var = &v;
  d.F_var = &C_v;
  d.name = name;

  f2d.push_back(d);
}

void Solver::constraint(Field3D &v, Field3D &C_v, const std::string name) {

  if (name.empty()) {
    throw BoutException("ERROR: Constraint requested for variable with empty name\n");
  }

  TRACE("Constrain 3D scalar: Solver::constraint(%s)", name.c_str());

#if CHECK > 0
  if (varAdded(name))
    throw BoutException("Variable '%s' already added to Solver", name.c_str());
#endif

  if (!has_constraints)
    throw BoutException("ERROR: This solver doesn't support constraints\n");

  if (initialised)
    throw BoutException("Error: Cannot add constraints to solver after initialisation\n");

  VarStr<Field3D> d;
  
  d.constraint = true;
  d.var = &v;
  d.F_var = &C_v;
  d.location = v.getLocation();
  d.name = name;
  
  f3d.push_back(d);
}

void Solver::constraint(Vector2D &v, Vector2D &C_v, const std::string name) {

  if (name.empty()) {
    throw BoutException("ERROR: Constraint requested for variable with empty name\n");
  }

  TRACE("Constrain 2D vector: Solver::constraint(%s)", name.c_str());

#if CHECK > 0  
  if (varAdded(name))
    throw BoutException("Variable '%s' already added to Solver", name.c_str());
#endif

  if (!has_constraints)
    throw BoutException("ERROR: This solver doesn't support constraints\n");

  if (initialised)
    throw BoutException("Error: Cannot add constraints to solver after initialisation\n");

  VarStr<Vector2D> d;
  
  d.constraint = true;
  d.var = &v;
  d.F_var = &C_v;
  d.covariant = v.covariant;
  d.name = name;
  
  v2d.push_back(d);

  // Add suffix, depending on co- /contravariance
  if (v.covariant) {
    constraint(v.x, C_v.x, d.name+"_x");
    constraint(v.y, C_v.y, d.name+"_x");
    constraint(v.z, C_v.z, d.name+"_x");
  } else {
    constraint(v.x, C_v.x, d.name+"x");
    constraint(v.y, C_v.y, d.name+"x");
    constraint(v.z, C_v.z, d.name+"x");
  }
}

void Solver::constraint(Vector3D &v, Vector3D &C_v, const std::string name) {

  if (name.empty()) {
    throw BoutException("ERROR: Constraint requested for variable with empty name\n");
  }

  TRACE("Constrain 3D vector: Solver::constraint(%s)", name.c_str());

#if CHECK > 0  
  if (varAdded(name))
    throw BoutException("Variable '%s' already added to Solver", name.c_str());
#endif

  if (!has_constraints)
    throw BoutException("ERROR: This solver doesn't support constraints\n");

  if (initialised)
    throw BoutException("Error: Cannot add constraints to solver after initialisation\n");

  VarStr<Vector3D> d;
  
  d.constraint = true;
  d.var = &v;
  d.F_var = &C_v;
  d.covariant = v.covariant;
  d.name = name;
  
  v3d.push_back(d);

  // Add suffix, depending on co- /contravariance
  if (v.covariant) {
    constraint(v.x, C_v.x, d.name+"_x");
    constraint(v.y, C_v.y, d.name+"_x");
    constraint(v.z, C_v.z, d.name+"_x");
  } else {
    constraint(v.x, C_v.x, d.name+"x");
    constraint(v.y, C_v.y, d.name+"x");
    constraint(v.z, C_v.z, d.name+"x");
  }
}

/**************************************************************************
 * Solver main loop: Initialise, run, and finish
 **************************************************************************/

int Solver::solve(int NOUT, BoutReal TIMESTEP) {
  
  Options *globaloptions = Options::getRoot(); // Default from global options
  
  if(NOUT < 0) {
    /// Get options
    OPTION(globaloptions, NOUT, 1);
    OPTION(globaloptions, TIMESTEP, 1.0);
    
    // Check specific solver options, which override global options
    OPTION(options, NOUT, NOUT);
    options->get("output_step", TIMESTEP, TIMESTEP);
  }

  /// syncronize timestep with those set to the monitors
  if (timestep > 0){
    if (!isMultiple(timestep,TIMESTEP)){
      throw BoutException("A monitor requested a timestep not compatible with the output_step!");
    }
    if (timestep < TIMESTEP*1.5){
      freqDefault=TIMESTEP/timestep+.5;
      NOUT*=freqDefault;
      TIMESTEP=timestep;
    } else {
      freqDefault = 1;
      // update old monitors
      int fac=timestep/TIMESTEP+.5;
      for (const auto &i: monitors){
        i->freq=i->freq*fac;
      }
    }
  }
  for (const auto &i: monitors){
    if (i->timestep < 0){
      i->timestep=timestep*freqDefault;
      i->freq=freqDefault;
    }
  }


  output_progress.write(_("Solver running for %d outputs with output timestep of %e\n"), NOUT, TIMESTEP);
  if (freqDefault > 1)
    output_progress.write(_("Solver running for %d outputs with monitor timestep of %e\n"),
                          NOUT/freqDefault, TIMESTEP*freqDefault);
  
  // Initialise
  if (init(NOUT, TIMESTEP)) {
    throw BoutException(_("Failed to initialise solver-> Aborting\n"));
  }
  initCalled=true;
  
  /// Run the solver
  output_info.write(_("Running simulation\n\n"));

  time_t start_time = time(nullptr);
  output_progress.write(_("\nRun started at  : %s\n"), toString(start_time).c_str());
  
  Timer timer("run"); // Start timer
  
  bool restart;
  OPTION(globaloptions, restart, false);
  bool append;
  OPTION(globaloptions, append, false);
  bool dump_on_restart;
  OPTION(globaloptions, dump_on_restart, !restart || !append);
  if ( dump_on_restart ) {
    /// Write initial state as time-point 0
    
    // Run RHS once to ensure all variables set
    if (run_rhs(simtime)) {
      throw BoutException("Physics RHS call failed\n");
    }
    
    // Call monitors so initial values are written to output dump files
    if (call_monitors(simtime, -1, NOUT)){
      throw BoutException("Initial monitor call failed!");
    }
  }
  
  int status;
  try {
    status = run();

    time_t end_time = time(nullptr);
    output_progress.write(_("\nRun finished at  : %s\n"), toString(end_time).c_str());
    output_progress.write(_("Run time : "));

    int dt = end_time - start_time;
    int i = static_cast<int>(dt / (60. * 60.));
    if (i > 0) {
      output_progress.write("%d h ", i);
      dt -= i*60*60;
    }
    i = static_cast<int>(dt / 60.);
    if (i > 0) {
      output_progress.write("%d m ", i);
      dt -= i*60;
    }
    output_progress.write("%d s\n", dt);
  } catch (BoutException &e) {
    output_error << "Error encountered in solver run\n";
    output_error << e.what() << endl;
    throw;
  }

  return status;
}


/**************************************************************************
 * Initialisation
 **************************************************************************/

int Solver::init(int UNUSED(nout), BoutReal UNUSED(tstep)) {
  
  TRACE("Solver::init()");

  if (initialised)
    throw BoutException(_("ERROR: Solver is already initialised\n"));

  output_progress.write(_("Initialising solver\n"));

  MPI_Comm_size(BoutComm::get(), &NPES);
  MPI_Comm_rank(BoutComm::get(), &MYPE);
  
  /// Mark as initialised. No more variables can be added
  initialised = true;

  return 0;
}

void Solver::outputVars(Datafile &outputfile, bool save_repeat) {
  /// Add basic variables to the file
  outputfile.addOnce(simtime,  "tt");
  outputfile.addOnce(iteration, "hist_hi");

  // Add 2D and 3D evolving fields to output file
  for(const auto& f : f2d) {
    // Add to dump file (appending)
    outputfile.add(*(f.var), f.name.c_str(), save_repeat);
  }  
  for(const auto& f : f3d) {
    // Add to dump file (appending)
    outputfile.add(*(f.var), f.name.c_str(), save_repeat);
    
    if(mms) {
      // Add an error variable
      outputfile.add(*(f.MMS_err), ("E_" + f.name).c_str(), save_repeat);
    }
  }
}

/////////////////////////////////////////////////////

/// Method to add a Monitor to the Solver
/// Note that behaviour changes if init() is called,
/// as the timestep cannot be changed afterwards
void Solver::addMonitor(Monitor * mon, MonitorPosition pos) {
  if (mon->timestep > 0){ // not default
    if (!initCalled && timestep < 0){
      timestep = mon->timestep;
    }
    if (!isMultiple(timestep,mon->timestep))
      throw BoutException(_("Couldn't add Monitor: %g is not a multiple of %g!")
                          ,timestep,mon->timestep);
    if (mon->timestep > timestep*1.5){
      mon->freq=(mon->timestep/timestep)+.5;
    } else { // mon.timestep is truly smaller
      if (initCalled)
        throw BoutException(_("Solver::addMonitor: Cannot reduce timestep \
(from %g to %g) after init is called!")
                            ,timestep,mon->timestep);
      int multi = timestep/mon->timestep+.5;
      timestep=mon->timestep;
      for (const auto &i: monitors){
        i->freq=i->freq*multi;
      }
      // update freqDefault so that monitors with no timestep are called at the
      // output frequency
      freqDefault *= multi;

      mon->freq=1;
    }
  } else {
    mon->freq = freqDefault;
  }
  mon->is_added = true; // Records that monitor has been added to solver so timestep should not be updated
  if(pos == Solver::FRONT) {
    monitors.push_front(mon);
  }else
    monitors.push_back(mon);
}

void Solver::removeMonitor(Monitor * f) {
  monitors.remove(f);
}

extern bool user_requested_exit;
int Solver::call_monitors(BoutReal simtime, int iter, int NOUT) {
  bool abort;
  MPI_Allreduce(&user_requested_exit,&abort,1,MPI_C_BOOL,MPI_LOR,MPI_COMM_WORLD);
  if(abort){
    NOUT=iter+1;
  }
  if(mms) {
    // Calculate MMS errors
    calculate_mms_error(simtime);
  }
  
  ++iter;
  try {
    // Call monitors
    for (const auto &it : monitors){
      if ((iter % it->freq)==0){
        // Call each monitor one by one
        int ret = it->call(this, simtime,iter/it->freq-1, NOUT/it->freq);
        if(ret)
          throw BoutException(_("Monitor signalled to quit"));
      }
    }
  } catch (BoutException &e) {
    for (const auto &it : monitors){
      it->cleanup();
    }
    output_error.write(_("Monitor signalled to quit\n"));
    throw;
  }

  if ( iter == NOUT ){
    for (const auto &it : monitors){
      it->cleanup();
    }
  }

  if (abort) {
    // restart file should be written by physics model
    output.write("User signalled to quit. Returning\n");
    return 1;
  }
  
  return 0;
}

int Solver::resetRHSCounter() {
  int t = rhs_ncalls;
  rhs_ncalls = 0;
  return t;
}

int Solver::resetRHSCounter_i() {
  int t = rhs_ncalls_i;
  rhs_ncalls_i = 0;
  return t;
}

int Solver::resetRHSCounter_e() {
  int t = rhs_ncalls_e;
  rhs_ncalls_e = 0;
  return t;
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
  
  int local_N = size(mesh->getRegion2D("RGN_NOBNDRY")) * (n2d + mesh->LocalNz*n3d);
  
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
  
  // Add the points which will be evolved in the boundaries
  local_N += size(mesh->getRegion2D("RGN_BNDRY")) * n2dbndry
      + size(mesh->getRegion3D("RGN_BNDRY")) * n3dbndry;
  
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

/// Perform an operation at a given Ind2D (jx,jy) location, moving data between BOUT++ and CVODE
void Solver::loop_vars_op(Ind2D i2d, BoutReal *udata, int &p, SOLVER_VAR_OP op, bool bndry) {
  int nz = mesh->LocalNz;
  
  switch(op) {
  case LOAD_VARS: {
    /// Load variables from IDA into BOUT++
    
    // Loop over 2D variables
    for(const auto& f : f2d) {
      if(bndry && !f.evolve_bndry)
        continue;
      (*f.var)[i2d] = udata[p];
      p++;
    }
    
    for (int jz=0; jz < nz; jz++) {
      
      // Loop over 3D variables
      for(const auto& f : f3d) {
        if(bndry && !f.evolve_bndry)
          continue;
        (*f.var)[f.var->getMesh()->ind2Dto3D(i2d, jz)] = udata[p];
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
      (*f.F_var)[i2d] = udata[p];
      p++;
    }
    
    for (int jz=0; jz < nz; jz++) {
      
      // Loop over 3D variables
      for(const auto& f : f3d) {
        if(bndry && !f.evolve_bndry)
          continue;
        (*f.F_var)[f.F_var->getMesh()->ind2Dto3D(i2d, jz)] = udata[p];
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
    
    for (int jz=0; jz < nz; jz++) {
      
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
      udata[p] = (*f.var)[i2d];
      p++;
    }
    
    for (int jz=0; jz < nz; jz++) {
      
      // Loop over 3D variables
      for(const auto& f : f3d) {
        if(bndry && !f.evolve_bndry)
          continue;
        udata[p] = (*f.var)[f.var->getMesh()->ind2Dto3D(i2d, jz)];
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
      udata[p] = (*f.F_var)[i2d];
      p++;
    }
    
    for (int jz=0; jz < nz; jz++) {
      
      // Loop over 3D variables
      for(const auto& f : f3d) {
        if(bndry && !f.evolve_bndry)
          continue;
        udata[p] = (*f.F_var)[f.F_var->getMesh()->ind2Dto3D(i2d, jz)];
        p++;
      }
    }
    break;
  }
  }
}

/// Loop over variables and domain. Used for all data operations for consistency
void Solver::loop_vars(BoutReal *udata, SOLVER_VAR_OP op) {
  int p = 0; // Counter for location in udata array
  
  // All boundaries
  for(const auto &i2d : mesh->getRegion2D("RGN_BNDRY")) {
    loop_vars_op(i2d, udata, p, op, true);
  }
  
  // Bulk of points
  for(const auto &i2d : mesh->getRegion2D("RGN_NOBNDRY")) {
    loop_vars_op(i2d, udata, p, op, false);
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
      throw BoutException(_("Variable '%s' not initialised"), f.name.c_str());

  for(const auto& f : f3d) 
    if(!f.var->isAllocated())
      throw BoutException(_("Variable '%s' not initialised"), f.name.c_str());
  
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
    if(f.var->getLocation() != (f.F_var)->getLocation()) {
      throw BoutException(_("Time derivative at wrong location - Field is at %s, derivative is at %s for field '%s'\n"),strLocation(f.var->getLocation()), strLocation(f.F_var->getLocation()),f.name.c_str());
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
  Field3D index(-1, mesh); // Set to -1, indicating out of domain

  int n2d = f2d.size();
  int n3d = f3d.size();

  int ind = localStart;

  int nz = mesh->LocalNz;
  
  // Find how many boundary cells are evolving
  int n2dbndry = 0;
  for (const auto &f : f2d) {
    if (f.evolve_bndry)
      ++n2dbndry;
  }
  int n3dbndry = 0;
  for (const auto &f : f3d) {
    if (f.evolve_bndry)
      ++n3dbndry;
  }
  
  if (n2dbndry + n3dbndry > 0) {
    // Some boundary points evolving

    for (const auto &i2d : mesh->getRegion2D("RGN_BNDRY")) {
      // Zero index contains 2D and 3D variables
      index[mesh->ind2Dto3D(i2d, 0)] = ind;
      ind += n2dbndry + n3dbndry;

      for (int jz = 1; jz < nz; jz++) {
        index[mesh->ind2Dto3D(i2d, jz)] = ind;
        ind += n3dbndry;
      }
    }
  }

  // Bulk of points
  for (const auto &i2d : mesh->getRegion2D("RGN_NOBNDRY")) {
    // Zero index contains 2D and 3D variables
    index[mesh->ind2Dto3D(i2d, 0)] = ind;
    ind += n2d + n3d;

    for (int jz = 1; jz < nz; jz++) {
      index[mesh->ind2Dto3D(i2d, jz)] = ind;
      ind += n3d;
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

void Solver::post_rhs(BoutReal UNUSED(t)) {
#if CHECK > 0
  for(const auto& f : f3d) {
    if(!f.F_var->isAllocated())
      throw BoutException(_("Time derivative for variable '%s' not set"), f.name.c_str());
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
    ASSERT1(f.var->getLocation() == f.F_var->getLocation());
    ASSERT1(f.var->getMesh() == f.F_var->getMesh());
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
  {
    TRACE("Solver checking time derivatives");
    for(const auto& f : f3d) {
      TRACE("Variable: %s", f.name.c_str());
      checkData(*f.F_var);
    }
  }
#endif
}

bool Solver::varAdded(const std::string &name) {
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

  return prefunc != nullptr;
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
