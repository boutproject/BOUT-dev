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

#include "bout/build_config.hxx"

#include "bout/array.hxx"
#include "bout/assert.hxx"
#include "bout/boutcomm.hxx"
#include "bout/boutexception.hxx"
#include "bout/field_factory.hxx"
#include "bout/initialprofiles.hxx"
#include "bout/interpolation.hxx"
#include "bout/msg_stack.hxx"
#include "bout/output.hxx"
#include "bout/region.hxx"
#include "bout/solver.hxx"
#include "bout/sys/timer.hxx"
#include "bout/sys/uuid.h"

#include <cmath>
#include <cstring>
#include <ctime>
#include <numeric>
#include <set>

// Implementations:
#include "impls/adams_bashforth/adams_bashforth.hxx"
#include "impls/arkode/arkode.hxx"
#include "impls/cvode/cvode.hxx"
#include "impls/euler/euler.hxx"
#include "impls/ida/ida.hxx"
#include "impls/imex-bdf2/imex-bdf2.hxx"
#include "impls/petsc/petsc.hxx"
#include "impls/power/power.hxx"
#include "impls/pvode/pvode.hxx"
#include "impls/rk3-ssp/rk3-ssp.hxx"
#include "impls/rk4/rk4.hxx"
#include "impls/rkgeneric/rkgeneric.hxx"
#include "impls/slepc/slepc.hxx"
#include "impls/snes/snes.hxx"
#include "impls/split-rk/split-rk.hxx"

// Static member variables

int* Solver::pargc = nullptr;
char*** Solver::pargv = nullptr;

/**************************************************************************
 * Constructor
 **************************************************************************/

Solver::Solver(Options* opts)
    : options(opts == nullptr ? &Options::root()["solver"] : opts),
      NPES(BoutComm::size()), MYPE(BoutComm::rank()),
      monitor_timestep((*options)["monitor_timestep"]
                           .doc("Call monitors on internal timesteps")
                           .withDefault(false)),
      save_repeat_run_id((*options)["save_repeat_run_id"]
                             .doc("Write run_id and run_restart_from at every output "
                                  "timestep, to make it easier to concatenate output "
                                  "data sets in time")
                             .withDefault(false)),
      is_nonsplit_model_diffusive(
          (*options)["is_nonsplit_model_diffusive"]
              .doc("If not a split operator, treat RHS as diffusive?")
              .withDefault(true)),
      mms((*options)["mms"]
              .doc("Use Method of Manufactured Solutions to track error scaling")
              .withDefault(false)),
      mms_initialise((*options)["mms_initialise"]
                         .doc("Use MMS solution for field initial conditions")
                         .withDefault(mms)),
      number_output_steps(
          (*options)["nout"]
              .doc("Number of output steps. Overrides global setting.")
              .withDefault(
                  Options::root()["nout"].doc("Number of output steps").withDefault(1))),
      output_timestep(
          (*options)["output_step"]
              .doc("Output time step size. Overrides global 'timestep' setting.")
              .withDefault(Options::root()["timestep"]
                               .doc("Output time step size")
                               .withDefault(1.0))) {}

/**************************************************************************
 * Add physics models
 **************************************************************************/

void Solver::setModel(PhysicsModel* m) {
  if (model) {
    throw BoutException("Solver can only evolve one model");
  }

  if (initialised) {
    throw BoutException("Solver already initialised");
  }

  // Initialise them model, which specifies which variables to evolve
  m->initialise(this);

  model = m;
}

/**************************************************************************
 * Add fields
 **************************************************************************/

void Solver::add(Field2D& v, const std::string& name, const std::string& description) {
  TRACE("Adding 2D field: Solver::add({:s})", name);

#if CHECK > 0
  if (varAdded(name)) {
    throw BoutException("Variable '{:s}' already added to Solver", name);
  }
#endif

  if (initialised) {
    throw BoutException("Error: Cannot add to solver after initialisation\n");
  }

  // Set boundary conditions
  v.setBoundary(name);
  ddt(v).copyBoundary(v); // Set boundary to be the same as v

  VarStr<Field2D> d;

  d.var = &v;
  d.F_var = &ddt(v);
  d.location = v.getLocation();
  d.name = name;
  d.description = description;

#if BOUT_USE_TRACK
  v.name = name;
#endif

  /// Generate initial perturbations.
  /// NOTE: This could be done in init, but this would prevent the user
  ///       from modifying the initial perturbation (e.g. to prevent unphysical situations)
  ///       before it's loaded into the solver. If restarting, this perturbation
  ///       will be over-written anyway
  if (mms_initialise) {
    // Load solution at t = 0

    FieldFactory* fact = FieldFactory::get();

    v = fact->create2D("solution", Options::getRoot()->getSection(name), v.getMesh());
  } else {
    initial_profile(name, v);
  }

  if (mms) {
    // Allocate storage for error variable
    d.MMS_err = bout::utils::make_unique<Field2D>(zeroFrom(v));
  }

  // Check if the boundary regions should be evolved
  // First get option from section "all"
  // then use that as default for specific section
  d.evolve_bndry = Options::root()["all"]["evolve_bndry"].withDefault(false);
  d.evolve_bndry = Options::root()[name]["evolve_bndry"].withDefault(d.evolve_bndry);

  v.applyBoundary(true);

  f2d.emplace_back(std::move(d));
}

void Solver::add(Field3D& v, const std::string& name, const std::string& description) {
  TRACE("Adding 3D field: Solver::add({:s})", name);

  Mesh* mesh = v.getMesh();

#if CHECK > 0
  if (varAdded(name)) {
    throw BoutException("Variable '{:s}' already added to Solver", name);
  }
#endif

  if (initialised) {
    throw BoutException("Error: Cannot add to solver after initialisation\n");
  }

  // Set boundary conditions
  v.setBoundary(name);
  ddt(v).copyBoundary(v); // Set boundary to be the same as v

  if (mesh->StaggerGrids && (v.getLocation() != CELL_CENTRE)) {
    output_info.write("\tVariable {:s} shifted to {:s}\n", name,
                      toString(v.getLocation()));
    ddt(v).setLocation(v.getLocation()); // Make sure both at the same location
  }

  VarStr<Field3D> d;

  d.var = &v;
  d.F_var = &ddt(v);
  d.location = v.getLocation();
  d.name = name;
  d.description = description;

#if BOUT_USE_TRACK
  v.name = name;
#endif

  if (mms_initialise) {
    // Load solution at t = 0
    FieldFactory* fact = FieldFactory::get();

    v = fact->create3D("solution", &Options::root()[name], mesh, v.getLocation());

  } else {
    initial_profile(name, v);
  }

  if (mms) {
    d.MMS_err = bout::utils::make_unique<Field3D>(zeroFrom(v));
  }

  // Check if the boundary regions should be evolved
  // First get option from section "all"
  // then use that as default for specific section
  d.evolve_bndry = Options::root()["all"]["evolve_bndry"].withDefault(false);
  d.evolve_bndry = Options::root()[name]["evolve_bndry"].withDefault(d.evolve_bndry);

  v.applyBoundary(true);     // Make sure initial profile obeys boundary conditions
  v.setLocation(d.location); // Restore location if changed

  f3d.emplace_back(std::move(d));
}

void Solver::add(Vector2D& v, const std::string& name, const std::string& description) {
  TRACE("Adding 2D vector: Solver::add({:s})", name);

  if (varAdded(name)) {
    throw BoutException("Variable '{:s}' already added to Solver", name);
  }

  if (initialised) {
    throw BoutException("Error: Cannot add to solver after initialisation\n");
  }

  // Set boundary conditions
  v.setBoundary(name);
  ddt(v).copyBoundary(v); // Set boundary to be the same as v

  VarStr<Vector2D> d;

  d.var = &v;
  d.F_var = &ddt(v);
  d.covariant = v.covariant;
  d.name = name;
  d.description = description;

  /// NOTE: No initial_profile call, because this will be done for each
  ///       component individually.

  /// Add suffix, depending on co- /contravariance
  if (v.covariant) {
    add(v.x, d.name + "_x");
    add(v.y, d.name + "_y");
    add(v.z, d.name + "_z");
  } else {
    add(v.x, d.name + "x");
    add(v.y, d.name + "y");
    add(v.z, d.name + "z");
  }

  /// Make sure initial profile obeys boundary conditions
  v.applyBoundary(true);
  v2d.emplace_back(std::move(d));
}

void Solver::add(Vector3D& v, const std::string& name, const std::string& description) {
  TRACE("Adding 3D vector: Solver::add({:s})", name);

  if (varAdded(name)) {
    throw BoutException("Variable '{:s}' already added to Solver", name);
  }

  if (initialised) {
    throw BoutException("Error: Cannot add to solver after initialisation\n");
  }

  // Set boundary conditions
  v.setBoundary(name);
  ddt(v).copyBoundary(v); // Set boundary to be the same as v

  VarStr<Vector3D> d;

  d.var = &v;
  d.F_var = &ddt(v);
  d.covariant = v.covariant;
  d.name = name;
  d.description = description;

  // Add suffix, depending on co- /contravariance
  if (v.covariant) {
    add(v.x, d.name + "_x");
    add(v.y, d.name + "_y");
    add(v.z, d.name + "_z");
  } else {
    add(v.x, d.name + "x");
    add(v.y, d.name + "y");
    add(v.z, d.name + "z");
  }

  v.applyBoundary(true);
  v3d.emplace_back(std::move(d));
}

/**************************************************************************
 * Constraints
 **************************************************************************/

void Solver::constraint(Field2D& v, Field2D& C_v, std::string name) {
  TRACE("Constrain 2D scalar: Solver::constraint({:s})", name);

  if (name.empty()) {
    throw BoutException("ERROR: Constraint requested for variable with empty name\n");
  }

#if CHECK > 0
  if (varAdded(name)) {
    throw BoutException("Variable '{:s}' already added to Solver", name);
  }
#endif

  if (!has_constraints) {
    throw BoutException("ERROR: This solver doesn't support constraints\n");
  }

  if (initialised) {
    throw BoutException("Error: Cannot add constraints to solver after initialisation\n");
  }

  VarStr<Field2D> d;

  d.constraint = true;
  d.var = &v;
  d.F_var = &C_v;
  d.name = std::move(name);

  f2d.emplace_back(std::move(d));
}

void Solver::constraint(Field3D& v, Field3D& C_v, std::string name) {
  TRACE("Constrain 3D scalar: Solver::constraint({:s})", name);

  if (name.empty()) {
    throw BoutException("ERROR: Constraint requested for variable with empty name\n");
  }

#if CHECK > 0
  if (varAdded(name)) {
    throw BoutException("Variable '{:s}' already added to Solver", name);
  }
#endif

  if (!has_constraints) {
    throw BoutException("ERROR: This solver doesn't support constraints\n");
  }

  if (initialised) {
    throw BoutException("Error: Cannot add constraints to solver after initialisation\n");
  }

  VarStr<Field3D> d;

  d.constraint = true;
  d.var = &v;
  d.F_var = &C_v;
  d.location = v.getLocation();
  d.name = std::move(name);

  f3d.emplace_back(std::move(d));
}

void Solver::constraint(Vector2D& v, Vector2D& C_v, std::string name) {
  TRACE("Constrain 2D vector: Solver::constraint({:s})", name);

  if (name.empty()) {
    throw BoutException("ERROR: Constraint requested for variable with empty name\n");
  }

#if CHECK > 0
  if (varAdded(name)) {
    throw BoutException("Variable '{:s}' already added to Solver", name);
  }
#endif

  if (!has_constraints) {
    throw BoutException("ERROR: This solver doesn't support constraints\n");
  }

  if (initialised) {
    throw BoutException("Error: Cannot add constraints to solver after initialisation\n");
  }

  // Add suffix, depending on co- /contravariance
  if (v.covariant) {
    constraint(v.x, C_v.x, name + "_x");
    constraint(v.y, C_v.y, name + "_y");
    constraint(v.z, C_v.z, name + "_z");
  } else {
    constraint(v.x, C_v.x, name + "x");
    constraint(v.y, C_v.y, name + "y");
    constraint(v.z, C_v.z, name + "z");
  }

  VarStr<Vector2D> d;

  d.constraint = true;
  d.var = &v;
  d.F_var = &C_v;
  d.covariant = v.covariant;
  d.name = std::move(name);

  v2d.emplace_back(std::move(d));
}

void Solver::constraint(Vector3D& v, Vector3D& C_v, std::string name) {
  TRACE("Constrain 3D vector: Solver::constraint({:s})", name);

  if (name.empty()) {
    throw BoutException("ERROR: Constraint requested for variable with empty name\n");
  }

#if CHECK > 0
  if (varAdded(name)) {
    throw BoutException("Variable '{:s}' already added to Solver", name);
  }
#endif

  if (!has_constraints) {
    throw BoutException("ERROR: This solver doesn't support constraints\n");
  }

  if (initialised) {
    throw BoutException("Error: Cannot add constraints to solver after initialisation\n");
  }

  // Add suffix, depending on co- /contravariance
  if (v.covariant) {
    constraint(v.x, C_v.x, name + "_x");
    constraint(v.y, C_v.y, name + "_y");
    constraint(v.z, C_v.z, name + "_z");
  } else {
    constraint(v.x, C_v.x, name + "x");
    constraint(v.y, C_v.y, name + "y");
    constraint(v.z, C_v.z, name + "z");
  }

  VarStr<Vector3D> d;

  d.constraint = true;
  d.var = &v;
  d.F_var = &C_v;
  d.covariant = v.covariant;
  d.name = std::move(name);

  v3d.emplace_back(std::move(d));
}

/**************************************************************************
 * Solver main loop: Initialise, run, and finish
 **************************************************************************/

int Solver::solve(int nout, BoutReal timestep) {

  Options& globaloptions = Options::root(); // Default from global options

  if (nout < 0) {
    nout = number_output_steps;
    timestep = output_timestep;
  } else {
    number_output_steps = nout;
    output_timestep = timestep;
  }

  finaliseMonitorPeriods(nout, timestep);

  output_progress.write(
      _("Solver running for {:d} outputs with output timestep of {:e}\n"), nout,
      timestep);
  if (default_monitor_period > 1) {
    output_progress.write(
        _("Solver running for {:d} outputs with monitor timestep of {:e}\n"),
        nout / default_monitor_period, timestep * default_monitor_period);
  }

  // Initialise
  if (init()) {
    throw BoutException(_("Failed to initialise solver-> Aborting\n"));
  }

  // Set the run ID
  run_restart_from = run_id; // Restarting from the previous run ID
  run_id = createRunID();

  // Put the run ID into the options tree
  // Forcing in case the value has been previously set
  Options::root()["run"]["run_id"].force(run_id, "Output");
  Options::root()["run"]["run_restart_from"].force(run_restart_from, "Output");

  /// Run the solver
  output_info.write(_("Running simulation\n\n"));
  output_info.write("Run ID: {:s}\n", run_id);
  if (run_restart_from != default_run_id) {
    output_info.write("Restarting from ID: {:s}\n", run_restart_from);
  }

  time_t start_time = time(nullptr);
  output_progress.write(_("\nRun started at  : {:s}\n"), toString(start_time));

  Timer timer("run"); // Start timer

  const bool restart =
      globaloptions["restart"].doc("Load state from restart files?").withDefault(false);

  const bool append =
      globaloptions["append"]
          .doc("Add new outputs to the end of existing files? If false, overwrite files.")
          .withDefault(false);
  const bool dump_on_restart = globaloptions["dump_on_restart"]
                                   .doc("Write initial state as time point 0?")
                                   .withDefault(!restart || !append);

  // Run RHS once to ensure all variables set
  if (run_rhs(simtime)) {
    throw BoutException("Physics RHS call failed\n");
  }

  // Check for unused/mistyped options
  const bool validate_input = globaloptions["input"]["validate"]
                                  .doc("Check for unused options and stop")
                                  .withDefault(false);
  bout::checkForUnusedOptions();
  if (validate_input) {
    return 0;
  }

  if (dump_on_restart) {

    /// Write initial state as time-point 0

    // Call monitors so initial values are written to output dump files
    if (call_monitors(simtime, -1, nout)) {
      throw BoutException("Initial monitor call failed!");
    }

    // Reset iteration counter to undo the increment from the initial call_monitors().
    // That call was either at t=0, or a repeat of the last output before restarting.
    resetIterationCounter(getIterationCounter() - 1);
  }

  int status;
  try {
    status = run();

    time_t end_time = time(nullptr);
    output_progress.write(_("\nRun finished at  : {:s}\n"), toString(end_time));
    output_progress.write(_("Run time : "));

    int dt = end_time - start_time;
    int i = static_cast<int>(dt / (60. * 60.));
    if (i > 0) {
      output_progress.write("{:d} h ", i);
      dt -= i * 60 * 60;
    }
    i = static_cast<int>(dt / 60.);
    if (i > 0) {
      output_progress.write("{:d} m ", i);
      dt -= i * 60;
    }
    output_progress.write("{:d} s\n", dt);
  } catch (BoutException& e) {
    output_error << "Error encountered in solver run\n";
    output_error << e.what() << endl;
    throw;
  }

  return status;
}

std::string Solver::createRunID() const {

  std::string result;
  result.resize(36);

  if (MYPE == 0) {
    // Generate a unique ID for this run
#if BOUT_HAS_UUID_SYSTEM_GENERATOR
    uuids::uuid_system_generator gen{};
#else
    std::random_device rd;
    auto seed_data = std::array<int, std::mt19937::state_size>{};
    std::generate(std::begin(seed_data), std::end(seed_data), std::ref(rd));
    std::seed_seq seq(std::begin(seed_data), std::end(seed_data));
    std::mt19937 generator(seq);
    uuids::uuid_random_generator gen{generator};
#endif

    result = uuids::to_string(gen());
  }

  // All ranks have same run_id
  // Standard representation of UUID is always 36 characters
  MPI_Bcast(&result[0], 36, MPI_CHAR, 0, BoutComm::get());

  return result;
}

std::string Solver::getRunID() const {
  AUTO_TRACE();
  if (run_id == default_run_id) {
    throw BoutException("run_id not set!");
  }
  return run_id;
}

std::string Solver::getRunRestartFrom() const {
  AUTO_TRACE();
  // Check against run_id, because this might not be a restarted run
  if (run_id == default_run_id) {
    throw BoutException("run_restart_from not set!");
  }
  return run_restart_from;
}

void Solver::writeToModelOutputFile(const Options& options) {
  if (model == nullptr) {
    return;
  }
  model->writeOutputFile(options);
}

/**************************************************************************
 * Initialisation
 **************************************************************************/

int Solver::init() {

  TRACE("Solver::init()");

  if (initialised) {
    throw BoutException(_("ERROR: Solver is already initialised\n"));
  }

  output_progress.write(_("Initialising solver\n"));

  /// Mark as initialised. No more variables can be added
  initialised = true;

  return 0;
}

void Solver::outputVars(Options& output_options, bool save_repeat) {
  Timer time("io");
  output_options["tt"].force(simtime, "Solver");
  output_options["hist_hi"].force(iteration, "Solver");

  output_options["run_id"]
      .doc("UUID for this simulation")
      .assignRepeat(run_id, "t", save_repeat and save_repeat_run_id, "Solver");
  output_options["run_restart_from"]
      .doc("run_id of the simulation this one was restarted from."
           "'zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz' means the run is not a restart, "
           "or the previous run did not have a run_id.")
      .assignRepeat(run_restart_from, "t", save_repeat and save_repeat_run_id, "Solver");

  // Add 2D and 3D evolving fields to output file
  for (const auto& f : f2d) {
    // Add to dump file (appending)
    output_options[f.name].assignRepeat(*(f.var), "t", save_repeat, "Solver");
    output_options[f.name].attributes["description"] = f.description;
  }
  for (const auto& f : f3d) {
    // Add to dump file (appending)
    output_options[f.name].assignRepeat(*(f.var), "t", save_repeat, "Solver");
    output_options[f.name].attributes["description"] = f.description;
    if (mms) {
      // Add an error variable
      output_options["E_" + f.name].assignRepeat(*(f.MMS_err), "t", save_repeat,
                                                 "Solver");
      output_options["E_" + f.name].attributes["description"] = f.description;
    }
  }

  if (save_repeat) {
    // Do not save if save_repeat=false so we avoid adding diagnostic variables to restart
    // files, otherwise they might cause errors if the solver type is changed before
    // restarting

    // Add solver diagnostics to output file
    for (const auto& d : diagnostic_int) {
      output_options[d.name].assignRepeat(*(d.var), "t", true, "Solver");
      output_options[d.name].attributes["description"] = d.description;
    }
    for (const auto& d : diagnostic_BoutReal) {
      output_options[d.name].assignRepeat(*(d.var), "t", true, "Solver");
      output_options[d.name].attributes["description"] = d.description;
    }
  }
}

void Solver::readEvolvingVariablesFromOptions(Options& options) {
  run_id = options["run_id"].withDefault(default_run_id);
  simtime = options["tt"].as<BoutReal>();
  iteration = options["hist_hi"].withDefault<int>(0);

  for (auto& f : f2d) {
    if (options.isSet(f.name)) {
      *(f.var) = options[f.name].as<Field2D>();
    } else {
      output_warn.write("Restart does not contain Field2D '{}' => Initialising", f.name);
    }
  }
  for (const auto& f : f3d) {
    if (options.isSet(f.name)) {
      *(f.var) = options[f.name].as<Field3D>();
      if (mms) {
        *(f.MMS_err) = options["E_" + f.name].as<Field3D>();
      }
    } else {
      output_warn.write("Restart does not contain Field3D '{}' => Initialising", f.name);
    }
  }
}

/////////////////////////////////////////////////////

BoutReal Solver::adjustMonitorPeriods(Monitor* new_monitor) {

  if (new_monitor->timestep < 0) {
    // The timestep will get adjusted when we call solve
    new_monitor->period = default_monitor_period;
    return internal_timestep;
  }

  if (!initialised && internal_timestep < 0) {
    // This is the first monitor to be added
    return new_monitor->timestep;
  }

  if (!isMultiple(internal_timestep, new_monitor->timestep)) {
    throw BoutException(_("Couldn't add Monitor: {:g} is not a multiple of {:g}!"),
                        internal_timestep, new_monitor->timestep);
  }

  if (new_monitor->timestep > internal_timestep * 1.5) {
    // Monitor has a larger timestep
    new_monitor->period =
        static_cast<int>(std::round(new_monitor->timestep / internal_timestep));
    return internal_timestep;
  }

  // Monitor timestep is smaller, so we need to adjust our timestep,
  // along with that of all of the other monitors

  if (initialised) {
    throw BoutException(
        _("Solver::addMonitor: Cannot reduce timestep (from {:g} to {:g}) "
          "after init is called!"),
        internal_timestep, new_monitor->timestep);
  }

  // This is the relative increase in timestep
  const auto multiplier =
      static_cast<int>(std::round(internal_timestep / new_monitor->timestep));
  for (const auto& monitor : monitors) {
    monitor.monitor->period *= multiplier;
  }

  // Update default_monitor_frequency so that monitors with no
  // timestep are called at the output frequency
  default_monitor_period *= multiplier;

  // This monitor is now the fastest monitor
  return new_monitor->timestep;
}

void Solver::finaliseMonitorPeriods(int& NOUT, BoutReal& output_timestep) {
  // Synchronise timestep with those of the monitors
  if (internal_timestep > 0) {
    if (!isMultiple(internal_timestep, output_timestep)) {
      throw BoutException(
          "A monitor requested a timestep not compatible with the output_step!");
    }
    if (internal_timestep < output_timestep * 1.5) {
      default_monitor_period =
          static_cast<int>(std::round(output_timestep / internal_timestep));
      NOUT *= default_monitor_period;
      output_timestep = internal_timestep;
    } else {
      default_monitor_period = 1;
      // update old monitors
      const auto multiplier =
          static_cast<int>(std::round(internal_timestep / output_timestep));
      for (const auto& monitor : monitors) {
        monitor.monitor->period *= multiplier;
      }
    }
  }
  int count = 0;
  // Now set any monitors which still have the default
  // timestep/period, and set the time_dimension for each monitor
  for (auto& monitor : monitors) {
    if (monitor.monitor->timestep < 0) {
      monitor.monitor->timestep = internal_timestep * default_monitor_period;
      monitor.monitor->period = default_monitor_period;
      monitor.time_dimension = "t";
    } else {
      monitor.time_dimension = fmt::format("t{}", ++count);
    }
  }
}

void Solver::addMonitor(Monitor* monitor, MonitorPosition pos) {

  internal_timestep = adjustMonitorPeriods(monitor);

  monitor->is_added = true;

  if (pos == MonitorPosition::FRONT) {
    monitors.push_front({monitor, ""});
  } else {
    monitors.push_back({monitor, ""});
  }
}

void Solver::removeMonitor(Monitor* f) {
  monitors.remove_if([&f](auto& monitor) { return monitor.monitor == f; });
}

extern bool user_requested_exit;
int Solver::call_monitors(BoutReal simtime, int iter, int NOUT) {
  bool abort;
  bout::globals::mpi->MPI_Allreduce(&user_requested_exit, &abort, 1, MPI_C_BOOL, MPI_LOR,
                                    BoutComm::get());
  if (abort) {
    NOUT = iter + 1;
  }
  if (mms) {
    // Calculate MMS errors
    calculate_mms_error(simtime);
  }

  ++iter;
  try {
    // We need to write each time dimension a maximum of once per
    // timestep. The set of unique time dimensions may be the same
    // size or smaller than the set of monitors, so we need to keep
    // track of the unique dimensions each timestep.
    std::set<std::string> seen_time_dimensions;

    // Call monitors
    for (const auto& monitor : monitors) {
      if ((iter % monitor.monitor->period) == 0) {
        // Call each monitor one by one
        const int ret =
            monitor.monitor->call(this, simtime, iter / monitor.monitor->period - 1,
                                  NOUT / monitor.monitor->period);
        if (ret != 0) {
          throw BoutException(_("Monitor signalled to quit (return code {})"), ret);
        }
        // Write the monitor's diagnostics to the main output file
        Options monitor_dump;
        monitor.monitor->outputVars(monitor_dump, monitor.time_dimension);
        model->writeOutputFile(monitor_dump, monitor.time_dimension);
        // This monitor's time dimension needs writing out
        seen_time_dimensions.insert(monitor.time_dimension);
      }
    }
    // Write all the unique time dimensions that were advanced this timestep
    for (const auto& time_dimension : seen_time_dimensions) {
      Options time_dump;
      time_dump[time_dimension].assignRepeat(simtime, time_dimension);
      model->writeOutputFile(time_dump, time_dimension);
    }

    model->finishOutputTimestep();
  } catch (const BoutException& e) {
    for (const auto& monitor : monitors) {
      monitor.monitor->cleanup();
    }
    output_error.write(_("Monitor signalled to quit (exception {})\n"), e.what());
    throw;
  }

  // Check if any of the monitors has asked to quit
  bout::globals::mpi->MPI_Allreduce(&user_requested_exit, &abort, 1, MPI_C_BOOL, MPI_LOR,
                                    BoutComm::get());

  if (iter == NOUT || abort) {
    for (const auto& monitor : monitors) {
      monitor.monitor->cleanup();
    }
  }

  if (abort) {
    // restart file should be written by physics model
    output.write(_("User signalled to quit. Returning\n"));
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

bool Solver::splitOperator() { return model->splitOperator(); }

/////////////////////////////////////////////////////

void Solver::addTimestepMonitor(TimestepMonitorFunc f) {
  timestep_monitors.push_front(f);
}

void Solver::removeTimestepMonitor(TimestepMonitorFunc f) { timestep_monitors.remove(f); }

int Solver::call_timestep_monitors(BoutReal simtime, BoutReal lastdt) {
  if (!monitor_timestep) {
    return 0;
  }

  for (const auto& monitor : timestep_monitors) {
    const int ret = monitor(this, simtime, lastdt);
    if (ret != 0) {
      return ret; // Return first time an error is encountered
    }
  }

  // Call physics model monitor
  return model->runTimestepMonitor(simtime, lastdt);
}

/**************************************************************************
 * Useful routines (protected)
 **************************************************************************/

int Solver::getLocalN() {

  // Cache the value, so this is not repeatedly called.
  // This value should not change after initialisation
  static int cacheLocalN{-1};
  if (cacheLocalN != -1) {
    return cacheLocalN;
  }

  // Must be initialised
  ASSERT0(initialised);

  // Return the number of points to evolve in f, plus the accumulator value.
  // If f.evolve_bndry, includes the boundary (NB: not guard!) points
  auto local_N_sum = [](int value, const auto& field) -> int {
    const auto boundary_size =
        field.evolve_bndry ? size(field.var->getRegion("RGN_BNDRY")) : 0;
    return value + boundary_size + size(field.var->getRegion("RGN_NOBNDRY"));
  };

  const auto local_N_2D = std::accumulate(begin(f2d), end(f2d), 0, local_N_sum);
  const auto local_N_3D = std::accumulate(begin(f3d), end(f3d), 0, local_N_sum);
  const auto local_N = local_N_2D + local_N_3D;

  cacheLocalN = local_N;

  return local_N;
}

std::unique_ptr<Solver> Solver::create(Options* opts) {
  return SolverFactory::getInstance().create(opts);
}

std::unique_ptr<Solver> Solver::create(const SolverType& type, Options* opts) {
  return SolverFactory::getInstance().create(type, opts);
}

/**************************************************************************
 * Looping over variables
 *
 * NOTE: This part is very inefficient, and should be replaced ASAP
 * Is the interleaving of variables needed or helpful to the solver?
 **************************************************************************/

/// Perform an operation at a given Ind2D (jx,jy) location, moving data between BOUT++ and CVODE
void Solver::loop_vars_op(Ind2D i2d, BoutReal* udata, int& p, SOLVER_VAR_OP op,
                          bool bndry) {
  // Use global mesh: FIX THIS!
  Mesh* mesh = bout::globals::mesh;

  int nz = mesh->LocalNz;

  switch (op) {
  case SOLVER_VAR_OP::LOAD_VARS: {
    /// Load variables from IDA into BOUT++

    // Loop over 2D variables
    for (const auto& f : f2d) {
      if (bndry && !f.evolve_bndry) {
        continue;
      }
      (*f.var)[i2d] = udata[p];
      p++;
    }

    for (int jz = 0; jz < nz; jz++) {

      // Loop over 3D variables
      for (const auto& f : f3d) {
        if (bndry && !f.evolve_bndry) {
          continue;
        }
        (*f.var)[f.var->getMesh()->ind2Dto3D(i2d, jz)] = udata[p];
        p++;
      }
    }
    break;
  }
  case SOLVER_VAR_OP::LOAD_DERIVS: {
    /// Load derivatives from IDA into BOUT++
    /// Used for preconditioner

    // Loop over 2D variables
    for (const auto& f : f2d) {
      if (bndry && !f.evolve_bndry) {
        continue;
      }
      (*f.F_var)[i2d] = udata[p];
      p++;
    }

    for (int jz = 0; jz < nz; jz++) {

      // Loop over 3D variables
      for (const auto& f : f3d) {
        if (bndry && !f.evolve_bndry) {
          continue;
        }
        (*f.F_var)[f.F_var->getMesh()->ind2Dto3D(i2d, jz)] = udata[p];
        p++;
      }
    }

    break;
  }
  case SOLVER_VAR_OP::SET_ID: {
    /// Set the type of equation (Differential or Algebraic)

    // Loop over 2D variables
    for (const auto& f : f2d) {
      if (bndry && !f.evolve_bndry) {
        continue;
      }
      if (f.constraint) {
        udata[p] = 0;
      } else {
        udata[p] = 1;
      }
      p++;
    }

    for (int jz = 0; jz < nz; jz++) {

      // Loop over 3D variables
      for (const auto& f : f3d) {
        if (bndry && !f.evolve_bndry) {
          continue;
        }
        if (f.constraint) {
          udata[p] = 0;
        } else {
          udata[p] = 1;
        }
        p++;
      }
    }

    break;
  }
  case SOLVER_VAR_OP::SAVE_VARS: {
    /// Save variables from BOUT++ into IDA (only used at start of simulation)

    // Loop over 2D variables
    for (const auto& f : f2d) {
      if (bndry && !f.evolve_bndry) {
        continue;
      }
      udata[p] = (*f.var)[i2d];
      p++;
    }

    for (int jz = 0; jz < nz; jz++) {

      // Loop over 3D variables
      for (const auto& f : f3d) {
        if (bndry && !f.evolve_bndry) {
          continue;
        }
        udata[p] = (*f.var)[f.var->getMesh()->ind2Dto3D(i2d, jz)];
        p++;
      }
    }
    break;
  }
    /// Save time-derivatives from BOUT++ into CVODE (returning RHS result)
  case SOLVER_VAR_OP::SAVE_DERIVS: {

    // Loop over 2D variables
    for (const auto& f : f2d) {
      if (bndry && !f.evolve_bndry) {
        continue;
      }
      udata[p] = (*f.F_var)[i2d];
      p++;
    }

    for (int jz = 0; jz < nz; jz++) {

      // Loop over 3D variables
      for (const auto& f : f3d) {
        if (bndry && !f.evolve_bndry) {
          continue;
        }
        udata[p] = (*f.F_var)[f.F_var->getMesh()->ind2Dto3D(i2d, jz)];
        p++;
      }
    }
    break;
  }
  }
}

/// Loop over variables and domain. Used for all data operations for consistency
void Solver::loop_vars(BoutReal* udata, SOLVER_VAR_OP op) {
  // Use global mesh: FIX THIS!
  Mesh* mesh = bout::globals::mesh;

  int p = 0; // Counter for location in udata array

  // All boundaries
  for (const auto& i2d : mesh->getRegion2D("RGN_BNDRY")) {
    loop_vars_op(i2d, udata, p, op, true);
  }

  // Bulk of points
  for (const auto& i2d : mesh->getRegion2D("RGN_NOBNDRY")) {
    loop_vars_op(i2d, udata, p, op, false);
  }
}

void Solver::load_vars(BoutReal* udata) {
  // Make sure data is allocated
  for (const auto& f : f2d) {
    f.var->allocate();
  }
  for (const auto& f : f3d) {
    f.var->allocate();
    f.var->setLocation(f.location);
  }

  loop_vars(udata, SOLVER_VAR_OP::LOAD_VARS);

  // Mark each vector as either co- or contra-variant

  for (const auto& v : v2d) {
    v.var->covariant = v.covariant;
  }
  for (const auto& v : v3d) {
    v.var->covariant = v.covariant;
  }
}

void Solver::load_derivs(BoutReal* udata) {
  // Make sure data is allocated
  for (const auto& f : f2d) {
    f.F_var->allocate();
  }
  for (const auto& f : f3d) {
    f.F_var->allocate();
    f.F_var->setLocation(f.location);
  }

  loop_vars(udata, SOLVER_VAR_OP::LOAD_DERIVS);

  // Mark each vector as either co- or contra-variant

  for (const auto& v : v2d) {
    v.F_var->covariant = v.covariant;
  }
  for (const auto& v : v3d) {
    v.F_var->covariant = v.covariant;
  }
}

// This function only called during initialisation
void Solver::save_vars(BoutReal* udata) {
  for (const auto& f : f2d) {
    if (!f.var->isAllocated()) {
      throw BoutException(_("Variable '{:s}' not initialised"), f.name);
    }
  }

  for (const auto& f : f3d) {
    if (!f.var->isAllocated()) {
      throw BoutException(_("Variable '{:s}' not initialised"), f.name);
    }
  }

  // Make sure vectors in correct basis
  for (const auto& v : v2d) {
    if (v.covariant) {
      v.var->toCovariant();
    } else {
      v.var->toContravariant();
    }
  }
  for (const auto& v : v3d) {
    if (v.covariant) {
      v.var->toCovariant();
    } else {
      v.var->toContravariant();
    }
  }

  loop_vars(udata, SOLVER_VAR_OP::SAVE_VARS);
}

void Solver::save_derivs(BoutReal* dudata) {
  // Make sure vectors in correct basis
  for (const auto& v : v2d) {
    if (v.covariant) {
      v.F_var->toCovariant();
    } else {
      v.F_var->toContravariant();
    }
  }
  for (const auto& v : v3d) {
    if (v.covariant) {
      v.F_var->toCovariant();
    } else {
      v.F_var->toContravariant();
    }
  }

  // Make sure 3D fields are at the correct cell location
  for (const auto& f : f3d) {
    if (f.var->getLocation() != (f.F_var)->getLocation()) {
      throw BoutException(_("Time derivative at wrong location - Field is at {:s}, "
                            "derivative is at {:s} for field '{:s}'\n"),
                          toString(f.var->getLocation()),
                          toString(f.F_var->getLocation()), f.name);
    }
  }

  loop_vars(dudata, SOLVER_VAR_OP::SAVE_DERIVS);
}

void Solver::set_id(BoutReal* udata) { loop_vars(udata, SOLVER_VAR_OP::SET_ID); }

Field3D Solver::globalIndex(int localStart) {
  // Use global mesh: FIX THIS!
  Mesh* mesh = bout::globals::mesh;

  Field3D index(-1, mesh); // Set to -1, indicating out of domain

  int n2d = f2d.size();
  int n3d = f3d.size();

  int ind = localStart;

  int nz = mesh->LocalNz;

  // Find how many boundary cells are evolving
  int n2dbndry = 0;
  for (const auto& f : f2d) {
    if (f.evolve_bndry) {
      ++n2dbndry;
    }
  }
  int n3dbndry = 0;
  for (const auto& f : f3d) {
    if (f.evolve_bndry) {
      ++n3dbndry;
    }
  }

  if (n2dbndry + n3dbndry > 0) {
    // Some boundary points evolving

    for (const auto& i2d : mesh->getRegion2D("RGN_BNDRY")) {
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
  for (const auto& i2d : mesh->getRegion2D("RGN_NOBNDRY")) {
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

int Solver::run_rhs(BoutReal t, bool linear) {
  int status;

  Timer timer("rhs");

  if (model->splitOperator()) {
    // Run both parts

    int nv = getLocalN();
    // Create two temporary arrays for system state
    Array<BoutReal> tmp(nv);
    Array<BoutReal> tmp2(nv);

    save_vars(tmp.begin()); // Copy variables into tmp
    pre_rhs(t);
    status = model->runConvective(t, linear);
    post_rhs(t); // Check variables, apply boundary conditions

    load_vars(tmp.begin());   // Reset variables
    save_derivs(tmp.begin()); // Save time derivatives
    pre_rhs(t);
    status = model->runDiffusive(t, linear);
    post_rhs(t);
    save_derivs(tmp2.begin()); // Save time derivatives
    for (BoutReal *t = tmp.begin(), *t2 = tmp2.begin(); t != tmp.end(); ++t, ++t2) {
      *t += *t2;
    }
    load_derivs(tmp.begin()); // Put back time-derivatives
  } else {
    pre_rhs(t);
    status = model->runRHS(t, linear);
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
int Solver::run_convective(BoutReal t, bool linear) {
  int status;

  Timer timer("rhs");
  pre_rhs(t);
  if (model->splitOperator()) {
    status = model->runConvective(t, linear);
  } else if (!is_nonsplit_model_diffusive) {
    // Return total
    status = model->runRHS(t, linear);
  } else {
    // Zero if not split
    for (const auto& f : f3d) {
      *(f.F_var) = 0.0;
    }
    for (const auto& f : f2d) {
      *(f.F_var) = 0.0;
    }
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
  if (model->splitOperator()) {

    status = model->runDiffusive(t, linear);
    post_rhs(t);
  } else if (is_nonsplit_model_diffusive) {
    // Return total
    status = model->runRHS(t, linear);
  } else {
    // Zero if not split
    for (const auto& f : f3d) {
      *(f.F_var) = 0.0;
    }
    for (const auto& f : f2d) {
      *(f.F_var) = 0.0;
    }
    status = 0;
  }
  rhs_ncalls_i++;
  return status;
}

void Solver::pre_rhs(BoutReal t) {

  // Apply boundary conditions to the values
  for (const auto& f : f2d) {
    if (!f.constraint) { // If it's not a constraint
      f.var->applyBoundary(t);
    }
  }

  for (const auto& f : f3d) {
    if (!f.constraint) {
      f.var->applyBoundary(t);
    }
  }
}

void Solver::post_rhs(BoutReal UNUSED(t)) {
#if CHECK > 0
  for (const auto& f : f3d) {
    if (!f.F_var->isAllocated()) {
      throw BoutException(_("Time derivative for variable '{:s}' not set"), f.name);
    }
  }
#endif
  // Make sure vectors in correct basis
  for (const auto& v : v2d) {
    if (v.covariant) {
      v.F_var->toCovariant();
    } else {
      v.F_var->toContravariant();
    }
  }
  for (const auto& v : v3d) {
    if (v.covariant) {
      v.F_var->toCovariant();
    } else {
      v.F_var->toContravariant();
    }
  }

  // Make sure 3D fields are at the correct cell location, etc.
  for ([[maybe_unused]] const auto& f : f3d) {
    ASSERT1_FIELDS_COMPATIBLE(*f.var, *f.F_var);
  }

  // Apply boundary conditions to the time-derivatives
  for (const auto& f : f2d) {
    if (!f.constraint
        && f.evolve_bndry) { // If it's not a constraint and if the boundary is evolving
      f.var->applyTDerivBoundary();
    }
  }

  for (const auto& f : f3d) {
    if (!f.constraint && f.evolve_bndry) {
      f.var->applyTDerivBoundary();
    }
  }
#if CHECK > 2
  {
    TRACE("Solver checking time derivatives");
    for (const auto& f : f3d) {
      TRACE("Variable: {:s}", f.name);
      checkData(*f.F_var);
    }
  }
#endif
}

bool Solver::varAdded(const std::string& name) {
  return contains(f2d, name) || contains(f3d, name) || contains(v2d, name)
         || contains(v3d, name);
}

bool Solver::hasPreconditioner() { return model->hasPrecon(); }

int Solver::runPreconditioner(BoutReal t, BoutReal gamma, BoutReal delta) {
  return model->runPrecon(t, gamma, delta);
}

bool Solver::hasJacobian() { return model->hasJacobian(); }
int Solver::runJacobian(BoutReal time) { return model->runJacobian(time); }

// Add source terms to time derivatives
void Solver::add_mms_sources(BoutReal t) {
  if (!mms) {
    return;
  }

  FieldFactory* fact = FieldFactory::get();

  // Iterate over 2D variables
  for (const auto& f : f2d) {
    *f.F_var += fact->create2D("source", Options::getRoot()->getSection(f.name),
                               f.var->getMesh(), (f.var)->getLocation(), t);
  }

  for (const auto& f : f3d) {
    *f.F_var += fact->create3D("source", Options::getRoot()->getSection(f.name),
                               f.var->getMesh(), (f.var)->getLocation(), t);
  }
}

// Calculate
void Solver::calculate_mms_error(BoutReal t) {
  FieldFactory* fact = FieldFactory::get();

  for (const auto& f : f3d) {
    Field3D solution = fact->create3D("solution", Options::getRoot()->getSection(f.name),
                                      f.var->getMesh(), (f.var)->getLocation(), t);

    *(f.MMS_err) = *(f.var) - solution;
  }
}
