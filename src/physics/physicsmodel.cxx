/**************************************************************************
 * Base class for Physics Models
 *
 * Changelog:
 * 
 * 2013-08 Ben Dudson <benjamin.dudson@york.ac.uk>
 *    * Initial version
 * 
 **************************************************************************
 * Copyright 2013 B.D.Dudson
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

#define BOUT_NO_USING_NAMESPACE_BOUTGLOBALS
#include <bout/physicsmodel.hxx>
#undef BOUT_NO_USING_NAMESPACE_BOUTGLOBALS

#include <bout/mesh.hxx>

#include <fmt/core.h>

#include <string>

namespace bout {
/// Name of the directory for restart files
std::string getRestartDirectoryName(Options& options) {
  if (options["restartdir"].isSet()) {
    // Solver-specific restart directory
    return options["restartdir"].withDefault<std::string>("data");
  }
  // Use the root data directory
  return options["datadir"].withDefault<std::string>("data");
}

std::string getRestartFilename(Options& options, int rank) {
  return fmt::format("{}/BOUT.restart.{}.nc", bout::getRestartDirectoryName(options),
                     rank);
}

std::string getOutputFilename(Options& options, int rank) {
  return fmt::format("{}/BOUT.dmp.experimental.{}.nc",
                     options["datadir"].withDefault<std::string>("data"), rank);
}
} // namespace bout

PhysicsModel::PhysicsModel()
    : mesh(bout::globals::mesh), dump(bout::globals::dump),
      restart_file(bout::getRestartFilename(Options::root(), BoutComm::rank())),
      output_file(bout::getOutputFilename(Options::root(), BoutComm::rank()),
                  Options::root()["append"]
                          .doc("Add output data to existing (dump) files?")
                          .withDefault(false)
                      ? bout::OptionsNetCDF::FileMode::append
                      : bout::OptionsNetCDF::FileMode::replace),
      output_enabled(Options::root()["output"]["enabled"]
                                             .doc("Write output files")
                                             .withDefault(true)),
      restart_enabled(Options::root()["restart_files"]["enabled"]
                          .doc("Write restart files")
                          .withDefault(true)) {}

void PhysicsModel::initialise(Solver* s) {
  if (initialised) {
    return; // Ignore second initialisation
  }
  initialised = true;

  // Set the protected variable, so user code can
  // call the solver functions
  solver = s;

  bout::experimental::addBuildFlagsToOptions(output_options);
  mesh->outputVars(output_options);

  // Restart option
  const bool restarting = Options::root()["restart"].withDefault(false);

  if (restarting) {
    restart_options = restart_file.read();
  }

  // Call user init code to specify evolving variables
  if (init(restarting)) {
    throw BoutException("Couldn't initialise physics model");
  }

  // Post-initialise, which reads restart files
  // This function can be overridden by the user
  if (postInit(restarting)) {
    throw BoutException("Couldn't restart physics model");
  }
}

int PhysicsModel::runRHS(BoutReal time) {
  return rhs(time);
}

bool PhysicsModel::splitOperator() {
  return splitop;
}

int PhysicsModel::runConvective(BoutReal time) {
  return convective(time);
}

int PhysicsModel::runDiffusive(BoutReal time, bool linear) {
  return diffusive(time, linear);
}

bool PhysicsModel::hasPrecon() { return (userprecon != nullptr); }

int PhysicsModel::runPrecon(BoutReal t, BoutReal gamma, BoutReal delta) {
  if(!userprecon)
    return 1;
  return (*this.*userprecon)(t, gamma, delta);
}

bool PhysicsModel::hasJacobian() { return (userjacobian != nullptr); }

int PhysicsModel::runJacobian(BoutReal t) {
  if (!userjacobian)
    return 1;
  return (*this.*userjacobian)(t);
}

void PhysicsModel::bout_solve(Field2D &var, const char *name,
                              const std::string& description) {
  // Add to solver
  solver->add(var, name, description);
}

void PhysicsModel::bout_solve(Field3D &var, const char *name,
                              const std::string& description) {
  solver->add(var, name, description);
}

void PhysicsModel::bout_solve(Vector2D &var, const char *name,
                              const std::string& description) {
  solver->add(var, name, description);
}

void PhysicsModel::bout_solve(Vector3D &var, const char *name,
                              const std::string& description) {
  solver->add(var, name, description);
}

int PhysicsModel::postInit(bool restarting) {
  TRACE("PhysicsModel::postInit");

  if (restarting) {
    solver->readEvolvingVariablesFromOptions(restart_options);
  }

  if (restart_enabled) {
    solver->outputVars(restart_options, false);
    bout::globals::mesh->outputVars(restart_options);

    restart_options["BOUT_VERSION"].force(bout::version::as_double, "PhysicsModel");

    // Write _everything_ to restart file
    restart_file.write(restart_options);
  }

  // Add monitor to the solver which calls restart.write() and
  // PhysicsModel::outputMonitor()
  solver->addMonitor(&modelMonitor);

  return 0;
}

void PhysicsModel::writeRestartFile() {
  if (restart_enabled) {
    restart_file.write(restart_options);
  }
}

void PhysicsModel::writeOutputFile() {
  writeOutputFile(output_options);
}

void PhysicsModel::writeOutputFile(const Options& options) {
  if (output_enabled) {
    output_file.write(options);
  }
}

int PhysicsModel::PhysicsModelMonitor::call(Solver* solver, BoutReal simtime,
                                            int iteration, int nout) {
  // Restart file variables
  solver->outputVars(model->restart_options, false);
  model->writeRestartFile();

  // Main output file variables
  // t_array for backwards compatibility? needed?
  model->output_options["t_array"].force(simtime);
  model->output_options["t_array"].attributes["time_dimension"] = "t";

  model->output_options["t"].force(simtime);
  model->output_options["t"].attributes["time_dimension"] = "t";

  model->output_options["iteration"].force(iteration);
  model->output_options["iteration"].attributes["time_dimension"] = "t";

  solver->outputVars(model->output_options, true);
  model->writeOutputFile();

  // Call user output monitor
  return model->outputMonitor(simtime, iteration, nout);
}
