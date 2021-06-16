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
void DataFileFacade::add(ValueType value, const std::string& name, bool save_repeat) {
  data.emplace_back(name, value, save_repeat);
}

bool DataFileFacade::write() {
  for (const auto& item : data) {
    bout::utils::visit(bout::OptionsConversionVisitor{Options::root(), item.name},
                       item.value);
    if (item.repeat) {
      Options::root()[item.name].attributes["time_dimension"] = "t";
    }
  }
  writeDefaultOutputFile();
  return true;
}
} // namespace bout

PhysicsModel::PhysicsModel()
    : mesh(bout::globals::mesh),
      output_file(bout::getOutputFilename(Options::root()),
                  Options::root()["append"]
                          .doc("Add output data to existing (dump) files?")
                          .withDefault(false)
                      ? bout::OptionsNetCDF::FileMode::append
                      : bout::OptionsNetCDF::FileMode::replace),
      output_enabled(Options::root()["output"]["enabled"]
                         .doc("Write output files")
                         .withDefault(true)),
      restart_file(bout::getRestartFilename(Options::root())),
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
    restartVars(restart_options);
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

void PhysicsModel::outputVars(Options& options) {
  for (const auto& item : dump.getData()) {
    bout::utils::visit(bout::OptionsConversionVisitor{options, item.name}, item.value);
    if (item.repeat) {
      options[item.name].attributes["time_dimension"] = "t";
    }
  }
}

void PhysicsModel::restartVars(Options& options) {
  for (const auto& item : restart.getData()) {
    bout::utils::visit(bout::OptionsConversionVisitor{options, item.name}, item.value);
    if (item.repeat) {
      options[item.name].attributes["time_dimension"] = "t";
    }
  }
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
    output_file.write(options, "t");
  }
}

void PhysicsModel::writeOutputFile(const Options& options, const std::string& time_dimension) {
  if (output_enabled) {
    output_file.write(options, time_dimension);
  }
}

void PhysicsModel::finishOutputTimestep() const {
  if (output_enabled) {
    output_file.verifyTimesteps();
  }
}

int PhysicsModel::PhysicsModelMonitor::call(Solver* solver, BoutReal simtime,
                                            int iteration, int nout) {
  // Restart file variables
  solver->outputVars(model->restart_options, false);
  model->restartVars(model->restart_options);
  model->writeRestartFile();

  // Main output file variables
  // t_array for backwards compatibility? needed?
  model->output_options["t_array"].assignRepeat(simtime);
  model->output_options["iteration"].assignRepeat(iteration);

  solver->outputVars(model->output_options, true);
  model->outputVars(model->output_options);
  model->writeOutputFile();

  // Call user output monitor
  return model->outputMonitor(simtime, iteration, nout);
}
