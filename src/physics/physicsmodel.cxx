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
} // namespace bout

PhysicsModel::PhysicsModel()
    : mesh(bout::globals::mesh), dump(bout::globals::dump),
      restart_file(bout::getRestartFilename(Options::root(), BoutComm::rank()),
                   bout::experimental::OptionsNetCDF::FileMode::append),
      modelMonitor(this) {}

void PhysicsModel::initialise(Solver* s) {
  if (initialised) {
    return; // Ignore second initialisation
  }
  initialised = true;

  // Set the protected variable, so user code can
  // call the solver functions
  solver = s;

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

  using namespace bout::experimental;

  if (restarting) {
    solver->readEvolvingVariablesFromOptions(restart_options);
  }

  const bool restart_enabled = Options::root()["restart"]["enabled"].withDefault(true);

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

int PhysicsModel::PhysicsModelMonitor::call(Solver* solver, BoutReal simtime, int iter,
                                            int nout) {
  solver->outputVars(model->restart_options, false);
  model->restart_file.write(model->restart_options);

  // Call user output monitor
  return model->outputMonitor(simtime, iter, nout);
}
