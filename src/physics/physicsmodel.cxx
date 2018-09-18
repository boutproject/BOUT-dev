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

#include <bout/physicsmodel.hxx>

PhysicsModel::PhysicsModel()
    : solver(nullptr), modelMonitor(this), splitop(false), userprecon(nullptr),
      userjacobian(nullptr), initialised(false) {

  // Set up restart file
  restart = Datafile(Options::getRoot()->getSection("restart"));
}

PhysicsModel::~PhysicsModel() {
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

void PhysicsModel::bout_solve(Field2D &var, const char *name) {
  // Add to solver
  solver->add(var, name);
}

void PhysicsModel::bout_solve(Field3D &var, const char *name) {
  solver->add(var, name);
}

void PhysicsModel::bout_solve(Vector2D &var, const char *name) {
  solver->add(var, name);
}

void PhysicsModel::bout_solve(Vector3D &var, const char *name) {
  solver->add(var, name);
}

int PhysicsModel::postInit(bool restarting) {
  TRACE("PhysicsModel::postInit");
  
  // Add the solver variables to the restart file
  // Second argument specifies no time history
  solver->outputVars(restart, false);

  string restart_dir;  ///< Directory for restart files
  string dump_ext, restart_ext;  ///< Dump, Restart file extension
  
  Options *options = Options::getRoot();
  if (options->isSet("restartdir")) {
    // Solver-specific restart directory
    options->get("restartdir", restart_dir, "data");
  } else {
    // Use the root data directory
    options->get("datadir", restart_dir, "data");
  }
  /// Get restart file extension
  options->get("dump_format", dump_ext, "nc");
  options->get("restart_format", restart_ext, dump_ext);

  string filename = restart_dir + "/BOUT.restart."+restart_ext;
  if (restarting) {
    output.write("Loading restart file: %s\n", filename.c_str());

    /// Load restart file
    if (!restart.openr(filename.c_str()))
      throw BoutException("Error: Could not open restart file\n");
    if (!restart.read())
      throw BoutException("Error: Could not read restart file\n");
    restart.close();
  }

  // Add mesh information to restart file
  // Note this is done after reading, so mesh variables
  // are not overwritten.
  mesh->outputVars(restart);
  // Version expected by collect routine
  restart.addOnce(const_cast<BoutReal &>(BOUT_VERSION), "BOUT_VERSION");

  /// Open the restart file for writing
  if (!restart.openw(filename.c_str()))
    throw BoutException("Error: Could not open restart file for writing\n");

  // Add monitor to the solver which calls restart.write() and
  // PhysicsModel::outputMonitor()
  solver->addMonitor(&modelMonitor);

  return 0;
}
