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

PhysicsModel::PhysicsModel()
    : mesh(bout::globals::mesh), dump(bout::globals::dump), modelMonitor(this) {

  // Set up restart file
  restart = Datafile(Options::getRoot()->getSection("restart"));
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
  
  // Add the solver variables to the restart file
  // Second argument specifies no time history
  solver->outputVars(restart, false);

  auto& options = Options::root();

  const std::string restart_dir = options["restartdir"]
                                      .doc("Directory for restart files")
                                      .withDefault(options["datadir"]);

  const std::string restart_ext = options["restart_format"]
                                      .doc("Restart file extension")
                                      .withDefault(options["dump_format"]);

  const std::string filename = restart_dir + "/BOUT.restart." + restart_ext;
  if (restarting) {
    output.write("Loading restart file: {:s}\n", filename);

    /// Load restart file
    if (!restart.openr(filename))
      throw BoutException("Error: Could not open restart file {:s}\n", filename);
    if (!restart.read())
      throw BoutException("Error: Could not read restart file {:s}\n", filename);
    restart.close();
  }

  // Add mesh information to restart file
  // Note this is done after reading, so mesh variables
  // are not overwritten.
  bout::globals::mesh->outputVars(restart);
  // Version expected by collect routine
  restart.addOnce(const_cast<BoutReal &>(bout::version::as_double), "BOUT_VERSION");

  /// Open the restart file for writing
  if (!restart.openw(filename))
    throw BoutException("Error: Could not open restart file for writing\n");

  if (restarting) {
    // Write variables to the restart files so that the initial data is not lost if there is
    // a crash before modelMonitor is called for the first time
    if (!restart.write()) {
      throw BoutException("Error: Failed to write initial data back to restart file");
    }
  }

  // Add monitor to the solver which calls restart.write() and
  // PhysicsModel::outputMonitor()
  solver->addMonitor(&modelMonitor);

  return 0;
}
