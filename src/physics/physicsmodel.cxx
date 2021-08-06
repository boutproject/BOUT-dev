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

PhysicsModel::PhysicsModel() : modelMonitor(this) {

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
  
  std::string restart_dir;  ///< Directory for restart files
  std::string dump_ext, restart_ext;  ///< Dump, Restart file extension
  
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

  std::string filename = restart_dir + "/BOUT.restart."+restart_ext;

  // Add the solver variables to the restart file
  // Second argument specifies no time history
  // Open and close the restart file first so that it knows it's filename - needed so
  // can_write_strings() works and we can skip writing run_id for HDF5 files.
  if (!restart.openr("%s",filename.c_str())) {
    throw BoutException("Error: Could not open restart file %s\n", filename.c_str());
  }
  restart.close();
  solver->outputVars(restart, false);

  if (restarting) {
    output.write("Loading restart file: %s\n", filename.c_str());

    /// Load restart file
    if (!restart.openr("%s",filename.c_str()))
      throw BoutException("Error: Could not open restart file %s\n", filename.c_str());
    if (!restart.read())
      throw BoutException("Error: Could not read restart file %s\n", filename.c_str());
    restart.close();
  }

  // Add mesh information to restart file
  // Note this is done after reading, so mesh variables
  // are not overwritten.
  bout::globals::mesh->outputVars(restart);
  // Version expected by collect routine
  restart.addOnce(const_cast<BoutReal &>(BOUT_VERSION), "BOUT_VERSION");

  /// Open the restart file for writing
  if (!restart.openw("%s",filename.c_str()))
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
