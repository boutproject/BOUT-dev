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

PhysicsModel::~PhysicsModel() {
}

int PhysicsModel::runRHS(BoutReal time) {
  if (!userrhs)
    return 1;
  return (*this.*userrhs)(time);
}

bool PhysicsModel::splitOperator() {
  return (userconv != 0);
}

int PhysicsModel::runConvective(BoutReal time) {
  if(!userconv)
    return 1;
  
  return (*this.*userconv)(time);
}

int PhysicsModel::runDiffusive(BoutReal time) {
  if(!userdiff)
    return 1;
  
  return (*this.*userdiff)(time);
}

bool PhysicsModel::hasPrecon() {
  return (userprecon != 0);
}

int PhysicsModel::runPrecon(BoutReal t, BoutReal gamma, BoutReal delta) {
  if(!userprecon)
    return 1;
  return (*this.*userprecon)(t, gamma, delta);
}

bool PhysicsModel::hasJacobian() {
  return (userjacobian != 0);
}

int PhysicsModel::runJacobian(BoutReal t) {
  if(!userjacobian)
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
