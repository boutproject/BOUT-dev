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

class PhysicsModel;

#ifndef __PHYSICS_MODEL_H__
#define __PHYSICS_MODEL_H__

#include <bout.hxx>
#include <options.hxx>
#include <msg_stack.hxx>
#include "solver.hxx"

class PhysicsModel {
public:
  typedef int (PhysicsModel::*rhsfunc)(BoutReal t);
  typedef int (PhysicsModel::*preconfunc)(BoutReal t, BoutReal gamma, BoutReal delta);
  typedef int (PhysicsModel::*jacobianfunc)(BoutReal t);
  
  PhysicsModel() : solver(0), userrhs(&PhysicsModel::rhs), userconv(0), userdiff(0), 
                   userprecon(0), userjacobian(0) {}
  ~PhysicsModel();
  
  int initialise(Solver *s, bool restarting) {
    solver = s;
    return init(restarting); // Call user code
  }
  
  int runRHS(BoutReal time);
  
  bool splitOperator();
  int runConvective(BoutReal time);
  int runDiffusive(BoutReal time);
  
  bool hasPrecon();
  int runPrecon(BoutReal t, BoutReal gamma, BoutReal delta);
  
  bool hasJacobian();
  int runJacobian(BoutReal t);
protected:

  // These two functions implemented by user code
  virtual int init(bool restarting) = 0;
  virtual int rhs(BoutReal t) {return 1;} 

  // Functions called by the user to set callback functions
  void setRHS(rhsfunc fset);
  void setSplitOperator(rhsfunc conv, rhsfunc diff);
  void setPrecon(preconfunc pset);
  void setJacobian(jacobianfunc pset);

  Solver *solver;
  void bout_solve(Field2D &var, const char *name);
  void bout_solve(Field3D &var, const char *name);
  void bout_solve(Vector2D &var, const char *name);
  void bout_solve(Vector3D &var, const char *name);
  
  bool bout_constrain(Field3D &var, Field3D &F_var, const char *name);
private:
  rhsfunc      userrhs;
  rhsfunc      userconv, userdiff; // Split operator functions
  preconfunc   userprecon;
  jacobianfunc userjacobian;
};

// Macro to define a simple main() which creates
// the given model and runs it.
#define BOUTMAIN(ModelClass)               \
  int main(int argc, char **argv) {        \
    BoutInitialise(argc, argv);            \
    ModelClass *model = new ModelClass();  \
    Solver *solver = Solver::create();     \
    solver->setModel(model);               \
    solver->addMonitor(bout_monitor);      \
    solver->solve();                       \
    delete model;                          \
    delete solver;                         \
    BoutFinalise();                        \
    return 0;                              \
  }



#endif // __PHYSICS_MODEL_H__

