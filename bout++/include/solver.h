/**************************************************************************
 * Base class for all solvers. Specifies required interface functions
 *
 * Changelog:
 * 
 * 2009-08 Ben Dudson, Sean Farley
 *    * Major overhaul, and changed API. Trying to make consistent
 *      interface to PETSc and SUNDIALS solvers
 * 
 **************************************************************************
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

#ifndef __SOLVER_H__
#define __SOLVER_H__

#include "field2d.h"
#include "field3d.h"
#include "vector2d.h"
#include "vector3d.h"

#include "globals.h"

#include <string>
using std::string;

#define SolverType const char*
#define SOLVERCVODE       "cvode"
#define SOLVERPVODE       "pvode"
#define SOLVERIDA         "ida"
#define SOLVERPETSC       "petsc"

enum SOLVER_VAR_OP {LOAD_VARS, LOAD_DERIVS, SET_ID, SAVE_VARS, SAVE_DERIVS};

/// RHS function pointer
typedef int (*rhsfunc)(BoutReal);

/// User-supplied preconditioner function
typedef int (*PhysicsPrecon)(BoutReal t, BoutReal gamma, BoutReal delta);

/// User-supplied Jacobian function
typedef int (*Jacobian)(BoutReal t);

/// Solution monitor, called each timestep
typedef int (*MonitorFunc)(BoutReal simtime, int iter, int NOUT);

class Solver {
 public:
  Solver();
  virtual ~Solver() { }
  
  // Routines to add variables. Solvers can just call these
  // (or leave them as-is)
  virtual void add(Field2D &v, Field2D &F_v, const char* name);
  virtual void add(Field3D &v, Field3D &F_v, const char* name);
  virtual void add(Vector2D &v, Vector2D &F_v, const char* name);
  virtual void add(Vector3D &v, Vector3D &F_v, const char* name);
  
  virtual bool constraints() {return has_constraints; } ///< Returns true if constraints available

  // Add constraint functions (optional)
  virtual void constraint(Field2D &v, Field2D &C_v, const char* name);
  virtual void constraint(Field3D &v, Field3D &C_v, const char* name);
  virtual void constraint(Vector2D &v, Vector2D &C_v, const char* name);
  virtual void constraint(Vector3D &v, Vector3D &C_v, const char* name);

  /// Specify a preconditioner (optional)
  virtual void setPrecon(PhysicsPrecon f) {}
  
  /// Specify a Jacobian (optional)
  virtual void setJacobian(Jacobian j) {}

  /// Initialise the solver, passing the RHS function
  /// NOTE: nout and tstep should be passed to run, not init.
  ///       Needed because of how the PETSc TS code works
  virtual int init(rhsfunc f, int argc, char **argv, bool restarting, int nout, BoutReal tstep);
  
  /// Run the solver, calling MonitorFunc nout times, at intervals of tstep
  virtual int run(MonitorFunc f) = 0;
  
  /// Clean-up code. Some solvers (PETSc) need things to be cleaned up early
  //virtual int free() {}; 
  
  // Solver status. Optional functions used to query the solver
  virtual int n2Dvars() const {return f2d.size();}  ///< Number of 2D variables. Vectors count as 3
  virtual int n3Dvars() const {return f3d.size();}  ///< Number of 3D variables. Vectors count as 3

  BoutReal rhs_wtime; ///< Wall time used in RHS
  int rhs_ncalls; ///< Number of calls to the RHS function

  void setRestartDir(const string &dir);
  void setRestartDir(const char* dir) {string s = string(dir); setRestartDir(s); }
  
  static Solver* Create();
  static Solver* Create(SolverType &type);
protected:

  int NPES, MYPE; ///< Number of processors and this processor's index
  
  /// Calculate the number of evolving variables on this processor
  int getLocalN();
  
  /// A structure to hold an evolving variable
  template <class T>
    struct VarStr {
      bool constraint;
      T *var;
      T *F_var;
      CELL_LOC location; // For fields and vector components
      bool covariant; // For vectors
      
      string name;    // Name of the variable
    };
  
  /// Vectors of variables to evolve
  vector< VarStr<Field2D> > f2d;
  vector< VarStr<Field3D> > f3d;
  vector< VarStr<Vector2D> > v2d;
  vector< VarStr<Vector3D> > v3d;
  
  Datafile restart; ///< Restart file object

  string restartdir;  ///< Directory for restart files
  string restartext;  ///< Restart file extension
  int archive_restart;

  bool has_constraints; ///< Can this solver handle constraints? Set to true if so.
  bool initialised; ///< Has init been called yet?

  BoutReal simtime;  ///< Current simulation time
  int iteration; ///< Current iteration (output time-step) number
};

#endif // __SOLVER_H__
