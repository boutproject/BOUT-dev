/**************************************************************************
 * Base class for all solvers. Specifies required interface functions
 *
 * Changelog:
 * 
 * 2009-08 Ben Dudson, Sean Farley
 *    * Major overhaul, and changed API. Trying to make consistent
 *      interface to PETSc and SUNDIALS solvers
 * 
 * 2013-08 Ben Dudson
 *    * Added OO-style API, to allow multiple physics models to coexist
 *      For now both APIs are supported
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

class Solver;

#include <bout_types.hxx>

///////////////////////////////////////////////////////////////////
// C function pointer types

/// RHS function pointer
typedef int (*rhsfunc)(BoutReal); // C-style function pointer

/// User-supplied preconditioner function
typedef int (*PhysicsPrecon)(BoutReal t, BoutReal gamma, BoutReal delta);

/// User-supplied Jacobian function
typedef int (*Jacobian)(BoutReal t);


/// Solution monitor, called each timestep
typedef int (*MonitorFunc)(Solver *solver, BoutReal simtime, int iter, int NOUT);

///////////////////////////////////////////////////////////////////

#ifndef __SOLVER_H__
#define __SOLVER_H__

#include "globals.hxx"
#include "field2d.hxx"
#include "field3d.hxx"
#include "vector2d.hxx"
#include "vector3d.hxx"

#include "physicsmodel.hxx"

#include <string>
#include <list>
using std::string;

#define SolverType const char*
#define SOLVERCVODE       "cvode"
#define SOLVERPVODE       "pvode"
#define SOLVERIDA         "ida"
#define SOLVERPETSC       "petsc"
#define SOLVERKARNIADAKIS "karniadakis"
#define SOLVERRK4         "rk4"
#define SOLVEREULER       "euler"

enum SOLVER_VAR_OP {LOAD_VARS, LOAD_DERIVS, SET_ID, SAVE_VARS, SAVE_DERIVS};

///////////////////////////////////////////////////////////////////

class Solver {
 public:
  Solver(Options *opts = NULL);
  virtual ~Solver() { }

  /////////////////////////////////////////////
  // New API
  
  void setModel(PhysicsModel *model); ///< Specify physics model to solve

  /////////////////////////////////////////////
  // Old API
  
  void setRHS(rhsfunc f) { phys_run = f; } ///< Set the RHS function
  virtual void setPrecon(PhysicsPrecon f) {} ///< Specify a preconditioner (optional)
  virtual void setJacobian(Jacobian j) {} ///< Specify a Jacobian (optional)
  virtual void setSplitOperator(rhsfunc fC, rhsfunc fD); ///< Split operator solves

  /////////////////////////////////////////////
  // Routines to add variables. Solvers can just call these
  // (or leave them as-is)
  virtual void add(Field2D &v, const char* name);
  virtual void add(Field3D &v, const char* name);
  virtual void add(Vector2D &v, const char* name);
  virtual void add(Vector3D &v, const char* name);
  
  virtual bool constraints() {return has_constraints; } ///< Returns true if constraints available

  // Add constraint functions (optional)
  virtual void constraint(Field2D &v, Field2D &C_v, const char* name);
  virtual void constraint(Field3D &v, Field3D &C_v, const char* name);
  virtual void constraint(Vector2D &v, Vector2D &C_v, const char* name);
  virtual void constraint(Vector3D &v, Vector3D &C_v, const char* name);
  
  /// Set a maximum internal timestep (only for explicit schemes)
  virtual void setMaxTimestep(BoutReal dt) {max_dt = dt;}
  /// Return the current internal timestep 
  virtual BoutReal getCurrentTimestep() {return 0.0;}
  
  int solve();

  /// Initialise the solver
  /// NOTE: nout and tstep should be passed to run, not init.
  ///       Needed because of how the PETSc TS code works
  virtual int init(bool restarting, int nout, BoutReal tstep);
  
  void addMonitor(MonitorFunc f);     ///< Add a monitor function to be called every output
  void removeMonitor(MonitorFunc f);  ///< Remove a monitor function previously added

  /// Run the solver, calling monitors nout times, at intervals of tstep
  virtual int run() = 0;
  
  // Solver status. Optional functions used to query the solver
  virtual int n2Dvars() const {return f2d.size();}  ///< Number of 2D variables. Vectors count as 3
  virtual int n3Dvars() const {return f3d.size();}  ///< Number of 3D variables. Vectors count as 3
  
  int rhs_ncalls; ///< Number of calls to the RHS function

  void setRestartDir(const string &dir);
  void setRestartDir(const char* dir) {string s = string(dir); setRestartDir(s); }
  
  static Solver* create(Options *opts = NULL);
  static Solver* create(SolverType &type, Options *opts = NULL);
  
  static void setArgs(int &c, char **&v) { pargc = &c; pargv = &v;}
protected:
  bool restarting;
  
  // Command-line arguments
  static int* pargc;
  static char*** pargv;

  // Settings to use during initialisation (set by constructor)
  Options *options;

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

  bool has_constraints; ///< Can this solver.hxxandle constraints? Set to true if so.
  bool initialised; ///< Has init been called yet?

  BoutReal simtime;  ///< Current simulation time
  int iteration; ///< Current iteration (output time-step) number

  int run_rhs(BoutReal t); ///< Run the user's RHS function
  int run_convective(BoutReal t); ///< Calculate only the convective parts
  int run_diffusive(BoutReal t); ///< Calculate only the diffusive parts
  
  int call_monitors(BoutReal simtime, int iter, int NOUT); ///< Calls all monitor functions

  // Loading data from BOUT++ to/from solver
  void load_vars(BoutReal *udata);
  void load_derivs(BoutReal *udata);
  int save_vars(BoutReal *udata);
  void save_derivs(BoutReal *dudata);
  
  BoutReal max_dt; ///< Maximum internal timestep
 private:
  PhysicsModel *model;    ///< physics model being evolved
  
  rhsfunc phys_run;       ///< The user's RHS function
  bool split_operator;
  rhsfunc phys_conv, phys_diff; ///< Convective and Diffusive parts (if split operator)
  
  std::list<MonitorFunc> monitors; ///< List of monitor functions

  void post_rhs(); // Should be run after user RHS is called
  
  // Loading data from BOUT++ to/from solver
  void loop_vars_op(int jx, int jy, BoutReal *udata, int &p, SOLVER_VAR_OP op);
  void loop_vars(BoutReal *udata, SOLVER_VAR_OP op);

  bool varAdded(const string &name); // Check if a variable has already been added
};

#endif // __SOLVER_H__
