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
#include <boutexception.hxx>
#include <unused.hxx>
#include "bout/monitor.hxx"

///////////////////////////////////////////////////////////////////
// C function pointer types

/// RHS function pointer
typedef int (*rhsfunc)(BoutReal); // C-style function pointer

/// User-supplied preconditioner function
typedef int (*PhysicsPrecon)(BoutReal t, BoutReal gamma, BoutReal delta);

/// User-supplied Jacobian function
typedef int (*Jacobian)(BoutReal t);


/// Solution monitor, called each timestep
typedef int (*TimestepMonitorFunc)(Solver *solver, BoutReal simtime, BoutReal lastdt);

///////////////////////////////////////////////////////////////////

#ifndef __SOLVER_H__
#define __SOLVER_H__

//#include "globals.hxx"
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
#define SOLVERSLEPC       "slepc"
#define SOLVERKARNIADAKIS "karniadakis"
#define SOLVERRK4         "rk4"
#define SOLVEREULER       "euler"
#define SOLVERRK3SSP      "rk3ssp"
#define SOLVERPOWER       "power"
#define SOLVERARKODE	  "arkode"
#define SOLVERIMEXBDF2    "imexbdf2"
#define SOLVERSNES        "snes"
#define SOLVERRKGENERIC   "rkgeneric"

enum SOLVER_VAR_OP {LOAD_VARS, LOAD_DERIVS, SET_ID, SAVE_VARS, SAVE_DERIVS};

///////////////////////////////////////////////////////////////////

/*!
 * Interface to integrators, mainly for time integration
 *
 * 
 * Creation
 * --------
 * 
 * Solver is a base class and can't be created directly:
 * 
 *     Solver *solver = Solver(); // Error
 * 
 * Instead, use the create() static function:
 * 
 *     Solver *solver = Solver::create(); // ok
 *
 * By default this will use the options in the "solver" section
 * of the options, equivalent to:
 * 
 *     Options *opts = Options::getRoot()->getSection("solver");
 *     Solver *solver = Solver::create(opts);
 * 
 * To use a different set of options, for example if there are
 * multiple solvers, use a different option section:
 *
 *     Options *opts = Options::getRoot()->getSection("anothersolver");
 *     Solver *anothersolver = Solver::create(opts);
 *
 * Problem specification
 * ---------------------
 * 
 * The equations to be solved are specified in a PhysicsModel object
 *
 *     class MyProblem : public PhysicsModel {
 *       protected:
 *         // This function called once at beginning
 *         int init(bool restarting) {
 *           SOLVE_FOR(f);   // Specify variables to solve
 *           // Set f to initial value if needed
 *           // otherwise read from input options
 *           return 0;
 *         }
 *
 *         // This function called to evaluate time derivatives
 *         int rhs(BoutReal t) {
 *            ddt(f) = 1.0; // Calculate time derivatives
 *            return 0;
 *         }
 *       private:
 *         Field3D f; // A variable to evolve
 *     }
 * 
 * The init() and rhs() functions must be defined, but there
 * are other functions which can be defined. See PhysicsModel
 * documentation for details.
 * 
 * Create an object, then add to the solver:
 * 
 *     MyProblem *prob = MyProblem();
 *     solver->setModel(prob);
 * 
 * Running simulation
 * ------------------
 * 
 * To run a calculation
 * 
 *     solver->solve();
 * 
 * This will use NOUT and TIMESTEP in the solver options 
 * (specified during creation). If these are not present
 * then the global options will be used.
 * 
 * To specify NOUT and TIMESTEP, pass the values to solve:
 *
 *     solver->solve(NOUT, TIMESTEP);
 */
class Solver {
 public:
  Solver(Options *opts = nullptr);
  virtual ~Solver();

  /////////////////////////////////////////////
  // New API
  
  /*!
   * Specify physics model to solve. Currently only one model
   * can be evolved by a Solver.
   */ 
  virtual void setModel(PhysicsModel *model);

  /////////////////////////////////////////////
  // Old API
  
  virtual void setRHS(rhsfunc f) { phys_run = f; } ///< Set the RHS function
  void setPrecon(PhysicsPrecon f) {prefunc = f;} ///< Specify a preconditioner (optional)
  virtual void setJacobian(Jacobian UNUSED(j)) {} ///< Specify a Jacobian (optional)
  virtual void setSplitOperator(rhsfunc fC, rhsfunc fD); ///< Split operator solves
  
  /////////////////////////////////////////////
  // Monitors
  
  enum MonitorPosition {BACK, FRONT}; ///< A type to set where in the list monitors are added

  /// Add a monitor to be called every output
  void addMonitor(Monitor * f, MonitorPosition pos=FRONT);
  void removeMonitor(Monitor * f);  ///< Remove a monitor function previously added

  void addTimestepMonitor(TimestepMonitorFunc f);    ///< Add a monitor function to be called every timestep
  void removeTimestepMonitor(TimestepMonitorFunc f); ///< Remove a previously added timestep monitor

  /////////////////////////////////////////////
  // Routines to add variables. Solvers can just call these
  // (or leave them as-is)
  
  /*!
   * Add a variable to be solved. This must be done
   * in the initialisation stage, before the simulation starts.
   */ 
  virtual void add(Field2D &v, const char* name);
  virtual void add(Field3D &v, const char* name);
  virtual void add(Vector2D &v, const char* name);
  virtual void add(Vector3D &v, const char* name);
  
  /*!
   * Returns true if constraints available
   */ 
  virtual bool constraints() {return has_constraints; }
  
  /*!
   * Add constraint functions (optional). These link a variable
   * v to a control parameter C_v such that v is adjusted 
   * to keep C_v = 0.
   */
  virtual void constraint(Field2D &v, Field2D &C_v, const char* name);
  virtual void constraint(Field3D &v, Field3D &C_v, const char* name);
  virtual void constraint(Vector2D &v, Vector2D &C_v, const char* name);
  virtual void constraint(Vector3D &v, Vector3D &C_v, const char* name);
  
  /// Set a maximum internal timestep (only for explicit schemes)
  virtual void setMaxTimestep(BoutReal dt) {max_dt = dt;}
  /// Return the current internal timestep 
  virtual BoutReal getCurrentTimestep() {return 0.0;}
  
  /*!
   * Start the solver. By default solve() uses options
   * to determine the number of steps and the output timestep.
   * If nout and dt are specified here then the options are not used
   * 
   * @param[in] nout   Number of output timesteps
   * @param[in] dt     The time between outputs
   */
  int solve(int nout=-1, BoutReal dt=0.0);

  /// Initialise the solver
  /// NOTE: nout and tstep should be passed to run, not init.
  ///       Needed because of how the PETSc TS code works
  virtual int init(int nout, BoutReal tstep);

  /*!
   * Run the solver, calling monitors nout times, at intervals of tstep 
   * This function is called by solve(), and is specific to each solver type
   * 
   * This should probably be protected, since it shouldn't be called
   * by users. 
   */
  virtual int run() = 0;

  //Should wipe out internal field vector and reset from current field object data
  virtual void resetInternalFields(){
    throw BoutException("resetInternalFields not supported by this Solver");}

  // Solver status. Optional functions used to query the solver
  virtual int n2Dvars() const {return f2d.size();}  ///< Number of 2D variables. Vectors count as 3
  virtual int n3Dvars() const {return f3d.size();}  ///< Number of 3D variables. Vectors count as 3
  
  int rhs_ncalls,rhs_ncalls_e,rhs_ncalls_i; ///< Number of calls to the RHS function
  
  /*!
   * Test if this solver supports split operators (e.g. implicit/explicit)
   */ 
  bool splitOperator() {return split_operator;}

  bool canReset;
  
  /// Add evolving variables to output (dump) file or restart file
  ///
  /// @param[inout] outputfile   The file to add variable to
  /// @param[in] save_repeat    If true, add variables with time dimension
  void outputVars(Datafile &outputfile, bool save_repeat=true);

  /*!
   * Create a Solver object. This uses the "type" option
   * in the given Option section to determine which solver
   * type to create.
   */ 
  static Solver* create(Options *opts = NULL);
  
  /*!
   * Create a Solver object, specifying the type
   */ 
  static Solver* create(SolverType &type, Options *opts = NULL);
  
  /*!
   * Pass the command-line arguments. This static function is
   * called by BoutInitialise, and puts references
   * into protected variables. These may then be used by Solvers
   * to control behavior
   * 
   */ 
  static void setArgs(int &c, char **&v) { pargc = &c; pargv = &v;}
  
protected:
  
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
      T *MMS_err;        // Error for MMS
      CELL_LOC location; // For fields and vector components
      bool covariant; // For vectors
      bool evolve_bndry; // Are the boundary regions being evolved?

      string name;    // Name of the variable
    };
  
  /// Vectors of variables to evolve
  vector< VarStr<Field2D> > f2d;
  vector< VarStr<Field3D> > f3d;
  vector< VarStr<Vector2D> > v2d;
  vector< VarStr<Vector3D> > v3d;
  
  bool has_constraints; ///< Can this solver.hxxandle constraints? Set to true if so.
  bool initialised; ///< Has init been called yet?

  BoutReal simtime;  ///< Current simulation time
  int iteration; ///< Current iteration (output time-step) number

  int run_rhs(BoutReal t); ///< Run the user's RHS function
  int run_convective(BoutReal t); ///< Calculate only the convective parts
  int run_diffusive(BoutReal t, bool linear=true); ///< Calculate only the diffusive parts
  
  int call_monitors(BoutReal simtime, int iter, int NOUT); ///< Calls all monitor functions
  
  bool monitor_timestep; ///< Should timesteps be monitored?
  int call_timestep_monitors(BoutReal simtime, BoutReal lastdt);

  bool have_user_precon(); // Do we have a user preconditioner?
  int run_precon(BoutReal t, BoutReal gamma, BoutReal delta);
  
  // Loading data from BOUT++ to/from solver
  void load_vars(BoutReal *udata);
  void load_derivs(BoutReal *udata);
  void save_vars(BoutReal *udata);
  void save_derivs(BoutReal *dudata);
  void set_id(BoutReal *udata);
  
  // 
  const Field3D globalIndex(int localStart);
  
  BoutReal max_dt; ///< Maximum internal timestep
  
private:
  bool initCalled=false; ///< Has the init function of the solver been called?
  int freqDefault=1;     ///< Default sampling rate at which to call monitors - same as output to screen
  BoutReal timestep=-1; ///< timestep - shouldn't be changed after init is called.
  PhysicsModel *model;    ///< physics model being evolved

  rhsfunc phys_run;       ///< The user's RHS function
  PhysicsPrecon prefunc;  // Preconditioner
  bool split_operator;
  rhsfunc phys_conv, phys_diff; ///< Convective and Diffusive parts (if split operator)

  bool mms; ///< Enable sources and solutions for Method of Manufactured Solutions
  bool mms_initialise; ///< Initialise variables to the manufactured solution

  void add_mms_sources(BoutReal t);
  void calculate_mms_error(BoutReal t);
  
  std::list<Monitor*> monitors; ///< List of monitor functions
  std::list<TimestepMonitorFunc> timestep_monitors; ///< List of timestep monitor functions

  void pre_rhs(BoutReal t); // Should be run before user RHS is called
  void post_rhs(BoutReal t); // Should be run after user RHS is called
  
  // Loading data from BOUT++ to/from solver
  void loop_vars_op(int jx, int jy, BoutReal *udata, int &p, SOLVER_VAR_OP op, bool bndry);
  void loop_vars(BoutReal *udata, SOLVER_VAR_OP op);

  bool varAdded(const string &name); // Check if a variable has already been added
};

#endif // __SOLVER_H__
