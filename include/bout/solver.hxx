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

#ifndef SOLVER_H
#define SOLVER_H

#include "bout/build_config.hxx"

#include "bout/bout_types.hxx"
#include "bout/boutexception.hxx"
#include "bout/monitor.hxx"
#include "bout/options.hxx"
#include "bout/unused.hxx"

#include <memory>

///////////////////////////////////////////////////////////////////
// C function pointer types

class Solver;

/// RHS function pointer
using rhsfunc = int (*)(BoutReal);

/// User-supplied preconditioner function
using PhysicsPrecon = int (*)(BoutReal t, BoutReal gamma, BoutReal delta);

/// User-supplied Jacobian function
using Jacobian = int (*)(BoutReal t);

/// Solution monitor, called each timestep
using TimestepMonitorFunc = int (*)(Solver* solver, BoutReal simtime, BoutReal lastdt);

#include "bout/field2d.hxx"
#include "bout/field3d.hxx"
#include "bout/generic_factory.hxx"
#include "bout/vector2d.hxx"
#include "bout/vector3d.hxx"

#define BOUT_NO_USING_NAMESPACE_BOUTGLOBALS
#include "physicsmodel.hxx"
#undef BOUT_NO_USING_NAMESPACE_BOUTGLOBALS

#include <list>
#include <string>

using SolverType = std::string;
constexpr auto SOLVERCVODE = "cvode";
constexpr auto SOLVERPVODE = "pvode";
constexpr auto SOLVERIDA = "ida";
constexpr auto SOLVERPETSC = "petsc";
constexpr auto SOLVERSLEPC = "slepc";
constexpr auto SOLVERRK4 = "rk4";
constexpr auto SOLVEREULER = "euler";
constexpr auto SOLVERRK3SSP = "rk3ssp";
constexpr auto SOLVERPOWER = "power";
constexpr auto SOLVERARKODE = "arkode";
constexpr auto SOLVERIMEXBDF2 = "imexbdf2";
constexpr auto SOLVERSNES = "snes";
constexpr auto SOLVERRKGENERIC = "rkgeneric";

enum class SOLVER_VAR_OP { LOAD_VARS, LOAD_DERIVS, SET_ID, SAVE_VARS, SAVE_DERIVS };

/// A type to set where in the list monitors are added
enum class MonitorPosition { BACK, FRONT };

class SolverFactory : public Factory<Solver, SolverFactory, Options*> {
public:
  static constexpr auto type_name = "Solver";
  static constexpr auto section_name = "solver";
  static constexpr auto option_name = "type";
  static constexpr auto default_type =
#if BOUT_HAS_CVODE
      SOLVERCVODE;
#elif BOUT_HAS_IDA
      SOLVERIDA;
#else
      SOLVERPVODE;
#endif
};

/// Simpler name for Factory registration helper class
///
/// Usage:
///
///     #include <bout/solverfactory.hxx>
///     namespace {
///     RegisterSolver<MySolver> registersolvermine("mysolver");
///     }
template <typename DerivedType>
using RegisterSolver = SolverFactory::RegisterInFactory<DerivedType>;

using RegisterUnavailableSolver = SolverFactory::RegisterUnavailableInFactory;

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
 *     auto solver = Solver::create(); // ok
 *
 * By default this will use the options in the "solver" section
 * of the options, equivalent to:
 *
 *     Options *opts = Options::getRoot()->getSection("solver");
 *     auto solver = Solver::create(opts);
 *
 * To use a different set of options, for example if there are
 * multiple solvers, use a different option section:
 *
 *     Options *opts = Options::getRoot()->getSection("anothersolver");
 *     auto anothersolver = Solver::create(opts);
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
  Solver(Options* opts = nullptr);
  virtual ~Solver() = default;

  /////////////////////////////////////////////
  // New API

  /// Specify physics model to solve. Currently only one model can be
  /// evolved by a Solver.
  virtual void setModel(PhysicsModel* model);

  /////////////////////////////////////////////
  // Monitors

  // Alternative names so that Solver::BACK and Solver::FRONT can be used as names for
  // MonitorPositions, for backward compatibility.
  static constexpr MonitorPosition BACK = MonitorPosition::BACK;
  static constexpr MonitorPosition FRONT = MonitorPosition::FRONT;

  /// Add a \p monitor to be called regularly
  ///
  /// The frequency at which the \p monitor is called is set by
  /// `Monitor::timestep`. By default this is every output
  /// timestep. When a new Monitor with a smaller timestep is added,
  /// the solver attemps to adjust its internal timestep to match, as
  /// well as adjusting the timesteps of the current set of Monitors
  /// to be multiples of the new timestep. If this is not possible,
  /// `addMonitor` will throw an exception.
  ///
  /// Adding new Monitors after the Solver has been initialised is
  /// only possible if their timestep is a multiple of the Solver's
  /// timestep. Smaller timesteps will throw an exception.
  void addMonitor(Monitor* monitor, MonitorPosition pos = MonitorPosition::FRONT);
  /// Remove a monitor function previously added
  void removeMonitor(Monitor* monitor);

  /// Add a monitor function to be called every timestep
  void addTimestepMonitor(TimestepMonitorFunc monitor);
  /// Remove a previously added timestep monitor
  void removeTimestepMonitor(TimestepMonitorFunc monitor);

  /////////////////////////////////////////////
  // Routines to add variables. Solvers can just call these
  // (or leave them as-is)

  /// Add a variable to be solved. This must be done in the
  /// initialisation stage, before the simulation starts.
  virtual void add(Field2D& v, const std::string& name,
                   const std::string& description = "");
  virtual void add(Field3D& v, const std::string& name,
                   const std::string& description = "");
  virtual void add(Vector2D& v, const std::string& name,
                   const std::string& description = "");
  virtual void add(Vector3D& v, const std::string& name,
                   const std::string& description = "");

  /// Returns true if constraints available
  virtual bool constraints() { return has_constraints; }

  /// Add constraint functions (optional). These link a variable v to
  /// a control parameter C_v such that v is adjusted to keep C_v = 0.
  virtual void constraint(Field2D& v, Field2D& C_v, std::string name);
  virtual void constraint(Field3D& v, Field3D& C_v, std::string name);
  virtual void constraint(Vector2D& v, Vector2D& C_v, std::string name);
  virtual void constraint(Vector3D& v, Vector3D& C_v, std::string name);

  /// Set a maximum internal timestep (only for explicit schemes)
  virtual void setMaxTimestep([[maybe_unused]] BoutReal dt) {}
  /// Return the current internal timestep
  virtual BoutReal getCurrentTimestep() { return 0.0; }

  /// Start the solver. By default solve() uses options
  /// to determine the number of steps and the output timestep.
  /// If nout and dt are specified here then the options are not used
  ///
  /// @param[in] nout     Number of output timesteps to run for
  /// @param[in] timestep The time between outputs
  int solve(int nout = -1, BoutReal timestep = 0.0);

  /// Initialise the solver
  virtual int init();

  /// Run the solver, calling monitors nout times, at intervals of
  /// tstep. This function is called by solve(), and is specific to
  /// each solver type
  ///
  /// This should probably be protected, since it shouldn't be called
  /// by users.
  virtual int run() = 0;

  /// Should wipe out internal field vector and reset from current field object data
  virtual void resetInternalFields() {
    throw BoutException("resetInternalFields not supported by this Solver");
  }

  // Solver status. Optional functions used to query the solver
  /// Number of 2D variables. Vectors count as 3
  virtual int n2Dvars() const { return static_cast<int>(f2d.size()); }
  /// Number of 3D variables. Vectors count as 3
  virtual int n3Dvars() const { return static_cast<int>(f3d.size()); }

  /// Get and reset the number of calls to the RHS function
  int resetRHSCounter();
  /// Same but for explicit timestep counter - for IMEX
  int resetRHSCounter_e();
  /// Same but fur implicit timestep counter - for IMEX
  int resetRHSCounter_i();

  /// Test if this solver supports split operators (e.g. implicit/explicit)
  bool splitOperator();

  bool canReset{false};

  /// Add evolving variables to output (dump) file or restart file
  ///
  /// @param[inout] outputfile  The `Options` to add variable to
  /// @param[in] save_repeat    If true, add variables with time dimension
  virtual void outputVars(Options& output_options, bool save_repeat = true);

  /// Copy evolving variables out of \p options
  virtual void readEvolvingVariablesFromOptions(Options& options);

  /// Create a Solver object. This uses the "type" option in the given
  /// Option section to determine which solver type to create.
  static std::unique_ptr<Solver> create(Options* opts = nullptr);

  /// Create a Solver object, specifying the type
  static std::unique_ptr<Solver> create(const SolverType& type, Options* opts = nullptr);

  /// Pass the command-line arguments. This static function is
  /// called by BoutInitialise, and puts references
  /// into protected variables. These may then be used by Solvers
  /// to control behavior
  static void setArgs(int& c, char**& v) {
    pargc = &c;
    pargv = &v;
  }

  /// A unique identifier for this run. Throws if the identifier hasn't been set yet.
  std::string getRunID() const;
  /// The run from which this was restarted. Throws if the identifier hasn't been set yet.
  std::string getRunRestartFrom() const;

  /// Get the number of completed output steps
  int getIterationCounter() const { return iteration; }

  /// Add one to the iteration count, used by BoutMonitor, but could be called by a
  // user-defined monitor (if `bout_run()` is not used)
  int incrementIterationCounter() { return ++iteration; }

  /// Write \p options to the model's output file
  void writeToModelOutputFile(const Options& options);

protected:
  /// Number of command-line arguments
  static int* pargc;
  /// Command-line arguments
  static char*** pargv;

  /// Settings to use during initialisation (set by constructor)
  Options* options{nullptr};

  /// Number of processors
  int NPES{1};
  /// This processor's index
  int MYPE{0};

  /// Calculate the number of evolving variables on this processor
  int getLocalN();

  /// A structure to hold an evolving variable
  template <class T>
  struct VarStr {
    bool constraint{false};              /// Does F_var represent a constraint?
    T* var{nullptr};                     /// The evolving variable
    T* F_var{nullptr};                   /// The time derivative or constraint on var
    std::unique_ptr<T> MMS_err{nullptr}; /// Error for MMS
    CELL_LOC location{CELL_DEFAULT};     /// For fields and vector components
    bool covariant{false};               /// For vectors
    bool evolve_bndry{false};            /// Are the boundary regions being evolved?
    std::string name;                    /// Name of the variable
    std::string description{""};         /// Description of what the variable is
  };

  /// Does \p var represent field \p name?
  template <class T>
  friend bool operator==(const VarStr<T>& var, const std::string& name) {
    return var.name == name;
  }

  /// Does \p vars contain a field with \p name?
  template <class T>
  bool contains(const std::vector<VarStr<T>>& vars, const std::string& name) {
    const auto in_vars = std::find(begin(vars), end(vars), name);
    return in_vars != end(vars);
  }

  /// Vectors of variables to evolve
  std::vector<VarStr<Field2D>> f2d;
  std::vector<VarStr<Field3D>> f3d;
  std::vector<VarStr<Vector2D>> v2d;
  std::vector<VarStr<Vector3D>> v3d;

  /// Vectors of diagnostic variables to save
  std::vector<VarStr<int>> diagnostic_int;
  std::vector<VarStr<BoutReal>> diagnostic_BoutReal;
  void add_int_diagnostic(int& i, const std::string& name,
                          const std::string& description = "") {
    VarStr<int> v;
    v.var = &i;
    v.name = name;
    v.description = description;
    diagnostic_int.emplace_back(std::move(v));
  };
  void add_BoutReal_diagnostic(BoutReal& r, const std::string& name,
                               const std::string& description = "") {
    VarStr<BoutReal> v;
    v.var = &r;
    v.name = name;
    v.description = description;
    diagnostic_BoutReal.emplace_back(std::move(v));
  };

  /// Can this solver handle constraints? Set to true if so.
  bool has_constraints{false};
  /// Has init been called yet?
  bool initialised{false};

  /// Current simulation time
  BoutReal simtime{0.0};

  /// Run the user's RHS function
  int run_rhs(BoutReal t, bool linear = false);
  /// Calculate only the convective parts
  int run_convective(BoutReal t, bool linear = false);
  /// Calculate only the diffusive parts
  int run_diffusive(BoutReal t, bool linear = false);

  /// Reset the iteration counter
  void resetIterationCounter(int value = 0) { iteration = value; }

  /// Calls all monitor functions
  ///
  /// There are two important things to note about how \p iter is
  /// passed along to each monitor:
  /// - The solvers all start their iteration numbering from zero, so the
  ///   initial state is calculated at \p iter = -1
  /// - Secondly, \p iter is passed along to each monitor *relative to
  ///   that monitor's period*
  ///
  /// In practice, this means that each monitor is called like:
  ///
  ///     monitor->call(solver, simulation_time,
  ///                   ((iter + 1) / monitor->period) - 1,
  ///                   NOUT / monitor->period);
  ///
  /// e.g. for a monitor with period 10, passing \p iter = 9 will
  /// result in it being called with a value of `(9 + 1)/10 - 1 == 0`
  int call_monitors(BoutReal simtime, int iter, int NOUT);

  /// Should timesteps be monitored?
  bool monitor_timestep{false};
  int call_timestep_monitors(BoutReal simtime, BoutReal lastdt);

  /// Do we have a user preconditioner?
  bool hasPreconditioner();
  /// Run the user preconditioner
  int runPreconditioner(BoutReal time, BoutReal gamma, BoutReal delta);

  /// Do we have a user Jacobian?
  bool hasJacobian();
  /// Run the user Jacobian
  int runJacobian(BoutReal time);

  // Loading data from BOUT++ to/from solver
  void load_vars(BoutReal* udata);
  void load_derivs(BoutReal* udata);
  void save_vars(BoutReal* udata);
  void save_derivs(BoutReal* dudata);
  void set_id(BoutReal* udata);

  /// Returns a Field3D containing the global indices
  Field3D globalIndex(int localStart);

  /// Maximum internal timestep
  BoutReal max_dt{-1.0};

  /// Helper struct for Monitor and the name of its time dimension
  struct MonitorInfo {
    /// Non-owning pointer to a monitor instance
    Monitor* monitor;
    /// Name of the time dimension for writing outputs to
    /// file. Monitors with the same period as the outputs have plain
    /// "t", all other monitors have an ascending integer suffix,
    /// starting with the front monitor. WARNING! This could change if
    /// there are different monitors added on subsequent restarts, or
    /// if they are added in a different order.
    std::string time_dimension;
  };

  /// Get the list of monitors
  auto getMonitors() const -> const std::list<MonitorInfo>& { return monitors; }

  /// Get the currently set number of output steps requested
  int getNumberOutputSteps() const { return number_output_steps; }
  /// Change the number of requested output steps
  void setNumberOutputSteps(int nout) { number_output_steps = nout; }
  /// Get the currently set output timestep
  BoutReal getOutputTimestep() const { return output_timestep; }

private:
  /// Generate a random UUID (version 4) and broadcast it to all processors
  std::string createRunID() const;

  /// Default value for `run_id`. Use 'z' because it is not a valid
  /// hex character, so this is an invalid UUID
  static constexpr auto default_run_id = "zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz";

  /// Randomly generated run ID
  /// Initialise with 36 characters so the allocated array is the right size
  std::string run_id = default_run_id;
  /// The run from which this was restarted.
  std::string run_restart_from = "yyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyy";
  /// Save `run_id` and `run_restart_from` every output
  bool save_repeat_run_id{false};

  /// Current iteration (output time-step) number
  int iteration{0};

  /// Number of calls to the RHS function
  int rhs_ncalls{0};
  /// Number of calls to the explicit (convective) RHS function
  int rhs_ncalls_e{0};
  /// Number of calls to the implicit (diffusive) RHS function
  int rhs_ncalls_i{0};
  /// Default sampling rate at which to call monitors - same as output to screen
  int default_monitor_period{1};
  /// timestep - shouldn't be changed after init is called.
  BoutReal internal_timestep{-1};
  /// Physics model being evolved
  PhysicsModel* model{nullptr};

  /// Should non-split physics models be treated as diffusive?
  bool is_nonsplit_model_diffusive{true};

  /// Enable sources and solutions for Method of Manufactured Solutions
  bool mms{false};
  /// Initialise variables to the manufactured solution
  bool mms_initialise{false};

  void add_mms_sources(BoutReal t);
  void calculate_mms_error(BoutReal t);

  /// List of monitor functions
  std::list<MonitorInfo> monitors;
  /// List of timestep monitor functions
  std::list<TimestepMonitorFunc> timestep_monitors;

  /// Should be run before user RHS is called
  void pre_rhs(BoutReal t);
  /// Should be run after user RHS is called
  void post_rhs(BoutReal t);

  /// Loading data from BOUT++ to/from solver
  void loop_vars_op(Ind2D i2d, BoutReal* udata, int& p, SOLVER_VAR_OP op, bool bndry);
  void loop_vars(BoutReal* udata, SOLVER_VAR_OP op);

  /// Check if a variable has already been added
  bool varAdded(const std::string& name);

  /// (Possibly) adjust the periods of \p monitor, and the `monitors`
  /// timesteps, returning the new Solver timestep
  BoutReal adjustMonitorPeriods(Monitor* monitor);

  /// Fix all the monitor periods based on \p output_timestep, as well
  /// as adjusting \p NOUT and \p output_timestep to be consistent
  void finaliseMonitorPeriods(int& NOUT, BoutReal& output_timestep);

  /// Number of requested output steps
  int number_output_steps;
  /// Requested timestep between outputs
  BoutReal output_timestep;
};

#endif // SOLVER_H
