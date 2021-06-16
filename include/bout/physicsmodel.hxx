/*!************************************************************************
 * \file physicsmodel.hxx
 * 
 * @brief Base class for Physics Models
 * 
 * 
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

#include "bout.hxx"
#include "msg_stack.hxx"
#include "options.hxx"
#include "options_netcdf.hxx"
#include "solver.hxx"
#include "unused.hxx"
#include "utils.hxx"
#include "bout/macro_for_each.hxx"
#include "bout/sys/variant.hxx"

#include <vector>
#include <type_traits>

class Mesh;

namespace bout {
/// Stop-gap shim for `DataFile`: allows PhysicsModels to still save
/// outputs using the old `DataFile` interface without upgrading to
/// the new `OptionsNetCDF`.
///
/// Stores pointers to objects in a vector, using a variant, which can
/// then be pushed into an `Options` and written to file using the
/// default implementaion of `PhysicsModel::outputVars`
class DataFileFacade {
public:
  /// This is `Options::ValueType` but using pointers rather than values
  using ValueType = bout::utils::variant<bool*, int*, BoutReal*, std::string*, Field2D*,
                                         Field3D*, FieldPerp*, Array<BoutReal>*,
                                         Matrix<BoutReal>*, Tensor<BoutReal>*>;

  template <class T>
  void addRepeat(T& value, std::string name) {
    add(&value, std::move(name), true);
  }
  template <class T>
  void addOnce(T& value, std::string name) {
    add(&value, std::move(name), false);
  }
  template <class T>
  void add(T& value, const std::string& name, bool save_repeat = false) {
    add(&value, std::move(name), save_repeat);
  }

  void addRepeat(ValueType value, std::string name) { add(value, std::move(name), true); }
  void addOnce(ValueType value, std::string name) { add(value, std::move(name), false); }
  void add(ValueType value, const std::string& name, bool save_repeat = false);

  /// Write stored data to file immediately
  bool write();

private:
  /// Helper struct to save enough information so that we can save an
  /// object to file later
  struct StoredValue {
    StoredValue(std::string name_, ValueType value_, bool repeat_)
        : name(std::move(name_)), value(value_), repeat(repeat_) {}
    std::string name;
    ValueType value;
    bool repeat;
  };
  std::vector<StoredValue> data;

public:
  const std::vector<StoredValue>& getData() { return data; }
};

struct OptionsConversionVisitor {
private:
  Options& options;
  std::string name;
public:
  OptionsConversionVisitor(Options& options_, std::string name_)
      : options(options_), name(std::move(name_)) {}
  template <class T>
  void operator()(T* value) {
    options[name] = *value;
  }
};
} // namespace bout

/*!
  Base class for physics models
 */
class PhysicsModel {
public:
  using preconfunc = int (PhysicsModel::*)(BoutReal t, BoutReal gamma, BoutReal delta);
  using jacobianfunc = int (PhysicsModel::*)(BoutReal t);

  template <class Model, typename = typename std::enable_if_t<
                             std::is_base_of<PhysicsModel, Model>::value>>
  using ModelPreconFunc = int (Model::*)(BoutReal t, BoutReal gamma, BoutReal delta);
  template <class Model, typename = typename std::enable_if_t<
                             std::is_base_of<PhysicsModel, Model>::value>>
  using ModelJacobianFunc = int (Model::*)(BoutReal t);

  PhysicsModel();
  PhysicsModel(Mesh* mesh_, const bout::OptionsNetCDF& output_, bool output_enabled_,
               const bout::OptionsNetCDF& restart_, bool restart_enabled_)
      : mesh(mesh_), output_file(output_), output_enabled(output_enabled_),
        restart_file(restart_), restart_enabled(restart_enabled_) {}

  virtual ~PhysicsModel() = default;
  
  Mesh* mesh{nullptr};
  bout::DataFileFacade dump{};
  bout::DataFileFacade restart{};

  /*!
   * Initialse the model, calling the init() and postInit() methods
   *
   * Note: this is usually only called by the Solver
   */
  void initialise(Solver* s);

  /*!
   * Run the RHS function, to calculate the time derivatives
   *
   * Input
   * -----
   *
   * @param[in] time  The simulation time
   *
   * The system state should be in the evolving variables
   *
   * Output
   * ------
   * 
   * The time derivatives will be put in the ddt() variables
   * 
   * Returns a flag: 0 indicates success, non-zero an error flag
   */
  int runRHS(BoutReal time);
  
  /*!
   * True if this model uses split operators
   */ 
  bool splitOperator();
  
  /*!
   * Run the convective (usually explicit) part of the model
   */
  int runConvective(BoutReal time);
  
  /*!
   * Run the diffusive (usually implicit) part of the model
   */
  int runDiffusive(BoutReal time, bool linear);
  
  /*!
   * True if a preconditioner has been defined
   */ 
  bool hasPrecon();
  
  /*!
   * Run the preconditioner. The system state should be in the 
   * evolving variables, and the vector to be solved in the ddt() variables.
   * The result will be put in the ddt() variables.
   *
   * Note: this is usually only called by the Solver
   *
   */
  int runPrecon(BoutReal t, BoutReal gamma, BoutReal delta);
  
  /*!
   * True if a Jacobian function has been defined
   */
  bool hasJacobian();
  
  /*!
   * Run the Jacobian-vector multiplication function
   * 
   * Note: this is usually only called by the Solver
   */ 
  int runJacobian(BoutReal t);

  int runTimestepMonitor(BoutReal simtime, BoutReal dt) {return timestepMonitor(simtime, dt);}
  
  /// Write \p options to `output_file`
  void writeOutputFile(const Options& options);
  /// Write variables with \p time_dimension from \p options to `output_file`
  void writeOutputFile(const Options& options, const std::string& time_dimension);

  /// Finish the output for this timestep, verifying all evolving
  /// variables have the correct length
  void finishOutputTimestep() const;
protected:
  
  // The init and rhs functions are implemented by user code to specify problem
  /*!
   * @brief This function is called once by the solver at the start of a simulation.
   * 
   * A valid PhysicsModel must implement this function
   * 
   * Variables should be read from the inputs, and the variables to 
   * be evolved should be specified.
   */
  virtual int init(bool restarting) = 0;
  
  /// Post-initialise. This reads the restart file
  ///
  /// @param[in] restarting   If true, will load state from restart file
  ///
  virtual int postInit(bool restarting);
  
  /*!
   * @brief This function is called by the time integration solver
   * at least once per time step
   * 
   * Variables being evolved will be set by the solver
   * before the call, and this function must calculate
   * and set the time-derivatives.
   *
   * By default this function just returns an error,
   * which will stop the simulation.
   */
  virtual int rhs(BoutReal UNUSED(t)) {return 1;}

  /// Output additional variables other than the evolving variables
  virtual void outputVars(Options& options);
  /// Add additional variables other than the evolving variables to the restart files
  virtual void restartVars(Options& options);

  /* 
     If split operator is set to true, then
     convective() and diffusive() are called instead of rhs()
     
     For implicit-explicit schemes, convective() will typically
     be treated explicitly, whilst diffusive() will be treated implicitly.
     For unsplit methods, both convective and diffusive will be called
     and the sum used to evolve the system:
     rhs() = convective() + diffusive()
   */
  virtual int convective(BoutReal UNUSED(t)) {return 1;}
  virtual int diffusive(BoutReal UNUSED(t)) {return 1;}
  virtual int diffusive(BoutReal t, bool UNUSED(linear)) { return diffusive(t); }
  
  /*!
   * Implemented by user code to monitor solution at output times
   */
  virtual int outputMonitor(BoutReal UNUSED(simtime), int UNUSED(iter), int UNUSED(NOUT)) {
    return 0;
  }
  
  /*!
   * Timestep monitor. If enabled by setting solver:monitor_timestep=true
   * then this function is called every internal timestep.
   */
  virtual int timestepMonitor(BoutReal UNUSED(simtime), BoutReal UNUSED(dt)) {return 0;}

  

  // Functions called by the user to set callback functions

  /// Specify that this model is split into a convective and diffusive part
  void setSplitOperator(bool split=true) {splitop = split;}

  /// Specify a preconditioner function
  void setPrecon(preconfunc pset) { userprecon = pset; }
  template <class Model>
  void setPrecon(ModelPreconFunc<Model> preconditioner) {
    userprecon = static_cast<preconfunc>(preconditioner);
  }

  /// Specify a Jacobian-vector multiply function
  void setJacobian(jacobianfunc jset) { userjacobian = jset; }
  template <class Model>
  void setJacobian(ModelJacobianFunc<Model> jacobian) {
    userjacobian = static_cast<jacobianfunc>(jacobian);
  }

  /// This is set by a call to initialise, and can be used by models to specify evolving variables
  Solver* solver{nullptr};

  /*!
   * Specify a variable for the solver to evolve
   *
   * @param[in] var  The variable to evolve
   * @param[in] name The name to use for variable initialisation and output
   * 
   * Note that the variable must not be destroyed (e.g. go out of scope)
   * after this call, since a pointer to \p var is stored in the solver.
   *
   * To evolve the state, the solver will set \p var, and the user-supplied
   * rhs() function should calculate ddt(var).
   */
  void bout_solve(Field2D &var, const char *name, const std::string& description="");
  void bout_solve(Field3D &var, const char *name, const std::string& description="");
  void bout_solve(Vector2D &var, const char *name, const std::string& description="");
  void bout_solve(Vector3D &var, const char *name, const std::string& description="");

  /// Helper function for reading from restart_options
  Options& readFromRestartFile(const std::string& name) {
    return restart_options[name];
  }

  /// Write the restart file to disk now
  void writeRestartFile();
  /// Write the output file to disk now
  void writeOutputFile();

  /*!
   * Specify a constrained variable \p var, which will be
   * adjusted to make \p F_var equal to zero.
   * If the solver does not support constraints then this will throw an exception
   * 
   * @param[in] var  The variable the solver should modify
   * @param[in] F_var  The control variable, which the user will set
   * @param[in] name   The name to use for initialisation and output
   * 
   */ 
  bool bout_constrain(Field3D &var, Field3D &F_var, const char *name);

  /*!
   * Monitor class for PhysicsModel
   */
  class PhysicsModelMonitor : public Monitor {
  public:
    PhysicsModelMonitor() = delete;
    PhysicsModelMonitor(PhysicsModel *model) : model(model) {}
    int call(Solver* solver, BoutReal simtime, int iter, int nout) override;

  private:
    PhysicsModel *model;
  };

private:
  /// State for outputs
  Options output_options;
  /// File to write the outputs to
  bout::OptionsNetCDF output_file;
  /// Should we write output files
  bool output_enabled{true};
  /// Stores the state for restarting
  Options restart_options;
  /// File to write the restart-state to
  bout::OptionsNetCDF restart_file;
  /// Should we write restart files
  bool restart_enabled{true};
  /// Split operator model?
  bool splitop{false};
  /// Pointer to user-supplied preconditioner function
  preconfunc userprecon{nullptr};
  /// Pointer to user-supplied Jacobian-vector multiply function
  jacobianfunc userjacobian{nullptr};
  /// True if model already initialised
  bool initialised{false};
  /// write restarts and pass outputMonitor method inside a Monitor subclass
  PhysicsModelMonitor modelMonitor{this};
};

/*!
 * Macro to define a simple main() which creates
 * the given model and runs it. This should be sufficient
 * for most use cases, but a user can define their own
 * main() function if needed.
 *
 * Example
 * -------
 *
 * class MyModel : public PhysicsModel {
 *   ..
 * };
 *
 * BOUTMAIN(MyModel);
 */
#define BOUTMAIN(ModelClass)                                       \
  int main(int argc, char** argv) {                                \
    int init_err = BoutInitialise(argc, argv);                     \
    if (init_err < 0) {                                            \
      return 0;                                                    \
    }                                                              \
    if (init_err > 0) {                                            \
      return init_err;                                             \
    }                                                              \
    try {                                                          \
      auto model = bout::utils::make_unique<ModelClass>();         \
      auto solver = Solver::create();                              \
      solver->setModel(model.get());                               \
      auto bout_monitor = bout::utils::make_unique<BoutMonitor>(); \
      solver->addMonitor(bout_monitor.get(), Solver::BACK);        \
      solver->solve();                                             \
    } catch (const BoutException& e) {                             \
      output << "Error encountered\n";                             \
      output << e.getBacktrace() << endl;                          \
      MPI_Abort(BoutComm::get(), 1);                               \
    }                                                              \
    BoutFinalise();                                                \
    return 0;                                                      \
  }

/// Macro to replace solver->add, passing variable name
#define SOLVE_FOR1(var) solver->add(var, #var);
#define SOLVE_FOR2(var1, var2) { \
  solver->add(var1, #var1);       \
  solver->add(var2, #var2);}
#define SOLVE_FOR3(var1, var2, var3) { \
  solver->add(var1, #var1);             \
  solver->add(var2, #var2);             \
  solver->add(var3, #var3);}
#define SOLVE_FOR4(var1, var2, var3, var4) { \
  solver->add(var1, #var1);             \
  solver->add(var2, #var2);             \
  solver->add(var3, #var3);             \
  solver->add(var4, #var4);}
#define SOLVE_FOR5(var1, var2, var3, var4, var5) { \
  solver->add(var1, #var1);             \
  solver->add(var2, #var2);             \
  solver->add(var3, #var3);             \
  solver->add(var4, #var4);             \
  solver->add(var5, #var5);}
#define SOLVE_FOR6(var1, var2, var3, var4, var5, var6) { \
  solver->add(var1, #var1);             \
  solver->add(var2, #var2);             \
  solver->add(var3, #var3);             \
  solver->add(var4, #var4);             \
  solver->add(var5, #var5);             \
  solver->add(var6, #var6);}

/// Add fields to the solver.
/// This should accept up to ten arguments
#define SOLVE_FOR(...)                  \
  { MACRO_FOR_EACH(SOLVE_FOR1, __VA_ARGS__) }

#endif // __PHYSICS_MODEL_H__

