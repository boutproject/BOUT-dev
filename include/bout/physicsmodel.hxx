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

#include <bout.hxx>
#include <options.hxx>
#include <msg_stack.hxx>
#include "solver.hxx"
#include "unused.hxx"

/*!
  Base class for physics models
 */
class PhysicsModel {
public:
  typedef int (PhysicsModel::*preconfunc)(BoutReal t, BoutReal gamma, BoutReal delta);
  typedef int (PhysicsModel::*jacobianfunc)(BoutReal t);
  
  PhysicsModel();
  
  virtual ~PhysicsModel();
  
  /*!
   * Initialse the model, calling the init() and postInit() methods
   *
   * Note: this is usually only called by the Solver
   */
  void initialise(Solver *s) {
    if(initialised)
      return; // Ignore second initialisation
    initialised = true;
    
    // Restart option
    bool restarting;
    Options::getRoot()->get("restart", restarting, false);
    
    // Set the protected variable, so user code can
    // call the solver functions
    solver = s;

    // Call user init code to specify evolving variables
    if ( init(restarting) ) {
      throw BoutException("Couldn't initialise physics model");
    }
    
    // Post-initialise, which reads restart files
    // This function can be overridden by the user
    if (postInit(restarting)) {
      throw BoutException("Couldn't restart physics model");
    }
  }
  
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
  void setPrecon(preconfunc pset) {userprecon = pset;}

  /// Specify a Jacobian-vector multiply function
  void setJacobian(jacobianfunc jset) {userjacobian = jset;}

  /// This is set by a call to initialise, and can be used by models to specify evolving variables
  Solver *solver;

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
  void bout_solve(Field2D &var, const char *name);
  void bout_solve(Field3D &var, const char *name);
  void bout_solve(Vector2D &var, const char *name);
  void bout_solve(Vector3D &var, const char *name);

  /// Stores the state for restarting
  Datafile restart; 

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
    int call(Solver* UNUSED(solver), BoutReal simtime, int iter, int nout) {
      // Save state to restart file
      model->restart.write();
      // Call user output monitor
      return model->outputMonitor(simtime, iter, nout);
    }
  private:
    PhysicsModel *model;
  };

  /// write restarts and pass outputMonitor method inside a Monitor subclass
  PhysicsModelMonitor modelMonitor;
private:
  bool splitop; ///< Split operator model?
  preconfunc   userprecon; ///< Pointer to user-supplied preconditioner function
  jacobianfunc userjacobian; ///< Pointer to user-supplied Jacobian-vector multiply function
  
  bool initialised; ///< True if model already initialised
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
#define BOUTMAIN(ModelClass)                          \
  int main(int argc, char **argv) {                   \
    int init_err = BoutInitialise(argc, argv);        \
    if (init_err < 0)				      \
      return 0;                                       \
    else if (init_err > 0) 			      \
      return init_err;				      \
    try {                                             \
      ModelClass *model = new ModelClass();           \
      Solver *solver = Solver::create();              \
      solver->setModel(model);                        \
      Monitor * bout_monitor = new BoutMonitor();     \
      solver->addMonitor(bout_monitor, Solver::BACK); \
      solver->outputVars(dump);                       \
      solver->solve();                                \
      delete model;                                   \
      delete solver;                                  \
      delete bout_monitor;                            \
    }catch (BoutException &e) {                       \
      output << "Error encountered\n";                \
      output << e.what() << endl;                     \
      MPI_Abort(BoutComm::get(), 1);                  \
    }                                                 \
    BoutFinalise();                                   \
    return 0;                                         \
  }

/// Macro to replace solver->add, passing variable name
#define SOLVE_FOR(var) solver->add(var, #var)
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

#endif // __PHYSICS_MODEL_H__

