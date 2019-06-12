/*
  Testing use of the solver class for custom integration
  
  Illustrates how to use the Solver class as an ODE integrator
 */

#include <bout/physicsmodel.hxx>  // Note: Need to use new API


// This class represents a sub-problem to be solved
class MyFunction : public PhysicsModel {
public:
  int init(bool UNUSED(restarting)) {
    solver->add(result, "result");
    return 0;
  }
  
  int rhs(BoutReal UNUSED(time)) {
    ddt(result) = 1.0;
    return 0;
  }
  
  int outputMonitor(BoutReal simtime, int UNUSED(iter), int UNUSED(NOUT)) {
    output.write("MyFunction: time = %e\n", simtime);
    return 0;
  }
  
  Field3D result;
private:
  
};


// This class represents the top-level model being solved
class TestIntegrate : public PhysicsModel {
public:
  ~TestIntegrate() {
    delete model;
    delete ode;
  }
  
  int init(bool UNUSED(restarting)) {
    
    // Create a model
    model = new MyFunction();
    
    // Create a solver, passing the options section "ode"
    ode = Solver::create(Options::getRoot()->getSection("ode"));
    
    // Specify the model to be solved
    ode->setModel(model);
    ode->solve(5, 0.1); // Number of outputs, step
    
    solver->add(f, "f");
    return 0;
  }
  
  int rhs(BoutReal UNUSED(time)) {
    ddt(f) = model->result;
    return 0;
  }
private:
  
  Field3D f; // Some variable being evolved
  
  Solver *ode; // Integration solver
  MyFunction *model;
};


BOUTMAIN(TestIntegrate);
