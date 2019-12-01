#include <cbout/bout.h>

int rhs(Options* ddt, Options* state, BoutReal time) {
  Field3D *f;
  Field3D_create(&f);

  /* Get values from the current state */
  Options_get_Field3D(&f, Options_get(state, "f"));

  /* Calculate time derivatives */
  Field3D_scale(f, -1.0);

  /* Set values in the time derivatives */
  Options_set_Field3D(Options_get(ddt, "f"), f);
  
  Field3D_delete(f);  
  return 0;
}

int main(int argc, char** argv) {
  /* Start BOUT++ */
  bout_initialise(argc, argv);

  /* Create a physics model */
  PhysicsModel* model;
  PhysicsModel_create(&model);
  PhysicsModel_evolve(model, "f");
  PhysicsModel_set_rhs(model, &rhs);

  /* Create a solver */
  Solver *solver;
  Solver_create(&solver);
  Solver_add_model(solver, model);

  /* Run simulation */
  Solver_solve(solver);
  
  /* Free objects */
  Solver_delete(solver);
  PhysicsModel_delete(model);
  
  /* Shut down BOUT++ */
  bout_finalise();
  return 0;
}
