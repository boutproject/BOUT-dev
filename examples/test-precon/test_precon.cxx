/*
 * Test simple implicit preconditioning
 * 
 *
 * PETSc flags for testing:
 * solver_type=petsc -ts_type theta -ts_theta_theta 0.5 -{ksp,snes,ts}_monitor
 */

#include <bout.hxx>
#include <boutmain.hxx>
#include <invert_parderiv.hxx>

int precon(BoutReal t, BoutReal cj, BoutReal delta); // Preconditioner
int jacobian(BoutReal t); // Jacobian-vector multiply

Field3D u, v; // Evolving variables

InvertPar *inv; // Parallel inversion class

int physics_init(bool restarting) {
  // Set variables to evolve
  SOLVE_FOR2(u,v);
  
  // Give the solver the preconditioner function
  solver->setPrecon(precon);
  
  // Set Jacobian
  solver->setJacobian(jacobian);
  
  // Initialise parallel inversion class
  inv = InvertPar::Create();
  inv->setCoefA(1.0);
  
  return 0;
}

int physics_run(BoutReal t) {
  mesh->communicate(u,v);
  
  ddt(u) = Grad_par(v);
  ddt(v) = Grad_par(u);
  
  return 0;
}

/*********************************************************
 * Preconditioner
 *
 * o System state in variables (as in rhs function)
 * o Values to be inverted in F_vars
 *
 * o Return values should be in vars (overwriting system state)
 * 
 *********************************************************/
int precon(BoutReal t, BoutReal gamma, BoutReal delta) {
  // Communicate vector to be inverted
  mesh->communicate(ddt(u), ddt(v));
  
  //output << "t = " << t << " Gamma = " << gamma << endl;
  
  // First matrix 
  //  | I   -U |
  //  | 0    I |
  u = ddt(u) + gamma*Grad_par(ddt(v));
  v = ddt(v);
  
  // Second matrix, containing Schur complement
  // | (1 - UL)^-1  0 | 
  // |   0          I |
  inv->setCoefB(-SQ(gamma));
  u = inv->solve(u);
  
  // Third matrix
  // |  I  0 |
  // | -L  I |
  mesh->communicate(u);
  v = gamma*Grad_par(u) + v;

  u.applyBoundary("dirichlet");
  v.applyBoundary("dirichlet");
  
  return 0;
}

/*********************************************************
 * Jacobian
 * 
 * o System state in variables (as in RHS function)
 * o Vector to be multiplied is in time derivatives
 *
 * o Output Jacobian-vector multiplied Jv should be in variables
 *
 * enable by setting solver / use_jacobian = true in BOUT.inp
 *********************************************************/

int jacobian(BoutReal t) {
  mesh->communicate(ddt(u), ddt(v));
  u = Grad_par(ddt(v));
  v = Grad_par(ddt(u));
  u.applyBoundary("dirichlet");
  v.applyBoundary("dirichlet");
  return 0;
}
