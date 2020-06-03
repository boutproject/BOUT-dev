/*
 * Test simple implicit preconditioning
 * 
 *
 * PETSc flags for testing:
 * solver:type=petsc -ts_type theta -ts_theta_theta 0.5 -{ksp,snes,ts}_monitor
 */

#include <bout.hxx>
#include <boutmain.hxx>
#include <invert_parderiv.hxx>

int precon(BoutReal t, BoutReal cj, BoutReal delta); // Preconditioner
int jacobian(BoutReal t); // Jacobian-vector multiply

Field3D u, v; // Evolving variables

std::unique_ptr<InvertPar> inv{nullptr}; // Parallel inversion class

int physics_init(bool UNUSED(restarting)) {
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

int physics_run(BoutReal UNUSED(t)) {
  u.getMesh()->communicate(u, v);

  ddt(u) = Grad_par(v);
  ddt(v) = Grad_par(u);
  
  return 0;
}

/*********************************************************
 * Preconditioner
 *
 * o System state in variables (as in rhs function)
 * o Values to be inverted in time derivatives
 *
 * o Return values should be in time derivatives
 * 
 *********************************************************/
int precon(BoutReal UNUSED(t), BoutReal gamma, BoutReal UNUSED(delta)) {
  auto* mesh = u.getMesh();

  // Communicate vector to be inverted
  mesh->communicate(ddt(u), ddt(v));
  
  //output << "t = " << t << " Gamma = " << gamma << endl;
  
  // First matrix 
  //  | I   -U |
  //  | 0    I |
  ddt(u) = ddt(u) + gamma*Grad_par(ddt(v));
  //ddt(v) = ddt(v);
  
  // Second matrix, containing Schur complement
  // | (1 - UL)^-1  0 | 
  // |   0          I |
  inv->setCoefB(-SQ(gamma));
  ddt(u) = inv->solve(ddt(u));
  
  // Third matrix
  // |  I  0 |
  // | -L  I |
  mesh->communicate(ddt(u));
  ddt(v) = gamma*Grad_par(ddt(u)) + ddt(v);

  ddt(u).applyBoundary("dirichlet");
  ddt(v).applyBoundary("dirichlet");
  
  return 0;
}

/*********************************************************
 * Jacobian
 * 
 * o System state in variables (as in RHS function)
 * o Vector to be multiplied is in time derivatives
 *
 * o Output Jacobian-vector multiplied Jv should be in time derivatives
 *
 * enable by setting solver / use_jacobian = true in BOUT.inp
 *********************************************************/

int jacobian(BoutReal UNUSED(t)) {
  auto* mesh = u.getMesh();

  mesh->communicate(ddt(u), ddt(v));
  Field3D utmp = Grad_par(ddt(v)); // Shouldn't overwrite ddt(u) before using it
  ddt(v) = Grad_par(ddt(u));
  ddt(u) = utmp;
  ddt(u).applyBoundary("dirichlet");
  ddt(v).applyBoundary("dirichlet");
  return 0;
}

