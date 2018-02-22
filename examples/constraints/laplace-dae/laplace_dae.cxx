/*************************************************************
 * 
 * 
 *************************************************************/

#include <bout.hxx>
#include <boutmain.hxx>
#include <initialprofiles.hxx>
#include <invert_laplace.hxx>
#include <boutexception.hxx>

Field3D U, Apar;   // Evolving variables

Field3D phi;  // Electrostatic potential: Delp2(U) = phi
Field3D jpar; // Parallel current: Delp2(Apar) = jpar
Field3D phibdry; // Used for calculating error in the boundary

bool constraint;

int flags;

// Preconditioner
int precon_phi(BoutReal t, BoutReal cj, BoutReal delta);
int jacobian(BoutReal t); // Jacobian-vector multiply
int jacobian_constrain(BoutReal t); // Jacobian-vector multiply

int physics_init(bool restarting) {
  // Give the solver two RHS functions
  
  // Get options
  Options *options = Options::getRoot();
  options = options->getSection("dae");
  OPTION(options, constraint, true);
  OPTION(options, flags, 0);
    
  // Just solving one variable, U
  SOLVE_FOR2(U, Apar);
  
  if(constraint) {
    phi = invert_laplace(U, flags);
    // Add phi equation as a constraint
    if(!bout_constrain(phi, ddt(phi), "phi"))
      throw BoutException("Solver does not support constraints");
    
    // Set preconditioner
    solver->setPrecon(precon_phi);
    
    // Set Jacobian
    solver->setJacobian(jacobian_constrain);
    
    phibdry.setBoundary("phi");
  }else {
    // Save phi to file every timestep
    SAVE_REPEAT(phi);
    phi.setBoundary("phi");
    
    // Set Jacobian
    solver->setJacobian(jacobian);
  }
  
  SAVE_REPEAT(jpar);
  jpar.setBoundary("jpar");

  return 0;
}

int physics_run(BoutReal time) {

  if(constraint) {
    mesh->communicate(Apar, phi);
    
    // phi is solved as a constraint (sparse Jacobian)
    // Calculate the error, and return in ddt(phi)
    ddt(phi) = Delp2(phi) - U;
    //mesh->communicate(ddt(phi));
    
    // Now the error in the boundary (quite inefficient)
    phibdry = phi;
    phibdry.applyBoundary();
    phibdry -= phi; // Contains error in the boundary
    
    ddt(phi).setBoundaryTo(phibdry);
    
  }else {
    mesh->communicate(U, Apar);
    
    // Solving for phi here (dense Jacobian)
    output << "U " << max(U) << endl;
    phi = invert_laplace(U, flags);
    phi.applyBoundary();
  }
  
  jpar = Delp2(Apar);
  jpar.applyBoundary();
  mesh->communicate(jpar, phi);
  
  output << "phi " << max(phi) << endl;
  
  for(int y=0;y<5;y++) {
    for(int x=0;x<5;x++)
      output << phi(x,y,64) << ", ";
    output << endl;
  }
  

  ddt(U) = Grad_par(jpar);
  ddt(Apar) = Grad_par(phi);
  
  return 0;
}

/*******************************************************************************
 * Preconditioner for when phi solved as a constraint
 * Currently only possible with the IDA solver
 *
 * o System state in variables (as in rhs function)
 * o Values to be inverted in time derivatives
 * 
 * o Return values should be in time derivatives
 *******************************************************************************/

int precon_phi(BoutReal t, BoutReal cj, BoutReal delta) {
  // Not preconditioning U or Apar equation
  
  ddt(phi) = invert_laplace(ddt(phi) - ddt(U), flags);
  
  return 0;
}

/*******************************************************************************
 * Jacobian-vector multiply
 *
 * Input
 *   System state is in variables
 *   Vector v is in time-derivatives
 * Output
 *   Jacobian-vector multiplied Jv should be in time derivatives
 * 
 *******************************************************************************/


/// Jacobian when solving phi in RHS
int jacobian(BoutReal t) {
  Field3D Jphi = invert_laplace(ddt(U), flags); // Inversion makes this dense
  mesh->communicate(Jphi, ddt(Apar));
  Field3D Jjpar = Delp2(ddt(Apar));
  mesh->communicate(Jjpar);
  
  ddt(U) = Grad_par(Jjpar);  // Dense matrix in evolving U
  ddt(Apar) = Grad_par(Jphi);
  
  return 0;
}

/// Jacobian when solving phi as a constraint.
/// No inversion, only sparse Delp2 and Grad_par operators 
int jacobian_constrain(BoutReal t) {
  
  mesh->communicate(ddt(Apar), ddt(phi));
  Field3D Jjpar = Delp2(ddt(Apar));
  mesh->communicate(Jjpar);
  
  U    = Grad_par(Jjpar);
  Apar = Grad_par(ddt(phi));
  
  phi  = Delp2(ddt(phi)) - ddt(U);
  
  phibdry = ddt(phi);
  phibdry.applyBoundary();
  phibdry -= ddt(phi); // Contains error in the boundary
  
  phi.setBoundaryTo(phibdry);
  
  ddt(phi) = phi;
  ddt(U) = U;
  ddt(Apar) = Apar;

  return 0;
}
