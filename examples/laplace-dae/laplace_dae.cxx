/*************************************************************
 * 
 * 
 *************************************************************/

#include <bout.hxx>
#include <boutmain.hxx>
#include <initialprofiles.hxx>
#include <invert_laplace.hxx>
#include <boutexception.hxx>

Field3D U;   // Evolving variable

Field3D phi; // Potential used for advection
Field3D phibdry; // Used for calculating error in the boundary

bool constraint;

int flags;

// Preconditioner
int precon_phi(BoutReal t, BoutReal cj, BoutReal delta);

int physics_init(bool restarting) {
  // Give the solver two RHS functions
  
  // Get options
  Options *options = Options::getRoot();
  options = options->getSection("dae");
  OPTION(options, constraint, true);
  OPTION(options, flags, 0);
    
  // Just solving one variable, U
  SOLVE_FOR(U);
  
  if(constraint) {
    phi = invert_laplace(U, flags);
    // Add phi equation as a constraint
    if(!bout_constrain(phi, ddt(phi), "phi"))
      throw BoutException("Solver does not support constraints");
    // Set preconditioner
    solver->setPrecon(precon_phi);
    
    phibdry.setBoundary("phi");
  }else {
    // Save phi to file every timestep
    SAVE_REPEAT(phi);
  }

  return 0;
}

int physics_run(BoutReal time) {

  if(constraint) {
    // Need communication
    mesh->communicate(U, phi);
    
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
    // Need communication
    mesh->communicate(U);
    
    // Solving for phi here (dense Jacobian)
    phi = invert_laplace(U, flags);
  }
  
  // Form of advection operator for reduced MHD type models
  ddt(U) = -bracket(phi, U, BRACKET_SIMPLE);
  
  return 0;
}

/*******************************************************************************
 * Preconditioner for when phi solved as a constraint
 * Currently only possible with the IDA solver
 *
 * o System state in variables (as in rhs function)
 * o Values to be inverted in F_vars
 * 
 * o Return values should be in vars (overwriting system state)
 *******************************************************************************/

int precon_phi(BoutReal t, BoutReal cj, BoutReal delta) {
  // Not preconditioning U equation
  U = ddt(U);
  
  phi = invert_laplace(ddt(phi) - ddt(U), flags);
  
  return 0;
}
