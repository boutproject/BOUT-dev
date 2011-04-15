/**************************************************************************
 * Similar set of equations to JOREK
 * 
 **************************************************************************/

#include "bout.h"

// Evolving quantities
Field3D rho, T, u, vpar, psi;
Field3D phi, Jpar;

Field2D rho0;
Field2D B0;
Vector2D B0vec;

Vector3D vExB, vD; // Velocities

// options

bool nonlinear;
bool full_bfield; // If true, use divergence-free expression for B

int physics_init(bool restarting) {
  
  Options *globalOptions = Options::getRoot();
  Options *options = globalOptions->getSection("jorek");
  
  OPTION(options, nonlinear,           false);
  OPTION(options, full_bfield,         false);
  
  // Set B field vector
  B0vec.covariant = false;
  B0vec.x = 0.;
  B0vec.y = Bpxy / hthe;
  B0vec.z = 0.;
  
  vExB.setBoundary("v");
  vD.setBoundary("v");
}

// Parallel gradient along perturbed field-line
const Field3D Grad_parP(const Field3D &f, CELL_LOC loc = CELL_DEFAULT) {
  // Derivative along equilibrium field-line
  Field3D result = Grad_par(f, loc);
  
  if(nonlinear) {
    if(full_bfield) {
      // Use full expression for perturbed B
      Vector3D Btilde = Curl(B0vec * psi);
      result += Btilde * Grad(f) / B0;
    }else {
      // Simplified expression
      result -= b0xGrad_dot_Grad(psi, f);
    }
  }
  return result;
}

int physics_run(BoutReal t) {
  
  // Invert laplacian for phi
  phi = invert_laplace(U, phi_flags, NULL);
  // Apply a boundary condition on phi for target plates
  phi.applyBoundary();
  
  // Communicate variables
  mesh->communicate(rho, T, u, phi, vpar, psi);
  
  // Get J from Psi
  Jpar = Delp2(psi);
  Jpar.applyBoundary();
  mesh->communicate(Jpar);
  
  ddt(u) = (B0^2) * b0xGrad_dot_Grad(psi, J0, CELL_CENTRE) 
    - (B0^2)*Grad_parP(Jpar, CELL_CENTRE)
    + b0xcv*Grad(P);
  
  if(nonlinear) {
    ddt(u) -= b0xGrad_dot_Grad(phi, u);    // Advection
  }
  
  Field3D rhot = rho0;
  Field3D Tt = T0;
  if(nonlinear) {
    rhot += rho;
    Tt += T;
  }

  if(flux_methods) {
    // ExB velocity
    vExB = (B0vec ^ Grad_perp(phi))/(B0*B0);
    vExB.applyBoundary();
    
    ////////// Density equation ////////////////
    
    // Diffusive flux (perpendicular)
    vD = Grad_perp(rho);
    vD.applyBoundary();
    
    if(nonlinear) {
      ddt(rho) = -Div(vExB + vD, rho0 + rho);
    }else {
      ddt(rho) = -Div(vExB + vD, rho0);
    }
    
    ////////// Temperature equation ////////////
  
    ddt(T) = -b0xGrad_dot_Grad(phi, T);
    
  }else {
    // Use analytic expressions, expand terms
    
    // Divergence of ExB velocity (neglecting parallel term)
    Field3D divExB = b0xcv*Grad(phi)/B0 - b0xGrad_dot_Grad(1./B0, phi);
    
    ddt(rho) = -b0xGrad_dot_Grad(phi, rhot) // Advection 
      - divExB*rhot; // Compression
    
    
    ddt(T) = -b0xGrad_dot_Grad(phi, Tt) 
      - (2./3.)*Tt*divExB;
  }
  
}
