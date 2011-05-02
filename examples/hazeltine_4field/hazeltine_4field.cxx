/**************************************************************
 * 4-field model by Hazeltine, Kotschenreuther and Morrison
 * Phys. Fluids vol 28 (8) 1985, p2472
 *
 * B.Dudson, University of York, Aug 2010 
 **************************************************************/

#include <bout.hxx>
#include <boutmain.hxx>

// Evolving fields
Field3D U, Psi, v, p;

// Auxilliary variables
Field3D phi, J;

// Resistivity
Field2D eta;

// Equilibrium quantities
Field2D J0;     // Parallel current
Field2D P0;     // Pressure
Vector2D b0xcv; // Curvature

// Factors
BoutReal epsilon, delta, beta, tau; // tau = Ti / Te



int physics_init(bool restarting)
{
  // Load 2D profiles
  mesh->get(J0, "Jpar0");    // A / m^2
  mesh->get(P0, "pressure"); // Pascals

  // Load curvature term
  b0xcv.covariant = false; // Read contravariant components
  mesh->get(b0xcv, "bxcv"); // mixed units x: T y: m^-2 z: m^-2
  
  // Load metrics
  if(mesg->get(Rxy,  "Rxy")) { // m
    output.write("Error: Cannot read Rxy from grid\n");
    return 1;
  }
  if(mesh->get(Bpxy, "Bpxy")) { // T
    output.write("Error: Cannot read Bpxy from grid\n");
    return 1;
  }
  mesh->get(Btxy, "Btxy"); // T
  mesh->get(B0,   "Bxy");  // T
  mesh->get(hthe, "hthe"); // m
  mesh->get(I,    "sinty");// m^-2 T^-1
  
  // Read options
  options.setSection("4field");
  
  
  if(ShiftXderivs) {
    if(IncIntShear) {
      // BOUT-06 style, using d/dx = d/dpsi + I * d/dz
      IntShiftTorsion = I;
      
    }else {
      // Dimits style, using local coordinate system
      b0xcv.z += I*b0xcv.x;
      I = 0.0;  // I disappears from metric
    }
  }

  if(mesh->get(Bbar, "bmag")) // Typical magnetic field
    Bbar = 1.0;
  
  
  
  
  Va = sqrt(Bbar*Bbar / (MU0*density*Mi));
  
  

  bout_solve(U, "U");
  bout_solve(Psi, "Psi");
  bout_solve(v, "v");
  bout_solve(p, "p");
  
  J.setBoundary("J");
  phi.setBoundary("phi");
  
  dump.add(J, "J", 1);
  dump.add(phi, "phi", 1);
  
  return 0;
}

const Field3D Grad_parP(const Field3D &f)
{
  return Grad_par(f) - bracket(Psi, f);
}

int physics_run(BoutReal t)
{
  // Invert vorticity to get phi
  phi = invert_laplace(U, phi_flags, NULL);
  phi.applyBoundary();
  
  // Communicate variables
  mesh->communicate(U, Psi, v, p, phi);
  
  // Calculate J from Psi
  J = Delp2(Psi);
  J.applyBoundary();
  mesh->communicate(J);
  
  // Vorticity
  
  // Note: Implementing [Delp phi ; Delp p]
  Field3D delp2p = Delp2(p);
  Field3D brackpphi = bracket(p, phi);
  mesh->communicate(delp2p, brackpphi);
  
  ddt(U) = 
    -bracket(phi + delta*tau*p, U)
    - \Grad_parP(J)
    + 0.5*delta*tau*( 
                     bracket(p, U)
                     - bracket(phi, delp2p)
                     - Delp2(brackpphi)
                      );
  
  // Parallel electric field
  ddt(Psi) = 
    -Grad_parP(phi)
    + eta * J
    + delta*Grad_parP(p);

  // Parallel velocity
  ddt(v) = 
    -bracket(phi, p) 
    - 0.5*(1.+tau)*Grad_parP(p)
    + delta*beta*tau*(
                      0.5*(1.+tau)*bracket(p,v)
                      );
  
  // Pressure
  ddt(p) = 
    -bracket(phi, p) 
    + beta*( 
            -Grad_parP(v + 2.*delta*J) 
            + 0.5*(1.+tau)*eta*Delp2(p)
             );
    
  return 0;
}
