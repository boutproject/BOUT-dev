/**************************************************************************
 * Similar set of equations to JOREK
 * 
 **************************************************************************/

#include <bout.h>
#include <invert_laplace.h>

// Evolving quantities
Field3D rho, T, u, vpar, psi;
// Derived quantities
Field3D Jpar, phi; // Parallel current, electric potential

Field2D D_perp, chi_perp, chi_par; // Particle and heat diffusion coefficients
Field2D eta0;  // Resistivity
Field3D eta;
BoutReal viscos_par, viscos_perp; // Viscosity coefficients

Field2D rho0, T0; // Equilibrium mass density and temperature
Field2D B0, J0, P0;
Vector2D b0xcv; // Curvature term
Vector2D B0vec; // B0 field vector

int phi_flags;

const BoutReal MU0 = 4.0e-7*PI;
const BoutReal Charge = 1.602e-19; // electron charge e (C)
const BoutReal Mi = 2.0*1.6726e-27; // Ion mass

// Normalisation factors
BoutReal Bnorm, Tnorm;
BoutReal wci, cs, rho_s;

// options

bool nonlinear;
bool full_bfield;   // If true, use divergence-free expression for B
bool flux_method;   // Use flux methods in rho and T equations
bool full_v_method; // Calculate full velocity equation

Vector3D vExB, vD; // Velocities

// Communication objects
FieldGroup comms;

int physics_init(bool restarting) {
  
  output.write("Solving JOREK-like reduced MHD equations\n");
  output.write("\tFile    : %s\n", __FILE__);
  output.write("\tCompiled: %s at %s\n", __DATE__, __TIME__);

  //////////////////////////////////////////////////////////////
  // Load data from the grid

  // Load 2D profiles
  mesh->get(J0, "Jpar0");    // A / m^2
  mesh->get(P0, "pressure"); // Pascals

  // Load curvature term
  b0xcv.covariant = false; // Read contravariant components
  mesh->get(b0xcv, "bxcv"); // mixed units x: T y: m^-2 z: m^-2

  // Load metrics
  // Metric coefficients
  Field2D Rxy, Bpxy, Btxy, hthe;
  Field2D I; // Shear factor
  
  if(mesh->get(Rxy,  "Rxy")) { // m
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

  Options *globalOptions = Options::getRoot();
  Options *options = globalOptions->getSection("jorek");
  
  OPTION(options, nonlinear,           false);
  OPTION(options, full_bfield,         false);
  OPTION(options, flux_method,         false);
  OPTION(options, full_v_method,       false);
  
  //////////////////////////////////////////////////////////////
  // SHIFTED RADIAL COORDINATES

  if(mesh->ShiftXderivs) {
    if(mesh->IncIntShear) {
      // BOUT-06 style, using d/dx = d/dpsi + I * d/dz
      mesh->IntShiftTorsion = I;
      
    }else {
      // Dimits style, using local coordinate system
      b0xcv.z += I*b0xcv.x;
      I = 0.0;  // I disappears from metric
    }
  }

  //////////////////////////////////////////////////////////////
  // NORMALISE QUANTITIES
  
  if(mesh->get(Bnorm, "bmag")) // Typical magnetic field (T)
    Bnorm = 1.0;
  
  wci = Charge * Bnorm / Mi; // Cyclotron angular frequency
  cs = sqrt(Charge * Tnorm / Mi); // Sound speed
  rho_s = cs / wci; // Larmor radius
  
  // Normalising times to 1/wci and length to rho_s
  
  
  
  //////////////////////////////////////////////////////////////
  // CALCULATE METRICS
  
  mesh->g11 = (Rxy*Bpxy)^2;
  mesh->g22 = 1.0 / (hthe^2);
  mesh->g33 = (I^2)*mesh->g11 + (B0^2)/mesh->g11;
  mesh->g12 = 0.0;
  mesh->g13 = -I*mesh->g11;
  mesh->g23 = -Btxy/(hthe*Bpxy*Rxy);
  
  mesh->J = hthe / Bpxy;
  mesh->Bxy = B0;
  
  mesh->g_11 = 1.0/mesh->g11 + ((I*Rxy)^2);
  mesh->g_22 = (B0*hthe/Bpxy)^2;
  mesh->g_33 = Rxy*Rxy;
  mesh->g_12 = Btxy*hthe*I*Rxy/Bpxy;
  mesh->g_13 = I*Rxy*Rxy;
  mesh->g_23 = Btxy*hthe*Rxy/Bpxy;
  
  mesh->geometry(); // Calculate quantities from metric tensor

  // Set B field vector
  B0vec.covariant = false;
  B0vec.x = 0.;
  B0vec.y = Bpxy / hthe;
  B0vec.z = 0.;
  
  vExB.setBoundary("v");
  vD.setBoundary("v");
  
  eta = eta0;

  // SET EVOLVING VARIABLES

  SOLVE_FOR5(rho, T, u, vpar, psi);
  
  comms.add(rho);
  comms.add(T);
  comms.add(u);
  comms.add(phi);
  comms.add(vpar);
  comms.add(psi);

  return 0;
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
  phi = invert_laplace(u, phi_flags, NULL);
  // Apply a boundary condition on phi for target plates
  phi.applyBoundary();
  
  // Communicate variables
  mesh->communicate(comms);
  
  // Get J from Psi
  Jpar = -B0*Delp2(psi);
  Jpar.applyBoundary();
  mesh->communicate(Jpar);
  
  Field3D rhot = rho0;
  Field3D Tt = T0;
  Field3D P = rho*T0 + T*rho0; // Perturbed pressure
  if(nonlinear) {
    rhot += rho;
    Tt += T;
    P += rho*T;
    
    eta = eta0*((Tt/Tnorm)^(-1.5)); // Update resistivity
  }

  if(flux_method) {
    // ExB velocity
    vExB = (B0vec ^ Grad_perp(phi))/(B0*B0);
    vExB.applyBoundary();
    
    ////////// Density equation ////////////////
    
    // Diffusive flux (perpendicular)
    vD = -D_perp * Grad_perp(rho);
    vD.applyBoundary();
    
    ddt(rho) = -Div(vExB + vD, rhot);
    
    ////////// Temperature equation ////////////
  
    vD = -chi_perp * Grad_perp(T);
    vD.applyBoundary();
    
    ddt(T) = 
      - b0xGrad_dot_Grad(phi, Tt)
      - (2./3.)*Tt*Div(vExB)
      - (Div(vD, T) - Div_par_K_Grad_par(chi_par, T))/rhot
      ;
  }else {
    // Use analytic expressions, expand terms
    
    // Divergence of ExB velocity (neglecting parallel term)
    Field3D divExB = b0xcv*Grad(phi)/B0 - b0xGrad_dot_Grad(1./B0, phi);
    
    ddt(rho) = -b0xGrad_dot_Grad(phi, rhot) // Advection 
      - divExB*rhot; // Compression
    
    ddt(T) = -b0xGrad_dot_Grad(phi, Tt) 
      - (2./3.)*Tt*divExB;
  }
  
  if(full_v_method) {
    vExB = (B0vec ^ Grad_perp(phi))/(B0*B0);
    
    ddt(vExB) = (-Grad(P))/rho;
    
    // Use this to calculate a vorticity and parallel velocity
    ddt(u) = B0vec * Curl(ddt(vExB));
    ddt(vpar) = B0vec * ddt(vExB);
  }else {
    // Split into vorticity and parallel velocity equations analytically
    
    ////////// Vorticity equation ////////////

    ddt(u) = (
	      (B0^2)*Grad_parP(Jpar/B0, CELL_CENTRE)
	      - (B0^2) * b0xGrad_dot_Grad(psi, J0/B0, CELL_CENTRE) 
	      + 2.*b0xcv*Grad(P)  // curvature term
	      ) / rhot;
    
    if(nonlinear) {
      ddt(u) -= b0xGrad_dot_Grad(phi, u);    // Advection
    }
    
    // Viscosity terms 
    if(viscos_par > 0.0)
      ddt(u) += viscos_par * Grad2_par2(u); // Parallel viscosity
    
    if(viscos_perp > 0.0)
      ddt(u) += viscos_perp * Delp2(u);     // Perpendicular viscosity
    
    ////////// Parallel velocity equation ////////////
    
    ddt(vpar) = -Grad_parP(P + rho0*T0, CELL_YLOW);
    if(nonlinear)
      ddt(vpar) -= b0xGrad_dot_Grad(phi, vpar); // Advection
  }

  ////////// Magnetic potential equation ////////////

  ddt(psi) = -Grad_parP(B0*phi, CELL_CENTRE) / B0 + eta*Jpar;
  
  return 0;
}
