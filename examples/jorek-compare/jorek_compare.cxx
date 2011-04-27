/**************************************************************************
 * Similar set of equations to JOREK
 * 
 **************************************************************************/

#include <bout.hxx>
#include <boutmain.hxx>

#include <invert_laplace.hxx>
#include <math.h>

// Evolving quantities
Field3D rho, Te, Ti, u, vpar, psi;
// Derived quantities
Field3D Jpar, phi; // Parallel current, electric potential

// Equilibrium quantities
Field2D rho0, Te0, Ti0; // Equilibrium mass density, electron and ion temperature
Field2D B0, J0, P0;
Vector2D b0xcv; // Curvature term
Vector2D B0vec; // B0 field vector

// Dissipation coefficients
Field2D D_perp; // Particle diffusion coefficient
Field2D chi_eperp, chi_epar; // Electron heat diffusion coefficients
Field2D chi_iperp, chi_ipar; // Ion heat diffusion coefficients
Field2D eta0;  // Resistivity
Field3D eta;
BoutReal viscos_par, viscos_perp; // Viscosity coefficients

int phi_flags;

// Constants
const BoutReal MU0 = 4.0e-7*PI;
const BoutReal Charge = 1.602e-19; // electron charge e (C)
const BoutReal Mi = 2.0*1.6726e-27; // Ion mass

// Normalisation factors
BoutReal Tnorm, rhonorm; // Partial normalisation to rho and MU0. Temperature normalised

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

  Options *globalOptions = Options::getRoot();
  Options *options = globalOptions->getSection("jorek");

  //////////////////////////////////////////////////////////////
  // Load data from the grid

  // Load 2D profiles
  if(mesh->get(J0, "Jpar0"));    // A / m^2
  
  if(mesh->get(rho0, "Ni0")) {
    output << "Warning: No density profile available\n";
    BoutReal d0;
    options->get("density", d0, 1e20);
    rho0 = d0;
  }
  
  // Read temperature
  mesh->get(Te0, "Te0");
  mesh->get(Ti0, "Ti0");

  // Try reading pressure profile (in Pascals)
  if(mesh->get(P0, "pressure")) {
    // Just calculate from Temp and density
    P0 = Charge * (Ti0 + Te0) * rho0;
  }else {
    // Make sure that density and temperature are consistent with pressure
    
    Field2D factor = P0 / (Charge * (Ti0 + Te0) * rho0);
    
    output.write("\tPressure factor %e -> %e\n", min(factor,true), max(factor, true));
    
    // Multiply temperatures by this factor
    Te0 *= factor;
    Ti0 *= factor;
  }
  rho0 *= Mi; // Convert density to mass density [kg / m^3]

  // Load dissipation coefficients, override in options file
  if(options->isSet("D_perp")) {
    BoutReal tmp;
    options->get("D_perp", tmp, 0.0);
    D_perp = tmp;
  }else mesh->get(D_perp, "D_perp");

  if(options->isSet("chi_eperp")) {
    BoutReal tmp;
    options->get("chi_eperp", tmp, 0.0);
    chi_eperp = tmp;
  }else mesh->get(chi_eperp, "chi_eperp");

  if(options->isSet("chi_iperp")) {
    BoutReal tmp;
    options->get("chi_iperp", tmp, 0.0);
    chi_iperp = tmp;
  }else mesh->get(chi_iperp, "chi_iperp");

  if(options->isSet("chi_epar")) {
    BoutReal tmp;
    options->get("chi_epar", tmp, 0.0);
    chi_epar = tmp;
  }else mesh->get(chi_epar, "chi_epar");

  if(options->isSet("chi_ipar")) {
    BoutReal tmp;
    options->get("chi_ipar", tmp, 0.0);
    chi_ipar = tmp;
  }else mesh->get(chi_ipar, "chi_ipar");

  if(options->isSet("eta")) {
    BoutReal tmp;
    options->get("eta", tmp, 0.0);
    eta0 = tmp;
  }else mesh->get(eta, "eta");

  if(options->isSet("viscos_perp")) {
    BoutReal tmp;
    options->get("viscos_perp", tmp, 0.0);
    viscos_perp = tmp;
  }else mesh->get(viscos_perp, "viscos_perp");

  if(options->isSet("viscos_par")) {
    BoutReal tmp;
    options->get("viscos_par", tmp, 0.0);
    viscos_par = tmp;
  }else mesh->get(viscos_par, "viscos_par");

  // Load curvature term
  b0xcv.covariant = false; // Read contravariant components
  mesh->get(b0xcv, "bxcv"); // mixed units x: T y: m^-2 z: m^-2
  
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
  
  rhonorm = max(rho0, true); // Maximum over all grid
  Tnorm = Mi / (MU0 * Charge * rhonorm); // Temperature normalisation

  SAVE_ONCE2(rhonorm, Tnorm); // Save normalisation factors to file

  // Normalise quantities
  
  P0 *= MU0;
  J0 *= MU0;
  rho0 /= rhonorm;
  Te0 /= Tnorm;
  Ti0 /= Tnorm;
  
  eta0        *= sqrt(rhonorm / MU0);
  viscos_perp *= sqrt(MU0 / rhonorm);
  viscos_par  *= sqrt(MU0 / rhonorm);
  D_perp      *= sqrt(MU0 * rhonorm);
  chi_eperp   *= sqrt(MU0 / rhonorm);
  chi_epar    *= sqrt(MU0 / rhonorm);
  chi_iperp   *= sqrt(MU0 / rhonorm);
  chi_ipar    *= sqrt(MU0 / rhonorm);
  
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

  SOLVE_FOR6(rho, Te, Ti, u, vpar, psi);
  
  comms.add(rho);
  comms.add(Te);
  comms.add(Ti);
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
  Field3D Tet = Te0;
  Field3D Tit = Ti0;
  Field3D P = rho*(Te0+Ti0) + (Te+Ti)*rho0; // Perturbed pressure
  if(nonlinear) {
    rhot += rho;
    Tet += Te;
    Tit += Ti;
    P += rho*(Te+Ti);
    
    eta = eta0*((Tet/Te0)^(-1.5)); // Update resistivity based on Te
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
    
    ////////// Temperature equations ////////////
  
    vD = -chi_eperp * Grad_perp(Te);
    vD.applyBoundary();
    
    ddt(Te) = 
      - b0xGrad_dot_Grad(phi, Tet)
      - (2./3.)*Tet*Div(vExB)
      - (Div(vD, Te) - Div_par_K_Grad_par(chi_epar, Te))/rhot
      ;
    
    vD = -chi_iperp * Grad_perp(Ti);
    vD.applyBoundary();
    
    ddt(Ti) = 
      - b0xGrad_dot_Grad(phi, Tit)
      - (2./3.)*Tit*Div(vExB)
      - (Div(vD, Ti) - Div_par_K_Grad_par(chi_ipar, Ti))/rhot
      ;
  }else {
    // Use analytic expressions, expand terms
    
    // Divergence of ExB velocity (neglecting parallel term)
    Field3D divExB = b0xcv*Grad(phi)/B0 - b0xGrad_dot_Grad(1./B0, phi);
    
    ddt(rho) = -b0xGrad_dot_Grad(phi, rhot) // Advection 
      - divExB*rhot; // Compression
    
    ddt(Te) = 
      -b0xGrad_dot_Grad(phi, Tet) 
      - (2./3.)*Tet*divExB
      ;

    ddt(Ti) = 
      -b0xGrad_dot_Grad(phi, Tit) 
      - (2./3.)*Tit*divExB
      ;
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
    
    ddt(vpar) = -Grad_parP(P + P0, CELL_YLOW);
    if(nonlinear)
      ddt(vpar) -= b0xGrad_dot_Grad(phi, vpar); // Advection
  }

  ////////// Magnetic potential equation ////////////

  ddt(psi) = -Grad_parP(B0*phi, CELL_CENTRE) / B0 + eta*Jpar;
  
  return 0;
}
