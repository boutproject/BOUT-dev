/*******************************************************************************
 * 2-fluid equations
 * Same as Maxim's version of BOUT - simplified 2-fluid for benchmarking
 *******************************************************************************/

#include <bout.hxx>
#include <boutmain.hxx>

#include <initialprofiles.hxx>
#include <derivs.hxx>
#include <interpolation.hxx>

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

// 2D initial profiles
Field2D Ni0, Ti0, Te0, Vi0, phi0, Ve0, rho0, Ajpar0;
// Staggered versions of initial profiles
Field2D Ni0_maybe_ylow, Te0_maybe_ylow;
Vector2D b0xcv; // for curvature terms

// 3D evolving fields
Field3D rho, Te, Ni, Ajpar, Vi, Ti;

// Derived 3D variables
Field3D phi, Apar, Ve, jpar;

// Non-linear coefficients
Field3D nu, mu_i, kapa_Te, kapa_Ti;

// 3D total values
Field3D Nit, Tit, Tet, Vit;

// pressures
Field3D pei, pe;
Field2D pei0, pe0;

// Metric coefficients
Field2D Rxy, Bpxy, Btxy, hthe;

// parameters
BoutReal Te_x, Ti_x, Ni_x, Vi_x, bmag, rho_s, fmei, AA, ZZ;
BoutReal lambda_ei, lambda_ii;
BoutReal nu_hat, mui_hat, wci, nueix, nuiix;
BoutReal beta_p;

// settings
bool estatic, ZeroElMass; // Switch for electrostatic operation (true = no Apar)
BoutReal zeff, nu_perp;
bool evolve_rho, evolve_te, evolve_ni, evolve_ajpar, evolve_vi, evolve_ti;
BoutReal ShearFactor;

int phi_flags, apar_flags; // Inversion flags

// Field routines
int solve_phi_tridag(Field3D &r, Field3D &p, int flags);
int solve_apar_tridag(Field3D &aj, Field3D &ap, int flags);


FieldGroup comms; // Group of variables for communications

Coordinates *coord; // Coordinate system

CELL_LOC maybe_ylow;

int physics_init(bool restarting) {
  Field2D I; // Shear factor 
  
  output.write("Solving 6-variable 2-fluid equations\n");
  
  /************* LOAD DATA FROM GRID FILE ****************/

  // Load 2D profiles (set to zero if not found)
  GRID_LOAD(Ni0);
  GRID_LOAD(Ti0);
  GRID_LOAD(Te0);
  GRID_LOAD(Vi0);
  GRID_LOAD(Ve0);
  GRID_LOAD(phi0);
  GRID_LOAD(rho0);
  GRID_LOAD(Ajpar0);

  // Load magnetic curvature term
  b0xcv.covariant = false; // Read contravariant components
  mesh->get(b0xcv, "bxcv"); // b0xkappa terms

  // Coordinate system
  coord = mesh->coordinates();

  // Load metrics
  GRID_LOAD(Rxy);
  GRID_LOAD(Bpxy);
  GRID_LOAD(Btxy);
  GRID_LOAD(hthe);
  mesh->get(coord->dx,   "dpsi");
  mesh->get(I,    "sinty");

  // Load normalisation values
  GRID_LOAD(Te_x);
  GRID_LOAD(Ti_x);
  GRID_LOAD(Ni_x);
  GRID_LOAD(bmag);

  Ni_x *= 1.0e14;
  bmag *= 1.0e4;

  /*************** READ OPTIONS *************************/

  // Read some parameters
  Options *globalOptions = Options::getRoot();
  Options *options = globalOptions->getSection("2fluid");
  OPTION(options, AA, 2.0);
  OPTION(options, ZZ, 1.0);

  OPTION(options, estatic,     false);
  OPTION(options, ZeroElMass,  false);
  OPTION(options, zeff,        1.0);
  OPTION(options, nu_perp,     0.0);
  OPTION(options, ShearFactor, 1.0);
  
  OPTION(options, phi_flags,   0);
  OPTION(options, apar_flags,  0);
  
  (globalOptions->getSection("Ni"))->get("evolve", evolve_ni,    true);
  (globalOptions->getSection("rho"))->get("evolve", evolve_rho,   true);
  (globalOptions->getSection("vi"))->get("evolve", evolve_vi,   true);
  (globalOptions->getSection("te"))->get("evolve", evolve_te,   true);
  (globalOptions->getSection("ti"))->get("evolve", evolve_ti,   true);
  (globalOptions->getSection("Ajpar"))->get("evolve", evolve_ajpar, true);
  
  if(ZeroElMass)
    evolve_ajpar = false; // Don't need ajpar - calculated from ohm's law

  /************* SHIFTED RADIAL COORDINATES ************/

  bool ShiftXderivs;
  globalOptions->get("shiftXderivs", ShiftXderivs, false); // Read global flag
  if(ShiftXderivs) {
    ShearFactor = 0.0;  // I disappears from metric
    b0xcv.z += I*b0xcv.x;
  }

  /************** CALCULATE PARAMETERS *****************/

  rho_s = 1.02*sqrt(AA*Te_x)/ZZ/bmag;
  fmei  = 1./1836.2/AA;

  lambda_ei = 24.-log(sqrt(Ni_x)/Te_x);
  lambda_ii = 23.-log(ZZ*ZZ*ZZ*sqrt(2.*Ni_x)/pow(Ti_x, 1.5));
  wci       = 9.58e3*ZZ*bmag/AA;
  nueix     = 2.91e-6*Ni_x*lambda_ei/pow(Te_x, 1.5);
  nuiix     = 4.78e-8*pow(ZZ,4.)*Ni_x*lambda_ii/pow(Ti_x, 1.5)/sqrt(AA);
  nu_hat    = zeff*nueix/wci;

  if(nu_perp < 1.e-10) {
    mui_hat      = (3./10.)*nuiix/wci;
  } else
    mui_hat      = nu_perp;

  if(estatic) {
    beta_p    = 1.e-29;
  }else
    beta_p    = 4.03e-11*Ni_x*Te_x/bmag/bmag;

  Vi_x = wci * rho_s;

  /************** PRINT Z INFORMATION ******************/
  
  BoutReal hthe0;
  if(mesh->get(hthe0, "hthe0") == 0) {
    output.write("    ****NOTE: input from BOUT, Z length needs to be divided by %e\n", hthe0/rho_s);
  }

  /************** NORMALISE QUANTITIES *****************/

  output.write("\tNormalising to rho_s = %e\n", rho_s);

  // Normalise profiles
  Ni0 /= Ni_x/1.0e14;
  Ti0 /= Te_x;
  Te0 /= Te_x;
  phi0 /= Te_x;
  Vi0 /= Vi_x;

  // Normalise curvature term
  b0xcv.x /= (bmag/1e4);
  b0xcv.y *= rho_s*rho_s;
  b0xcv.z *= rho_s*rho_s;
  
  // Normalise geometry 
  Rxy /= rho_s;
  hthe /= rho_s;
  I *= rho_s*rho_s*(bmag/1e4)*ShearFactor;
  coord->dx /= rho_s*rho_s*(bmag/1e4);

  // Normalise magnetic field
  Bpxy /= (bmag/1.e4);
  Btxy /= (bmag/1.e4);
  coord->Bxy  /= (bmag/1.e4);

  // calculate pressures
  pei0 = (Ti0 + Te0)*Ni0;
  pe0 = Te0*Ni0;

  /**************** CALCULATE METRICS ******************/

  coord->g11 = SQ(Rxy*Bpxy);
  coord->g22 = 1.0 / SQ(hthe);
  coord->g33 = SQ(I)*coord->g11 + SQ(coord->Bxy)/coord->g11;
  coord->g12 = 0.0;
  coord->g13 = -I*coord->g11;
  coord->g23 = -Btxy/(hthe*Bpxy*Rxy);
  
  coord->J = hthe / Bpxy;
  
  coord->g_11 = 1.0/coord->g11 + SQ(I*Rxy);
  coord->g_22 = SQ(coord->Bxy*hthe/Bpxy);
  coord->g_33 = Rxy*Rxy;
  coord->g_12 = Btxy*hthe*I*Rxy/Bpxy;
  coord->g_13 = I*Rxy*Rxy;
  coord->g_23 = Btxy*hthe*Rxy/Bpxy;

  coord->geometry();
  
  /**************** SET EVOLVING VARIABLES *************/

  // Tell BOUT++ which variables to evolve
  // add evolving variables to the communication object
  if(evolve_rho) {
    bout_solve(rho, "rho");
    comms.add(rho);
    output.write("rho\n");
  }else
    initial_profile("rho", rho);

  if(evolve_ni) {
    bout_solve(Ni, "Ni");
    comms.add(Ni);
    output.write("ni\n");
  }else
    initial_profile("Ni", Ni);

  if(evolve_te) {
    bout_solve(Te, "Te");
    comms.add(Te);
    output.write("te\n");
  }else
    initial_profile("Te", Te);

  if(evolve_ajpar) {
    bout_solve(Ajpar, "Ajpar");
    comms.add(Ajpar);
    output.write("ajpar\n");
  }else {
    initial_profile("Ajpar", Ajpar);
    if(ZeroElMass)
      dump.add(Ajpar, "Ajpar", 1); // output calculated Ajpar
  }

  if(evolve_vi) {
    bout_solve(Vi, "Vi");
    comms.add(Vi);
    output.write("vi\n");
  }else
    initial_profile("Vi", Vi);

  if(evolve_ti) {
    bout_solve(Ti, "Ti");
    comms.add(Ti);
    output.write("ti\n");
  }else
    initial_profile("Ti", Ti);
  
  // Set boundary conditions
  jpar.setBoundary("jpar");

  /************** SETUP COMMUNICATIONS **************/

  // add extra variables to communication
  comms.add(phi);
  comms.add(Apar);

  // Add any other variables to be dumped to file
  dump.add(phi,  "phi",  1);
  dump.add(Apar, "Apar", 1);
  dump.add(jpar, "jpar", 1);

  dump.add(Ni0, "Ni0", 0);
  dump.add(Te0, "Te0", 0);
  dump.add(Ti0, "Ti0", 0);

  dump.add(Te_x,  "Te_x", 0);
  dump.add(Ti_x,  "Ti_x", 0);
  dump.add(Ni_x,  "Ni_x", 0);
  dump.add(rho_s, "rho_s", 0);
  dump.add(wci,   "wci", 0);

  if (mesh->StaggerGrids) {
    maybe_ylow = CELL_YLOW;
  } else {
    maybe_ylow = CELL_CENTRE;
  }
  Vi = interp_to(Vi,maybe_ylow);
  Ni0_maybe_ylow = interp_to(Ni0, maybe_ylow, RGN_NOBNDRY);
  Te0_maybe_ylow = interp_to(Te0, maybe_ylow, RGN_NOBNDRY);

  return(0);
}

// just define a macro for V_E dot Grad
#define vE_Grad(f, p) ( b0xGrad_dot_Grad(p, f) / coord->Bxy )

int physics_run(BoutReal t) {
  // Solve EM fields

  solve_phi_tridag(rho, phi, phi_flags);

  if(estatic || ZeroElMass) {
    // Electrostatic operation
    Apar = 0.0;
  }else {
    solve_apar_tridag(Ajpar, Apar, apar_flags); // Linear Apar solver
  }

  // Communicate variables
  mesh->communicate(comms);

  // Update profiles
  Nit = Ni0;  //+ Ni.DC();
  Tit = Ti0; // + Ti.DC();
  Tet = Te0; // + Te.DC();
  Vit = Vi0; // + Vi;

  // Update non-linear coefficients on the mesh
  nu      = nu_hat * Nit / pow(Tet,1.5);
  mu_i    = mui_hat * Nit / sqrt(Tit);
  kapa_Te = 3.2*(1./fmei)*(wci/nueix)*pow(Tet,2.5);
  kapa_Ti = 3.9*(wci/nuiix)*pow(Tit,2.5);
  
  // note: nonlinear terms are not here
  pei = (Te0+Ti0)*Ni + (Te + Ti)*Ni0;
  pe  = Te0*Ni + Te*Ni0;
  
  if(ZeroElMass) {
    // Set jpar,Ve,Ajpar neglecting the electron inertia term
    jpar = ((Te0_maybe_ylow*Grad_par(Ni, maybe_ylow)) - (Ni0_maybe_ylow*Grad_par(phi, maybe_ylow)))/interp_to(fmei*0.51*nu,maybe_ylow);
    
    // Set boundary conditions on jpar (in BOUT.inp)
    jpar.applyBoundary();
    
    // Need to communicate jpar
    mesh->communicate(jpar);

    Ve = Vi - jpar/Ni0_maybe_ylow;
    Ajpar = Ve;
  }else {
    
    Ve = Ajpar + Apar;
    jpar = Ni0_maybe_ylow*(Vi - Ve);
  }

  // DENSITY EQUATION

  ddt(Ni) = 0.0;
  if(evolve_ni) {
    ddt(Ni) -= vE_Grad(Ni0, phi);
  }

  // ION VELOCITY

  ddt(Vi) = 0.0;
  if(evolve_vi) {
    ddt(Vi) -= vE_Grad(Vi0, phi) + vE_Grad(Vi, phi0) + vE_Grad(Vi, phi);
    ddt(Vi) -= Vpar_Grad_par(Vi0, Vi) + Vpar_Grad_par(Vi, Vi0) + Vpar_Grad_par(Vi, Vi);
    ddt(Vi) -= Grad_par(pei)/Ni0_maybe_ylow;
  }

  // ELECTRON TEMPERATURE

  ddt(Te) = 0.0;
  if(evolve_te) {
    ddt(Te) -= vE_Grad(Te0, phi) + vE_Grad(Te, phi0) + vE_Grad(Te, phi);
    ddt(Te) -= Vpar_Grad_par(Ve, Te0) + Vpar_Grad_par(Ve0, Te) + Vpar_Grad_par(Ve, Te);
    ddt(Te) += 1.333*Te0*( V_dot_Grad(b0xcv, pe)/Ni0 - V_dot_Grad(b0xcv, phi) );
    ddt(Te) += 3.333*Te0*V_dot_Grad(b0xcv, Te);
    ddt(Te) += (0.6666667/Ni0)*Div_par_K_Grad_par(kapa_Te, Te);
  }

  // ION TEMPERATURE

  ddt(Ti) = 0.0;
  if(evolve_ti) {
    ddt(Ti) -= vE_Grad(Ti0, phi) + vE_Grad(Ti, phi0) + vE_Grad(Ti, phi);
    ddt(Ti) -= Vpar_Grad_par(Vi, Ti0) + Vpar_Grad_par(Vi0, Ti) + Vpar_Grad_par(Vi, Ti);
    ddt(Ti) += 1.333*( Ti0*V_dot_Grad(b0xcv, pe)/Ni0 - Ti*V_dot_Grad(b0xcv, phi) );
    ddt(Ti) -= 3.333*Ti0*V_dot_Grad(b0xcv, Ti);
    ddt(Ti) += (0.6666667/Ni0)*Div_par_K_Grad_par(kapa_Ti, Ti);
  }

  // VORTICITY

  ddt(rho) = 0.0;
  if(evolve_rho) {

    ddt(rho) += SQ(coord->Bxy)*Div_par(jpar, CELL_CENTRE);
  }
  

  // AJPAR
  

  ddt(Ajpar) = 0.0;
  if(evolve_ajpar) {
    ddt(Ajpar) += (1./fmei)*Grad_par(phi, maybe_ylow);
    ddt(Ajpar) -= (1./fmei)*(Te0_maybe_ylow/Ni0_maybe_ylow)*Grad_par(Ni, maybe_ylow);
    //ddt(Ajpar) -= (1./fmei)*1.71*Grad_par(Te);
    ddt(Ajpar) += 0.51*interp_to(nu, maybe_ylow)*jpar/Ni0_maybe_ylow;
  }

  return(0);
}

/*******************************************************************************
 *                       FAST LINEAR FIELD SOLVERS
 *******************************************************************************/

#include <invert_laplace.hxx>

// Performs inversion of rho (r) to get phi (p)
int solve_phi_tridag(Field3D &r, Field3D &p, int flags) {
  
  //output.write("Solving phi: %e, %e -> %e\n", max(abs(r)), min(Ni0), max(abs(r/Ni0)));

  if(invert_laplace(r/Ni0, p, flags, NULL)) {
    return 1;
  }

  //Field3D pertPi = Ti*Ni0 + Ni*Ti0;
  //p -= pertPi/Ni0;
  return(0);
}

int solve_apar_tridag(Field3D &aj, Field3D &ap, int flags) {
  static Field2D a;
  static int set = 0;

  if(set == 0) {
    // calculate a
    a = (-0.5*beta_p/fmei)*Ni0;
    set = 1;
  }

  if(invert_laplace(a*(Vi - aj), ap, flags, &a))
    return 1;

  return(0);
}
