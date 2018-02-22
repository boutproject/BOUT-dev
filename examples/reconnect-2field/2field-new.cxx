/*****************************************************************************
 * 2 field (Apar, vorticity) model for benchmarking 
 * simple slab reconnection model
 *****************************************************************************/

#include <bout.hxx>
#include <bout/boutmain.hxx>

#include <bout/invert_laplace.hxx>
#include <bout/initialprofiles.hxx>

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

// 2D initial profiles
Field2D Jpar0, Te0, Ni0;

// 3D evolving fields
Field3D U, Apar;

// Derived 3D variables
Field3D phi, jpar;

// External coil field
Field3D Apar_ext, Jpar_ext, Phi0_ext, Upar0_ext;

// Metric coefficients
Field2D Rxy, Bpxy, Btxy, hthe, Bxy;

// Constants
const BoutReal MU0 = 4.0e-7*PI;
const BoutReal Charge = 1.60217646e-19; // electron charge e (C)
const BoutReal Mi = 2.0*1.67262158e-27; // Ion mass
const BoutReal Me = 9.1093816e-31;  // Electron mass
const BoutReal Me_Mi = Me / Mi; // Electron mass / Ion mass

// normalisation parameters
BoutReal Tenorm, Nenorm, Bnorm;
BoutReal Cs, rho_s, wci, beta_hat;

BoutReal eta, mu_perp, mu_par;
BoutReal eta_hyper_par, eta_hyper_perp;
BoutReal mu_hyper_par,  mu_hyper_perp;

// Poisson brackets: b0 x Grad(f) dot Grad(g) / B = [f, g]
// Method to use: BRACKET_ARAKAWA, BRACKET_STD or BRACKET_SIMPLE
BRACKET_METHOD bm; // Bracket method for advection terms

int phi_flags; // Inversion flags

bool nonlinear;
bool parallel_lc;
bool include_jpar0;
int  jpar_bndry, jpar_ext_bndry;
bool test_advection;
bool include_apar_ext, include_phi0_ext;

void smooth_bndry(const Field3D &jpar, int jpar_bndry);

int physics_init(bool restarting) {

  // Load 2D profiles
  GRID_LOAD3(Jpar0, Te0, Ni0);
  Ni0 *= 1e20; // To m^-3
  
  // Load metrics
  GRID_LOAD(Rxy);
  GRID_LOAD(Bpxy);
  GRID_LOAD(Btxy);
  GRID_LOAD(hthe);
  mesh->get(mesh->Bxy,  "Bxy");
  
  // Read some parameters
  Options *globalOptions = Options::getRoot();
  Options *options = globalOptions->getSection("2field");
  
  // normalisation values
  OPTION(options, test_advection, false);
  OPTION(options, nonlinear, false);
  OPTION(options, parallel_lc, true);
  OPTION(options, include_jpar0, true);
  OPTION(options, jpar_bndry, 0);

  OPTION(options, eta,      1.e-3); // Normalised resistivity
  OPTION(options, mu_perp,  1.e-7); // Normalised viscosity
  OPTION(options, mu_par,   1.e-7); // Normalised viscosity
  OPTION(options, eta_hyper_perp, 0.e-3); // Normalised resistivity
  OPTION(options, eta_hyper_par,  0.e-3); // Normalised resistivity
  OPTION(options, mu_hyper_perp,  0.e-3); // Normalised viscosity
  OPTION(options, mu_hyper_par,   0.e-3); // Normalised viscosity
  
  OPTION(options, phi_flags,   0);
  
  int bracket_method;
  OPTION(options, bracket_method, 0);
  switch(bracket_method) {
  case 0: {
    bm = BRACKET_STD; 
    output << "\tBrackets: default differencing\n";
    break;
  }
  case 1: {
    bm = BRACKET_SIMPLE; 
    output << "\tBrackets: simplified operator\n";
    break;
  }
  case 2: {
    bm = BRACKET_ARAKAWA; 
    output << "\tBrackets: Arakawa scheme\n";
    break;
  }
  case 3: {
    bm = BRACKET_CTU; 
    output << "\tBrackets: Corner Transport Upwind method\n";
    break;
  }
  default:
    output << "ERROR: Invalid choice of bracket method. Must be 0 - 3\n";
    return 1;
  }

  ///////////////////////////////////////////////////
  // Normalisation
  
  Tenorm = max(Te0, true);
  if(Tenorm < 1)
    Tenorm = 1000;
  Nenorm = max(Ni0, true);
  if(Nenorm < 1)
    Nenorm = 1.e19;
  Bnorm  = max(mesh->Bxy, true);
  
  // Sound speed in m/s
  Cs = sqrt(Charge*Tenorm / Mi);

  // drift scale
  rho_s = Cs * Mi / (Charge * Bnorm);
  
  // Ion cyclotron frequency
  wci = Charge * Bnorm / Mi;
  
  beta_hat = MU0 * Charge*Tenorm * Nenorm / (Bnorm*Bnorm);
  
  output << "\tNormalisations:" << endl;
  output << "\tCs       = " << Cs << endl;
  output << "\trho_s    = " << rho_s << endl;
  output << "\twci      = " << wci << endl; 
  output << "\tbeta_hat = " << beta_hat << endl; 
  
  SAVE_ONCE3(Tenorm, Nenorm, Bnorm);
  SAVE_ONCE4(Cs, rho_s, wci, beta_hat);
  
  // Normalise geometry 
  Rxy  /= rho_s;
  hthe /= rho_s;
  mesh->dx /= rho_s*rho_s*Bnorm;

  // Normalise magnetic field
  Bpxy /= Bnorm;
  Btxy /= Bnorm;
  mesh->Bxy /= Bnorm;
  Bxy = mesh->Bxy;

  // Plasma quantities
  Jpar0 /= Nenorm*Charge*Cs;

  // CALCULATE METRICS

  mesh->g11 = (Rxy*Bpxy)^2;
  mesh->g22 = 1.0 / (hthe^2);
  mesh->g33 = Bxy^2/mesh->g11;
  mesh->g12 = 0.0;
  mesh->g13 = 0.;
  mesh->g23 = -Btxy/(hthe*Bpxy*Rxy);
  
  mesh->J = hthe / Bpxy;
  
  mesh->g_11 = 1.0/mesh->g11;
  mesh->g_22 = (Bxy*hthe/Bpxy)^2;
  mesh->g_33 = Rxy*Rxy;
  mesh->g_12 = 0.;
  mesh->g_13 = 0.;
  mesh->g_23 = Btxy*hthe*Rxy/Bpxy;

  mesh->geometry();

  // Tell BOUT++ which variables to evolve
  SOLVE_FOR2(U, Apar);
  
  // Set boundary conditions
  jpar.setBoundary("jpar");
  phi.setBoundary("phi");
  
  // Add any other variables to be dumped to file
  SAVE_REPEAT2(phi, jpar);
  SAVE_ONCE(Jpar0);
  
  // Generate external field
  Options *Apar_ext_options = globalOptions->getSection("Apar_ext");
  OPTION(Apar_ext_options, include_apar_ext,    false);
  OPTION(Apar_ext_options, jpar_ext_bndry, 0);

  Options *Phi0_ext_options = globalOptions->getSection("Phi0_ext");
  OPTION(Phi0_ext_options, include_phi0_ext,   false);  

  if(include_apar_ext) {
    output << "\tInitializing Apar_ext\n";
    initial_profile("Apar_ext", Apar_ext);
    Jpar_ext = - Delp2(Apar_ext);
    Jpar_ext.applyBoundary();
    smooth_bndry(Jpar_ext,jpar_ext_bndry);
    SAVE_ONCE2(Apar_ext,Jpar_ext);
  }
  if(include_phi0_ext) {
    output << "\tInitializing Phi0_ext\n";
    initial_profile("Phi0_ext", Phi0_ext); 
    Upar0_ext = Delp2(Phi0_ext)/Bxy;
    Upar0_ext.applyBoundary();
    SAVE_ONCE2(Phi0_ext,Upar0_ext);
  }

  return 0;
}

const Field3D Grad_par0_LtoC(const Field3D &f) {
  Field3D result;
  if(parallel_lc) {
    result = Grad_par_LtoC(f);
  }else{
     result = Grad_par(f);
  }
  return result;
}

const Field3D Grad_par0_CtoL(const Field3D &f) {
  Field3D result;
  if(parallel_lc) {
    result = Grad_par_CtoL(f);
  }else{
     result = Grad_par(f);
  }
  return result;
}

const Field3D Grad_par1(const Field3D &f) {
  Field3D result=0.;  

  result -= beta_hat * bracket(Apar, f, BRACKET_ARAKAWA);
  if (include_apar_ext) {
    result -= beta_hat * bracket(Apar_ext, f, BRACKET_ARAKAWA);
  }
  return result;
}

void smooth_bndry(const Field3D &f, int bndry)
{
    // Boundary in jpar
    if(mesh->firstX()) {
      for(int i=bndry;i>=0;i--)
	for(int j=0;j<mesh->LocalNy;j++)
	  for(int k=0;k<mesh->LocalNz;k++) {
	    f[i][j][k] = f[i+1][j][k];
	  }
    }
    if(mesh->lastX()) {
      for(int i=mesh->LocalNx-bndry-1;i<mesh->LocalNx;i++)
	for(int j=0;j<mesh->LocalNy;j++)
	  for(int k=0;k<mesh->LocalNz;k++) {
	    f[i][j][k] = f[i-1][j][k];
	  }
    }
}

int physics_run(BoutReal t) {
  // Solve EM fields

  // U = (1/B) * Delp2(phi)
  phi = invert_laplace(Bxy*U, phi_flags);
  phi.applyBoundary(); // For target plates only
  
  mesh->communicate(U, phi, Apar);
  
  jpar = -Delp2(Apar);
  jpar.applyBoundary();
  mesh->communicate(jpar);
  
  if(jpar_bndry > 0) {
    smooth_bndry(jpar,jpar_bndry);
  }

  // VORTICITY
  ddt(U)=0.;

  if(!test_advection){
    ddt(U) += SQ(Bxy)*Grad_par0_LtoC(jpar/Bxy);
    if(include_apar_ext)
      ddt(U) += SQ(Bxy)*Grad_par0_LtoC(Jpar_ext/Bxy);  
  }
  if(include_jpar0){
    // Grad_par0(Jpar0) should vanish
    ddt(U) += SQ(Bxy)*Grad_par1(Jpar0/Bxy);
  }

  if(include_phi0_ext) 
    ddt(U) -= bracket(Phi0_ext, U, bm);

  if(nonlinear) {
    ddt(U) -= bracket(phi, U, bm); // ExB advection
    ddt(U) -= SQ(Bxy)*Grad_par1(jpar/Bxy);
    if(include_apar_ext)
      ddt(U) += SQ(Bxy)*Grad_par1(Jpar_ext/Bxy);
  }
  
  if(mu_perp > 0.)
    ddt(U) += mu_perp*Delp2(U);
  if(mu_par > 0.)
    ddt(U) += mu_par*Grad2_par2(U);
  if(mu_hyper_perp > 0.)
    ddt(U) -= mu_hyper_perp*Delp2(Delp2(U));
  if(mu_hyper_par > 0.)
    ddt(U) += mu_hyper_par*Grad2_par2(Grad2_par2(U));

  // APAR
  ddt(Apar)=0.; 

  if(!test_advection)
    ddt(Apar) -= Grad_par0_CtoL(phi) / beta_hat;

  if (include_phi0_ext) {
    ddt(Apar) -= Grad_par0_CtoL(Phi0_ext) / beta_hat;
    ddt(Apar) -= Grad_par1(Phi0_ext) / beta_hat;
  }
  if (nonlinear) {
    ddt(Apar) -= Grad_par1(phi) / beta_hat;
  }

  if(eta > 0.) {
    ddt(Apar) -= eta*jpar / beta_hat;
    if (include_apar_ext) 
      ddt(Apar) -= eta*Jpar_ext / beta_hat;
  }
  if(eta_hyper_perp > 0.) {
    ddt(Apar) += eta_hyper_perp*Delp2(jpar) / beta_hat;
    if (include_apar_ext) 
      ddt(Apar) += eta_hyper_perp*Delp2(Jpar_ext) / beta_hat;
  }
  if(eta_hyper_par > 0.)  { 
    ddt(Apar) += eta_hyper_par*Grad_par0_LtoC(Bxy*Grad_par0_CtoL(jpar/Bxy)) / beta_hat;
    if (include_apar_ext)
      ddt(Apar) += eta_hyper_par*Grad_par0_LtoC(Bxy*Grad_par0_CtoL(Jpar_ext/Bxy)) / beta_hat;
  } 
  return 0;
}

