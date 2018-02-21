/*****************************************************************************
 * 2 field (Apar, vorticity) model for benchmarking 
 * simple slab reconnection model
 *****************************************************************************/

#include <bout.hxx>
#include <bout/boutmain.hxx>

#include <bout/invert_laplace.hxx>
#include <bout/invert_parderiv.hxx>
#include <bout/initialprofiles.hxx>
#include <bout/constants.hxx>

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

// 2D initial profiles
Field2D Jpar0, Te0, Ni0;

// 3D evolving fields
Field3D Upar, Apar;

// Derived 3D variables
Field3D Phi, Jpar;

// External fields
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

BoutReal eta, mu;

// Poisson brackets: b0 x Grad(f) dot Grad(g) / B = [f, g]
// Method to use: BRACKET_ARAKAWA, BRACKET_STD or BRACKET_SIMPLE
BRACKET_METHOD bm; // Bracket method for advection terms

int phi_flags; // Inversion flags

bool nonlinear;
bool parallel_lc;
bool include_jpar0;
int jpar_bndry;

void smooth_bndry(Field3D f, int bndry = 2);
void set_bndry(Field3D f, BoutReal val=0.0, int bndry=2);
int precon(BoutReal t, BoutReal cj, BoutReal delta); // Preconditioner
InvertPar *inv; // Parallel inversion class used in preconditioner

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
  Bxy=mesh->Bxy;
  
  // Read some parameters
  Options *globalOptions = Options::getRoot();
  Options *options = globalOptions->getSection("2field");
  
  // normalisation values
  OPTION(options, nonlinear, false);
  OPTION(options, parallel_lc, true);
  OPTION(options, include_jpar0, true);
  OPTION(options, jpar_bndry, 0);

  OPTION(options, eta, 1e-3); // Normalised resistivity
  OPTION(options, mu, 1.e-3); // Normalised vorticity
  
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
  Bxy /= Bnorm;
  
  // Plasma quantities
  Jpar0 /= Nenorm*Charge*Cs;

  // CALCULATE METRICS

  mesh->g11 = (Rxy*Bpxy)^2;
  mesh->g22 = 1.0 / (hthe^2);
  mesh->g33 = (mesh->Bxy^2)/mesh->g11;
  mesh->g12 = 0.0;
  mesh->g13 = 0.;
  mesh->g23 = -Btxy/(hthe*Bpxy*Rxy);
  
  mesh->J = hthe / Bpxy;
  
  mesh->g_11 = 1.0/mesh->g11;
  mesh->g_22 = (mesh->Bxy*hthe/Bpxy)^2;
  mesh->g_33 = Rxy*Rxy;
  mesh->g_12 = 0.;
  mesh->g_13 = 0.;
  mesh->g_23 = Btxy*hthe*Rxy/Bpxy;

  mesh->geometry();

  // Tell BOUT++ which variables to evolve
  SOLVE_FOR2(Upar, Apar);
  
  // Set boundary conditions
  Jpar.setBoundary("Jpar");
  Phi.setBoundary("Phi");
  
  // Add any other variables to be dumped to file
  SAVE_REPEAT2(Phi, Jpar);
  SAVE_ONCE(Jpar0);
  
  // Generate external fields
  Apar_ext=0;
  Jpar_ext=0;
  initial_profile("Apar_ext", Apar_ext);
  Jpar_ext = -Delp2(Apar_ext);
  set_bndry(Jpar_ext,0.0,2);
  // Use vacuum field if requested
  mesh->communicate(Apar_ext,Jpar_ext);
  SAVE_ONCE2(Apar_ext,Jpar_ext);

  Phi0_ext=0;
  Upar0_ext=0;
  initial_profile("Phi0_ext", Phi0_ext);
  Upar0_ext = -Delp2(Phi0_ext)/Bxy;
  set_bndry(Upar0_ext,0.0,2);
  mesh->communicate(Phi0_ext,Upar0_ext);
  SAVE_ONCE2(Phi0_ext,Upar0_ext);

  // Give the solver the preconditioner function
  solver->setPrecon(precon);
  // Initialise parallel inversion class
  inv = InvertPar::Create();
  inv->setCoefA(1.0);
  Upar.setBoundary("Upar");
  Apar.setBoundary("Apar");
  
  return 0;
}

const Field3D Grad_parP_LtoC(const Field3D &f) {
  Field3D result;
  if(parallel_lc) {
    result = Grad_par_LtoC(f);
    if(nonlinear) {
      result -= beta_hat * bracket(Apar_ext+Apar, f, BRACKET_ARAKAWA);
    }else
      result -= beta_hat * bracket(Apar_ext, f, BRACKET_ARAKAWA);
  }else {
    if(nonlinear) {
      result = Grad_parP((Apar+Apar_ext)*beta_hat, f);
    }else {
      result = Grad_parP(Apar_ext*beta_hat, f);
    }
  }
  return result;
}

const Field3D Grad_parP_CtoL(const Field3D &f) {
  Field3D result;
  if(parallel_lc) {
    result = Grad_par_CtoL(f);
    if(nonlinear) {
      result -= beta_hat * bracket(Apar + Apar_ext, f, BRACKET_ARAKAWA);
    }else {
      result -= beta_hat * bracket(Apar_ext, f, BRACKET_ARAKAWA);
    }
  }else {
    if(nonlinear) {
      result = Grad_parP((Apar+Apar_ext)*beta_hat, f);
    }else {
      result = Grad_parP(Apar_ext*beta_hat, f);
    }
  }
  return result;
}

int physics_run(BoutReal t) {
  // Solve EM fields

  // Upar = (1/B) * Delp2(Phi)
  Phi = invert_laplace(Bxy*Upar, phi_flags);
  Phi.applyBoundary(); // For target plates only
  
  mesh->communicate(Upar, Phi, Apar);
  
  Jpar = -Delp2(Apar+Apar_ext); // Jpar includes external current
  Jpar.applyBoundary();
  mesh->communicate(Jpar);
  
///  if(jpar_bndry > 0) 
///    smooth_bndry(Jpar,jpar_bndry);

  // VORTICITY
  ddt(Upar) = SQ(Bxy)*Grad_parP_LtoC(Jpar/Bxy);

  if(include_jpar0) {
   ddt(Upar) -= SQ(Bxy)*beta_hat * bracket(Apar+Apar_ext, Jpar0/Bxy, BRACKET_ARAKAWA);
  }

    //ExB advection
  ddt(Upar) -= bracket(Phi0_ext, Upar, bm);   
 // ddt(Upar) -= bracket(Phi, Upar0_ext, bm);  
  if(nonlinear) {
    ddt(Upar) -= bracket(Phi, Upar, bm);  
  }
    //Viscosity
  if(mu > 0.)
    ddt(Upar) += mu*Delp2(Upar);

  // APAR
 //   ddt(Apar) = -Grad_parP_CtoL(Phi) / beta_hat;
    ddt(Apar) = -Grad_parP_CtoL(Phi+Phi0_ext) / beta_hat;

  if(eta > 0.)
    ddt(Apar) -= eta*Jpar / beta_hat;

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
int precon(BoutReal t, BoutReal gamma, BoutReal delta) {
  mesh->communicate(ddt(Apar));
  Field3D Jp;
  
  Jp = -Delp2(ddt(Apar));
  mesh->communicate(Jp);

///  if(jpar_bndry > 0)
///    smooth_bndry(Jp,jpar_bndry);
  
  Field3D Upar1 = ddt(Upar) + gamma*SQ(Bxy)*Grad_par_LtoC(Jp/Bxy);
  
  inv->setCoefB(-SQ(gamma*Bxy)/beta_hat);
  ddt(Upar) = inv->solve(Upar1);
  ddt(Upar).applyBoundary();
  
  Field3D Phip = invert_laplace(Bxy*ddt(Upar), phi_flags);
  mesh->communicate(Phip);
  
  ddt(Apar) = ddt(Apar) - (gamma / beta_hat)*Grad_par_CtoL(Phip);
  ddt(Apar).applyBoundary();

  return 0;
}

void smooth_bndry(Field3D f, int bndry)
{
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

void set_bndry(Field3D f, BoutReal val, int bndry)
{
    if(mesh->firstX()) {
      for(int i=bndry;i>=0;i--)
	  for(int j=0;j<mesh->LocalNy;j++)
	  for(int k=0;k<mesh->LocalNz;k++) {
	    f[i][j][k] = val;
	  }
    }
    if(mesh->lastX()) {
      for(int i=mesh->LocalNx-bndry-1;i<mesh->LocalNx;i++)
	  for(int j=0;j<mesh->LocalNy;j++)
	  for(int k=0;k<mesh->LocalNz;k++) {
	    f[i][j][k] = val;
	  }
    }
}

