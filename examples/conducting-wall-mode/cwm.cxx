/*******************************************************************************
 * Written by Brett Friedman: friedman@physics.ucla.edu
 * Linear Conducting Wall Mode Instability
 * Mode discoverd by H.L. Berk et. al. 1993
 * Model version in the code created by M. Umansky and J. Myra.
 *******************************************************************************/
#include <bout.hxx>
#include <boutmain.hxx>

#include <initialprofiles.hxx>
#include <derivs.hxx>
#include <interpolation.hxx>
#include <invert_laplace.hxx>

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

// 2D initial profiles
Field2D Ni0, Ti0, Te0, Vi0, phi0, Ve0, rho0, Ajpar0;
Vector2D b0xcv; // for curvature terms

// 3D evolving fields
Field3D rho, te, ajpar;

// Derived 3D variables
Field3D phi, jpar;

// e-i Collision frequency
Field3D nu;

// Phi boundary conditions
Field3D dphi_bc_ydown, dphi_bc_yup;

// Metric coefficients
Field2D Rxy, Bpxy, Btxy, hthe,Zxy;

// parameters
BoutReal Te_x, Ti_x, Ni_x, Vi_x, bmag, rho_s, fmei, AA, ZZ;
BoutReal lambda_ei, lambda_ii;
BoutReal nu_hat, wci, nueix;

bool bout_exb;  // Use BOUT-06 expression for ExB velocity

BoutReal zeff, nu_perp;
BoutReal ShearFactor;
BoutReal time_step;

bool filter_z;
int filter_z_mode;

int phi_flags; // Inversion flags

// Coefficients for linear sheath problem
Field2D LAMBDA1, LAMBDA2;

// My ixseps variables
int my_ixseps;

// Communication object
FieldGroup comms;

int physics_init(bool restarting)
{
  Field2D I; // Shear factor 
  
  output.write("Solving 6-variable 2-fluid equations\n");

  /************* LOAD DATA FROM GRID FILE ****************/

  // Load 2D profiles (set to zero if not found)
  mesh->get(Ni0,    "Ni0");
  mesh->get(Ti0,    "Ti0");
  mesh->get(Te0,    "Te0");
  mesh->get(Vi0,    "Vi0");
  mesh->get(Ve0,    "Ve0");
  mesh->get(phi0,   "phi0");
  mesh->get(rho0,   "rho0");
  mesh->get(Ajpar0, "Ajpar0");


  // Load magnetic curvature term
  b0xcv.covariant = false; // Read contravariant components
  mesh->get(b0xcv, "bxcv"); // b0xkappa terms

  // Load metrics
  mesh->get(Rxy,  "Rxy");
  mesh->get(Zxy,  "Zxy");
  mesh->get(Bpxy, "Bpxy");
  mesh->get(Btxy, "Btxy");
  mesh->get(hthe, "hthe");
  mesh->get(mesh->dx,   "dpsi");
  mesh->get(I,    "sinty");
  mesh->get(mesh->zShift, "qinty");

  // Load normalisation values
  mesh->get(Te_x, "Te_x");
  mesh->get(Ti_x, "Ti_x");
  mesh->get(Ni_x, "Ni_x");
  mesh->get(bmag, "bmag");

  // Get separatrix location
  mesh->get(my_ixseps, "ixseps1");


  Ni_x *= 1.0e14;
  bmag *= 1.0e4;


  /*************** READ OPTIONS *************************/
  // Read some parameters

  Options *globalOptions = Options::getRoot();

  globalOptions->get("TIMESTEP", time_step, 1.0);

  Options *options = globalOptions->getSection("2fluid");
  OPTION(options, AA, 4.0); // <=> options.get("AA", AA, 1.0);
  OPTION(options, ZZ, 1.0);

  OPTION(options, zeff,        1.0);
  OPTION(options, nu_perp,     0.0); 
  OPTION(options, ShearFactor, 1.0); 
  OPTION(options, bout_exb,    false);
  
  OPTION(options, phi_flags,   0);

  // Toroidal filtering
  OPTION(options, filter_z,          false);  // Filter a single n
  OPTION(options, filter_z_mode,     1);


  /************* SHIFTED RADIAL COORDINATES ************/

  if(mesh->ShiftXderivs) {
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
  nu_hat    = zeff*nueix/wci;
  
  Vi_x = wci * rho_s;

  output.write("Collisions: nueix = %e, nu_hat = %e\n", nueix, nu_hat);

  /************** PRINT Z INFORMATION ******************/
  
  BoutReal hthe0;
  if(mesh->get(hthe0, "hthe0") == 0) {
    output.write("    ****NOTE: input from BOUT, Z length needs to be divided by %e\n", hthe0/rho_s);
  }


  /************** NORMALISE QUANTITIES *****************/

  output.write("\tNormalising to rho_s = %e\n", rho_s);

  // Normalise profiles
  Ni0  /= Ni_x/1.0e14;
  Ti0  /= Te_x;
  Te0  /= Te_x;
  phi0 /= Te_x;
  Vi0  /= Vi_x;

  // Normalise curvature term
  b0xcv.x /= (bmag/1e4);
  b0xcv.y *= rho_s*rho_s;
  b0xcv.z *= rho_s*rho_s;
  
  // Normalise geometry 
  Rxy /= rho_s;
  hthe /= rho_s;
  I *= rho_s*rho_s*(bmag/1e4)*ShearFactor;
  mesh->dx /= rho_s*rho_s*(bmag/1e4);

  // Normalise magnetic field
  Bpxy /= (bmag/1.e4);
  Btxy /= (bmag/1.e4);
  mesh->Bxy  /= (bmag/1.e4);

  // Set nu
  nu = nu_hat * Ni0 / (Te0^1.5);

  /**************** CALCULATE METRICS ******************/

  mesh->g11 = (Rxy*Bpxy)^2;
  mesh->g22 = 1.0 / (hthe^2);
  mesh->g33 = (I^2)*mesh->g11 + (mesh->Bxy^2)/mesh->g11;
  mesh->g12 = 0.0;
  mesh->g13 = -I*mesh->g11;
  mesh->g23 = -Btxy/(hthe*Bpxy*Rxy);
  
  mesh->J = hthe / Bpxy;
  
  mesh->g_11 = 1.0/mesh->g11 + ((I*Rxy)^2);
  mesh->g_22 = (mesh->Bxy*hthe/Bpxy)^2;
  mesh->g_33 = Rxy*Rxy;
  mesh->g_12 = Btxy*hthe*I*Rxy/Bpxy;
  mesh->g_13 = I*Rxy*Rxy;
  mesh->g_23 = Btxy*hthe*Rxy/Bpxy;

  mesh->geometry();

  /**************** SET EVOLVING VARIABLES *************/

  // Tell BOUT++ which variables to evolve
  // add evolving variables to the communication object
  bout_solve(rho, "rho");
  comms.add(rho);
  
  bout_solve(te, "te");
  comms.add(te);
   
  /************** SETUP COMMUNICATIONS **************/

  // add extra variables to communication
  comms.add(phi);

  /*************** DUMP VARIABLES TO OUTPUT**********/
  dump.add(phi,  "phi",  1);

  SAVE_ONCE4(Ni0,Te0,phi0,rho0);
  SAVE_ONCE5(Rxy,Bpxy,Btxy,Zxy,hthe);
  dump.add(mesh->Bxy, "Bxy", 0);  dump.add(my_ixseps, "ixseps", 0);

  SAVE_ONCE3(Te_x,Ti_x,Ni_x);
  SAVE_ONCE6(AA,ZZ,zeff,rho_s,wci,bmag);
  dump.add(mesh->local_nx, "ngx", 0);
  dump.add(mesh->local_ny, "ngy", 0);
  dump.add(mesh->local_nz, "ngz", 0);
  SAVE_ONCE3(nu_hat,time_step,hthe0);
  
  return(0);
}
// End of physics_init()
//////////////////////////////////////////////////////////////////

// Routines for ExB terms (end of this file)
const Field3D vE_Grad(const Field2D &f, const Field3D &p);
const Field3D vE_Grad(const Field3D &f, const Field2D &p);
const Field3D vE_Grad(const Field3D &f, const Field3D &p);

// Boundary Functions
void phi_sheath_bndryconds();
void bndry_yup_Grad_par(Field3D &var, const Field3D &value);
void bndry_ydown_Grad_par(Field3D &var, const Field3D &value);

///////////////////////////////////////////////////////////////////
// Function called at each time step
// Time derivatives calculated here
int physics_run(BoutReal t)
{

  // Invert vorticity to get phi
  
  // Solves \nabla^2_\perp x + (1./c)*\nabla_perp c\cdot\nabla_\perp x + a x = b
  // Arguments are:   (b,   bit-field, a,    c)
  // Passing NULL -> missing term
  phi = invert_laplace(rho/Ni0, phi_flags, NULL, NULL);

  // Communicate variables
  mesh->communicate(comms);
 
  phi_sheath_bndryconds();
  
  // Evolve rho and te
  ddt(rho) = -Ni0/(fmei*0.51*nu)*Grad2_par2(phi);

  ddt(te) = -vE_Grad(Te0, phi);

  // Z filtering
  if(filter_z) {
    // Filter out all except filter_z_mode
    
    ddt(rho) = filter(ddt(rho), filter_z_mode);
    ddt(te) = filter(ddt(te), filter_z_mode);
  }


  return 0;
}
//End of physics_run
/////////////////////////////////////////////////////////////////



/****************BOUNDARY FUNCTIONS*****************************/
// Sheath Boundary Conditions on Phi
// Linearized
void phi_sheath_bndryconds()
{
  LAMBDA1 = 0.0;
  LAMBDA2 = 1.0;
  //LAMBDA1 = 1.0;
  //LAMBDA2 = log(2.0*sqrt(PI*fmei));

  
  dphi_bc_ydown = fmei*0.51*nu*(LAMBDA1*phi + LAMBDA2*te);
  dphi_bc_yup = -fmei*0.51*nu*(LAMBDA1*phi + LAMBDA2*te);

  bndry_ydown_Grad_par(phi,dphi_bc_ydown);
  bndry_yup_Grad_par(phi,dphi_bc_yup);
    
}

// Boundary gradient to specified Field3D object
void bndry_yup_Grad_par(Field3D &var, const Field3D &value)
{

  RangeIterator xrup = mesh->iterateBndryUpperY();

  for(xrup.first(); !xrup.isDone(); xrup.next())
    for(int jy=mesh->yend+1; jy<mesh->local_ny; jy++)
      for(int jz=0; jz<mesh->local_nz; jz++) {

	var[xrup.ind][jy][jz] = var[xrup.ind][jy-1][jz] + mesh->dy[xrup.ind][jy]*sqrt(mesh->g_22[xrup.ind][jy])*value[xrup.ind][jy][jz];

    }

}

void bndry_ydown_Grad_par(Field3D &var, const Field3D &value)
{

  RangeIterator xrdn = mesh->iterateBndryLowerY();

  for(xrdn.first(); !xrdn.isDone(); xrdn.next())
    for(int jy=mesh->ystart-1; jy>=0; jy--)
      for(int jz=0; jz<mesh->local_nz; jz++) {

	var[xrdn.ind][jy][jz] = var[xrdn.ind][jy+1][jz] - mesh->dy[xrdn.ind][jy]*sqrt(mesh->g_22[xrdn.ind][jy])*value[xrdn.ind][jy][jz];

    }

}

/////////////////////////////////////////////////////////////////
// ExB terms. These routines allow comparisons with BOUT-06
// if bout_exb=true is set in BOUT.inp
/////////////////////////////////////////////////////////////////
const Field3D vE_Grad(const Field2D &f, const Field3D &p)
{
  Field3D result;
  if(bout_exb) {
    // Use a subset of terms for comparison to BOUT-06
    result = VDDX(DDZ(p), f);
  }else {
    // Use full expression with all terms
    result = b0xGrad_dot_Grad(p, f) / mesh->Bxy;
  }
  return result;
}

const Field3D vE_Grad(const Field3D &f, const Field2D &p)
{
  Field3D result;
  if(bout_exb) {
    // Use a subset of terms for comparison to BOUT-06
    result = VDDZ(-DDX(p), f);
  }else {
    // Use full expression with all terms
    result = b0xGrad_dot_Grad(p, f) / mesh->Bxy;
  }
  return result;
}

const Field3D vE_Grad(const Field3D &f, const Field3D &p)
{
  Field3D result;
  if(bout_exb) {
    // Use a subset of terms for comparison to BOUT-06
    result = VDDX(DDZ(p), f) + VDDZ(-DDX(p), f);
  }else {
    // Use full expression with all terms
    result = b0xGrad_dot_Grad(p, f) / mesh->Bxy;
  }
  return result;
}
