/*******************************************************************************
 * 2-fluid equations
 * Same as Maxim's version of BOUT - simplified 2-fluid for benchmarking
 *******************************************************************************/

#include "bout.h"
#include "initialprofiles.h"
#include "derivs.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

// 2D initial profiles
Field2D Ni0, Ti0, Te0, Vi0;

// 3D evolving fields
Field3D  Te, Ni, Vi, Ti;

// 3D time-derivatives
Field3D  F_Te, F_Ni,  F_Vi, F_Ti;

// Non-linear coefficients
Field3D kapa_Te, kapa_Ti;

// 3D total values
Field3D Nit, Tit, Tet, Vit;

// pressures
Field3D peit, pe;
Field2D pei0, pe0;

// Metric coefficients
Field2D Rxy, Bpxy, Btxy, hthe;

// parameters
real Te_x, Ti_x, Ni_x, Vi_x, bmag, rho_s, fmei, AA, ZZ;
real lambda_ei, lambda_ii;
real nu_hat, mui_hat, wci, nueix, nuiix;

real chi_perp, D_perp, mu_perp;
real lambda_relax;

int physics_init()
{
  Field2D I; // Shear factor 
  
  output.write("Solving 6-variable 2-fluid equations\n");

  /************* LOAD DATA FROM GRID FILE ****************/

  // Load 2D profiles (set to zero if not found)
  GRID_LOAD(Ni0);
  GRID_LOAD(Ti0);
  GRID_LOAD(Te0);
  GRID_LOAD(Vi0);

  // Load metrics
  GRID_LOAD(Rxy);
  GRID_LOAD(Bpxy);
  GRID_LOAD(Btxy);
  GRID_LOAD(hthe);
  mesh->get(mesh->dx,   "dpsi");

  // Load normalisation values
  GRID_LOAD(Te_x);
  GRID_LOAD(Ti_x);
  GRID_LOAD(Ni_x);
  GRID_LOAD(bmag);

  Ni_x *= 1.0e14;
  bmag *= 1.0e4;

  /*************** READ OPTIONS *************************/

  // Read some parameters
  options.setSection("2fluid");
  OPTION(AA, 2.0);
  OPTION(ZZ, 1.0);

  OPTION(chi_perp,  0.6); // Read in m^2 / s 
  OPTION(D_perp,    0.6);
  OPTION(mu_perp,   0.6);

  OPTION(lambda_relax, 10.0);
  
  /************** CALCULATE PARAMETERS *****************/

  rho_s = 1.02*sqrt(AA*Te_x)/ZZ/bmag;
  fmei  = 1./1836.2/AA;

  lambda_ei = 24.-log(sqrt(Ni_x)/Te_x);
  lambda_ii = 23.-log(ZZ*ZZ*ZZ*sqrt(2.*Ni_x)/pow(Ti_x, 1.5));
  wci       = 9.58e3*ZZ*bmag/AA;
  nueix     = 2.91e-6*Ni_x*lambda_ei/pow(Te_x, 1.5);
  nuiix     = 4.78e-8*pow(ZZ,4.)*Ni_x*lambda_ii/pow(Ti_x, 1.5)/sqrt(AA);

  Vi_x = wci * rho_s;

  /************** PRINT Z INFORMATION ******************/
  
  real hthe0;
  if(GRID_LOAD(hthe0) == 0) {
    output.write("    ****NOTE: input from BOUT, Z length needs to be divided by %e\n", hthe0/rho_s);
  }

  /************** NORMALISE QUANTITIES *****************/

  output.write("\tNormalising to rho_s = %e\n", rho_s);

  // Normalise profiles
  Ni0 /= Ni_x/1.0e14;
  Ti0 /= Te_x;
  Te0 /= Te_x;
  Vi0 /= Vi_x;

   // Normalise geometry 
  Rxy /= rho_s;
  hthe /= rho_s;
  mesh->dx /= rho_s*rho_s*(bmag/1e4);

  // Normalise magnetic field
  Bpxy /= (bmag/1e4);
  Btxy /= (bmag/1e4);
  mesh->Bxy  /= (bmag/1e4);

  // calculate pressures
  pei0 = (Ti0 + Te0)*Ni0;
  pe0 = Te0*Ni0;

  // Normalise coefficients
  chi_perp /= rho_s*rho_s*wci;
  D_perp   /= rho_s*rho_s*wci;
  mu_perp  /= rho_s*rho_s*wci;

  chi_perp = 0.1;
  D_perp = 0.1;
  mu_perp = 0.1;
  
  output.write("Diffusion coefficients: chi %e D %e Mu %e\n",
	       chi_perp, D_perp, mu_perp);

  /**************** CALCULATE METRICS ******************/

  mesh->g11 = (Rxy*Bpxy)^2;
  mesh->g22 = 1.0 / (hthe^2);
  mesh->g33 = (mesh->Bxy^2)/mesh->g11;
  mesh->g12 = 0.0;
  mesh->g13 = 0.0;
  mesh->g23 = -Btxy/(hthe*Bpxy*Rxy);
  
  mesh->J = hthe / Bpxy;
  
  mesh->g_11 = 1.0/mesh->g11;
  mesh->g_22 = (mesh->Bxy*hthe/Bpxy)^2;
  mesh->g_33 = Rxy*Rxy;
  mesh->g_12 = 0.0;
  mesh->g_13 = 0.0;
  mesh->g_23 = Btxy*hthe*Rxy/Bpxy;


  /**************** SET EVOLVING VARIABLES *************/

  // Tell BOUT++ which variables to evolve
  // add evolving variables to the communication object

  Ni = Vi = Te = Ti = 0.0;
  bout_solve(Ni,    F_Ni,    "Ni");
  bout_solve(Vi,    F_Vi,    "Vi");
  bout_solve(Te,    F_Te,    "Te");
  bout_solve(Ti,    F_Ti,    "Ti");

  /************** SETUP COMMUNICATIONS **************/
  
  // Add any other variables to be dumped to file
  dump.add(Ni0, "Ni0", 0);
  dump.add(Te0, "Te0", 0);
  dump.add(Ti0, "Ti0", 0);

  dump.add(Te_x,  "Te_x", 0);
  dump.add(Ti_x,  "Ti_x", 0);
  dump.add(Ni_x,  "Ni_x", 0);
  dump.add(rho_s, "rho_s", 0);
  dump.add(wci,   "wci", 0);

  return(0);
}

// Operator for radial diffusive flux
/*
Field3D Div_X_K_Grad_X(const Field3D &difVi, const Field3D &Vi)
{
  Field2D sg = 1./sqrt(mesh->g_11);
  return difVi * D2DX2(Vi)/mesh->g_11
    + DDX( difVi * sg ) * DDX(Vi) * sg;
}
*/

// This version the same as in BOUT-06. Note the R's moved,and hthe added
Field3D Div_X_K_Grad_X(const Field3D &difVi, const Field3D &Vi)
{
  return difVi * ( (Rxy*Bpxy)^2 ) * D2DX2(Vi) 
    + (Bpxy / hthe) * DDX( difVi * Rxy * Rxy * Bpxy * hthe ) * DDX(Vi);
}

int physics_run(real t)
{
  //real bmk_t = MPI_Wtime();
  
  // Communicate variables
  mesh->communicate(Ni, Vi, Te, Ti);
  
  // Update profiles
  Nit = Ni0  + Ni.DC();
  Tit = Ti0  + Ti.DC();
  Tet = Te0  + Te.DC();
  Vit = Vi0  + Vi.DC();

  // Update non-linear coefficients on the mesh
  kapa_Te = 3.2*(1./fmei)*(wci/nueix)*(Tet^2.5);
  kapa_Ti = 3.9*(wci/nuiix)*(Tit^2.5);

  // note: nonlinear terms are not here
  peit = (Tet+Tit)*Nit;

  // DENSITY EQUATION
  F_Ni = -Vpar_Grad_par(Vit, Nit) 
    -Nit*Div_par(Vit) 
    +Div_X_K_Grad_X(D_perp*(Nit*0.0+1.0), Nit)
    ;


  // ION VELOCITY
  //F_Vi = -Grad_par(peit)/Nit -Vpar_Grad_par(Vit, Vit) + mu_perp*Delp2(Nit*Vit)/Nit;
  F_Vi = (
	  -Grad_par(peit) 
	  +Div_X_K_Grad_X(mu_perp*Nit, Vit)
	  )/Nit 
    -Vpar_Grad_par(Vit, Nit*Vit)/Nit 
    - F_Ni*Vit/Nit
    ;


  // ELECTRON TEMPERATURE
  F_Te = (Div_par_K_Grad_par(kapa_Te, Tet) 
	  +Div_X_K_Grad_X(chi_perp*Nit, Tet)
	  )/(1.5*Nit) 
    - F_Ni*Tet/Nit;


  // ION TEMPERATURE
  F_Ti = (Div_par_K_Grad_par(kapa_Ti, Tit) 
	  +Div_X_K_Grad_X(chi_perp*Nit, Tit)
	  )/(1.5*Nit)
    - F_Ni*Tit/Nit;


  // INNER TARGET PLATE
  
  bndry_ydown_flat(F_Ni); // Zero-gradient Ni
  bndry_ydown_relax_val(F_Vi, Vit, -3.095e4/Vi_x);
  bndry_ydown_relax_val(F_Te, Tet, 10./Te_x);
  bndry_ydown_relax_val(F_Ti, Tit, 10./Te_x);

  // OUTER TARGET PLATE

  bndry_yup_flat(F_Ni);
  bndry_yup_relax_val(F_Vi, Vit, 3.095e4/Vi_x);
  bndry_yup_relax_val(F_Te, Tet, 10./Te_x);
  bndry_yup_relax_val(F_Ti, Tit, 10./Te_x);
  
  // CORE BOUNDARY

  bndry_core_relax_val(F_Ni, Nit, 1e13/Ni_x);
  bndry_core_flat(F_Vi);
  bndry_core_relax_val(F_Te, Tet, 100./Te_x, lambda_relax);
  bndry_core_relax_val(F_Ti, Tit, 100./Te_x, lambda_relax);
  
  // PF BOUNDARY

  bndry_pf_relax_val(F_Ni, Nit, 1e12/Ni_x);
  bndry_pf_flat(F_Vi);
  bndry_pf_relax_val(F_Te, Tet, 10./Te_x);
  bndry_pf_relax_val(F_Ti, Tit, 10./Te_x);
  
  // OUTER BOUNDARY

  bndry_sol_relax_val(F_Ni, Nit, 1e12/Ni_x);
  bndry_sol_flat(F_Vi);
  bndry_sol_relax_val(F_Te, Tet, 10./Te_x);
  bndry_sol_relax_val(F_Ti, Tit, 10./Te_x);

  //output.write("TIMING: %e\n", MPI_Wtime() - bmk_t);

  return(0);
}
