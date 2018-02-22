/*******************************************************************************
 * UEDGE benchmark case
 *
 * Solves equations for 
 *  density Ni
 *  parallel ion velocity Vi
 *  electron and ion temperatures Te, Ti
 * 
 * Intended to be run for NZ=1 (i.e. X and Y only) for comparison with UEDGE
 *
 *******************************************************************************/

#include <bout.hxx>
#include <bout/boutmain.hxx>
#include <bout/derivs.hxx>

#include <cmath>

// 2D initial profiles
Field2D Ni0, Ti0, Te0, Vi0;

// 3D evolving fields
Field3D  Te, Ni, Vi, Ti;

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
BoutReal Te_x, Ti_x, Ni_x, Vi_x, bmag, rho_s, fmei, AA, ZZ;
BoutReal lambda_ei, lambda_ii;
BoutReal nu_hat, mui_hat, wci, nueix, nuiix;

BoutReal chi_perp, D_perp, mu_perp;

int physics_init(bool restarting)
{
  Field2D I; // Shear factor 
  
  output.write("Solving transport equations for Ni, Vi, Ti, Te\n");

  /////////////// LOAD DATA FROM GRID FILE //////////////

  // Load 2D profiles (set to zero if not found)
  GRID_LOAD(Ni0);
  GRID_LOAD(Ti0);
  GRID_LOAD(Te0);
  GRID_LOAD(Vi0);

  // Load metrics
  GRID_LOAD(Rxy);         // Major radius [m]
  GRID_LOAD2(Bpxy, Btxy); // Poloidal, Toroidal B field [T]
  GRID_LOAD(hthe);        // Poloidal arc length [m / radian]
  mesh->get(mesh->dx,   "dpsi");

  // Load normalisation values
  GRID_LOAD(Te_x);
  GRID_LOAD(Ti_x);
  GRID_LOAD(Ni_x);
  GRID_LOAD(bmag);

  Ni_x *= 1.0e14;
  bmag *= 1.0e4;

  /////////////// READ OPTIONS //////////////////////////

  // Read some parameters
  Options *globalOptions = Options::getRoot();
  Options *options = globalOptions->getSection("uedge");
  OPTION(options, AA, 2.0);
  OPTION(options, ZZ, 1.0);

  OPTION(options, chi_perp,  0.6); // Read in m^2 / s 
  OPTION(options, D_perp,    0.6);
  OPTION(options, mu_perp,   0.6);
  
  ////////////// CALCULATE PARAMETERS ///////////////////

  rho_s = 1.02*sqrt(AA*Te_x)/ZZ/bmag;
  fmei  = 1./1836.2/AA;

  lambda_ei = 24.-log(sqrt(Ni_x)/Te_x);
  lambda_ii = 23.-log(ZZ*ZZ*ZZ*sqrt(2.*Ni_x)/pow(Ti_x, 1.5));
  wci       = 9.58e3*ZZ*bmag/AA;
  nueix     = 2.91e-6*Ni_x*lambda_ei/pow(Te_x, 1.5);
  nuiix     = 4.78e-8*pow(ZZ,4.)*Ni_x*lambda_ii/pow(Ti_x, 1.5)/sqrt(AA);

  Vi_x = wci * rho_s;

  ///////////// PRINT Z INFORMATION /////////////////////
  
  BoutReal hthe0;
  if(GRID_LOAD(hthe0) == 0) {
    output.write("    ****NOTE: input from BOUT, Z length needs to be divided by %e\n", hthe0/rho_s);
  }

  ///////////// NORMALISE QUANTITIES ////////////////////

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

  /////////////// CALCULATE METRICS /////////////////

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
  
  mesh->geometry(); // Calculate other metrics
  
  //////////////// BOUNDARIES ///////////////////////
  // 
  // We want to apply the relaxing boundries to total density,
  // temperature etc.

  Ni0.applyBoundary("neumann");
  Te0.applyBoundary("neumann");
  Ti0.applyBoundary("neumann");

  Ni.setBackground(Ni0);
  Te.setBackground(Te0);
  Ti.setBackground(Ti0);
  
  ///////////// SET EVOLVING VARIABLES //////////////
  //
  // Tell BOUT++ which variables to evolve
  // add evolving variables to the communication object

  Ni = Vi = Te = Ti = 0.0;
  SOLVE_FOR4(Ni, Vi, Te, Ti);
  
  ///////////// ADD OUTPUT VARIABLES ////////////////
  //
  // Add any other variables to be dumped to file
  
  SAVE_ONCE3(Ni0, Te0, Ti0);                // Background quantities 
  SAVE_ONCE5(Te_x, Ti_x, Ni_x, rho_s, wci); // Normalisation factors

  /*
  dump.add(ddt(Ni), "ddt_ni", 1);
  dump.add(ddt(Ti), "ddt_ti", 1);
  dump.add(ddt(Te), "ddt_te", 1);
  dump.add(ddt(Vi), "ddt_vi", 1);
  */

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
Field3D Div_X_K_Grad_X(const Field3D &difFi, const Field3D &Fi)
{
  Field3D result;
  
  result = difFi * ( (Rxy*Bpxy)^2 ) * D2DX2(Fi) 
    + (Bpxy / hthe) * DDX( difFi * Rxy * Rxy * Bpxy * hthe ) * DDX(Fi);

  return result;
}

int physics_run(BoutReal t)
{
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
  
  peit = (Tet+Tit)*Nit;

  // DENSITY EQUATION
  ddt(Ni) = -Vpar_Grad_par(Vit, Nit) 
    -Nit*Div_par(Vit) 
    +Div_X_K_Grad_X(D_perp*(Nit*0.0+1.0), Nit)
    ;
  
  // ION VELOCITY
  ddt(Vi) = (
	     -Grad_par(peit) 
	     +Div_X_K_Grad_X(mu_perp*Nit, Vit)
	     )/Nit 
    -Vpar_Grad_par(Vit, Nit*Vit)/Nit 
    - ddt(Ni)*Vit/Nit
    ;
  
  // ELECTRON TEMPERATURE
  ddt(Te) = (Div_par_K_Grad_par(kapa_Te, Tet) 
	     +Div_X_K_Grad_X(chi_perp*Nit, Tet)
	     )/(1.5*Nit) 
    - ddt(Ni)*Tet/Nit;
  
  // ION TEMPERATURE
  ddt(Ti) = (Div_par_K_Grad_par(kapa_Ti, Tit) 
	     +Div_X_K_Grad_X(chi_perp*Nit, Tit)
	     )/(1.5*Nit)
    - ddt(Ni)*Tit/Nit;

  return(0);
}
