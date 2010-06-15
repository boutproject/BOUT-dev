/*******************************************************************************
 * 2-fluid equations
 * Same as Maxim's version of BOUT - simplified 2-fluid for benchmarking
 *******************************************************************************/

#include "bout.h"
#include "initialprofiles.h"
#include "derivs.h"
#include "interpolation.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

// 2D initial profiles
Field2D Ni0, Ti0, Te0, Vi0, phi0, Ve0, rho0, Ajpar0;
Vector2D b0xcv; // for curvature terms

// 3D evolving fields
Field3D rho, Te, Ni, Ajpar, Vi, Ti;

// 3D time-derivatives
Field3D F_rho, F_Te, F_Ni, F_Ajpar, F_Vi, F_Ti;

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
real Te_x, Ti_x, Ni_x, Vi_x, bmag, rho_s, fmei, AA, ZZ;
real lambda_ei, lambda_ii;
real nu_hat, mui_hat, wci, nueix, nuiix;
real beta_p;

// settings
bool estatic, ZeroElMass; // Switch for electrostatic operation (true = no Apar)
real zeff, nu_perp;
bool evolve_rho, evolve_te, evolve_ni, evolve_ajpar, evolve_vi, evolve_ti;
real ShearFactor;

int phi_flags, apar_flags; // Inversion flags

// Communication object
Communicator comms;
Communicator com_jp;

// Field routines
int solve_phi_tridag(Field3D &r, Field3D &p, int flags);
int solve_apar_tridag(Field3D &aj, Field3D &ap, int flags);

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
  GRID_LOAD(Ve0);
  GRID_LOAD(phi0);
  GRID_LOAD(rho0);
  GRID_LOAD(Ajpar0);

  // Load magnetic curvature term
  b0xcv.covariant = false; // Read contravariant components
  grid.get(b0xcv, "bxcv"); // b0xkappa terms

  // Load metrics
  GRID_LOAD(Rxy);
  GRID_LOAD(Bpxy);
  GRID_LOAD(Btxy);
  GRID_LOAD(hthe);
  grid.get(dx,   "dpsi");
  grid.get(I,    "sinty");
  grid.get(zShift, "qinty");

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

  OPTION(estatic,     false);
  OPTION(ZeroElMass,  false);
  OPTION(zeff,        1.0);
  OPTION(nu_perp,     0.0);
  OPTION(ShearFactor, 1.0);
  
  OPTION(phi_flags,   0);
  OPTION(apar_flags,  0);
  
  options.get("rho",   "evolve", evolve_rho,   true);
  options.get("Te",    "evolve", evolve_te,    true);
  options.get("Ni",    "evolve", evolve_ni,    true);
  options.get("Ajpar", "evolve", evolve_ajpar, true);
  options.get("Vi",    "evolve", evolve_vi,    true);
  options.get("Ti",    "evolve", evolve_ti,    true);

  if(ZeroElMass)
    evolve_ajpar = 0; // Don't need ajpar - calculated from ohm's law

  /************* SHIFTED RADIAL COORDINATES ************/

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
  
  real hthe0;
  if(mesh->get(hthe0, "hthe0") == 0) {
    output.write("    ****NOTE: input from BOUT, Z length needs to be divided by %e\n", hthe0/rho_s);
  }

  /************** SHIFTED GRIDS LOCATION ***************/

  // Velocities defined on cell boundaries
  Vi.setLocation(CELL_YLOW);
  Ajpar.setLocation(CELL_YLOW);

  // Apar and jpar too
  Apar.setLocation(CELL_YLOW); 
  jpar.setLocation(CELL_YLOW);

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
  dx /= rho_s*rho_s*(bmag/1e4);

  // Normalise magnetic field
  Bpxy /= (bmag/1.e4);
  Btxy /= (bmag/1.e4);
  Bxy  /= (bmag/1.e4);

  // calculate pressures
  pei0 = (Ti0 + Te0)*Ni0;
  pe0 = Te0*Ni0;

  /**************** CALCULATE METRICS ******************/

  g11 = (Rxy*Bpxy)^2;
  g22 = 1.0 / (hthe^2);
  g33 = (I^2)*g11 + (Bxy^2)/g11;
  g12 = 0.0;
  g13 = -I*g11;
  g23 = -Btxy/(hthe*Bpxy*Rxy);
  
  J = hthe / Bpxy;
  
  g_11 = 1.0/g11 + ((I*Rxy)^2);
  g_22 = (Bxy*hthe/Bpxy)^2;
  g_33 = Rxy*Rxy;
  g_12 = Btxy*hthe*I*Rxy/Bpxy;
  g_13 = I*Rxy*Rxy;
  g_23 = Btxy*hthe*Rxy/Bpxy;

  // Twist-shift. NOTE: Should really use qsafe rather than qinty (small correction)

  if((jyseps2_2 / MYSUB) == MYPE) {
    for(int i=0;i<ngx;i++)
      ShiftAngle[i] = zShift[i][MYSUB]; // MYSUB+MYG-1
  }
  if(NYPE > 1)
    MPI_Bcast(ShiftAngle, ngx, PVEC_REAL_MPI_TYPE,jyseps2_2/MYSUB, MPI_COMM_WORLD);

  /**************** SET EVOLVING VARIABLES *************/

  // Tell BOUT++ which variables to evolve
  // add evolving variables to the communication object
  if(evolve_rho) {
    bout_solve(rho,   F_rho,   "rho");
    comms.add(rho);
    output.write("rho\n");
  }else
    initial_profile("rho", rho);

  if(evolve_ni) {
    bout_solve(Ni,    F_Ni,    "Ni");
    comms.add(Ni);
    output.write("ni\n");
  }else
    initial_profile("Ni", Ni);

  if(evolve_te) {
    bout_solve(Te,    F_Te,    "Te");
    comms.add(Te);
    
    output.write("te\n");
  }else
    initial_profile("Te", Te);

  if(evolve_ajpar) {
    bout_solve(Ajpar, F_Ajpar, "Ajpar");
    comms.add(Ajpar);
    output.write("ajpar\n");
  }else {
    initial_profile("Ajpar", Ajpar);
    if(ZeroElMass)
      dump.add(Ajpar, "Ajpar", 1); // output calculated Ajpar
  }

  if(evolve_vi) {
    bout_solve(Vi,    F_Vi,    "Vi");
    comms.add(Vi);
    output.write("vi\n");
  }else
    initial_profile("Vi", Vi);

  if(evolve_ti) {
    bout_solve(Ti,    F_Ti,    "Ti");
    comms.add(Ti);
    output.write("ti\n");
  }else
    initial_profile("Ti", Ti);
  
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

  
  //dump.add(F_Ni, "F_Ni", 1);
  //dump.add(F_rho, "F_rho", 1);
  //dump.add(F_Ajpar, "F_Ajpar", 1);

  com_jp.add(jpar);
  
  return(0);
}

// just define a macro for V_E dot Grad
#define vE_Grad(f, p) ( b0xGrad_dot_Grad(p, f) / Bxy )

int physics_run(real t)
{
  // Solve EM fields

  solve_phi_tridag(rho, phi, phi_flags);

  if(estatic || ZeroElMass) {
    // Electrostatic operation
    Apar = 0.0;
  }else {
    solve_apar_tridag(Ajpar, Apar, apar_flags); // Linear Apar solver
  }

  // Communicate variables
  comms.run();

  // zero-gradient Y boundaries (temporary! for interchange test)

  // Update profiles
  Nit = Ni0;  //+ Ni.DC();
  Tit = Ti0; // + Ti.DC();
  Tet = Te0; // + Te.DC();
  Vit = Vi0; // + Vi;

  // Update non-linear coefficients on the mesh
  nu      = nu_hat * Nit / (Tet^1.5);
  mu_i    = mui_hat * Nit / (Tit^0.5);
  kapa_Te = 3.2*(1./fmei)*(wci/nueix)*(Tet^2.5);
  kapa_Ti = 3.9*(wci/nuiix)*(Tit^2.5);
  
  // note: nonlinear terms are not here
  pei = (Te0+Ti0)*Ni + (Te + Ti)*Ni0;
  pe  = Te0*Ni + Te*Ni0;
  
  if(ZeroElMass) {
    // Set jpar,Ve,Ajpar neglecting the electron inertia term
    jpar = ((Te0*Grad_par(Ni, CELL_YLOW)) - (Ni0*Grad_par(phi, CELL_YLOW)))/(fmei*0.51*nu);

    /*
    for(int jx=MXG;jx<ngx-MXG;jx++) {
      for(int jy=MYG;jy<ngy-MYG;jy++) {
	for(int jz=0;jz<ngz;jz++) {
	  jpar[jx][jy][jz] = ( (Te0[jx][jy] * (Ni[jx][jy+1][jz] - Ni[jx][jy][jz]))
			       - (Ni0[jx][jy] * (phi[jx][jy+1][jz] - phi[jx][jy][jz])) )
	    / (fmei * 0.51 * nu[jx][jy][jz] * dy[jx][jy] * sqrt(g_22[jx][jy]));
			       
	}
      }
    }
    */

    // Set radial boundary condition on jpar
    bndry_inner_flat(jpar);
    bndry_sol_flat(jpar);
    bndry_toroidal(jpar);
    
    // Need to communicate jpar
    com_jp.run();

    Ve = Vi - jpar/Ni0;
    Ajpar = Ve;
  }else {
    
    Ve = Ajpar + Apar;
    jpar = Ni0*(Vi - Ve);
  }

  // DENSITY EQUATION

  F_Ni = 0.0;
  if(evolve_ni) {
    F_Ni -= vE_Grad(Ni0, phi);

    /*
      F_Ni -= vE_Grad(Ni, phi0) + vE_Grad(Ni0, phi) + vE_Grad(Ni, phi);
      F_Ni -= Vpar_Grad_par(Vi, Ni0) + Vpar_Grad_par(Vi0, Ni) + Vpar_Grad_par(Vi, Ni);
      F_Ni -= Ni0*Div_par(Vi) + Ni*Div_par(Vi0) + Ni*Div_par(Vi);
      F_Ni += Div_par(jpar);
      F_Ni += 2.0*V_dot_Grad(b0xcv, pe);
      F_Ni -= 2.0*(Ni0*V_dot_Grad(b0xcv, phi) + Ni*V_dot_Grad(b0xcv, phi0) + Ni*V_dot_Grad(b0xcv, phi));
    */
  }

  // ION VELOCITY

  F_Vi = 0.0;
  if(evolve_vi) {
    F_Vi -= vE_Grad(Vi0, phi) + vE_Grad(Vi, phi0) + vE_Grad(Vi, phi);
    F_Vi -= Vpar_Grad_par(Vi0, Vi) + Vpar_Grad_par(Vi, Vi0) + Vpar_Grad_par(Vi, Vi);
    F_Vi -= Grad_par(pei)/Ni0;
  }

  // ELECTRON TEMPERATURE

  F_Te = 0.0;
  if(evolve_te) {
    F_Te -= vE_Grad(Te0, phi) + vE_Grad(Te, phi0) + vE_Grad(Te, phi);
    F_Te -= Vpar_Grad_par(Ve, Te0) + Vpar_Grad_par(Ve0, Te) + Vpar_Grad_par(Ve, Te);
    F_Te += 1.333*Te0*( V_dot_Grad(b0xcv, pe)/Ni0 - V_dot_Grad(b0xcv, phi) );
    F_Te += 3.333*Te0*V_dot_Grad(b0xcv, Te);
    F_Te += (0.6666667/Ni0)*Div_par_K_Grad_par(kapa_Te, Te);
  }

  // ION TEMPERATURE

  F_Ti = 0.0;
  if(evolve_ti) {
    F_Ti -= vE_Grad(Ti0, phi) + vE_Grad(Ti, phi0) + vE_Grad(Ti, phi);
    F_Ti -= Vpar_Grad_par(Vi, Ti0) + Vpar_Grad_par(Vi0, Ti) + Vpar_Grad_par(Vi, Ti);
    F_Ti += 1.333*( Ti0*V_dot_Grad(b0xcv, pe)/Ni0 - Ti*V_dot_Grad(b0xcv, phi) );
    F_Ti -= 3.333*Ti0*V_dot_Grad(b0xcv, Ti);
    F_Ti += (0.6666667/Ni0)*Div_par_K_Grad_par(kapa_Ti, Ti);
  }

  // VORTICITY

  F_rho = 0.0;
  if(evolve_rho) {
    /*
    F_rho -= vE_Grad(rho0, phi) + vE_Grad(rho, phi0) + vE_Grad(rho, phi);
    F_rho -= Vpar_Grad_par(Vi, rho0) + Vpar_Grad_par(Vi0, rho) + Vpar_Grad_par(Vi, rho);
    */
    
    //F_rho += 2.0*Bnorm*V_dot_Grad(b0xcv, pei);

    F_rho += Bxy*Bxy*Div_par(jpar, CELL_CENTRE);

    /*
    for(int jx=MXG;jx<ngx-MXG;jx++) {
      for(int jy=MYG;jy<ngy-MYG;jy++) {
	for(int jz=0;jz<ngz;jz++) {
	  F_rho[jx][jy][jz] = Bxy[jx][jy]*Bxy[jx][jy] * (jpar[jx][jy+1][jz] - jpar[jx][jy][jz]) / (dy[jx][jy] * sqrt(g_22[jx][jy]));
	}
      }
    }
    */

    //F_rho += 1e-2 * mu_i * Laplacian(rho);
  }
  

  // AJPAR
  

  F_Ajpar = 0.0;
  if(evolve_ajpar) {
    //F_Ajpar -= vE_Grad(Ajpar0, phi) + vE_Grad(Ajpar, phi0) + vE_Grad(Ajpar, phi);

    /*
    for(int jx=MXG;jx<ngx-MXG;jx++) {
      for(int jy=MYG;jy<ngy-MYG;jy++) {
	for(int jz=0;jz<ngz;jz++) {
	  F_Ajpar[jx][jy][jz] += (1./fmei) * (phi[jx][jy][jz] - phi[jx][jy-1][jz]) / (dy[jx][jy] * sqrt(g_22[jx][jy]));
	  F_Ajpar[jx][jy][jz] -= (1./fmei)*(Te0[jx][jy]/Ni0[jx][jy])*(Ni[jx][jy][jz] - Ni[jx][jy-1][jz]) / (dy[jx][jy] * sqrt(g_22[jx][jy]));
	}
      }
    }
    */

    F_Ajpar += (1./fmei)*Grad_par(phi, CELL_YLOW);
    F_Ajpar -= (1./fmei)*(Te0/Ni0)*Grad_par(Ni, CELL_YLOW);
    //F_Ajpar -= (1./fmei)*1.71*Grad_par(Te);
    F_Ajpar += 0.51*interp_to(nu, CELL_YLOW)*jpar/Ni0;
  }

  // RADIAL BOUNDARY CONDITIONS

  apply_boundary(F_rho, "rho");
  apply_boundary(F_Te, "Te");
  apply_boundary(F_Ni, "Ni");
  apply_boundary(F_Ajpar, "Ajpar");
  apply_boundary(F_Vi, "Vi");
  apply_boundary(F_Ti, "Ti");

  return(0);
}

/*******************************************************************************
 *                       FAST LINEAR FIELD SOLVERS
 *******************************************************************************/

#include "invert_laplace.h"

// Performs inversion of rho (r) to get phi (p)
int solve_phi_tridag(Field3D &r, Field3D &p, int flags)
{
  //output.write("Solving phi: %e, %e -> %e\n", max(abs(r)), min(Ni0), max(abs(r/Ni0)));

  if(invert_laplace(r/Ni0, p, flags, NULL)) {
    return 1;
  }

  //Field3D pertPi = Ti*Ni0 + Ni*Ti0;
  //p -= pertPi/Ni0;
  return(0);
}

int solve_apar_tridag(Field3D &aj, Field3D &ap, int flags)
{
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
