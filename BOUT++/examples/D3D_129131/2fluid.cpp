/*******************************************************************************
 * 2-fluid turbulence model
 * This version intended to have inputs as similar to BOUT-06 as possible
 * for cross-benchmarking etc.
 *******************************************************************************/

#include "bout.h"
#include "initialprofiles.h"
#include "derivs.h"
#include "interpolation.h"
#include "invert_laplace.h"
#include "topology.h"

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
real nuIonNeutral; // Ion-neutral collision rate (normalised by wci)

// settings
bool estatic, ZeroElMass; // Switch for electrostatic operation (true = no Apar)
real zeff, nu_perp;
bool evolve_rho, evolve_te, evolve_ni, evolve_ajpar, evolve_vi, evolve_ti;
real ShearFactor;

bool curv_upwind; // Use upwinding methods for curvature terms

bool laplace_extra_rho_term; // An extra first-order term in the vorticity inversion
bool vort_include_pi;    // Include Pi in vorticity

bool bout_exb;  // Use BOUT-06 expression for ExB velocity
bool bout_jpar; // Use BOUT-06 method for Jpar
bool OhmPe;     // Include the Pe term in Ohm's law

int bkgd;   // Profile options for coefficients (same options as BOUT-06)
int iTe_dc; // Profile evolution options

real lambda; // Boundary condition relaxation rate (negative)

bool stagger; // Use CtoL and LtoC for parallel derivs

// Switches for the equation terms
bool ni_ni1_phi0, ni_ni0_phi1, ni_ni1_phi1, ni_nit_phit;
bool ni_vi1_ni0,  ni_vi0_ni1,  ni_vi1_ni1,  ni_vit_nit;
bool ni_jpar1, ni_pe1, ni_ni1;
bool ni_ni0_curv_phi1, ni_ni1_curv_phi0, ni_ni1_curv_phi1, ni_nit_curv_phit;

bool rho_rho0_phi1, rho_rho1_phi0, rho_rho1_phi1;
bool rho_vi1_rho0, rho_vi0_rho1, rho_vi1_rho1;
bool rho_pei1, rho_jpar1, rho_rho1;

bool vi_vi0_phi1, vi_vi1_phi0, vi_vi1_phi1, vi_vit_phit;
bool vi_vi1_vi0, vi_vi0_vi1, vi_vi1_vi1, vi_vit_vit;
bool vi_pei1, vi_peit, vi_vi1;

bool te_te1_phi0, te_te0_phi1, te_te1_phi1;

bool ti_ti1_phi0, ti_ti0_phi1, ti_ti1_phi1;

bool relax_flat_bndry; // Use relaxing boundary conditions

int low_pass_z; // Low-pass filter result

int phi_flags, apar_flags; // Inversion flags

// Communication object
Communicator comms;
Communicator com_jp;

// BOUT-06 L1
const Field3D Div_par_CtoL(const Field3D &var)
{
  return Bxy * Grad_par_CtoL(var / Bxy);
}

int physics_init()
{
  Field2D I; // Shear factor 
  
  output.write("Solving 6-variable 2-fluid equations\n");

  ////////////////////////////////////////////////////////
  // LOAD DATA FROM GRID FILE

  // Load 2D profiles (set to zero if not found)
  grid_load2d(Ni0,    "Ni0");
  grid_load2d(Ti0,    "Ti0");
  grid_load2d(Te0,    "Te0");
  grid_load2d(Vi0,    "Vi0");
  grid_load2d(Ve0,    "Ve0");
  grid_load2d(phi0,   "phi0");
  grid_load2d(rho0,   "rho0");
  grid_load2d(Ajpar0, "Ajpar0");

  // Load magnetic curvature term
  b0xcv.covariant = false; // Read contravariant components
  grid_load2d(b0xcv, "bxcv"); // b0xkappa terms

  // Load metrics
  grid_load2d(Rxy,  "Rxy");
  grid_load2d(Bpxy, "Bpxy");
  grid_load2d(Btxy, "Btxy");
  grid_load2d(hthe, "hthe");
  grid_load2d(dx,   "dpsi");
  grid_load2d(I,    "sinty");
  grid_load2d(mesh->zShift, "qinty");

  // Load normalisation values
  grid_load(Te_x, "Te_x");
  grid_load(Ti_x, "Ti_x");
  grid_load(Ni_x, "Ni_x");
  grid_load(bmag, "bmag");

  Ni_x *= 1.0e14;
  bmag *= 1.0e4;

  ////////////////////////////////////////////////////////
  // READ OPTIONS

  // Read some parameters
  options.setSection("2fluid");
  options.get("AA", AA, 2.0);
  options.get("ZZ", ZZ, 1.0);

  OPTION(estatic,     false);
  OPTION(ZeroElMass,  false);
  OPTION(zeff,        1.0);
  OPTION(nu_perp,     0.0);
  OPTION(ShearFactor, 1.0);
  OPTION(OhmPe,       true);
  OPTION(bout_jpar,   false);
  OPTION(bout_exb,    false);
  OPTION(curv_upwind, false);

  OPTION(nuIonNeutral, -1.); 

  OPTION(bkgd,      2);
  OPTION(iTe_dc,    2);

  OPTION(stagger, false);

  OPTION(lambda,   -10.);
  if(lambda > 0.) {
    output.write("WARNING: lambda should be < 0. Reversing sign\n");
    lambda *= -1.0;
  }

  OPTION(relax_flat_bndry, true);

  OPTION(laplace_extra_rho_term, false);
  OPTION(vort_include_pi, false);

  options.get("low_pass_z",  low_pass_z,  -1);

  options.get("phi_flags",   phi_flags,   0);
  options.get("apar_flags",  apar_flags,  0);

  options.get("rho",   "evolve", evolve_rho,   true);
  options.get("Te",    "evolve", evolve_te,    true);
  options.get("Ni",    "evolve", evolve_ni,    true);
  options.get("Ajpar", "evolve", evolve_ajpar, true);
  options.get("Vi",    "evolve", evolve_vi,    true);
  options.get("Ti",    "evolve", evolve_ti,    true);
  
  if(ZeroElMass)
    evolve_ajpar = 0; // Don't need ajpar - calculated from ohm's law

  ////////////////////////////////////////////////////////
  // Equation terms

  if(evolve_ni) {
    options.setSection("Ni");
    options.get("ni1_phi0", ni_ni1_phi0, false);
    options.get("ni0_phi1", ni_ni0_phi1, false);
    options.get("ni1_phi1", ni_ni1_phi1, false);
    options.get("nit_phit", ni_nit_phit, false);
    options.get("vi1_ni0",  ni_vi1_ni0, false);
    options.get("vi0_ni1",  ni_vi0_ni1, false);
    options.get("vi1_ni1",  ni_vi1_ni1, false);
    options.get("vit_nit",  ni_vit_nit, false);
    options.get("jpar1",    ni_jpar1,  false);
    options.get("pe1",      ni_pe1,    false);
    options.get("ni0_curv_phi1", ni_ni0_curv_phi1, false);
    options.get("ni1_curv_phi0", ni_ni1_curv_phi0, false);
    options.get("ni1_curv_phi1", ni_ni1_curv_phi1, false);
    options.get("nit_curv_phit", ni_nit_curv_phit, false);
  }    

  if(evolve_rho) {
    options.setSection("rho");
    options.get("rho0_phi1", rho_rho0_phi1, false);
    options.get("rho1_phi0", rho_rho1_phi0, false);
    options.get("rho1_phi1", rho_rho1_phi1, false);
    options.get("vi1_rho0",  rho_vi1_rho0, false);
    options.get("vi0_rho1",  rho_vi0_rho1, false);
    options.get("vi1_rho1",  rho_vi1_rho1, false);
    options.get("pei1",   rho_pei1, false);
    options.get("jpar1",  rho_jpar1, false);
    options.get("rho1",   rho_rho1, false);
  }
  
  if(evolve_vi) {
    options.setSection("vi");
    options.get("vi0_phi1", vi_vi0_phi1, false);
    options.get("vi1_phi0", vi_vi1_phi0, false);
    options.get("vi1_phi1", vi_vi1_phi1, false);
    options.get("vit_phit", vi_vit_phit, false);
    options.get("vi1_vi0", vi_vi1_vi0, false);
    options.get("vi0_vi1", vi_vi0_vi1, false);
    options.get("vi1_vi1", vi_vi1_vi1, false);
    options.get("vit_vit", vi_vit_vit, false);
    options.get("pei1", vi_pei1, false);
    options.get("peit", vi_peit, false);
    options.get("vi1", vi_vi1, false);
  }

  if(evolve_te) {
    options.setSection("te");
    options.get("te1_phi0", te_te1_phi0, false);
    options.get("te0_phi1", te_te0_phi1, false);
    options.get("te1_phi1", te_te1_phi1, false);
  }

  if(evolve_ti) {
    options.setSection("ti");
    options.get("ti1_phi0", ti_ti1_phi0, false);
    options.get("ti0_phi1", ti_ti0_phi1, false);
    options.get("ti1_phi1", ti_ti1_phi1, false);
  }

  ////////////////////////////////////////////////////////
  // SHIFTED RADIAL COORDINATES

  if(mesh->ShiftXderivs) {
    ShearFactor = 0.0;  // I disappears from metric
    b0xcv.z += I*b0xcv.x;
  }

  ////////////////////////////////////////////////////////
  // CALCULATE PARAMETERS

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

  output.write("Collisions: nueix = %e, nu_hat = %e\n", nueix, nu_hat);

  ////////////////////////////////////////////////////////
  // PRINT Z INFORMATION
  
  real hthe0;
  if(grid_load(hthe0, "hthe0") == 0) {
    output.write("    ****NOTE: input from BOUT, Z length needs to be divided by %e\n", hthe0/rho_s);
  }

  ////////////////////////////////////////////////////////
  // SHIFTED GRIDS LOCATION

  // Velocities defined on cell boundaries
  Vi.setLocation(CELL_YLOW);
  Ajpar.setLocation(CELL_YLOW);

  // Apar and jpar too
  Apar.setLocation(CELL_YLOW); 
  jpar.setLocation(CELL_YLOW);

  ////////////////////////////////////////////////////////
  // NORMALISE QUANTITIES

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

  ////////////////////////////////////////////////////////
  // CALCULATE METRICS

  mesh->g11 = (Rxy*Bpxy)^2;
  mesh->g22 = 1.0 / (hthe^2);
  mesh->g33 = (I^2)*mesh->g11 + (Bxy^2)/mesh->g11;
  mesh->g12 = 0.0;
  mesh->g13 = -I*mesh->g11;
  mesh->g23 = -Btxy/(hthe*Bpxy*Rxy);
  
  J = hthe / Bpxy;
  
  mesh->g_11 = 1.0/mesh->g11 + ((I*Rxy)^2);
  mesh->g_22 = (Bxy*hthe/Bpxy)^2;
  mesh->g_33 = Rxy*Rxy;
  mesh->g_12 = Btxy*hthe*I*Rxy/Bpxy;
  mesh->g_13 = I*Rxy*Rxy;
  mesh->g_23 = Btxy*hthe*Rxy/Bpxy;
  
  ////////////////////////////////////////////////////////
  // Check twist-shift

  // Core
  if(YPROC(jyseps2_2) == PE_YIND) {
    for(int i=0;i<mesh->ngx;i++) {
      ShiftAngle[i] = mesh->zShift[i][MYG+MYSUB-1] - mesh->zShift[i][MYG+MYSUB]; // Jump across boundary
      //output.write("%d: %e\n", i, ShiftAngle[i]);
    }
  }else if(YPROC(jyseps1_1+1) == PE_YIND) {
    for(int i=0;i<mesh->ngx;i++) {
      ShiftAngle[i] = mesh->zShift[i][MYG-1] - mesh->zShift[i][MYG]; // Jump across boundary
      //output.write("%d: %e\n", i, ShiftAngle[i]);
    }
  }
  
  // Lower PF. Note by default no Twist-Shift used here, so need to switch on
  if(YPROC(jyseps1_1) == PE_YIND) {
    for(int i=0;i<mesh->ngx;i++) {
      ShiftAngle[i] = mesh->zShift[i][MYG+MYSUB-1] - mesh->zShift[i][MYG+MYSUB]; // Jump across boundary
      //output.write("%d: %e\n", i, ShiftAngle[i]);
    }
    TS_up_in = true; // Switch on twist-shift
    
  }else if(YPROC(jyseps2_2+1) == PE_YIND) {
    for(int i=0;i<mesh->ngx;i++) {
      ShiftAngle[i] = mesh->zShift[i][MYG-1] - mesh->zShift[i][MYG]; // Jump across boundary
      //output.write("%d: %e\n", i, ShiftAngle[i]);
    }
    TS_down_in = true;
  }

  // In the core, need to set ShiftAngle everywhere for ballooning initial condition
  MPI_Group group_world;
  MPI_Comm_group(MPI_COMM_WORLD, &group_world); // Group of all processors
  
  int *ranks = new int[NYPE];
  int npcore = 0;
  for(int p = YPROC(jyseps1_1+1); p <= YPROC(jyseps2_2);p++) {
    ranks[npcore] = PROC_NUM(PE_XIND, p);
    //output.write("%d: %d, %d\n", npcore, p, ranks[npcore]);
    npcore++;
  }
  
  MPI_Group grp;
  int ierr = MPI_Group_incl(group_world, npcore, ranks, &grp); // Create group
  //output.write("MPI_Group_incl: %d\n", ierr);
  
  MPI_Comm core_comm;
  MPI_Comm_create(MPI_COMM_WORLD, grp, &core_comm); // Create communicator
  
  delete[] ranks;

  if(MYPE_IN_CORE) {
    MPI_Bcast(ShiftAngle, mesh->ngx, PVEC_REAL_MPI_TYPE, npcore-1, core_comm);
  }

  //for(int i=0; i<mesh->ngx;i++)
  //  output.write("%d -> %e\n", i, ShiftAngle[i]);

  //MPI_Comm_free(&core_comm); // crashes
  //MPI_Group_free(&grp);

  ////////////////////////////////////////////////////////
  // SET EVOLVING VARIABLES

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

  if(!restarting) {
    // Smooth the initial perturbation a few times (includes comms)
    
    Ni = smooth_y(Ni);
    Ni = smooth_x(Ni);
    Ni = smooth_y(Ni);
    Ni = smooth_x(Ni);
    Ni = smooth_y(Ni);
    
    rho = smooth_y(rho);
    rho = smooth_x(rho);
    rho = smooth_y(rho);
    rho = smooth_x(rho);
    rho = smooth_y(rho);
    
    // Make sure initial perturbation obeys boundary condition
    
    if(relax_flat_bndry) {
      // Set all flat
      bndry_inner_flat(rho);
      bndry_sol_flat(rho);
      
      bndry_inner_flat(Te);
      bndry_sol_flat(Te);
      
      bndry_inner_flat(Ti);
      bndry_sol_flat(Ti);
      
      bndry_inner_flat(Ni);
      bndry_sol_flat(Ni);
      
      bndry_inner_flat(Ajpar);
      bndry_sol_flat(Ajpar);
      
      bndry_inner_flat(Vi);
      bndry_sol_flat(Vi);
    }else {
      apply_boundary(rho, "rho");
      apply_boundary(Te, "Te");
      apply_boundary(Ni, "Ni");
      apply_boundary(Ajpar, "Ajpar");
      apply_boundary(Vi, "Vi");
      apply_boundary(Ti, "Ti");
    }
  }
  
  ////////////////////////////////////////////////////////
  // SETUP COMMUNICATIONS

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
  
  com_jp.add(jpar);
  
  return(0);
}

// Routines for ExB terms (end of this file)
const Field2D vE_Grad(const Field2D &f, const Field2D &p);
const Field3D vE_Grad(const Field2D &f, const Field3D &p);
const Field3D vE_Grad(const Field3D &f, const Field2D &p);
const Field3D vE_Grad(const Field3D &f, const Field3D &p);

int physics_run(real t)
{
  //output.write("t = %e\n", t);
  
  ////////////////////////////////////////////////////////
  // Invert vorticity to get phi
  //
  // Solves \nabla^2_\perp x + \nabla_perp c\cdot\nabla_\perp x + a x = b
  // Arguments are:   (b,   bit-field, a,    c)
  // Passing NULL -> missing term
  
  if(laplace_extra_rho_term) {
    // Include the first order term Grad_perp Ni dot Grad_perp phi
    phi = invert_laplace(rho/Ni0, phi_flags, NULL, &Ni0);
  }else
    phi = invert_laplace(rho/Ni0, phi_flags);

  if(vort_include_pi) {
    // Include Pi term in vorticity
    phi -= (Ti*Ni0 + Ni*Te0) / Ni0;
  }

  ////////////////////////////////////////////////////////
  // Invert Ajpar to get Apar
  
  if(estatic || ZeroElMass) {
    // Electrostatic operation
    Apar = 0.0;
  }else {
    static Field2D acoeff;
    static bool aset = false;
    
    if(!aset) // calculate Apar coefficient
      acoeff = (-0.5*beta_p/fmei)*Ni0;
    aset = true;
  
    Apar = invert_laplace(acoeff*(Vi - Ajpar), apar_flags, &acoeff);
  }
  
  ////////////////////////////////////////////////////////
  // Communicate variables
  comms.run();

  ////////////////////////////////////////////////////////
  // Update profiles for calculating nu, mu_i, kapa_Te,i
  switch(bkgd) {
  case 0: { // Toroidal averages
    Nit = Ni0 + Ni.DC();
    Tit = Ti0 + Ti.DC();
    Tet = Te0 + Te.DC();
    Vit = Vi0 + Vi.DC();
    break;
  }
  case 1: { // Full perturbed values
    Nit = Ni0 + Ni;
    Tit = Ti0 + Ti;
    Tet = Te0 + Te;
    Vit = Vi0 + Vi;
    break; 
  }
  case 2: { // Unperturbed values
    Nit = Ni0;
    Tit = Ti0;
    Tet = Te0;
    Vit = Vi0;
    break;
  }
  default: {
    bout_error("ERROR: Invalid bkgd option\n");
  }
  }

  ////////////////////////////////////////////////////////
  // Update non-linear coefficients on the mesh
  nu      = nu_hat * Nit / (Tet^1.5);
  mu_i    = mui_hat * Nit / (Tit^0.5);
  kapa_Te = 3.2*(1./fmei)*(wci/nueix)*(Tet^2.5);
  kapa_Ti = 3.9*(wci/nuiix)*(Tit^2.5);
  
  // note: nonlinear terms are not here
  pei = (Te0+Ti0)*Ni + (Te + Ti)*Ni0;
  pe  = Te0*Ni + Te*Ni0;
  
  ////////////////////////////////////////////////////////
  if(ZeroElMass) {
    // Set jpar,Ve,Ajpar neglecting the electron inertia term

    if(!bout_jpar) {
      // Calculate Jpar, communicating across processors
      if(!stagger) {
	jpar = -(Ni0*Grad_par(phi, CELL_YLOW)) / (fmei*0.51*nu);
	
	if(OhmPe)
	  jpar += (Te0*Grad_par(Ni, CELL_YLOW)) / (fmei*0.51*nu);
      }else {
	jpar = -(Ni0*Grad_par_LtoC(phi))/(fmei*0.51*nu);
	
	if(OhmPe)
	  jpar += (Te0*Grad_par_LtoC(Ni)) / (fmei*0.51*nu);
      }
      
      // Set toroidal  boundary condition on jpar
      bndry_toroidal(jpar);
      
      // Need to communicate jpar
      com_jp.run();
      
    }else {
      // Use BOUT-06 method, no communications
      
      jpar.Allocate();
      for(int jx=0;jx<mesh->ngx;jx++)
	for(int jy=0;jy<mesh->ngy;jy++)
	  for(int jz=0;jz<mesh->ngz;jz++) {
	    real dNi_dpar, dPhi_dpar;
	  
	    // parallel derivs at left guard point
	    if (jy<jstart){
	      dNi_dpar=-1.5*Ni[jx][jy][jz] + 2.*Ni[jx][jy+1][jz] - 0.5*Ni[jx][jy+2][jz];
	      dPhi_dpar=-1.5*phi[jx][jy][jz] + 2.*phi[jx][jy+1][jz] - 0.5*phi[jx][jy+2][jz];
	    }else if (jy>jend){
	      // parallel derivs at right guard point
	      dNi_dpar = 1.5*Ni[jx][jy][jz] - 2.*Ni[jx][jy-1][jz] + 0.5*Ni[jx][jy-2][jz];
	      dPhi_dpar = 1.5*phi[jx][jy][jz] - 2.*phi[jx][jy-1][jz] + 0.5*phi[jx][jy-2][jz];
	    }else {
	      dNi_dpar = 0.5*Ni[jx][jy+1][jz] - 0.5*Ni[jx][jy-1][jz];
	      dPhi_dpar = 0.5*phi[jx][jy+1][jz] - 0.5*phi[jx][jy-1][jz];
	    }
	    
	    real c0=((Bpxy[jx][jy]/Bxy[jx][jy])/hthe[jx][jy])/dy[jx][jy];
	    dNi_dpar = dNi_dpar*c0;
	    dPhi_dpar = dPhi_dpar*c0;
	    
	    jpar[jx][jy][jz] = -(1./fmei)*(1./(0.51*nu[jx][jy][jz]))*Ni0[jx][jy]*dPhi_dpar;
	    if(OhmPe)
	      jpar[jx][jy][jz] += (1./fmei)*(1./(0.51*nu[jx][jy][jz]))*Te0[jx][jy]*dNi_dpar;
	      
	  }
    }
    
    apply_boundary(jpar, "jpar");
    
    Ve = Vi - jpar/Ni0;
    Ajpar = Ve;
  }else {
    
    Ve = Ajpar + Apar;
    jpar = Ni0*(Vi - Ve);
  }

  ////////////////////////////////////////////////////////
  // DENSITY EQUATION

  F_Ni = 0.0;
  if(evolve_ni) {
    
    if(ni_ni1_phi0)
      F_Ni -= vE_Grad(Ni, phi0);
    
    if(ni_ni0_phi1) 
      F_Ni -= vE_Grad(Ni0, phi);
   
    if(ni_ni1_phi1)
      F_Ni -= vE_Grad(Ni, phi);
    
    if(ni_nit_phit)
      F_Ni -= vE_Grad(Nit, phi0 + phi) - vE_Grad(Ni0, phi0);

    if(ni_vi1_ni0)
      F_Ni -= Vpar_Grad_par(Vi, Ni0);
    
    if(ni_vi0_ni1)
      F_Ni -= Vpar_Grad_par(Vi0, Ni);
    
    if(ni_vi1_ni1)
      F_Ni -= Vpar_Grad_par(Vi, Ni);

    if(ni_vit_nit)
      F_Ni -= Vpar_Grad_par(Vit, Nit) - Vpar_Grad_par(Vi0, Ni0);

    if(ni_jpar1) {
      if(stagger) {
	F_Ni += Div_par_CtoL(jpar);
      }else
	F_Ni += Div_par(jpar);
    }

    if(ni_pe1)
      F_Ni += 2.0*V_dot_Grad(b0xcv, pe);
    
    if(ni_ni0_curv_phi1)
      F_Ni -= 2.0*Ni0*V_dot_Grad(b0xcv, phi);
    
    if(ni_ni1_curv_phi0)
      F_Ni -= 2.0*Ni*V_dot_Grad(b0xcv, phi0);
    
    if(ni_ni1_curv_phi1)
      F_Ni -= 2.0*Ni*V_dot_Grad(b0xcv, phi);

    if(ni_nit_curv_phit)
      F_Ni -= 2.0*Nit*V_dot_Grad(b0xcv, phi+phi0) - 2.0*Ni0*V_dot_Grad(b0xcv, phi0);
    
    if(ni_ni1)
      F_Ni += mu_i * Delp2(Ni);
    
    //F_Ni -= Ni0*Div_par(Vi) + Ni*Div_par(Vi0) + Ni*Div_par(Vi);

    if(low_pass_z > 0)
      F_Ni = low_pass(F_Ni, low_pass_z);
  }

  ////////////////////////////////////////////////////////
  // ION VELOCITY

  F_Vi = 0.0;
  if(evolve_vi) {
    if(vi_vi0_phi1)
      F_Vi -= vE_Grad(Vi0, phi);

    if(vi_vi1_phi0)
      F_Vi -= vE_Grad(Vi, phi0);

    if(vi_vi1_phi1)
      F_Vi -= vE_Grad(Vi, phi);
    
    if(vi_vit_phit)
      F_Vi -= vE_Grad(Vit, phi+phi0) - vE_Grad(Vi0, phi+phi0);
    
    if(vi_vi1_vi0)
      F_Vi -= Vpar_Grad_par(Vi0, Vi);

    if(vi_vi0_vi1)
      F_Vi -= Vpar_Grad_par(Vi, Vi0);
    
    if(vi_vi1_vi1)
      F_Vi -= Vpar_Grad_par(Vi, Vi);

    if(vi_vit_vit)
      F_Vi -= Vpar_Grad_par(Vit, Vit) - Vpar_Grad_par(Vi0, Vi0);

    if(vi_pei1)
      F_Vi -= Grad_par(pei)/Ni0;

    if(vi_peit)
      F_Vi -= Grad_par(pei)/Nit;
    
    if(vi_vi1)
      F_Vi -= mu_i*Delp2(Vi);

    if(low_pass_z > 0)
      F_Vi = low_pass(F_Vi, low_pass_z);
  }

  ////////////////////////////////////////////////////////
  // ELECTRON TEMPERATURE

  F_Te = 0.0;
  if(evolve_te) {
    if(te_te1_phi0)
      F_Te -= vE_Grad(Te, phi0);
    if(te_te0_phi1)
      F_Te -= vE_Grad(Te0, phi);
    if(te_te1_phi1)
      F_Te -= vE_Grad(Te, phi);
    
    /*
    F_Te -= vE_Grad(Te0, phi) + vE_Grad(Te, phi0) + vE_Grad(Te, phi);
    F_Te -= Vpar_Grad_par(Ve, Te0) + Vpar_Grad_par(Ve0, Te) + Vpar_Grad_par(Ve, Te);
    F_Te += 1.333*Te0*( V_dot_Grad(b0xcv, pe)/Ni0 - V_dot_Grad(b0xcv, phi) );
    F_Te += 3.333*Te0*V_dot_Grad(b0xcv, Te);
    F_Te += (0.6666667/Ni0)*Div_par_K_Grad_par(kapa_Te, Te);

    */
    if(low_pass_z > 0)
      F_Te = low_pass(F_Te, low_pass_z);
  }

  ////////////////////////////////////////////////////////
  // ION TEMPERATURE

  F_Ti = 0.0;
  if(evolve_ti) {
    if(ti_ti1_phi0)
      F_Ti -= vE_Grad(Ti, phi0);
    if(ti_ti0_phi1)
      F_Ti -= vE_Grad(Ti0, phi);
    if(ti_ti1_phi1)
      F_Ti -= vE_Grad(Ti, phi);
    
    /*
    F_Ti -= vE_Grad(Ti0, phi) + vE_Grad(Ti, phi0) + vE_Grad(Ti, phi);
    F_Ti -= Vpar_Grad_par(Vi, Ti0) + Vpar_Grad_par(Vi0, Ti) + Vpar_Grad_par(Vi, Ti);
    F_Ti += 1.333*( Ti0*V_dot_Grad(b0xcv, pe)/Ni0 - Ti*V_dot_Grad(b0xcv, phi) );
    F_Ti -= 3.333*Ti0*V_dot_Grad(b0xcv, Ti);
    F_Ti += (0.6666667/Ni0)*Div_par_K_Grad_par(kapa_Ti, Ti);
    */

    if(low_pass_z > 0)
      F_Ti = low_pass(F_Ti, low_pass_z);
  }

  ////////////////////////////////////////////////////////
  // VORTICITY

  F_rho = 0.0;
  if(evolve_rho) {
    
    if(rho_rho0_phi1)
      F_rho -= vE_Grad(rho0, phi);

    if(rho_rho1_phi0)
      F_rho -= vE_Grad(rho, phi0);

    if(rho_rho1_phi1)
      F_rho -= vE_Grad(rho, phi);

    if(rho_vi1_rho0)
      F_rho -= Vpar_Grad_par(Vi, rho0);
    
    if(rho_vi0_rho1)
      F_rho -= Vpar_Grad_par(Vi0, rho);
    
    if(rho_vi1_rho1)
      F_rho -= Vpar_Grad_par(Vi, rho);
    
    if(rho_pei1) {
      if(curv_upwind) {
	F_rho += 2.0*Bxy*V_dot_Grad(b0xcv, pei);  // Use upwinding
      }else
	F_rho += 2.0*Bxy*b0xcv*Grad(pei);     // Use central differencing
    }    

    if(rho_jpar1) {
      if(stagger) {
	F_rho += Bxy*Bxy*Div_par_CtoL(jpar);
      }else 
	F_rho += Bxy*Bxy*Div_par(jpar, CELL_CENTRE);
    }

    if(rho_rho1)
      F_rho += mu_i * Delp2(rho);

    if(low_pass_z > 0)
      F_rho = low_pass(F_rho, low_pass_z);
  }
  
  ////////////////////////////////////////////////////////
  // AJPAR
  
  F_Ajpar = 0.0;
  if(evolve_ajpar) {
    //F_Ajpar -= vE_Grad(Ajpar0, phi) + vE_Grad(Ajpar, phi0) + vE_Grad(Ajpar, phi);
    //F_Ajpar -= (1./fmei)*1.71*Grad_par(Te);
    
    if(stagger) {
      F_Ajpar += (1./fmei)*Grad_par_LtoC(phi); // Right-hand differencing
    }else
      F_Ajpar += (1./fmei)*Grad_par(phi, CELL_YLOW);
    
    if(OhmPe) {
      if(stagger) {
	F_Ajpar -= (1./fmei)*(Tet/Nit)*Grad_par_LtoC(Ni);
      }else 
	F_Ajpar -= (1./fmei)*(Te0/Ni0)*Grad_par(Ni, CELL_YLOW);
    }
    
    F_Ajpar += 0.51*interp_to(nu, CELL_YLOW)*jpar/Ni0;

    if(low_pass_z > 0)
      F_Ajpar = low_pass(F_Ajpar, low_pass_z);
  }

  ////////////////////////////////////////////////////////
  // Profile evolution options

  switch(iTe_dc) {
  case 1: { // subtacting out toroidal averages for all fields
    if(evolve_ni)
      F_Ni -= F_Ni.DC();
    if(evolve_rho)
      F_rho -= F_rho.DC();
    if(evolve_te)
      F_Te -= F_Te.DC();
    if(evolve_ti)
      F_Ti -= F_Ti.DC();
    if(evolve_ajpar)
      F_Ajpar -= F_Ajpar.DC();
    break;
  }
  case 2: { // not subtacting out toroidal averages for any field
    break;
  }
  case 4: { // using toroidal averages in right-hand sides, e.g., axisymmetric mode
    if(evolve_ni)
      F_Ni = F_Ni.DC();
    if(evolve_rho)
      F_rho = F_rho.DC();
    if(evolve_te)
      F_Te = F_Te.DC();
    if(evolve_ti)
      F_Ti = F_Ti.DC();
    if(evolve_ajpar)
      F_Ajpar = F_Ajpar.DC();
    break;
  }
  default: {
    bout_error("ERROR: invalid option for iTe_dc\n");
  }
  }
  
  ////////////////////////////////////////////////////////
  // RADIAL BOUNDARY CONDITIONS

  if(relax_flat_bndry) {
    // BOUT-06 style relaxing boundary conditions
    
    // Zero-gradient at target plates 
    
    if(evolve_rho) {
      bndry_inner_relax_flat(F_rho, rho, lambda);
      bndry_sol_relax_flat(F_rho, rho, lambda);
      
      bndry_ydown_flat(F_rho);
      bndry_yup_flat(F_rho);
    }
      
    if(evolve_ni) {
      bndry_inner_relax_flat(F_Ni, Ni, lambda);
      bndry_sol_relax_flat(F_Ni, Ni, lambda);
      
      bndry_ydown_flat(F_Ni);
      bndry_yup_flat(F_Ni);
    }

    if(evolve_te) {
      bndry_inner_relax_flat(F_Te, Te, lambda);
      bndry_sol_relax_flat(F_Te, Te, lambda);

      bndry_ydown_flat(F_Te);
      bndry_yup_flat(F_Te);
    }
    
    if(evolve_ti) {
      bndry_inner_relax_flat(F_Ti, Ti, lambda);
      bndry_sol_relax_flat(F_Ti, Ti, lambda);

      bndry_ydown_flat(F_Ti);
      bndry_yup_flat(F_Ti);
    }

    if(evolve_ajpar) {
      bndry_inner_relax_flat(F_Ajpar, Ajpar, lambda);
      bndry_sol_relax_flat(F_Ajpar, Ajpar, lambda);
      
      bndry_ydown_flat(F_Ajpar);
      bndry_yup_flat(F_Ajpar);
    }
    
  }else {
    // Use the boundary condition specified in BOUT.inp
    
    apply_boundary(F_rho, "rho");
    apply_boundary(F_Te, "Te");
    apply_boundary(F_Ni, "Ni");
    apply_boundary(F_Ajpar, "Ajpar");
    apply_boundary(F_Vi, "Vi");
    apply_boundary(F_Ti, "Ti");
  }

  return(0);
}

/////////////////////////////////////////////////////////////////
// ExB terms. These routines allow comparisons with BOUT-06
// if bout_exb=true is set in BOUT.inp
/////////////////////////////////////////////////////////////////

const Field2D vE_Grad(const Field2D &f, const Field2D &p)
{
  Field2D result;
  if(bout_exb) {
    // Use a subset of terms for comparison to BOUT-06
    result = 0.0;
  }else {
    // Use full expression with all terms
    result = b0xGrad_dot_Grad(p, f) / Bxy;
  }
  return result;
}

const Field3D vE_Grad(const Field2D &f, const Field3D &p)
{
  Field3D result;
  if(bout_exb) {
    // Use a subset of terms for comparison to BOUT-06
    result = VDDX(DDZ(p), f);
    //result = DDX(DDZ(p)*f);
  }else {
    // Use full expression with all terms
    result = b0xGrad_dot_Grad(p, f) / Bxy;
  }
  return result;
}

const Field3D vE_Grad(const Field3D &f, const Field2D &p)
{
  Field3D result;
  if(bout_exb) {
    // Use a subset of terms for comparison to BOUT-06
    result = VDDZ(-DDX(p), f);
    //result = DDZ(-DDX(p) * f);
  }else {
    // Use full expression with all terms
    result = b0xGrad_dot_Grad(p, f) / Bxy;
  }
  return result;
}

const Field3D vE_Grad(const Field3D &f, const Field3D &p)
{
  Field3D result;
  if(bout_exb) {
    // Use a subset of terms for comparison to BOUT-06
    result = VDDX(DDZ(p), f) + VDDZ(-DDX(p), f);
    //result = DDX(DDZ(p) * f) + DDZ(-DDX(p) * f);
  }else {
    // Use full expression with all terms
    result = b0xGrad_dot_Grad(p, f) / Bxy;
  }
  return result;
}
