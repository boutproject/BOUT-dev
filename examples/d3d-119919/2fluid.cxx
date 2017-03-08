/*******************************************************************************
 * 2-fluid turbulence model
 * This version intended to have inputs as similar to BOUT-06 as possible
 * for cross-benchmarking etc.
 *******************************************************************************/

#include <bout.hxx>
#include <boutmain.hxx>

#include <initialprofiles.hxx>
#include <derivs.hxx>
#include <interpolation.hxx>
#include <invert_laplace.hxx>
#include <msg_stack.hxx>

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

// 2D initial profiles
Field2D Ni0, Ti0, Te0, Vi0, phi0, Ve0, rho0, Ajpar0;
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
BoutReal nuIonNeutral; // Ion-neutral collision rate (normalised by wci)

// settings
bool estatic, ZeroElMass; // Switch for electrostatic operation (true = no Apar)
BoutReal zeff, nu_perp;
bool evolve_rho, evolve_te, evolve_ni, evolve_ajpar, evolve_vi, evolve_ti;
BoutReal ShearFactor;

bool curv_upwind; // Use upwinding methods for curvature terms

bool laplace_extra_rho_term; // An extra first-order term in the vorticity inversion
bool vort_include_pi;    // Include Pi in vorticity

bool bout_jpar; // Use BOUT-06 method for Jpar
bool OhmPe;     // Include the Pe term in Ohm's law

int bkgd;   // Profile options for coefficients (same options as BOUT-06)
int iTe_dc; // Profile evolution options

bool stagger; // Use CtoL and LtoC for parallel derivs

// Poisson brackets: b0 x Grad(f) dot Grad(g) / B = [f, g]
// Method to use: BRACKET_ARAKAWA, BRACKET_STD or BRACKET_SIMPLE
BRACKET_METHOD bm; // Bracket method for advection terms

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

int lowPass_z; // Low-pass filter result

int phi_flags, apar_flags; // Inversion flags

// Group of objects for communications
FieldGroup comms;

int physics_init(bool restarting) {
  Field2D I; // Shear factor 
  
  output.write("Solving 6-variable 2-fluid equations\n");

  ////////////////////////////////////////////////////////
  // LOAD DATA FROM GRID FILE

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
  mesh->get(Bpxy, "Bpxy");
  mesh->get(Btxy, "Btxy");
  mesh->get(hthe, "hthe");
  
  Field2D dx;
  if(!mesh->get(dx,   "dpsi")) {
    output << "\tUsing dpsi as the x grid spacing\n";
    mesh->dx = dx; // Only use dpsi if found
  }else {
    // dx will have been read already from the grid
    output << "\tUsing dx as the x grid spacing\n";
  }
  mesh->get(I,    "sinty");
  Field2D qinty;
  if(!mesh->get(qinty, "qinty")) {
    output << "\tUsing qinty as the Z shift\n";
    mesh->zShift = qinty;
  }else {
    // Keep zShift
    output << "\tUsing zShift as the Z shift\n";
  }

  // Load normalisation values
  mesh->get(Te_x, "Te_x");
  mesh->get(Ti_x, "Ti_x");
  mesh->get(Ni_x, "Ni_x");
  mesh->get(bmag, "bmag");

  Ni_x *= 1.0e14;
  bmag *= 1.0e4;

  ////////////////////////////////////////////////////////
  // READ OPTIONS

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
  OPTION(options, OhmPe,       true);
  OPTION(options, bout_jpar,   false);
  OPTION(options, curv_upwind, false);

  // Choose method to use for Poisson bracket advection terms
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

  OPTION(options, nuIonNeutral, -1.); 

  OPTION(options, bkgd,      2);
  OPTION(options, iTe_dc,    2);

  OPTION(options, stagger, false);

  OPTION(options, laplace_extra_rho_term, false);
  OPTION(options, vort_include_pi, false);

  options->get("lowPass_z",  lowPass_z,  -1);

  options->get("phi_flags",   phi_flags,   0);
  options->get("apar_flags",  apar_flags,  0);
  
  (globalOptions->getSection("Ni"))->get("evolve", evolve_ni,    true);
  (globalOptions->getSection("rho"))->get("evolve", evolve_rho,   true);
  (globalOptions->getSection("vi"))->get("evolve", evolve_vi,   true);
  (globalOptions->getSection("te"))->get("evolve", evolve_te,   true);
  (globalOptions->getSection("ti"))->get("evolve", evolve_ti,   true);
  (globalOptions->getSection("Ajpar"))->get("evolve", evolve_ajpar, true);
  
  if(ZeroElMass)
    evolve_ajpar = 0; // Don't need ajpar - calculated from ohm's law

  ////////////////////////////////////////////////////////
  // Equation terms

  if(evolve_ni) {
    options = globalOptions->getSection("Ni");
    options->get("ni1_phi0", ni_ni1_phi0, false);
    options->get("ni0_phi1", ni_ni0_phi1, false);
    options->get("ni1_phi1", ni_ni1_phi1, false);
    options->get("nit_phit", ni_nit_phit, false);
    options->get("vi1_ni0",  ni_vi1_ni0, false);
    options->get("vi0_ni1",  ni_vi0_ni1, false);
    options->get("vi1_ni1",  ni_vi1_ni1, false);
    options->get("vit_nit",  ni_vit_nit, false);
    options->get("jpar1",    ni_jpar1,  false);
    options->get("pe1",      ni_pe1,    false);
    options->get("ni0_curv_phi1", ni_ni0_curv_phi1, false);
    options->get("ni1_curv_phi0", ni_ni1_curv_phi0, false);
    options->get("ni1_curv_phi1", ni_ni1_curv_phi1, false);
    options->get("nit_curv_phit", ni_nit_curv_phit, false);
  }    

  if(evolve_rho) {
    options = globalOptions->getSection("rho");
    options->get("rho0_phi1", rho_rho0_phi1, false);
    options->get("rho1_phi0", rho_rho1_phi0, false);
    options->get("rho1_phi1", rho_rho1_phi1, false);
    options->get("vi1_rho0",  rho_vi1_rho0, false);
    options->get("vi0_rho1",  rho_vi0_rho1, false);
    options->get("vi1_rho1",  rho_vi1_rho1, false);
    options->get("pei1",   rho_pei1, false);
    options->get("jpar1",  rho_jpar1, false);
    options->get("rho1",   rho_rho1, false);
  }
  
  if(evolve_vi) {
    options = globalOptions->getSection("vi");
    options->get("vi0_phi1", vi_vi0_phi1, false);
    options->get("vi1_phi0", vi_vi1_phi0, false);
    options->get("vi1_phi1", vi_vi1_phi1, false);
    options->get("vit_phit", vi_vit_phit, false);
    options->get("vi1_vi0", vi_vi1_vi0, false);
    options->get("vi0_vi1", vi_vi0_vi1, false);
    options->get("vi1_vi1", vi_vi1_vi1, false);
    options->get("vit_vit", vi_vit_vit, false);
    options->get("pei1", vi_pei1, false);
    options->get("peit", vi_peit, false);
    options->get("vi1", vi_vi1, false);
  }

  if(evolve_te) {
    options = globalOptions->getSection("te");
    options->get("te1_phi0", te_te1_phi0, false);
    options->get("te0_phi1", te_te0_phi1, false);
    options->get("te1_phi1", te_te1_phi1, false);
  }

  if(evolve_ti) {
    options = globalOptions->getSection("ti");
    options->get("ti1_phi0", ti_ti1_phi0, false);
    options->get("ti0_phi1", ti_ti0_phi1, false);
    options->get("ti1_phi1", ti_ti1_phi1, false);
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
  
  BoutReal hthe0;
  if(mesh->get(hthe0, "hthe0") == 0) {
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
  mesh->dx /= rho_s*rho_s*(bmag/1e4);

  // Normalise magnetic field
  Bpxy /= (bmag/1.e4);
  Btxy /= (bmag/1.e4);
  mesh->Bxy  /= (bmag/1.e4);

  // calculate pressures
  pei0 = (Ti0 + Te0)*Ni0;
  pe0 = Te0*Ni0;

  ////////////////////////////////////////////////////////
  // CALCULATE METRICS

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

  ////////////////////////////////////////////////////////
  // SET EVOLVING VARIABLES

  // Tell BOUT++ which variables to evolve
  // add evolving variables to the communication object
  if(evolve_rho) {
    SOLVE_FOR(rho);
    comms.add(rho);
    output.write("rho\n");
  }else {
    initial_profile("rho", rho);
    rho.setBoundary("rho");
    rho.applyBoundary();
  }

  if(evolve_ni) {
    SOLVE_FOR(Ni);
    comms.add(Ni);
    output.write("ni\n");
  }else {
    initial_profile("Ni", Ni);
    Ni.setBoundary("Ni");
    Ni.applyBoundary();
  }

  if(evolve_te) {
    SOLVE_FOR(Te);
    comms.add(Te);
    
    output.write("te\n");
  }else {
    initial_profile("Te", Te);
    Te.setBoundary("Te");
    Te.applyBoundary();
  }

  if(evolve_ajpar) {
    SOLVE_FOR(Ajpar);
    comms.add(Ajpar);
    output.write("ajpar\n");
  }else {
    initial_profile("Ajpar", Ajpar);
    if(ZeroElMass)
      dump.add(Ajpar, "Ajpar", 1); // output calculated Ajpar
    
    Ajpar.setBoundary("Ajpar");
    Ajpar.applyBoundary();
  }

  if(evolve_vi) {
    SOLVE_FOR(Vi);
    comms.add(Vi);
    output.write("vi\n");
  }else {
    initial_profile("Vi", Vi);
    Vi.setBoundary("Vi");
    Vi.applyBoundary();
  }

  if(evolve_ti) {
    SOLVE_FOR(Ti);
    comms.add(Ti);
    output.write("ti\n");
  }else {
    initial_profile("Ti", Ti);
    Ti.setBoundary("Ti");
    Ti.applyBoundary();
  }

  jpar.setBoundary("jpar");
  
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
  
  return(0);
}

// ExB terms using Poisson bracket
#define vE_Grad(f, p) ( bracket(p, f, bm) )

int physics_run(BoutReal t) {
  
  ////////////////////////////////////////////////////////
  // Invert vorticity to get phi
  //
  // Solves \nabla^2_\perp x + \nabla_perp c\cdot\nabla_\perp x + a x = b
  // Arguments are:   (b,   bit-field, a,    c)
  // Passing NULL -> missing term
  
  msg_stack.push("Solving for phi");
  if(laplace_extra_rho_term) {
    // Include the first order term Grad_perp Ni dot Grad_perp phi
    phi = invert_laplace(rho/Ni0, phi_flags, NULL, &Ni0);
  }else
    phi = invert_laplace(rho/Ni0, phi_flags);
  
  if(vort_include_pi) {
    // Include Pi term in vorticity
    phi -= (Ti*Ni0 + Ni*Te0) / Ni0;
  }
  msg_stack.pop();
  

  ////////////////////////////////////////////////////////
  // Invert Ajpar to get Apar
  
  msg_stack.push("Solving for Apar");
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
  msg_stack.pop();
  
  ////////////////////////////////////////////////////////
  // Communicate variables
  mesh->communicate(comms);

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
    throw BoutException("ERROR: Invalid bkgd option\n");
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
      
      // Need to communicate jpar
      mesh->communicate(jpar);
      
    }else {
      // Use BOUT-06 method, no communications
      
      jpar.allocate();
      for(int jx=0;jx<mesh->localNx;jx++)
	for(int jy=0;jy<mesh->localNy;jy++)
	  for(int jz=0;jz<mesh->localNz;jz++) {
	    BoutReal dNi_dpar, dPhi_dpar;
	  
	    // parallel derivs at left guard point
	    if (jy<mesh->ystart){
	      dNi_dpar=-1.5*Ni[jx][jy][jz] + 2.*Ni[jx][jy+1][jz] - 0.5*Ni[jx][jy+2][jz];
	      dPhi_dpar=-1.5*phi[jx][jy][jz] + 2.*phi[jx][jy+1][jz] - 0.5*phi[jx][jy+2][jz];
	    }else if (jy>mesh->yend){
	      // parallel derivs at right guard point
	      dNi_dpar = 1.5*Ni[jx][jy][jz] - 2.*Ni[jx][jy-1][jz] + 0.5*Ni[jx][jy-2][jz];
	      dPhi_dpar = 1.5*phi[jx][jy][jz] - 2.*phi[jx][jy-1][jz] + 0.5*phi[jx][jy-2][jz];
	    }else {
	      dNi_dpar = 0.5*Ni[jx][jy+1][jz] - 0.5*Ni[jx][jy-1][jz];
	      dPhi_dpar = 0.5*phi[jx][jy+1][jz] - 0.5*phi[jx][jy-1][jz];
	    }
	    
	    BoutReal c0=((Bpxy[jx][jy]/mesh->Bxy[jx][jy])/hthe[jx][jy])/mesh->dy[jx][jy];
	    dNi_dpar = dNi_dpar*c0;
	    dPhi_dpar = dPhi_dpar*c0;
	    
	    jpar[jx][jy][jz] = -(1./fmei)*(1./(0.51*nu[jx][jy][jz]))*Ni0[jx][jy]*dPhi_dpar;
	    if(OhmPe)
	      jpar[jx][jy][jz] += (1./fmei)*(1./(0.51*nu[jx][jy][jz]))*Te0[jx][jy]*dNi_dpar;
	      
	  }
    }
    
    jpar.applyBoundary();
    
    Ve = Vi - jpar/Ni0;
    Ajpar = Ve;
  }else {
    
    Ve = Ajpar + Apar;
    jpar = Ni0*(Vi - Ve);
  }

  ////////////////////////////////////////////////////////
  // DENSITY EQUATION

  msg_stack.push("Density equation");
  ddt(Ni) = 0.0;
  if(evolve_ni) {
    
    if(ni_ni1_phi0)
      ddt(Ni) -= vE_Grad(Ni, phi0);
    
    if(ni_ni0_phi1) 
      ddt(Ni) -= vE_Grad(Ni0, phi);
   
    if(ni_ni1_phi1)
      ddt(Ni) -= vE_Grad(Ni, phi);
    
    if(ni_nit_phit)
      ddt(Ni) -= vE_Grad(Nit, phi0 + phi) - vE_Grad(Ni0, phi0);

    if(ni_vi1_ni0)
      ddt(Ni) -= Vpar_Grad_par(Vi, Ni0);
    
    if(ni_vi0_ni1)
      ddt(Ni) -= Vpar_Grad_par(Vi0, Ni);
    
    if(ni_vi1_ni1)
      ddt(Ni) -= Vpar_Grad_par(Vi, Ni);

    if(ni_vit_nit)
      ddt(Ni) -= Vpar_Grad_par(Vit, Nit) - Vpar_Grad_par(Vi0, Ni0);

    if(ni_jpar1) {
      if(stagger) {
	ddt(Ni) += Div_par_CtoL(jpar);
      }else
	ddt(Ni) += Div_par(jpar);
    }

    if(ni_pe1)
      ddt(Ni) += 2.0*V_dot_Grad(b0xcv, pe);
    
    if(ni_ni0_curv_phi1)
      ddt(Ni) -= 2.0*Ni0*V_dot_Grad(b0xcv, phi);
    
    if(ni_ni1_curv_phi0)
      ddt(Ni) -= 2.0*Ni*V_dot_Grad(b0xcv, phi0);
    
    if(ni_ni1_curv_phi1)
      ddt(Ni) -= 2.0*Ni*V_dot_Grad(b0xcv, phi);

    if(ni_nit_curv_phit)
      ddt(Ni) -= 2.0*Nit*V_dot_Grad(b0xcv, phi+phi0) - 2.0*Ni0*V_dot_Grad(b0xcv, phi0);
    
    if(ni_ni1)
      ddt(Ni) += mu_i * Delp2(Ni);
    
    //ddt(Ni) -= Ni0*Div_par(Vi) + Ni*Div_par(Vi0) + Ni*Div_par(Vi);

    if(lowPass_z > 0)
      ddt(Ni) = lowPass(ddt(Ni), lowPass_z);
  }
  msg_stack.pop();

  ////////////////////////////////////////////////////////
  // ION VELOCITY
  
  msg_stack.push("Ion velocity equation");
  ddt(Vi) = 0.0;
  if(evolve_vi) {
    if(vi_vi0_phi1)
      ddt(Vi) -= vE_Grad(Vi0, phi);

    if(vi_vi1_phi0)
      ddt(Vi) -= vE_Grad(Vi, phi0);

    if(vi_vi1_phi1)
      ddt(Vi) -= vE_Grad(Vi, phi);
    
    if(vi_vit_phit)
      ddt(Vi) -= vE_Grad(Vit, phi+phi0) - vE_Grad(Vi0, phi+phi0);
    
    if(vi_vi1_vi0)
      ddt(Vi) -= Vpar_Grad_par(Vi0, Vi);

    if(vi_vi0_vi1)
      ddt(Vi) -= Vpar_Grad_par(Vi, Vi0);
    
    if(vi_vi1_vi1)
      ddt(Vi) -= Vpar_Grad_par(Vi, Vi);

    if(vi_vit_vit)
      ddt(Vi) -= Vpar_Grad_par(Vit, Vit) - Vpar_Grad_par(Vi0, Vi0);

    if(vi_pei1)
      ddt(Vi) -= Grad_par(pei)/Ni0;

    if(vi_peit)
      ddt(Vi) -= Grad_par(pei)/Nit;
    
    if(vi_vi1)
      ddt(Vi) -= mu_i*Delp2(Vi);

    if(lowPass_z > 0)
      ddt(Vi) = lowPass(ddt(Vi), lowPass_z);
  }
  msg_stack.pop();
  
  ////////////////////////////////////////////////////////
  // ELECTRON TEMPERATURE

  msg_stack.push("Electron temperature equation");
  ddt(Te) = 0.0;
  if(evolve_te) {
    if(te_te1_phi0)
      ddt(Te) -= vE_Grad(Te, phi0);
    if(te_te0_phi1)
      ddt(Te) -= vE_Grad(Te0, phi);
    if(te_te1_phi1)
      ddt(Te) -= vE_Grad(Te, phi);
    
    /*
    ddt(Te) -= vE_Grad(Te0, phi) + vE_Grad(Te, phi0) + vE_Grad(Te, phi);
    ddt(Te) -= Vpar_Grad_par(Ve, Te0) + Vpar_Grad_par(Ve0, Te) + Vpar_Grad_par(Ve, Te);
    ddt(Te) += 1.333*Te0*( V_dot_Grad(b0xcv, pe)/Ni0 - V_dot_Grad(b0xcv, phi) );
    ddt(Te) += 3.333*Te0*V_dot_Grad(b0xcv, Te);
    ddt(Te) += (0.6666667/Ni0)*Div_par_K_Grad_par(kapa_Te, Te);

    */
    if(lowPass_z > 0)
      ddt(Te) = lowPass(ddt(Te), lowPass_z);
  }
  msg_stack.pop();
  
  ////////////////////////////////////////////////////////
  // ION TEMPERATURE

  msg_stack.push("Ion temperature equation");
  ddt(Ti) = 0.0;
  if(evolve_ti) {
    if(ti_ti1_phi0)
      ddt(Ti) -= vE_Grad(Ti, phi0);
    if(ti_ti0_phi1)
      ddt(Ti) -= vE_Grad(Ti0, phi);
    if(ti_ti1_phi1)
      ddt(Ti) -= vE_Grad(Ti, phi);
    
    /*
    ddt(Ti) -= vE_Grad(Ti0, phi) + vE_Grad(Ti, phi0) + vE_Grad(Ti, phi);
    ddt(Ti) -= Vpar_Grad_par(Vi, Ti0) + Vpar_Grad_par(Vi0, Ti) + Vpar_Grad_par(Vi, Ti);
    ddt(Ti) += 1.333*( Ti0*V_dot_Grad(b0xcv, pe)/Ni0 - Ti*V_dot_Grad(b0xcv, phi) );
    ddt(Ti) -= 3.333*Ti0*V_dot_Grad(b0xcv, Ti);
    ddt(Ti) += (0.6666667/Ni0)*Div_par_K_Grad_par(kapa_Ti, Ti);
    */

    if(lowPass_z > 0)
      ddt(Ti) = lowPass(ddt(Ti), lowPass_z);
  }
  msg_stack.pop();
  
  ////////////////////////////////////////////////////////
  // VORTICITY

  msg_stack.push("Vorticity equation");
  ddt(rho) = 0.0;
  if(evolve_rho) {
    
    if(rho_rho0_phi1)
      ddt(rho) -= vE_Grad(rho0, phi);

    if(rho_rho1_phi0)
      ddt(rho) -= vE_Grad(rho, phi0);

    if(rho_rho1_phi1)
      ddt(rho) -= vE_Grad(rho, phi);

    if(rho_vi1_rho0)
      ddt(rho) -= Vpar_Grad_par(Vi, rho0);
    
    if(rho_vi0_rho1)
      ddt(rho) -= Vpar_Grad_par(Vi0, rho);
    
    if(rho_vi1_rho1)
      ddt(rho) -= Vpar_Grad_par(Vi, rho);
    
    if(rho_pei1) {
      if(curv_upwind) {
	ddt(rho) += 2.0*mesh->Bxy*V_dot_Grad(b0xcv, pei);  // Use upwinding
      }else
	ddt(rho) += 2.0*mesh->Bxy*b0xcv*Grad(pei);     // Use central differencing
    }    

    if(rho_jpar1) {
      if(stagger) {
	ddt(rho) += mesh->Bxy*mesh->Bxy*Div_par_CtoL(jpar);
      }else 
	ddt(rho) += mesh->Bxy*mesh->Bxy*Div_par(jpar, CELL_CENTRE);
    }

    if(rho_rho1)
      ddt(rho) += mu_i * Delp2(rho);

    if(lowPass_z > 0)
      ddt(rho) = lowPass(ddt(rho), lowPass_z);
  }
  msg_stack.pop();
  
  ////////////////////////////////////////////////////////
  // AJPAR
  
  msg_stack.push("Ajpar equation");
  ddt(Ajpar) = 0.0;
  if(evolve_ajpar) {
    //ddt(Ajpar) -= vE_Grad(Ajpar0, phi) + vE_Grad(Ajpar, phi0) + vE_Grad(Ajpar, phi);
    //ddt(Ajpar) -= (1./fmei)*1.71*Grad_par(Te);
    
    if(stagger) {
      ddt(Ajpar) += (1./fmei)*Grad_par_LtoC(phi); // Right-hand differencing
    }else
      ddt(Ajpar) += (1./fmei)*Grad_par(phi, CELL_YLOW);
    
    if(OhmPe) {
      if(stagger) {
	ddt(Ajpar) -= (1./fmei)*(Tet/Nit)*Grad_par_LtoC(Ni);
      }else 
	ddt(Ajpar) -= (1./fmei)*(Te0/Ni0)*Grad_par(Ni, CELL_YLOW);
    }
    
    ddt(Ajpar) += 0.51*interp_to(nu, CELL_YLOW)*jpar/Ni0;

    if(lowPass_z > 0)
      ddt(Ajpar) = lowPass(ddt(Ajpar), lowPass_z);
  }
  msg_stack.pop();
  
  ////////////////////////////////////////////////////////
  // Profile evolution options

  switch(iTe_dc) {
  case 1: { // subtacting out toroidal averages for all fields
    if(evolve_ni)
      ddt(Ni) -= ddt(Ni).DC();
    if(evolve_rho)
      ddt(rho) -= ddt(rho).DC();
    if(evolve_te)
      ddt(Te) -= ddt(Te).DC();
    if(evolve_ti)
      ddt(Ti) -= ddt(Ti).DC();
    if(evolve_ajpar)
      ddt(Ajpar) -= ddt(Ajpar).DC();
    break;
  }
  case 2: { // not subtacting out toroidal averages for any field
    break;
  }
  case 4: { // using toroidal averages in right-hand sides, e.g., axisymmetric mode
    if(evolve_ni)
      ddt(Ni) = ddt(Ni).DC();
    if(evolve_rho)
      ddt(rho) = ddt(rho).DC();
    if(evolve_te)
      ddt(Te) = ddt(Te).DC();
    if(evolve_ti)
      ddt(Ti) = ddt(Ti).DC();
    if(evolve_ajpar)
      ddt(Ajpar) = ddt(Ajpar).DC();
    break;
  }
  default: {
    throw BoutException("ERROR: invalid option for iTe_dc\n");
  }
  }
  
  return(0);
}

