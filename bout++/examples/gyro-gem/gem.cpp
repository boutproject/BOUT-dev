/****************************************************************
 * GEM Gyro-fluid model
 * 
 * 6 moments for each species
 *
 * This version uses global parameters for collisionality etc.
 ****************************************************************/

#include "bout.h"
#include "gyro_average.h"
#include "invert_laplace.h"

#include <cmath>

//////////////////////////////////////////
// Evolving quantities

// Ion species
Field3D Ni;   // Gyro-center density
Field3D ApUi; // beta_e*Apar + mu_i * Ui
Field3D Tipar, Tiperp; // Parallel or perpendicular temp
Field3D qipar, qiperp; // Parallel and perpendicular heat flux

// Electron species
Field3D Ne;
Field3D ApUe; // beta_e*Apar + mu_e * Ue
Field3D Tepar, Teperp;
Field3D qepar, qeperp;

//////////////////////////////////////////
// Derived quantities

Field3D phi;    // Electrostatic potential
Field3D Apar;   // Parallel vector potential

Field3D Ui, Ue; // Ion and electron parallel velocity

Field3D Jpar;   // Parallel current

//////////////////////////////////////////
// Equilibrium

Vector3D B0vec; // Equilibrium B field vector
Field2D Grad_par_logB; // Grad_par(log(B))

Field2D Ni0, Ne0; // Gyro-center densities
Field2D Ti0, Te0; // Starting isotropic temperatures

BoutReal tau_e, tau_i; // T_z / (Te * Z)
BoutReal mu_e, mu_i;   // M_z / (M_i * Z)

BoutReal beta_e; // Electron dynamical beta

BoutReal rho_e, rho_i; // Electron, ion gyroradius

// Collisional transport coefficients
const BoutReal eta     = 0.51;

const BoutReal alpha_e = 0.71;
const BoutReal kappa_e = 3.2;
const BoutReal pi_e    = 0.73;

// alpha_i = 0
const BoutReal kappa_i = 3.9;
const BoutReal pi_i    = 0.73;

Field3D Rei;

//////////////////////////////////////////
// Options

bool adiabatic_electrons;  // Solve adiabatic electrons
bool small_rho_e;          // Neglect electron gyro-radius
bool include_grad_par_B;   // Include terms like Grad_par(log(B))

BoutReal Landau; // Multiplier for Landau damping terms

BoutReal nu_e, nu_i; // Collisional dissipation

BoutReal nu_perp, nu_par;  // Artificial 

const BRACKET_METHOD bm = BRACKET_STD; // Method to use for brackets

int phi_flags, apar_flags; // Inversion flags

//////////////////////////////////////////
// Normalisation factors

BoutReal Lbar;   // Perpendicular scale length
BoutReal Tenorm; // Typical value of Te for normalisation
BoutReal Bbar;   // Magnetic field  

BoutReal Cs;    // Sound speed sqrt(Tenorm / Mi)
BoutReal Tbar;  // Timescale Lbar / Cs

FieldGroup comms; // Communications

////////////////////////////////////////////////////////////////////////
// Initialisation

int physics_init(bool restarting)
{
  //////////////////////////////////
  // Read options

  options.setSection("gem");
  
  OPTION(adiabatic_electrons, false);
  OPTION(small_rho_e, true);
  OPTION(include_grad_par_B, true);
  
  OPTION(Landau, 1.0);
  
  OPTION(nu_perp, 0.01); // Artificial perpendicular dissipation
  OPTION(nu_par, 3e-3);  // Artificial parallel dissipation
  
  OPTION(phi_flags,  0);
  OPTION(apar_flags, 0);

  //////////////////////////////////
  // Read profiles

  // Mesh
  Field2D Rxy, Bpxy, Btxy, Bxy, hthe;
  GRID_LOAD(Rxy);    // Major radius [m]
  GRID_LOAD(Bpxy);   // Poloidal B field [T]
  GRID_LOAD(Btxy);   // Toroidal B field [T]
  GRID_LOAD(Bxy);    // Total B field [T]
  GRID_LOAD(hthe);   // Poloidal arc length [m / radian]
  
  GRID_LOAD(Te0); // Electron temperature in eV
  GRID_LOAD(Ne0); // Electron number density in m^-3

  Field2D p_e = 1.602e-19 * Te0 * Ne0; // Electron pressure in Pascals
 
  //////////////////////////////////
  // Pick normalisation factors
  
  if(mesh->get(Lbar, "Lbar")) // Try to read from grid file
    Lbar = 1.0;
  OPTION(Lbar, Lbar); // Override in options file
  SAVE_ONCE(Lbar);   // Save in output file

  BoutReal AA; // Ion atomic mass
  BoutReal ZZ; // Ion charge
  OPTION(AA, 2.0); // Deuterium by default
  OPTION(ZZ, 1.0); 

  Tenorm = max(Te0,true); SAVE_ONCE(Tenorm); // Maximum value over the grid
  Cs = sqrt(1.602e-19*Tenorm / (AA*1.67262158e-27)); SAVE_ONCE(Cs); // Sound speed in m/s
  Tbar = Lbar / Cs; SAVE_ONCE(Tbar); // Timescale in seconds
  Bbar = max(Bxy, true); SAVE_ONCE(Bbar);
  
  beta_e =  4.e-7*PI * max(p_e,true) / (Bbar*Bbar); SAVE_ONCE(beta_e); 
  
  // Mass to charge ratios
  mu_i = 1. / ZZ;
  mu_e = -1. / (AA * 1860.);
  
  BoutReal t_e, t_i; // Braginskii collision times
  
  nu_e = Lbar / (Cs*t_e); SAVE_ONCE(nu_e);
  nu_i = Lbar / (Cs*t_i); SAVE_ONCE(nu_i);
  
  //////////////////////////////////
  // Metric tensor components
  
  // Normalise
  hthe /= Lbar;
  Bpxy /= Bbar;
  Btxy /= Bbar;
  Bxy  /= Bbar;
  hthe /= Lbar;
  mesh->dx   /= Lbar*Lbar*Bbar;
  
  // Metric components
  
  mesh->g11 = (Rxy*Bpxy)^2;
  mesh->g22 = 1.0 / (hthe^2);
  mesh->g33 = (Bxy^2)/mesh->g11;
  mesh->g12 = 0.0;
  mesh->g13 = 0.;
  mesh->g23 = -Btxy/(hthe*Bpxy*Rxy);
  
  mesh->J = hthe / Bpxy;
  mesh->Bxy = Bxy;
  
  mesh->g_11 = 1.0/mesh->g11;
  mesh->g_22 = (Bxy*hthe/Bpxy)^2;
  mesh->g_33 = Rxy*Rxy;
  mesh->g_12 = 0.;
  mesh->g_13 = 0.;
  mesh->g_23 = Btxy*hthe*Rxy/Bpxy;
  
  mesh->geometry();
  
  // Set B field vector
  
  B0vec.covariant = false;
  B0vec.x = 0.;
  B0vec.y = Bpxy / hthe;
  B0vec.z = 0.;

  // Precompute this for use in RHS
  if(include_grad_par_B) {
    Grad_par_logB = Grad_par(log(mesh->Bxy));
  }else
    Grad_par_logB = 0.;
  
  //////////////////////////////////
  
  // Add ion equations
  SOLVE_FOR6(Ni, ApUi, Tipar, Tiperp, qipar, qiperp);
  comms.add(Ni, ApUi, Tipar, Tiperp, qipar, qiperp);
  
  if(adiabatic_electrons) {
    // Solving with adiabatic electrons
    
  }else {
    // Add electron equations
    
    SOLVE_FOR6(Ne, ApUe, Tepar, Teperp, qepar, qeperp);
    comms.add(Ne, ApUe, Tepar, Teperp, qepar, qeperp);
  }
  
  //////////////////////////////////
  
  if(!restarting) {
    // Initial current
    
    Field2D Jpar0;
    if(mesh->get(Jpar0, "Jpar0") == 0) {
      // Initial current specified. Set parallel electron velocity
      
    }
    
    // Initial potential
    
    Field2D phi0;
    if(mesh->get(phi0, "phi0") == 0) {
      
    }
    
  }
 
  return 0;
}

////////////////////////////////////////////////////////////////////////
// Prototypes

const Field3D curvature(const Field3D &f);

const Field3D UE_Grad(const Field3D &f, const Field3D &phi);
const Field3D WE_Grad(const Field3D &f, const Field3D &Phi);

const Field3D Grad_parP(const Field3D &f);
const Field3D Div_parP(const Field3D &f);

////////////////////////////////////////////////////////////////////////
// RHS function

int physics_run(BoutReal time)
{
  // Quantities which depend on species
  Field3D phi_G, Phi_G; // Gyro-reduced potential
  Field3D S_D, K_par, K_perp, K_D; // Collisional dissipation terms
  
  
  ////////////////////////////////////////////
  // Adiabatic electrons
  
  if(adiabatic_electrons) {
    // Solve adiabatic electrons using surface-averaged phi
    
    Field2D phi_zonal = mesh->averageY(phi.DC()); // Average over Y and Z
    Ne = phi - phi_zonal;
    
    // Need to solve with polarisation!
  }

  ////////////////////////////////////////////
  // Polarisation equation (quasi-neutrality)
  
  if(small_rho_e) {
    // Neglect electron Larmor radius
    
    Field3D dn = Ne - gyroPade1(Ni, rho_i) - gyroPade2(Tiperp, rho_i);
    phi = invert_laplace(tau_i * dn / (rho_i * rho_i), phi_flags);
    phi -= tau_i * dn;
  }else {
    Field3D dn = gyroPade1(Ne, rho_e) + gyroPade2(Teperp, rho_e)
      - gyroPade1(Ni, rho_i) - gyroPade2(Tiperp, rho_i);
    
    // Neglect electron gyroscreening
    phi = invert_laplace(tau_i * dn / (rho_i * rho_i), phi_flags);
    phi -= tau_i * dn;
  }
  
  ////////////////////////////////////////////
  // Helmholtz equation for Apar
  
  Field2D a = beta_e * (1./mu_e - 1./mu_i);
  Apar = invert_laplace(ApUe/mu_e - ApUi/mu_i, apar_flags, &a);
  
  Ui = (ApUi - beta_e*Apar) / mu_i;
  Ue = (ApUe - beta_e*Apar) / mu_e;
  
  Jpar = Ui - Ue;

  ////////////////////////////////////////////
  // Resistivity
  
  Rei = mu_e*nu_e*(eta*Jpar + 
                   (alpha_e/kappa_e)*(qepar + qeperp + alpha_e*Jpar));
  
  ////////////////////////////////////////////
  // Electron equations
  
  if(!adiabatic_electrons) {
    // Electron equations
    
    if(small_rho_e) {
      // No gyro-averaging for small rho_e
      phi_G = phi;
      Phi_G = 0.0;
    }else {
      phi_G = gyroPade1(phi, rho_e);
      Phi_G = gyroPade2(phi, rho_e);
    }
    
    // Collisional dissipation
    S_D = (nu_e / (3.*pi_e)) * (Tepar - Teperp);
    K_par = mu_e*tau_e*nu_e*((5./2.)/kappa_e)*(qepar + 0.6*alpha_e*Jpar);
    K_perp = mu_e*tau_e*nu_e*((5./2.)/kappa_e)*(qeperp + 0.4*alpha_e*Jpar);
    K_D = 1.28*mu_e*tau_e*nu_e*((5./2.)/kappa_e)*(qepar - 1.5*qeperp);
    
    ddt(Ne) = -UE_Grad(Ne0 + Ne, phi_G) 
      - WE_Grad(Te0 + Teperp, Phi_G)
      - Div_parP(Ue)
      + curvature(phi_G + tau_e*Ne + 0.5*(tau_e*Tepar + tau_e*Teperp + Phi_G));
    
    ddt(ApUe) = -mu_e*UE_Grad(Ue, phi_G)
      - mu_e*WE_Grad(qeperp, Phi_G)
      - Grad_parP(phi_G + tau_e*(Ne0 + Te0 + Ne + Tepar))
      + mu_e * tau_e * curvature(2.*Ue + qepar + 0.5*qeperp)
      - tau_e * (Phi_G + tau_e*Teperp - tau_e*Tepar)*Grad_par_logB
      + Rei;
    
    ddt(Tepar) = - UE_Grad(Te0 + Tepar, phi_G)
      - 2.*Div_parP(Ue + qepar)
      + curvature(phi_G + tau_e*(Ne+Tepar) + 2.*tau_e*Tepar)
      - (Ue + qeperp)*Grad_par_logB
      - 2.*S_D;
    
    ddt(Teperp) = - UE_Grad(Te0 + Teperp, phi_G)
      - WE_Grad(Ne0 + Ne + 2.*(Te0 + Teperp), Phi_G)
      - Div_parP(qeperp)
      + 0.5*curvature(phi_G + Phi_G + tau_e*(Ne + Teperp) 
                      + 3.*(Phi_G + tau_e*Teperp))
      + (Ue + qeperp)*Grad_par_logB
      + S_D;
    
    ddt(qepar) = - UE_Grad(qepar, phi_G)
      - 1.5*(1./mu_e)*Grad_parP(tau_e*(Te0 + Tepar))
      + 0.5*mu_e*tau_e*curvature(3.*Ue + 8.*qepar)
      - Landau*(tau_e/mu_e)*(1. - 0.125*Grad2_par2(qepar))
      - (1./mu_e)*K_par
      - (1./mu_e)*K_D;
    
    ddt(qeperp) = - UE_Grad(qeperp, phi_G)
      - WE_Grad(Ue + 2.*qeperp, Phi_G)
      - (1./mu_e)*Grad_parP(Phi_G + tau_e*(Te0 + Teperp))
      + 0.5*tau_e*curvature(Ue + 6.*qeperp)
      - (tau_e/mu_e)*(Phi_G + tau_e*Teperp - tau_e*Tepar)*Grad_par_logB
      - (1./mu_e)*K_perp
      + (1./mu_e)*K_D;
    
  }
  
  ////////////////////////////////////////////
  // Ion equations
  
  // Calculate gyroreduced potentials
  phi_G = gyroPade1(phi, rho_i);
  Phi_G = gyroPade2(phi, rho_i);
  
  // Collisional dissipation
  S_D = (nu_i / (3.*pi_i)) * (Tipar - Tiperp);
  K_par = mu_i*tau_i*nu_i*((5./2.)/kappa_i)*qipar;
  K_perp = mu_i*tau_i*nu_i*((5./2.)/kappa_i)*qiperp;
  K_D = 1.28*mu_i*tau_i*nu_i*((5./2.)/kappa_i)*(qipar - 1.5*qiperp);

  ddt(Ni) = -UE_Grad(Ni0 + Ne, phi_G) 
    - WE_Grad(Ti0 + Tiperp, Phi_G)
    - Div_parP(Ui)
    + curvature(phi_G + tau_i*Ni + 0.5*(tau_i*Tipar + tau_e*Tiperp + Phi_G));
  
  ddt(ApUi) = -mu_i*UE_Grad(Ui, phi_G)
    - mu_i*WE_Grad(qiperp, Phi_G)
    - Grad_parP(phi_G + tau_i*(Ni0 + Ti0 + Ni + Tipar))
    + mu_i * tau_i * curvature(2.*Ui + qipar + 0.5*qiperp)
    - tau_i * (Phi_G + tau_i*Tiperp - tau_i*Tipar)*Grad_par_logB
    + Rei;
  
  ddt(Tipar) = - UE_Grad(Ti0 + Tipar, phi_G)
    - 2.*Div_parP(Ui + qipar)
    + curvature(phi_G + tau_i*(Ni+Tipar) + 2.*tau_i*Tipar)
    - (Ui + qiperp)*Grad_par_logB
    - 2.*S_D;
  
  ddt(Tiperp) = - UE_Grad(Ti0 + Tiperp, phi_G)
    - WE_Grad(Ni0 + Ni + 2.*(Ti0 + Tiperp), Phi_G)
    - Div_parP(qiperp)
    + 0.5*curvature(phi_G + Phi_G + tau_i*(Ni + Tiperp) 
                    + 3.*(Phi_G + tau_i*Tiperp))
    + (Ui + qiperp)*Grad_par_logB
    + S_D;
  
  ddt(qipar) = - UE_Grad(qipar, phi_G)
    - 1.5*(1./mu_i)*Grad_parP(tau_i*(Ti0 + Tipar))
    + 0.5*tau_i*curvature(3.*Ui + 8.*qipar)
    - (1./mu_e)*K_par
    - (1./mu_e)*K_D;
  
  ddt(qiperp) = - UE_Grad(qiperp, phi_G)
    - WE_Grad(Ui + 2.*qiperp, Phi_G)
    - (1./mu_i)*Grad_parP(Phi_G + tau_i*(Ti0 + Tiperp))
    + 0.5*tau_i*curvature(Ui + 6.*qiperp)
    - (tau_i/mu_i)*(Phi_G + tau_i*Tiperp - tau_i*Tipar)*Grad_par_logB
    - (1./mu_e)*K_perp
    + (1./mu_e)*K_D;

  return 0;
}

////////////////////////////////////////////////////////////////////////
// Curvature operator

// K(f) = Div((c/B^2) B x Grad(f))
// Simple implementation. Could be improved to eliminate the communication
const Field3D curvature(const Field3D &f)
{
  /*
  Vector3D gradf = Grad(f);
  // Set boundaries to zero-gradient
  gradf.x.applyBoundary("dirichlet");
  gradf.y.applyBoundary("dirichlet");
  gradf.z.applyBoundary("dirichlet");
  // Communicate
  mesh->communicate(gradf);
  
  return Div(B0vec ^ gradf / (mesh->Bxy*mesh->Bxy));
  */
  return -bracket(log(mesh->Bxy^2), f, bm);
}

////////////////////////////////////////////////////////////////////////
// Advection terms

const Field3D UE_Grad(const Field3D &f, const Field3D &phi)
{
  Field3D delp2 = Delp2(f);
  delp2.applyBoundary("dirichlet");
  mesh->communicate(delp2);
  
  return bracket(phi, f, bm)
    + nu_perp*Delp2( delp2 * ( (1./mesh->Bxy)^4 ) )
    - nu_par*Grad2_par2(f); // NB: This should be changed for variable B
}

const Field3D WE_Grad(const Field3D &f, const Field3D &Phi)
{
  return bracket(Phi, f, bm);
}

////////////////////////////////////////////////////////////////////////
// Parallel derivative

const Field3D Grad_parP(const Field3D &f)
{
  return Grad_par(f) - beta_e*bracket(Apar, f, bm);
}

const Field3D Div_parP(const Field3D &f)
{
  return mesh->Bxy*Grad_parP(f/mesh->Bxy);
}

