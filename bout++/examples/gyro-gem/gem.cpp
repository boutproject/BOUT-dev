/****************************************************************
 * GEM Gyro-fluid model
 * 
 * 6 moments for each species
 ****************************************************************/

#include "bout.h"

// Ion species
Field3D Ni;   // Gyro-center density
Field3D ApVi; // Combination of Apar and Vi
Field3D Ti_Tipar, Ti_Tiperp; // Temp + Parallel or perpendicular temp
Field3D qipar, qiperp; // Parallel and perpendicular heat flux

// Electron species
Field3D Ne;
Field3D ApVe;
Field3D Te_Tepar, Te_Teperp;
Field3D qepar, qeperp;

Field3D phi;    // Electrostatic potential

Vector3D B0vec; // Equilibrium B field vector

//////////////////////////////////////////
// Options

bool adiabatic_electrons;  // Solve adiabatic electrons
bool small_rho_e;          // Neglect electron gyro-radius

//////////////////////////////////////////
// Normalisation factors

BoutReal Lperp;  // Perpendicular scale length
BoutReal Tenorm; // Typical value of Te for normalisation

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

  //////////////////////////////////
  // Read profiles
  
  Field2D Te0;
  GRID_LOAD(Te0);

  Field2D pressure;
 
  //////////////////////////////////
  // Pick normalisation factors
  
  if(mesh->get(Lperp, "Lperp")) { // Try to read from input
    // Get the maximum gradient in the equilibrium pressure
    Lperp = 1. / max(abs(Grad(pressure)),true);
  }

  Tenorm = max(Te0,true); // Maximum value over the grid
  
  // Add factors to output file
  SAVE_ONCE(Lperp);
  SAVE_ONCE(Tenorm);
  
  //////////////////////////////////
  
  // Add ion equations
  SOLVE_FOR(Ni);
  SOLVE_FOR(ApVi);
  SOLVE_FOR(Ti_Tipar);
  SOLVE_FOR(Ti_Tiperp);
  SOLVE_FOR(qipar);
  SOLVE_FOR(qiperp);
  comms.add(Ni, ApVi, Ti_Tipar, Ti_Tiperp, qipar, qiperp);
  
  if(!adiabatic_electrons) {
    // Add electron equations
    
    SOLVE_FOR(Ne);
    SOLVE_FOR(ApVe);
    SOLVE_FOR(Te_Tepar);
    SOLVE_FOR(Te_Teperp);
    SOLVE_FOR(qepar);
    SOLVE_FOR(qeperp);
  }else {
    // Solving with adiabatic electrons
    
    
  }
}

////////////////////////////////////////////////////////////////////////
// RHS function

int physics_run(BoutReal time)
{
  ////////////////////////////////////////////
  // Polarisation equation (quasi-neutrality)
  


  if(adiabatic_electrons) {
    // Solve adiabatic electrons using surface-averaged phi
    
    Field2D phi_zonal = mesh->averageY(phi.DC());
    Ne = phi - phi_zonal;
  }else {
    // Electron equations
  }
  
  ////////////////////////////////////////////
  // Ion equations
  
  // Calculate gyroreduced potentials
  Field3D phi_G = gyroPade1(phi);
  Field3D Phi_G = gyroPade2(phi);
  
}

////////////////////////////////////////////////////////////////////////
// Curvature operator

// K(f) = Div((c/B^2) B x Grad(f))
// Simple implementation. Could be improved to eliminate the communication
const Field3D curvature(const Field3D &f)
{
  Vector3D gradf = Grad(f);
  // Set boundaries to zero-gradient
  gradf.x.applyBoundary("dirichlet");
  gradf.y.applyBoundary("dirichlet");
  gradf.z.applyBoundary("dirichlet");
  // Communicate
  mesh->communicate(gradf);
  
  return Div(B0vec ^ gradf / (mesh->Bxy*mesh->Bxy));
}

////////////////////////////////////////////////////////////////////////
// Density equation

const Field3D density_rhs()
{
  
}
