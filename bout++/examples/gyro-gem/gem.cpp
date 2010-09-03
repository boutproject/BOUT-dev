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

//////////////////////////////////////////
// Options

bool adiabatic_electrons;  // Solve adiabatic electrons
bool small_rho_e;          // Neglect electron gyro-radius

FieldGroup comms; // Communications

int physics_init(bool restarting)
{
  
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
  }
}

int physics_run(BoutReal time)
{

}
