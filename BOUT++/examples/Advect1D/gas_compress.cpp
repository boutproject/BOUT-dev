/*******************************************************************
 * Compressible gas-dynamics
 *
 * B.Dudson, December 2007
 *******************************************************************/

#include "bout.h"

// Evolving variables 
Field3D N, P; // Density, Pressure
Vector3D V;   // velocity

// Time-derivatives
Field3D F_N, F_P;
Vector3D F_V;

// parameters
real gamma_ratio;   // Ratio of specific heats
real nu;      // Viscosity
bool include_viscosity;

Vector2D g; // Acceleration

// Parallel communication object
Communicator comms;

int physics_init()
{
  // 2D initial profiles
  Field2D N0, P0;
  Vector2D V0;
  real v0_multiply;

  // Read initial conditions

  grid_load2d(N0, "density");
  grid_load2d(P0, "pressure");
  V0.covariant = false; // Read contravariant components
  V.covariant = false; // Evolve contravariant components
  grid_load2d(V0, "v");
  g.covariant = false;
  grid_load2d(g, "g");
  
  // read options

  if(options.getReal("gas", "gamma", gamma_ratio))
    gamma_ratio = 5.0 / 3.0;
  if(options.getReal("gas", "viscosity", nu))
    nu = 0.1;
  if(options.getBool("gas", "include_viscosity", include_viscosity))
    include_viscosity = false;
  if(options.getReal("gas", "v0_multiply", v0_multiply))
    v0_multiply = 1.0;

  V0 *= v0_multiply;

  // Set evolving variables
  
  bout_solve(N, F_N, "density");
  bout_solve(P, F_P, "pressure");
  bout_solve(V, F_V, "v");

  if(!restarting) {
    // Set variables to these values (+ the initial perturbation)
    // NOTE: This must be after the calls to bout_solve
    N += N0;
    P += P0;
    V += V0;
  }

  // set communications
  comms.add(N);
  comms.add(P);
  comms.add(V);
  
  return 0;
}

int physics_run(real t)
{
  // Run communications
  comms.run();

  // Density
  
  F_N = -V_dot_Grad(V, N) - N*Div(V);
  
  // Velocity 
  
  F_V = -V_dot_Grad(V, V) - Grad(P)/N + g;

  if(include_viscosity) {
    // Add viscosity
    
    F_V.y += nu*Laplacian(V.y);
    F_V.z += nu*Laplacian(V.z);
  }
  
  // Pressure

  F_P = -V_dot_Grad(V, P) - gamma_ratio*P*Div(V);

  // Set boundary conditions
  apply_boundary(F_N, "density");
  apply_boundary(F_P, "pressure");
  F_V.to_contravariant();
  apply_boundary(F_V, "v");

  return 0;
}


