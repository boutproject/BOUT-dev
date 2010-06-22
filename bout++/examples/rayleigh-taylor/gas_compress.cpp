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
bool sub_initial;

// 2D initial profiles
Field2D N0, P0;
Vector2D V0;

Vector2D g; // Acceleration

int physics_init(bool restarting)
{
  real v0_multiply;

  // Read initial conditions

  mesh->get(N0, "density");
  mesh->get(P0, "pressure");
  V0.covariant = false; // Read contravariant components
  V.covariant = false; // Evolve contravariant components
  mesh->get(V0, "v");
  g.covariant = false;
  mesh->get(g, "g");
  
  // read options
  
  options.setSection("gas");
  
  options.get("gamma",  gamma_ratio, 5.0/3.0);
  options.get("viscosity", nu, 0.1);
  options.get("include_viscosity", include_viscosity, false);
  options.get("v0_multiply", v0_multiply, 1.0);
  options.get("sub_initial", sub_initial, false);

  V0 *= v0_multiply;

  // Set evolving variables
  
  bout_solve(N, F_N, "density");
  bout_solve(P, F_P, "pressure");
  bout_solve(V, F_V, "v");

  if(!restarting) {
    // Apply boundary conditions
    apply_boundary(N, "density");
    apply_boundary(P, "pressure");
    V.to_contravariant();
    apply_boundary(V, "v");

    // Set variables to these values (+ the initial perturbation)
    // NOTE: This must be after the calls to bout_solve
    N += N0;
    P += P0;
    V += V0;

  }
  
  return 0;
}

int physics_run(real t)
{
  //output.write("Running %e\n", t);
  // Communicate variables
  mesh->communicate(N,P,V);

  // Density
  
  F_N = -V_dot_Grad(V, N) - N*Div(V);
 
  //output.write("N ");
 
  // Velocity 
  
  F_V = -V_dot_Grad(V, V) - Grad(P)/N + g;

  if(sub_initial) {
    F_V += Grad(P0)/N0 - g;
  }

  //output.write("V ");

  if(include_viscosity) {
    // Add viscosity
    
    F_V.y += nu*Laplacian(V.y);
    F_V.z += nu*Laplacian(V.z);

    //output.write("nu ");
  }
  
  // Pressure

  F_P = -V_dot_Grad(V, P) - gamma_ratio*P*Div(V);

  //output.write("P\n");

  // Set boundary conditions
  apply_boundary(F_N, "density");
  apply_boundary(F_P, "pressure");
  F_V.to_contravariant();
  apply_boundary(F_V, "v");

  //output.write("finished\n");

  return 0;
}


