/*******************************************************************
 * Compressible gas-dynamics
 *
 * B.Dudson, December 2007
 *******************************************************************/

#include "bout.h"

// Evolving variables 
Field3D N, P; // Density, Pressure
Vector3D V;   // velocity

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
  
  bout_solve(N, "density");
  bout_solve(P, "pressure");
  bout_solve(V, "v");

  if(!restarting) {
    // Apply boundary conditions
    apply_boundary(N, "density");
    apply_boundary(P, "pressure");
    V.toContravariant();
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
  // Communicate variables
  mesh->communicate(N,P,V);

  // Density
  
  ddt(N) = -V_dot_Grad(V, N) - N*Div(V);
 
  // Velocity 
  
  ddt(V) = -V_dot_Grad(V, V) - Grad(P)/N + g;

  if(sub_initial) {
    ddt(V) += Grad(P0)/N0 - g;
  }

  if(include_viscosity) {
    // Add viscosity
    
    ddt(V).y += nu*Laplacian(V.y);
    ddt(V).z += nu*Laplacian(V.z);
  }
  
  // Pressure

  ddt(P) = -V_dot_Grad(V, P) - gamma_ratio*P*Div(V);

  // Set boundary conditions
  apply_boundary(ddt(N), "density");
  apply_boundary(ddt(P), "pressure");
  ddt(V).toContravariant();
  apply_boundary(ddt(V), "v");

  return 0;
}


