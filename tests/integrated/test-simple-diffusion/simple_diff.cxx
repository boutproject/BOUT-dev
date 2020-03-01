/*******************************************************************
 * Compressible gas-dynamics
 *
 * B.Dudson, December 2007
 *******************************************************************/

#include <bout.hxx>
#include <boutmain.hxx>

// Evolving variables 
Field3D N, P; // Density, Pressure
Vector3D V;   // velocity

// parameters
BoutReal gamma_ratio;   // Ratio of specific heats
BoutReal nu;      // Viscosity
bool include_viscosity;

Vector2D g; // Acceleration

int physics_init(bool restarting)
{
  // 2D initial profiles
  Field2D N0, P0;
  Vector2D V0;
  BoutReal v0_multiply;

  // Read initial conditions

  mesh->get(N0, "density");
  //mesh->get(P0, "pressure");
  V0.covariant = false; // Read contravariant components
  V.covariant = false; // Evolve contravariant components
  mesh->get(V0, "v");
  g.covariant = false;
  mesh->get(g, "g");
  
  // read options
  
  Options *options = Options::getRoot();
  options = options->getSection("gas");
  options->get("gamma", gamma_ratio, 5./3.);
  options->get("viscosity", nu, 0.1);
  options->get("include_viscosity", include_viscosity, false);
  options->get("v0_multiply", v0_multiply, 1.0);
  
  V0 *= v0_multiply;

  // Set evolving variables
  
  bout_solve(N, "density");
  //bout_solve(P, "pressure");
  //bout_solve(V, "v");
  
  if(!restarting) {
    // Set variables to these values (+ the initial perturbation)
    // NOTE: This must be after the calls to bout_solve
    N += N0;
    //P += P0;
    //V += V0;
  }
  
  return 0;
}

int physics_run(BoutReal UNUSED(t))
{
  // Run communications
  //mesh->communicate(N,P,V);
  mesh->communicate(N);
  // Density
  
  ddt(N) = Laplace(N);
  
  // Velocity 
  
  //ddt(V) = 0;//-V_dot_Grad(V, V) - Grad(P)/N + g;

  if(include_viscosity) {
    // Add viscosity
    
    //ddt(V).y += nu*Laplace(V.y);
    //ddt(V).z += nu*Laplace(V.z);
  }
  
  // Pressure

  //ddt(P) = 0; //-V_dot_Grad(V, P) - gamma_ratio*P*Div(V);
  
  return 0;
}


