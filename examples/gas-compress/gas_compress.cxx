
#include "gas_compress.hxx"

int GasCompress::init(bool restarting) {

  // Read initial conditions

  mesh->get(N0, "density");
  mesh->get(P0, "pressure");
  V0.covariant = false; // Read contravariant components
  V.covariant = false; // Evolve contravariant components
  mesh->get(V0, "v");
  g.covariant = false;
  mesh->get(g, "g");
  
  // read options
  
  auto& options = Options::root()["gas"];
  gamma_ratio = options["gamma"].withDefault(0.1);
  include_viscosity = options["include_viscosity"].withDefault(false);
  
  BoutReal v0_multiply = options["v0_multiply"].withDefault(1.0);
  V0 *= v0_multiply;
  
  sub_initial = options["sub_initial"].withDefault(false);
  
  V.y.setLocation(CELL_YLOW); // Stagger
  
  // Set evolving variables
  
  bout_solve(N, "density");
  bout_solve(P, "pressure");
  bout_solve(V, "v");
  
  if(!restarting) {
    // Set variables to these values (+ the initial perturbation)
    // NOTE: This must be after the calls to bout_solve
    N += N0;
    P += P0;
    V += V0;
  }
  
  return 0;
}

int GasCompress::rhs(BoutReal UNUSED(time)) {
  // Run communications
  mesh->communicate(N,P,V);
  
  // Density
  ddt(N) = -Div(V, N);
  
  // Velocity

  if (sub_initial) {
    // Subtract force balance of initial profiles
    ddt(V) = -V_dot_Grad(V, V) - Grad(P - P0, CELL_VSHIFT) / N;
  } else {
    ddt(V) = -V_dot_Grad(V, V) - Grad(P, CELL_VSHIFT) / N + g;
  }

  if(include_viscosity) {
    // Add viscosity
    
    ddt(V).y += nu*Laplace(V.y);
    ddt(V).z += nu*Laplace(V.z);
  }
  
  // Pressure
  ddt(P) = -Div(V, P) - (gamma_ratio-1.)*P*Div(V);
  
  return 0;
}

// Create a simple main() function
BOUTMAIN(GasCompress);
