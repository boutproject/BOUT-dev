/**************************************************************************
 * Ideal MHD physics module for BOUT++
 * This version evolves the entire quantity (initial + perturbed)
 **************************************************************************/

#include <bout.hxx>
#include <boutmain.hxx>
#include <msg_stack.hxx>

// 3D evolving variables
Field3D rho, p; // density, pressure
Vector3D v, B;  // velocity, magnetic field

Field3D divB; // Divergence of B (for monitoring)

// parameters
BoutReal gamma;
bool include_viscos;
BoutReal viscos;

int physics_init(bool restarting) {
  // 2D initial profiles
  Field2D rho0, p0;
  Vector2D v0, B0;

  // read options
  Options *globalOptions = Options::getRoot();
  Options *options = globalOptions->getSection("mhd");
  OPTION(options, gamma,          5.0/3.0);
  OPTION(options, include_viscos, false);
  OPTION(options, viscos,         0.1);
  
  // Read 2D initial profiles
  GRID_LOAD(rho0);
  GRID_LOAD(p0);
  v0.covariant = true; // Read covariant components of v0
  GRID_LOAD(v0);
  B0.covariant = false; // Read contravariant components of B0
  GRID_LOAD(B0);

  // tell BOUT which variables to evolve
  
  bout_solve(rho, "density");
  bout_solve(p, "pressure");
  v.covariant = true; // evolve covariant components
  bout_solve(v, "v");
  B.covariant = false; // evolve contravariant components
  bout_solve(B, "B");

  output.write("dx[0,0] = %e, dy[0,0] = %e, dz = %e\n", 
	       mesh->dx[0][0], mesh->dy[0][0], mesh->dz);

  dump.add(divB, "divB", 1);
  
  divB.setBoundary("DivB"); // Set boundary conditions from options
  
  if(!restarting) {
    // Set variables to these values (+ the initial perturbation)
    // NOTE: This must be after the calls to bout_solve
    rho += rho0;
    p += p0;
    v += v0;
    B += B0;
    
    // Added this for modifying the Orszag-Tang vortex problem
    BoutReal v_fact;
    OPTION(options, v_fact, 1.0);
    v *= v_fact;
  }

  return 0;
}

int physics_run(BoutReal t) {
  // Communicate variables
  mesh->communicate(v, B, p, rho);

  msg_stack.push("F_rho");
  
  ddt(rho) = -V_dot_Grad(v, rho) - rho*Div(v);

  msg_stack.pop(); msg_stack.push("F_p");

  ddt(p) = -V_dot_Grad(v, p) - gamma*p*Div(v);
  
  msg_stack.pop(); msg_stack.push("F_v");
  
  ddt(v) = -V_dot_Grad(v, v) + ((Curl(B)^B) - Grad(p))/rho;

  if(include_viscos) {
    ddt(v).x += viscos * Laplacian(v.x);
    ddt(v).y += viscos * Laplacian(v.y);
    ddt(v).z += viscos * Laplacian(v.z);
  }
  
  msg_stack.pop(); msg_stack.push("F_B");
  
  ddt(B) = Curl(v^B);

  msg_stack.pop(); msg_stack.push("DivB");
  

  return 0;
}
