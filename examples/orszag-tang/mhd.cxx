/**************************************************************************
 * Ideal MHD physics module for BOUT++
 * This version evolves the entire quantity (initial + perturbed)
 **************************************************************************/

#include <bout/physicsmodel.hxx>

class MHD : public PhysicsModel {
private:
  // 3D evolving variables
  Field3D rho, p; // density, pressure
  Vector3D v, B;  // velocity, magnetic field

  Field3D divB; // Divergence of B (for monitoring)

  // parameters
  BoutReal g;
  bool include_viscos;
  BoutReal viscos;

  int init(bool restarting) override {
    // 2D initial profiles
    Field2D rho0, p0;
    Vector2D v0, B0;

    // read options
    auto globalOptions = Options::root();
    auto options = globalOptions["mhd"];
    OPTION(options, g, 5.0 / 3.0);
    OPTION(options, include_viscos, false);
    OPTION(options, viscos, 0.1);

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

    Coordinates *coord = mesh->coordinates();
    output.write("dx[0,0] = %e, dy[0,0] = %e, dz = %e\n", coord->dx(0, 0),
                 coord->dy(0, 0), coord->dz);

    SAVE_REPEAT(divB);

    divB.setBoundary("DivB"); // Set boundary conditions from options

    if (!restarting) {
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

  /// This function is called every output, before
  /// the data is written to file. It can therefore be used
  /// to calculate diagnostics
  ///
  /// @param[in] simtime   Simulation time
  /// @param[in] iter      Output step
  /// @param[in] NOUT      Total number of outputs requested
  int outputMonitor(BoutReal UNUSED(simtime), int UNUSED(iter),
                    int UNUSED(NOUT)) override {
    // Calculate divergence of magnetic field
    divB = Div(B);
    return 0;
  }

  int rhs(BoutReal t) override {
    // Communicate variables

    mesh->communicate(v, B, p, rho);

    {
      TRACE("ddt(rho)");
      ddt(rho) = -V_dot_Grad(v, rho) - rho * Div(v);
    }
    {
      TRACE("ddt(p)");
      ddt(p) = -V_dot_Grad(v, p) - g * p * Div(v);
    }
    {
      TRACE("ddt(v)");
      ddt(v) = -V_dot_Grad(v, v) + (cross(Curl(B), B) - Grad(p)) / rho;

      if (include_viscos) {
        ddt(v).x += viscos * Laplace(v.x);
        ddt(v).y += viscos * Laplace(v.y);
        ddt(v).z += viscos * Laplace(v.z);
      }
    }
    {
      TRACE("ddt(B)");

      ddt(B) = Curl(cross(v, B));
    }

    return 0;
  }
};

BOUTMAIN(MHD);
