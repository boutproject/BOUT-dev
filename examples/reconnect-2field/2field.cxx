/*****************************************************************************
 * 2 field (Apar, vorticity) model for benchmarking
 * simple slab reconnection model
 *****************************************************************************/

#include <bout/physicsmodel.hxx>

#include <interpolation.hxx>
#include <invert_laplace.hxx>
#include <invert_parderiv.hxx>
#include <initialprofiles.hxx>
#include <bout/constants.hxx>

class TwoField : public PhysicsModel {
private:
  // 2D initial profiles
  Field2D Jpar0, Te0, Ni0;

  // 3D evolving fields
  Field3D U, Apar;

  // Derived 3D variables
  Field3D phi, jpar;

  // External coil field
  Field3D Apar_ext, Jpar_ext, Phi0_ext, U0_ext;

  // Metric coefficients
  Field2D Rxy, Bpxy, Btxy, hthe;

  // Constants
  const BoutReal MU0 = 4.0e-7 * PI;
  const BoutReal Charge = 1.60217646e-19;   // electron charge e (C)
  const BoutReal Mi = 2.0 * 1.67262158e-27; // Ion mass
  const BoutReal Me = 9.1093816e-31;        // Electron mass
  const BoutReal Me_Mi = Me / Mi;           // Electron mass / Ion mass

  // normalisation parameters
  BoutReal Tenorm, Nenorm, Bnorm;
  BoutReal Cs, rho_s, wci, beta_hat;

  BoutReal eta, mu;

  // Poisson brackets: b0 x Grad(f) dot Grad(g) / B = [f, g]
  // Method to use: BRACKET_ARAKAWA, BRACKET_STD or BRACKET_SIMPLE
  BRACKET_METHOD bm; // Bracket method for advection terms


  bool nonlinear;
  bool include_jpar0;
  int jpar_bndry;

  std::unique_ptr<InvertPar> inv{nullptr}; // Parallel inversion class used in preconditioner

  // Coordinate system metric
  Coordinates *coord, *coord_ylow;

  // Inverts a Laplacian to get potential
  std::unique_ptr<Laplacian> phiSolver{nullptr};
  
protected:
  int init(bool UNUSED(restarting)) override {

    // Load 2D profiles
    GRID_LOAD(Jpar0, Te0, Ni0);
    Ni0 *= 1e20; // To m^-3
    
    // Coordinate system
    coord = mesh->getCoordinates();
    coord = mesh->getCoordinates(CELL_YLOW);

    // Load metrics
    GRID_LOAD(Rxy, Bpxy, Btxy, hthe);
    mesh->get(coord->Bxy, "Bxy");

    // Set locations of staggered fields
    Apar.setLocation(CELL_YLOW);

    // Read some parameters
    auto& options = Options::root()["2field"];

    // normalisation values
    nonlinear = options["nonlinear"].withDefault(false);
    include_jpar0 = options["include_jpar0"].withDefault(true);
    jpar_bndry = options["jpar_bndry"].withDefault(0);

    eta = options["eta"].doc("Normalised resistivity").withDefault(1e-3);
    mu = options["mu"].doc("Normalised vorticity").withDefault(1.e-3);

    switch (options["bracket_method"].withDefault(0)) {
    case 0: {
      bm = BRACKET_STD;
      output << "\tBrackets: default differencing\n";
      break;
    }
    case 1: {
      bm = BRACKET_SIMPLE;
      output << "\tBrackets: simplified operator\n";
      break;
    }
    case 2: {
      bm = BRACKET_ARAKAWA;
      output << "\tBrackets: Arakawa scheme\n";
      break;
    }
    case 3: {
      bm = BRACKET_CTU;
      output << "\tBrackets: Corner Transport Upwind method\n";
      break;
    }
    default:
      output << "ERROR: Invalid choice of bracket method. Must be 0 - 3\n";
      return 1;
    }

    ///////////////////////////////////////////////////
    // Normalisation

    Tenorm = max(Te0, true);
    if (Tenorm < 1) {
      Tenorm = 1000;
    }
    
    Nenorm = max(Ni0, true);
    if (Nenorm < 1) {
      Nenorm = 1.e19;
    }
    
    Bnorm = max(coord->Bxy, true);

    // Sound speed in m/s
    Cs = sqrt(Charge * Tenorm / Mi);

    // drift scale
    rho_s = Cs * Mi / (Charge * Bnorm);

    // Ion cyclotron frequency
    wci = Charge * Bnorm / Mi;

    beta_hat = MU0 * Charge * Tenorm * Nenorm / (Bnorm * Bnorm);

    output << "\tNormalisations:" << endl;
    output << "\tCs       = " << Cs << endl;
    output << "\trho_s    = " << rho_s << endl;
    output << "\twci      = " << wci << endl;
    output << "\tbeta_hat = " << beta_hat << endl;

    SAVE_ONCE(Tenorm, Nenorm, Bnorm);
    SAVE_ONCE(Cs, rho_s, wci, beta_hat);

    // Normalise geometry
    Rxy /= rho_s;
    hthe /= rho_s;
    coord->dx /= rho_s * rho_s * Bnorm;

    // Normalise magnetic field
    Bpxy /= Bnorm;
    Btxy /= Bnorm;
    coord->Bxy /= Bnorm;

    // Plasma quantities
    Jpar0 /= Nenorm * Charge * Cs;

    // CALCULATE METRICS

    coord->g11 = SQ(Rxy * Bpxy);
    coord->g22 = 1.0 / SQ(hthe);
    coord->g33 = SQ(coord->Bxy) / coord->g11;
    coord->g12 = 0.0;
    coord->g13 = 0.;
    coord->g23 = -Btxy / (hthe * Bpxy * Rxy);

    coord->J = hthe / Bpxy;

    coord->g_11 = 1.0 / coord->g11;
    coord->g_22 = SQ(coord->Bxy * hthe / Bpxy);
    coord->g_33 = Rxy * Rxy;
    coord->g_12 = 0.;
    coord->g_13 = 0.;
    coord->g_23 = Btxy * hthe * Rxy / Bpxy;

    coord->geometry();

    // Tell BOUT++ which variables to evolve
    SOLVE_FOR(U, Apar);

    // Set boundary conditions
    jpar.setBoundary("jpar");
    phi.setBoundary("phi");

    // Add any other variables to be dumped to file
    SAVE_REPEAT(phi, jpar);
    SAVE_ONCE(Jpar0);

    // Generate external field

    initial_profile("Apar_ext", Apar_ext);
    Jpar_ext = -Delp2(Apar_ext);
    SAVE_ONCE(Apar_ext, Jpar_ext);

    initial_profile("Phi0_ext", Phi0_ext);
    U0_ext = -Delp2(Phi0_ext) / coord->Bxy;
    SAVE_ONCE(Phi0_ext, U0_ext);

    // Give the solver the preconditioner function
    setPrecon((preconfunc)&TwoField::precon);
    
    // Initialise parallel inversion class
    inv = InvertPar::create();
    inv->setCoefA(1.0);
    U.setBoundary("U");
    Apar.setBoundary("Apar");

    // Create a solver for the Laplacian
    phiSolver = Laplacian::create();
    
    return 0;
  }

  const Field3D Grad_parP(const Field3D &f, CELL_LOC loc = CELL_DEFAULT) {
    Field3D result;
    if (nonlinear) {
      result = ::Grad_parP((Apar + Apar_ext) * beta_hat, f);
    } else {
      result = ::Grad_parP(Apar_ext * beta_hat, f);
    }
    if (mesh->StaggerGrids) {
      mesh->communicate(result);
    }
    return interp_to(result, loc);
  }

  int rhs(BoutReal UNUSED(time)) override {
    // Solve EM fields

    // U = (1/B) * Delp2(phi)
    phi = phiSolver->solve(coord->Bxy * U);
    phi.applyBoundary(); // For target plates only

    mesh->communicate(U, phi, Apar);

    jpar = -Delp2(Apar + Apar_ext); // total Apar
    jpar.applyBoundary();
    mesh->communicate(jpar);

    if (jpar_bndry > 0) {
      // Boundary in jpar
      if (mesh->firstX()) {
        for (int i = jpar_bndry; i >= 0; i--)
          for (int j = 0; j < mesh->LocalNy; j++)
            for (int k = 0; k < mesh->LocalNz; k++) {
              jpar(i, j, k) = jpar(i + 1, j, k);
            }
      }
      if (mesh->lastX()) {
        for (int i = mesh->LocalNx - jpar_bndry - 1; i < mesh->LocalNx; i++)
          for (int j = 0; j < mesh->LocalNy; j++)
            for (int k = 0; k < mesh->LocalNz; k++) {
              jpar(i, j, k) = jpar(i - 1, j, k);
            }
      }
    }

    // VORTICITY
    ddt(U) = SQ(coord->Bxy) * Grad_parP(jpar / coord_ylow->Bxy, CELL_CENTRE);

    if (include_jpar0) {
      ddt(U) -= SQ(coord->Bxy) * beta_hat *
                interp_to(
                    bracket(Apar + Apar_ext, Jpar0 / coord_ylow->Bxy, BRACKET_ARAKAWA),
                    CELL_CENTRE
                );
    }

    ddt(U) -= bracket(Phi0_ext, U, bm); // ExB advection
    // ddt(U) -= bracket(phi, U0_ext, bm); // ExB advection
    if (nonlinear) {
      ddt(U) -= bracket(phi, U, bm); // ExB advection
    }

    if (mu > 0.)
      ddt(U) += mu * Delp2(U);

    // APAR

    ddt(Apar) = -Grad_parP(phi, CELL_YLOW) / beta_hat;
    ddt(Apar) += -Grad_parP(Phi0_ext, CELL_YLOW) / beta_hat;

    if (eta > 0.)
      ddt(Apar) -= eta * jpar / beta_hat;

    return 0;
  }

public:
  /*********************************************************
   * Preconditioner
   *
   * o System state in variables (as in rhs function)
   * o Values to be inverted in time derivatives
   *
   * o Return values should be in time derivatives
   *
   *********************************************************/
  int precon(BoutReal UNUSED(t), BoutReal gamma, BoutReal UNUSED(delta)) {
    mesh->communicate(ddt(Apar));
    Field3D Jp = -Delp2(ddt(Apar));
    mesh->communicate(Jp);

    if (jpar_bndry > 0) {
      // Boundary in jpar
      if (mesh->firstX()) {
        for (int i = jpar_bndry; i >= 0; i--)
          for (int j = 0; j < mesh->LocalNy; j++)
            for (int k = 0; k < mesh->LocalNz; k++) {
              Jp(i,j,k) = Jp(i + 1,j,k);
            }
      }
      if (mesh->lastX()) {
        for (int i = mesh->LocalNx - jpar_bndry - 1; i < mesh->LocalNx; i++)
          for (int j = 0; j < mesh->LocalNy; j++)
            for (int k = 0; k < mesh->LocalNz; k++) {
              Jp(i,j,k) = Jp(i - 1,j,k);
            }
      }
    }

    Field3D U1 = ddt(U) + gamma * SQ(coord->Bxy)
                          * Grad_par(Jp / coord_ylow->Bxy, CELL_CENTRE);

    inv->setCoefB(-SQ(gamma * coord->Bxy) / beta_hat);
    ddt(U) = inv->solve(U1);
    ddt(U).applyBoundary();

    Field3D phip = phiSolver->solve(coord->Bxy * ddt(U));
    mesh->communicate(phip);

    ddt(Apar) = ddt(Apar) - (gamma / beta_hat) * Grad_par(phip, CELL_YLOW);
    ddt(Apar).applyBoundary();

    return 0;
  }
};

BOUTMAIN(TwoField);
