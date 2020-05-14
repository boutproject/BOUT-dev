#include <bout/physicsmodel.hxx>
#include <derivs.hxx>
#include <initialprofiles.hxx>
#include <invert_laplace.hxx>
#include <smoothing.hxx>

class ESEL : public PhysicsModel {
private:
  Field3D n, vort; // Evolving density, temp and vorticity
  Field3D N;       // ln(n)
  Field3D phi;
  Field2D B;                            // Magnetic field
  BoutReal D, mu;                       // Diffusion coefficients
  Field2D sigma_n, sigma_T, sigma_vort; // dissipation terms
  BoutReal zeta;                        // rho/R0
  BRACKET_METHOD bm;                    // Bracket method for advection terms
  std::unique_ptr<Laplacian> phiSolver{nullptr};                 // Laplacian solver for vort -> phi
  bool test_laplacian; // If true, compute the error on the Laplacian inversion and abort
  Field3D vort_error;

protected:
  int init(bool UNUSED(restart)) {

    auto& options = Options::root()["esel"];

    zeta = options["zeta"].withDefault(2.15e-3);
    D = options["D"].withDefault(1.97e-3);
    mu = options["mu"].withDefault(3.88e-2);
    test_laplacian = options["test_laplacian"].withDefault(false);

    // Set sources and sinks from input profile
    initial_profile("sigma_n", sigma_n);
    initial_profile("sigma_T", sigma_T);
    initial_profile("sigma_vort", sigma_vort);
    initial_profile("B", B);

    SAVE_ONCE(sigma_n);

    // Poisson brackets: b_hat x Grad(f) dot Grad(g) / B = [f, g]
    // Method to use: BRACKET_ARAKAWA, BRACKET_STD or BRACKET_SIMPLE
    // Choose method to use for Poisson bracket advection terms

    switch (options["bracket"].withDefault(2)) {
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
      output << "ERROR: Invalid choice of bracket method. Must be 0-3\n";
      return 1;
    }

    Coordinates* coord = sigma_n.getCoordinates();

    // generate coordinate system
    coord->Bxy = 1;

    coord->g11 = 1.0;
    coord->g22 = 1.0;
    coord->g33 = 1.0;
    coord->g12 = 0.0;
    coord->g13 = 0.0;
    coord->g23 = 0.0;

    coord->g_11 = 1.0;
    coord->g_22 = 1.0;
    coord->g_33 = 1.0;
    coord->g_12 = 0.0;
    coord->g_13 = 0.0;
    coord->g_23 = 0.0;

    coord->geometry();

    SOLVE_FOR(N, vort);
    SAVE_REPEAT(phi);
    if (test_laplacian) {
      SAVE_REPEAT(vort_error);
    }
    phiSolver = Laplacian::create();
    phi = 0.0; // Starting phi

    return 0;
  }

  Field3D C(const Field3D& f) { return zeta * DDZ(f); }

  int rhs(BoutReal UNUSED(time)) {
    N.getMesh()->communicate(N, vort);

    phiSolver->setCoefC2(N);
    phi = phiSolver->solve(vort, phi);

    phi.getMesh()->communicate(phi);

    if (test_laplacian) {

      Field3D vort2 = D2DX2(phi) + D2DZ2(phi) + DDX(N) * DDX(phi) + DDZ(N) * DDZ(phi);
      vort_error = (vort - vort2);

      bout::globals::dump.write();

      MPI_Barrier(BoutComm::get());

      return 1;
    }

    // Continuity equation:
    ddt(N) = bracket(phi, N, bm) + C(phi) - C(N) + D * Delp2(N) - sigma_n;

    // Vorticity equation:
    ddt(vort) = bracket(phi, vort, bm) - C(exp(N)) + mu * Delp2(vort) - sigma_vort * vort;

    // n.b bracket terms do not have minus sign before them because
    // B is pointing in -ve y direction in BOUT coordinates.
    //// This may be wrong, but it is probably consistently wrong

    return 0;
  }
};

BOUTMAIN(ESEL);
