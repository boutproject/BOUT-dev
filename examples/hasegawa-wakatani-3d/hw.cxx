/// 3D simulations of HW
///
/// This version uses indexed operators
/// which reduce the number of loops over the domain
///
/// GPU processing is enabled if BOUT_ENABLE_CUDA is defined
/// Profiling markers and ranges are set if USE_NVTX is defined
/// Based on Ben Dudson, Steven Glenn code, Yining Qin update 0521-2020

#include <iostream>

#include <bout/physicsmodel.hxx>
#include <bout/single_index_ops.hxx>
#include <invert_laplace.hxx>

#define DISABLE_RAJA 0
#include <bout/rajalib.hxx>

class HW3D : public PhysicsModel {
private:
  Field3D n, vort; // Evolving density and vorticity
  Field3D phi;     // Electrostatic potential

  // Model parameters
  BoutReal alpha;                       // Adiabaticity (~conductivity)
  BoutReal kappa;                       // Density gradient drive
  BoutReal Dvort, Dn;                   // Diffusion
  std::unique_ptr<Laplacian> phiSolver; // Laplacian solver for vort -> phi

public:
  // Note: rhs() function must be public, so that RAJA can use CUDA

  int init(bool UNUSED(restart)) override {

    auto& options = Options::root()["hw"];
    alpha = options["alpha"].withDefault(1.0);
    kappa = options["kappa"].withDefault(0.1);
    Dvort = options["Dvort"].doc("Vorticity diffusion (normalised)").withDefault(1e-2);
    Dn = options["Dn"].doc("Density diffusion (normalised)").withDefault(1e-2);

    SOLVE_FOR(n, vort);
    SAVE_REPEAT(phi);
    phiSolver = Laplacian::create(nullptr, CELL_CENTRE, mesh, solver, &dump);

    phi = 0.; // Starting phi

    return 0;
  }

  int rhs(BoutReal UNUSED(time)) override {
    // Solve for potential
    phi = phiSolver->solve(vort, phi);
    Field3D phi_minus_n = phi - n;
    // Communicate variables
    mesh->communicate(n, vort, phi, phi_minus_n);

    // Create accessors which enable fast access
    auto n_acc = FieldAccessor<>(n);
    auto vort_acc = FieldAccessor<>(vort);
    auto phi_acc = FieldAccessor<>(phi);
    auto phi_minus_n_acc = FieldAccessor<>(phi_minus_n);

    // Note: Capture class members, otherwise if accessed via `this` pointer
    //       an illegal memory access error may occur on a GPU
    BOUT_FOR_RAJA(i, n.getRegion("RGN_NOBNDRY"), CAPTURE(alpha, kappa, Dn, Dvort)) {

      BoutReal div_current = alpha * Div_par_Grad_par(phi_minus_n_acc, i);

      ddt(n_acc)[i] = -bracket(phi_acc, n_acc, i) - div_current - kappa * DDZ(phi_acc, i)
                  + Dn * Delp2(n_acc, i);

      ddt(vort_acc)[i] =
          -bracket(phi_acc, vort_acc, i) - div_current + Dvort * Delp2(vort_acc, i);

    };

    return 0;
  }
};

// Define a main() function
BOUTMAIN(HW3D);
