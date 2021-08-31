/*!******************************************************************
 * \file blob2d.cxx
 *
 *       2D simulations
 *
 *        NR Walkden, B Dudson  20 January 2012
 *******************************************************************/

#include <bout/physicsmodel.hxx> // Commonly used BOUT++ components
#include <derivs.hxx>            // To use DDZ()
#include <invert_laplace.hxx>    // Laplacian inversion

#include "bout/single_index_ops.hxx" // Operators at a single index

#define RUN_WITH_RAJA 0

#ifdef BOUT_HAS_RAJA
#include "RAJA/RAJA.hpp" // using RAJA lib
#endif

#if defined(BOUT_USE_CUDA) && defined(__CUDACC__)
#include <cuda_profiler_api.h>
#endif

/// 2D drift-reduced model, mainly used for blob studies
///
///
class Blob2D : public PhysicsModel {
public:
  // Evolving variables
  Field3D n, vort; ///< Density and vorticity

  // Auxilliary variables
  Field3D phi; ///< Electrostatic potential

  // Parameters
  BoutReal rho_s;   ///< Bohm gyro radius
  BoutReal Omega_i; ///< Ion cyclotron frequency
  BoutReal c_s;     ///< Bohm sound speed
  BoutReal n0;      ///< Reference density

  // Constants to calculate the parameters
  BoutReal Te0; ///< Isothermal temperature [eV]
  BoutReal B0;  ///< Constant magnetic field [T]
  BoutReal m_i; ///< Ion mass [kg]
  BoutReal m_e; ///< Electron mass [kg]
  BoutReal e;   ///< Electron charge [C]

  BoutReal D_n, D_vort; ///< Diffusion coefficients
  BoutReal R_c;         ///< Radius of curvature
  BoutReal L_par;       ///< Parallel connection length

  // Model options
  bool boussinesq;   ///< Use the Boussinesq approximation in vorticity
  bool compressible; ///< If allow inclusion of n grad phi term in density evolution
  bool sheath;       ///< Sheath connected?

  std::unique_ptr<Laplacian> phiSolver{
      nullptr}; ///< Performs Laplacian inversions to calculate phi

  int init(bool UNUSED(restarting)) override {

    /******************Reading options *****************/

    auto& globalOptions = Options::root();
    auto& options = globalOptions["model"];

    // Load system parameters
    Te0 = options["Te0"].doc("Temperature in eV").withDefault(30.0);

    e = options["e"].withDefault(1.602e-19);
    m_i = options["m_i"].withDefault(2 * 1.667e-27);
    m_e = options["m_e"].withDefault(9.11e-31);

    n0 = options["n0"].doc("Background density in cubic m").withDefault(1e19);
    D_vort = options["D_vort"].doc("Viscous diffusion coefficient").withDefault(0.0);
    D_n = options["D_n"].doc("Density diffusion coefficient").withDefault(0.0);

    R_c = options["R_c"].doc("Radius of curvature").withDefault(1.5);
    L_par = options["L_par"].doc("Parallel connection length").withDefault(10.0);

    B0 = options["B0"].doc("Value of magnetic field strength").withDefault(0.35);

    // System option switches

    compressible = options["compressible"]
                       .doc("Compressible ExB term in density equation")
                       .withDefault(false);
    boussinesq = options["boussinesq"]
                     .doc("Use Boussinesq approximation in vorticity")
                     .withDefault(true);
    sheath = options["sheath"].doc("Sheath closure").withDefault(true);

    /***************Calculate the Parameters **********/

    Omega_i = e * B0 / m_i;    // Cyclotron Frequency
    c_s = sqrt(e * Te0 / m_i); // Bohm sound speed
    rho_s = c_s / Omega_i;     // Bohm gyro-radius

    output.write(
        "\n\n\t----------Parameters: ------------ \n\tOmega_i = {:e} /s,\n\tc_s = "
        "{:e} m/s,\n\trho_s = {:e} m\n",
        Omega_i, c_s, rho_s);

    // Calculate delta_*, blob size scaling
    output.write("\tdelta_* = rho_s * (dn/n) * {:e} ",
                 pow(L_par * L_par / (R_c * rho_s), 1. / 5));

    /************ Create a solver for potential ********/

    if (boussinesq) {
      // BOUT.inp section "phiBoussinesq"
      phiSolver = Laplacian::create(Options::getRoot()->getSection("phiBoussinesq"));
    } else {
      // BOUT.inp section "phiSolver"
      phiSolver = Laplacian::create(Options::getRoot()->getSection("phiSolver"));
    }
    phi = 0.0; // Starting guess for first solve (if iterative)

    /************ Tell BOUT++ what to solve ************/

    SOLVE_FOR(n, vort);

    // Output phi
    SAVE_REPEAT(phi);
    SAVE_ONCE(rho_s, c_s, Omega_i);

    return 0;
  }

  int rhs(BoutReal UNUSED(t)) override {

    // Run communications
    ////////////////////////////////////////////////////////////////////////////
    mesh->communicate(n, vort);

    // Invert div(n grad(phi)) = grad(n) grad(phi) + n Delp_perp^2(phi) = omega
    ////////////////////////////////////////////////////////////////////////////

    if (!boussinesq) {
      // Including full density in vorticit inversion
      phiSolver->setCoefC(n); // Update the 'C' coefficient. See invert_laplace.hxx
      phi = phiSolver->solve(vort / n, phi); // Use previous solution as guess
    } else {
      // Background density only (1 in normalised units)
      phi = phiSolver->solve(vort);
    }

    mesh->communicate(phi);

    // Create data accessors for fast inner loop
    auto n_acc = FieldAccessor<>(n);
    auto vort_acc = FieldAccessor<>(vort);
    auto phi_acc = FieldAccessor<>(phi);

    const auto& region = n.getRegion("RGN_NOBNDRY"); // Region object

#if RUN_WITH_RAJA
    auto indices = region.getIndices(); // A std::vector of Ind3D objects
    Ind3D* ob_i = &(indices)[0];

    RAJA::forall<EXEC_POL>(RAJA::RangeSegment(0, indices.size()), [=] RAJA_DEVICE(int id) {
      int i = ob_i[id].ind;
#else
    BOUT_FOR(i, region) {
#endif
      ddt(n_acc)[i] = -bracket(phi_acc, n_acc, i) - 2 * DDZ(n_acc, i) * (rho_s / R_c)
                      + D_n * Delp2(n_acc, i);

      ddt(vort_acc)[i] = -bracket(phi_acc, vort_acc, i)
                         + 2 * DDZ(n_acc, i) * (rho_s / R_c)
                         + D_vort * Delp2(vort_acc, i);

#if RUN_WITH_RAJA
    });
#else
    }
#endif

    if (compressible) {
#if RUN_WITH_RAJA
      RAJA::forall<EXEC_POL>(RAJA::RangeSegment(0, indices.size()), [=] RAJA_DEVICE(int id) {
        int i = ob_i[id].ind;
#else
      BOUT_FOR(i, region) {
#endif
        ddt(n_acc)[i] -= 2 * n_acc[i] * DDZ(phi_acc, i) * (rho_s / R_c); // ExB Compression term
#if RUN_WITH_RAJA
      });
#else
      }
#endif
    }
    
    if (sheath) {
      // Sheath closure
#if RUN_WITH_RAJA
      RAJA::forall<EXEC_POL>(RAJA::RangeSegment(0, indices.size()), [=] RAJA_DEVICE(int id) {
        int i = ob_i[id].ind;
#else
      BOUT_FOR(i, region) {
#endif
        ddt(n_acc)[i] += n_acc[i] * phi_acc[i] * (rho_s / L_par);
        ddt(vort_acc)[i] += phi_acc[i] * (rho_s / L_par);
#if RUN_WITH_RAJA
      });
#else
      }
#endif
    }

    return 0;
  }
};

// Define a standard main() function
BOUTMAIN(Blob2D);
