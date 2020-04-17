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

/// 2D drift-reduced model, mainly used for blob studies
///
///
class Blob2D : public PhysicsModel {
private:
  // Evolving variables
  Field3D n, omega; ///< Density and vorticity

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

  std::unique_ptr<Laplacian> phiSolver{nullptr}; ///< Performs Laplacian inversions to calculate phi

protected:
  int init(bool UNUSED(restarting)) {

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

    output.write("\n\n\t----------Parameters: ------------ \n\tOmega_i = {:e} /s,\n\tc_s = "
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

    SOLVE_FOR(n, omega);

    // Output phi
    SAVE_REPEAT(phi);
    SAVE_ONCE(rho_s, c_s, Omega_i);

    return 0;
  }

  int rhs(BoutReal UNUSED(t)) {

    // Run communications
    ////////////////////////////////////////////////////////////////////////////
    mesh->communicate(n, omega);

    // Invert div(n grad(phi)) = grad(n) grad(phi) + n Delp_perp^2(phi) = omega
    ////////////////////////////////////////////////////////////////////////////

    if (!boussinesq) {
      // Including full density in vorticit inversion
      phiSolver->setCoefC(n); // Update the 'C' coefficient. See invert_laplace.hxx
      phi = phiSolver->solve(omega / n, phi); // Use previous solution as guess
    } else {
      // Background density only (1 in normalised units)
      phi = phiSolver->solve(omega);
    }

    mesh->communicate(phi);

    // Make sure fields have Coordinates
    // This sets the Field::fast_coords member to a Coordinate*
    // Not a long-term solution, but here until a better solution is found.
    n.fast_coords = n.getCoordinates();
    omega.fast_coords = omega.getCoordinates();
    phi.fast_coords = phi.getCoordinates();
    
    // Allocate arrays to store the time derivatives
    ddt(n).allocate();
    ddt(omega).allocate();
    // Iterate over the mesh except boundaries because derivatives are needed
    BOUT_FOR(i, n.getRegion("RGN_NOBNDRY")) {

      // Density Evolution
      /////////////////////////////////////////////////////////////////////////////

      ddt(n)[i] = -bracket(phi, n, i)             // ExB term
                  + 2 * DDZ(n, i) * (rho_s / R_c) // Curvature term
                  + D_n * Delp2(n, i);            // Diffusion term

      if (compressible) {
        ddt(n)[i] -= 2 * n[i] * DDZ(phi, i) * (rho_s / R_c); // ExB Compression term
      }

      if (sheath) {
        // Sheath closure
        ddt(n)[i] += n[i] * phi[i] * (rho_s / L_par);
      }

      // Vorticity evolution
      /////////////////////////////////////////////////////////////////////////////

      ddt(omega)[i] = -bracket(phi, omega, i)                // ExB term
                   + 2 * DDZ(n, i) * (rho_s / R_c) / n[i] // Curvature term
                   + D_vort * Delp2(omega, i) / n[i]      // Viscous diffusion term
          ;

      if (sheath) {
        ddt(omega)[i] += phi[i] * (rho_s / L_par);
      }
    }
    return 0;
  }
};

// Define a standard main() function
BOUTMAIN(Blob2D);
