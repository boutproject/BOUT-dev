/*******************************************************************
        2D simulations

         NR Walkden, B Dudson  20 January 2012
 *******************************************************************/

#include "bout.hxx"                  // Commonly used BOUT++ components
#include "bout/invert/laplacexz.hxx" // Laplacian inversion
#include "derivs.hxx"                // To use DDZ()

class Blob2D : public PhysicsModel {
private:
  // Evolving variables
  Field3D n, omega; // Density and vorticity

  // Auxilliary variables
  Field3D phi; // Electrostatic potential

  // Parameters
  // Bohm gyro radius, Ion cyclotron frequency, Bohm sound speed
  BoutReal rho_s, Omega_i, c_s, n0; 

  // Constants to calculate the parameters
  BoutReal Te0, e, B0, D_n, D_vort, m_i, m_e;
  BoutReal R_c;   // Radius of curvature
  BoutReal L_par; // Parallel connection length

  // Model options
  bool boussinesq;   // Use the Boussinesq approximation in vorticity
  bool compressible; // If allow inclusion of n grad phi term in density evolution
  bool sheath;       // Sheath connected?

  LaplaceXZ *phiSolver;

  int boussinesq_reuse; // Determines how long between updates of the density in the
                        // vorticity
  int boussinesq_used;  // How many times has it been reused

protected:
  int init(bool UNUSED(restarting)) {

    /******************Reading options *****************/

    auto globalOptions = Options::root();
    auto options = globalOptions["model"];

    // Load system parameters
    Te0 = options["Te0"].withDefault(30); // Temp in eV
    e = options["e"].withDefault(1.602e-19);
    m_i = options["m_i"].withDefault(2 * 1.667e-27);
    m_e = options["m_e"].withDefault(9.11e-31);

    n0 = options["n0"].withDefault(1e19);      // Background density in cubic m
    D_vort = options["D_vort"].withDefault(0); // Viscous diffusion coefficient
    D_n = options["D_n"].withDefault(0);       // Density diffusion coefficient

    R_c = options["R_c"].withDefault(1.5);    // Radius of curvature
    L_par = options["L_par"].withDefault(10); // Parallel connection length
    B0 = options["B0"].withDefault(0.35);     // Value of magnetic field strength

    // System option switches

    compressible = options["compressible"].withDefault(
        false); // Include compressible ExB term in density equation
    boussinesq = options["boussinesq"].withDefault(
        true); // Use Boussinesq approximation in vorticity
    sheath = options["sheath"].withDefault(true); // Sheath closure

    boussinesq_reuse = options["boussinesq_reuse"].withDefault(
        0);                                 // How many times to reuse n in vorticity?
    boussinesq_used = boussinesq_reuse + 1; // Ensure updated first time

    /***************Calculate the Parameters **********/

    Omega_i = e * B0 / m_i;    // Cyclotron Frequency
    c_s = sqrt(e * Te0 / m_i); // Bohm sound speed
    rho_s = c_s / Omega_i;     // Bohm gyro-radius

    output.write("\n\n\t----------Parameters: ------------ \n\tOmega_i = %e /s,\n\tc_s = "
                 "%e m/s,\n\trho_s = %e m\n",
                 Omega_i, c_s, rho_s);

    // Calculate delta_*, blob size scaling
    output.write("\tdelta_* = rho_s * (dn/n) * %e\n",
                 pow(L_par * L_par / (R_c * rho_s), 1. / 5));

    /************ Create a solver for potential ********/

    if (boussinesq) {
      // Use options in BOUT.inp section "phiBoussinesq"
      phiSolver = LaplaceXZ::create(mesh, &Options::root()["phiBoussinesq"]);
      // Set the coefficients once here
      phiSolver->setCoefs(Field2D(1.0), Field2D(0.0));
    } else {
      // Use options in BOUT.inp section "phiSolver"
      phiSolver = LaplaceXZ::create(mesh, &Options::root()["phiSolver"]);
      // Coefficients will be set every RHS call
    }
    phi = 0.0; // Starting guess for first solve (if iterative)

    /************ Tell BOUT++ what to solve ************/

    SOLVE_FOR2(n, omega);

    // Output phi
    SAVE_REPEAT(phi);
    SAVE_ONCE3(rho_s, c_s, Omega_i);

    return 0;
  }

  int rhs(BoutReal UNUSED(t)) {
    //Coordinates *coord = mesh->getCoordinates();

    // Run communications
    ////////////////////////////////////////////////////////////////////////////
    mesh->communicate(n, omega);

    // Invert div(n grad(phi)) = grad(n) grad(phi) + n Delp_perp^2(phi) = omega
    ////////////////////////////////////////////////////////////////////////////

    if (!boussinesq) {
      // Including full density in vorticity inversion
      boussinesq_used++;
      if (boussinesq_used > boussinesq_reuse) {
        // Update density
        phiSolver->setCoefs(n, 0.0);
        boussinesq_used = 0;
      }
      phi = phiSolver->solve(omega, phi); // Use previous solution as guess
    } else {
      // Background density only (1 in normalised units)
      // Coefficients already set in setup
      phi = phiSolver->solve(omega, phi);
    }

    mesh->communicate(phi);

    // Density Evolution
    /////////////////////////////////////////////////////////////////////////////

    ddt(n) = -bracket(phi, n, BRACKET_SIMPLE) // ExB term
             + 2 * DDZ(n) * (rho_s / R_c);     // Curvature term
    ddt(n) += D_n*Delp2(n, CELL_DEFAULT, false);

    // if(coord->is3D()){
    //   ddt(n) += Div_Perp_Lap_FV(D_n, n);   // Diffusion term
    // }else{
    //   ddt(n) += D_n*Delp2(n);
    // }
    
    if (compressible) {
      ddt(n) -= 2 * n * DDZ(phi) * (rho_s / R_c); // ExB Compression term
    }

    if (sheath) {
      ddt(n) +=
          n * phi * (rho_s / L_par); // - (n - 1)*(rho_s/L_par);      // Sheath closure
    }

    // Vorticity evolution
    /////////////////////////////////////////////////////////////////////////////

    ddt(omega) = -bracket(phi, omega, BRACKET_SIMPLE) // ExB term
                 + 2 * DDZ(n) * (rho_s / R_c) / n;     // Curvature term
    ddt(omega) += D_vort * Delp2(omega, CELL_DEFAULT, false)/n; // Viscous diffusion term

    // if(coord->is3D()){
    //   ddt(omega) += Div_Perp_Lap_FV(D_vort ,omega) / n ;  // Viscous diffusion term
    // }else{
    //   ddt(omega) += D_vort * Delp2(omega)/n; // Viscous diffusion term
    // }

    if (sheath) {
      ddt(omega) += phi * (rho_s / L_par);
    }

    return 0;
  }
};

BOUTMAIN(Blob2D);
