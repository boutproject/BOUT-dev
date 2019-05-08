/**************************************************************************
 * Similar set of equations to JOREK
 *
 **************************************************************************/

#include <bout/physicsmodel.hxx>

#include <bout/constants.hxx>
#include <invert_laplace.hxx>

class Jorek : public PhysicsModel {
private:
  // Evolving quantities
  Field3D rho, Te, Ti, U, Vpar, Apar;
  // Derived quantities
  Field3D Jpar, phi; // Parallel current, electric potential

  // Equilibrium quantities
  Field2D rho0, Te0, Ti0; // Equilibrium mass density, electron and ion temperature
  Field2D B0, J0, P0;
  Vector2D b0xcv; // Curvature term
  Vector2D B0vec; // B0 field vector

  // Dissipation coefficients
  Field2D D_perp;              // Particle diffusion coefficient
  Field2D chi_eperp, chi_epar; // Electron heat diffusion coefficients
  Field2D chi_iperp, chi_ipar; // Ion heat diffusion coefficients

  // Collisional terms
  BoutReal tau_enorm;
  Field3D tau_e; // electron collision time

  Field2D eta0; // Resistivity
  Field3D eta;

  BoutReal viscos_par, viscos_perp, viscos_coll; // Viscosity coefficients
  BoutReal hyperresist;                          // Hyper-resistivity coefficient

  int phi_flags;

  // Constants
  const BoutReal MU0 = 4.0e-7 * PI;
  const BoutReal Charge = 1.60217646e-19;   // electron charge e (C)
  const BoutReal Mi = 2.0 * 1.67262158e-27; // Ion mass
  const BoutReal Me = 9.1093816e-31;        // Electron mass
  const BoutReal Me_Mi = Me / Mi;           // Electron mass / Ion mass

  // Normalisation factors
  BoutReal Tnorm, rhonorm; // Partial normalisation to rho and MU0. Temperature normalised

  // options

  bool nonlinear;
  bool full_bfield;     // If true, use divergence-free expression for B
  bool flux_method;     // Use flux methods in rho and T equations
  int jpar_bndry_width; // Set jpar = 0 in a boundary region

  bool electron_density; // Solve Ne rather than Ni (adds Jpar term to density)

  bool vorticity_momentum; // Vorticity is curl of momentum, rather than velocity

  bool include_profiles; // Include zero-order equilibrium terms

  bool parallel_lc; // Use CtoL and LtoC differencing

  int low_pass_z; // Toroidal (Z) filtering of all variables

  Vector3D vExB, vD; // Velocities
  Field3D divExB;    // Divergence of ExB flow

  BoutReal Wei; // Factor for the electron-ion collision term
  bool ohmic_heating;

  // Poisson brackets: b0 x Grad(f) dot Grad(g) / B = [f, g]
  // Method to use: BRACKET_ARAKAWA, BRACKET_STD or BRACKET_SIMPLE
  BRACKET_METHOD bm; // Bracket method for advection terms

  // Communication objects
  FieldGroup comms;

  // Coordinate system
  Coordinates* coord;

  // Inverts a Laplacian to get potential
  Laplacian* phiSolver;

  int init(bool UNUSED(restarting)) override {

    output.write("Solving JOREK-like reduced MHD equations\n");
    output.write("\tFile    : %s\n", __FILE__);
    output.write("\tCompiled: %s at %s\n", __DATE__, __TIME__);

    auto globalOptions = Options::root();
    auto options = globalOptions["jorek"];

    //////////////////////////////////////////////////////////////
    // Load data from the grid

    // Load 2D profiles
    mesh->get(J0, "Jpar0"); // A / m^2

    if (mesh->get(rho0, "Ni0")) {
      output << "Warning: No density profile available\n";
      rho0 = options["density"].withDefault<BoutReal>(1.0);
    }
    rho0 *= 1e20; // Convert to m^[-3]

    // Read temperature
    mesh->get(Te0, "Te0");
    mesh->get(Ti0, "Ti0");

    // Try reading pressure profile (in Pascals)
    if (mesh->get(P0, "pressure")) {
      // Just calculate from Temp and density
      P0 = Charge * (Ti0 + Te0) * rho0;
    } else {
      // Make sure that density and temperature are consistent with pressure

      Field2D factor = P0 / (Charge * (Ti0 + Te0) * rho0);

      output.write("\tPressure factor %e -> %e\n", min(factor, true), max(factor, true));

      // Multiply temperatures by this factor
      Te0 *= factor;
      Ti0 *= factor;
    }
    rho0 *= Mi; // Convert density to mass density [kg / m^3]

    // Load dissipation coefficients, override in options file
    if (options["D_perp"].isSet()) {
      D_perp = options["D_perp"].withDefault<BoutReal>(0.0);
    } else
      mesh->get(D_perp, "D_perp");

    if (options["chi_eperp"].isSet()) {
      chi_eperp = options["chi_eperp"].withDefault<BoutReal>(0.0);
    } else
      mesh->get(chi_eperp, "chi_eperp");

    if (options["chi_iperp"].isSet()) {
      chi_iperp = options["chi_iperp"].withDefault<BoutReal>(0.0);
    } else
      mesh->get(chi_iperp, "chi_iperp");

    if (options["chi_epar"].isSet()) {
      chi_epar = options["chi_epar"].withDefault<BoutReal>(0.0);
    } else
      mesh->get(chi_epar, "chi_epar");

    if (options["chi_ipar"].isSet()) {
      chi_ipar = options["chi_ipar"].withDefault<BoutReal>(0.0);
    } else
      mesh->get(chi_ipar, "chi_ipar");

    if (options["viscos_perp"].isSet()) {
      viscos_perp = options["viscos_perp"].withDefault<BoutReal>(-1.0);
    } else
      mesh->get(viscos_perp, "viscos_perp");

    if (options["viscos_par"].isSet()) {
      viscos_par = options["viscos_par"].withDefault<BoutReal>(-1.0);
    } else
      mesh->get(viscos_par, "viscos_par");

    viscos_coll = options["viscos_coll"].withDefault(-1.0);

    // Load curvature term
    b0xcv.covariant = false;  // Read contravariant components
    mesh->get(b0xcv, "bxcv"); // mixed units x: T y: m^-2 z: m^-2

    // Metric coefficients
    Field2D Rxy, Bpxy, Btxy, hthe;
    Field2D I; // Shear factor

    coord = mesh->getCoordinates();

    if (mesh->get(Rxy, "Rxy")) { // m
      output_error.write("Error: Cannot read Rxy from grid\n");
      return 1;
    }
    if (mesh->get(Bpxy, "Bpxy")) { // T
      output_error.write("Error: Cannot read Bpxy from grid\n");
      return 1;
    }
    mesh->get(Btxy, "Btxy"); // T
    mesh->get(B0, "Bxy");    // T
    mesh->get(hthe, "hthe"); // m
    mesh->get(I, "sinty");   // m^-2 T^-1

    nonlinear = options["nonlinear"].withDefault(false);
    full_bfield = options["full_bfield"].withDefault(false);
    flux_method = options["flux_method"].withDefault(false);

    jpar_bndry_width = options["jpar_bndry_width"].withDefault(-1);

    hyperresist = options["hyperresist"].withDefault(-1);

    electron_density = options["electron_density"].withDefault(false);
    vorticity_momentum = options["vorticity_momentum"].withDefault(false);
    include_profiles = options["include_profiles"].withDefault(false);
    parallel_lc = options["parallel_lc"].withDefault(true);

    phi_flags = options["phi_flags"].withDefault(0);

    low_pass_z = options["low_pass_z"].withDefault(-1); // Default is no filtering

    Wei = options["Wei"].withDefault(1.0);

    ohmic_heating = options["ohmic_heating"].withDefault(true);

    switch (options["bracket_method"].withDefault<int>(0)) {
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

    //////////////////////////////////////////////////////////////
    // SHIFTED RADIAL COORDINATES

    // Check type of parallel transform
    std::string ptstr =
        Options::root()["mesh"]["paralleltransform"].withDefault<std::string>("identity");

    if (lowercase(ptstr) == "shifted") {
      // Dimits style, using local coordinate system
      b0xcv.z += I * b0xcv.x;
      I = 0.0; // I disappears from metric
    }

    //////////////////////////////////////////////////////////////
    // NORMALISE QUANTITIES

    rhonorm = max(rho0, true);             // Maximum over all grid
    BoutReal Temax = max(Te0, true);       // Maximum Te value
    Tnorm = Mi / (MU0 * Charge * rhonorm); // Temperature normalisation

    SAVE_ONCE(rhonorm, Tnorm); // Save normalisation factors to file

    // Normalise quantities

    P0 *= MU0;
    J0 *= MU0;
    rho0 /= rhonorm;
    Te0 /= Tnorm;
    Ti0 /= Tnorm;

    viscos_perp *= sqrt(MU0 / rhonorm);
    viscos_par *= sqrt(MU0 / rhonorm);
    D_perp *= sqrt(MU0 * rhonorm);
    chi_eperp *= sqrt(MU0 / rhonorm);
    chi_epar *= sqrt(MU0 / rhonorm);
    chi_iperp *= sqrt(MU0 / rhonorm);
    chi_ipar *= sqrt(MU0 / rhonorm);

    // Coulomb logarithm
    BoutReal CoulombLog = 6.6 - 0.5 * log(rhonorm / (Mi * 1e20)) + 1.5 * log(Temax);
    output << "\tCoulomb logarithm = " << CoulombLog << endl;

    // Factor in front of tau_e expression
    // tau_e = tau_enorm * Tet^1.5 / rhot
    tau_enorm = 3.44e11 * (Mi / rhonorm) * Tnorm * sqrt(Tnorm) / CoulombLog;
    output << "\ttau_enorm = " << tau_enorm;
    tau_enorm /= sqrt(MU0 * rhonorm); // Normalise
    output << "\tNormalised tau_enorm = " << tau_enorm << endl;

    // Calculate or read in the resistivity
    if (options["eta"].isSet()) {
      BoutReal etafactor = options["eta"].withDefault(0.0);
      // Calculate in normalised units
      eta0 = etafactor * Me * Mi
             / (1.96 * MU0 * rhonorm * Charge * Charge * tau_enorm * rho0);
    } else {
      mesh->get(eta0, "eta0");     // Read in SI units
      eta0 *= sqrt(rhonorm / MU0); // Normalise
    }

    //////////////////////////////////////////////////////////////
    // CALCULATE METRICS

    coord->g11 = SQ(Rxy * Bpxy);
    coord->g22 = 1.0 / SQ(hthe);
    coord->g33 = SQ(I) * coord->g11 + SQ(B0) / coord->g11;
    coord->g12 = 0.0;
    coord->g13 = -I * coord->g11;
    coord->g23 = -Btxy / (hthe * Bpxy * Rxy);

    coord->J = hthe / Bpxy;
    coord->Bxy = B0;

    coord->g_11 = 1.0 / coord->g11 + SQ(I * Rxy);
    coord->g_22 = SQ(B0 * hthe / Bpxy);
    coord->g_33 = Rxy * Rxy;
    coord->g_12 = Btxy * hthe * I * Rxy / Bpxy;
    coord->g_13 = I * Rxy * Rxy;
    coord->g_23 = Btxy * hthe * Rxy / Bpxy;

    coord->geometry(); // Calculate quantities from metric tensor

    // Set B field vector
    B0vec.covariant = false;
    B0vec.x = 0.;
    B0vec.y = Bpxy / hthe;
    B0vec.z = 0.;

    vExB.setBoundary("v");
    vD.setBoundary("v");

    Jpar.setBoundary("Jpar");

    phi.setBoundary("phi");

    // Set starting dissipation terms
    eta = eta0;
    tau_e = tau_enorm * pow(Te0, 1.5) / rho0;

    output.write("\tNormalised tau_e = %e -> %e\n", min(tau_e, true), max(tau_e, true));

    // Set locations for staggered grids
    vD.setLocation(CELL_VSHIFT);

    // SET EVOLVING VARIABLES

    SOLVE_FOR(rho, Te, Ti, U, Vpar, Apar);

    comms.add(rho, Te, Ti, U, Vpar, Apar);
    comms.add(phi);

    SAVE_ONCE(P0, J0, rho0, Te0, Ti0); // Save normalised profiles

    if (nonlinear) {
      SAVE_REPEAT(eta);
    } else {
      SAVE_ONCE(eta);
    }

    SAVE_REPEAT(phi, Jpar); // Save each timestep
    SAVE_REPEAT(divExB);
    
    // Create a solver for the Laplacian
    phiSolver = Laplacian::create();
    phiSolver->setFlags(phi_flags);
    if (vorticity_momentum) {
      phiSolver->setCoefC(rho0);
    }
    return 0;
  }

  // Parallel gradient along perturbed field-line
  const Field3D Grad_parP(const Field3D& f, CELL_LOC loc = CELL_DEFAULT) {
    // Derivative along equilibrium field-line
    Field3D result;

    if (parallel_lc) {
      if (loc == CELL_YLOW) {
        result = Grad_par_CtoL(f);
      } else
        result = Grad_par_LtoC(f);
    } else
      result = Grad_par(f, loc);

    if (nonlinear) {
      if (full_bfield) {
        // Use full expression for perturbed B
        Vector3D Btilde = Curl(B0vec * Apar / B0);
        result += Btilde * Grad(f) / B0;
      } else {
        // Simplified expression
        result -= bracket(Apar, f, BRACKET_ARAKAWA);
      }
    }
    return result;
  }

  const Field3D Div_parP(const Field3D& f, CELL_LOC loc = CELL_DEFAULT) {
    return B0 * Grad_parP(f / B0, loc);
  }

  int rhs(BoutReal t) override {
    TRACE("Started physics_run(%e)", t);

    // Invert laplacian for phi
    if (vorticity_momentum) {
      // Vorticity is b dot curl(rho * v)
      Field2D rprof = rho0;
      if (nonlinear) {
        rprof += DC(rho); // Axisymmetric rho only
        phiSolver->setCoefC(rprof);
      }
      phi = phiSolver->solve(B0 * U / rprof);
    } else {
      // Vorticity is b dot curl(v)
      phi = phiSolver->solve(B0 * U);
    }
    // Apply a boundary condition on phi for target plates
    phi.applyBoundary();

    // Communicate variables
    mesh->communicate(comms);

    // Get J from Psi
    Jpar = -Delp2(Apar);
    Jpar.applyBoundary();

    if (jpar_bndry_width > 0) {
      // Boundary in jpar
      if (mesh->firstX()) {
        for (int i = jpar_bndry_width; i >= 0; i--)
          for (int j = 0; j < mesh->LocalNy; j++)
            for (int k = 0; k < mesh->LocalNz; k++) {
              Jpar(i, j, k) = 0.5 * Jpar(i + 1, j, k);
            }
      }
      if (mesh->lastX()) {
        for (int i = mesh->LocalNx - jpar_bndry_width - 1; i < mesh->LocalNx; i++)
          for (int j = 0; j < mesh->LocalNy; j++)
            for (int k = 0; k < mesh->LocalNz; k++) {
              Jpar(i, j, k) = 0.5 * Jpar(i - 1, j, k);
            }
      }
    }

    mesh->communicate(Jpar);

    // Jpar = smooth_x(Jpar); // Smooth in x direction

    Field3D rhot = rho0;
    Field3D Tet = Te0;
    Field3D Tit = Ti0;
    Field3D P = rho * (Te0 + Ti0) + (Te + Ti) * rho0; // Perturbed pressure

    if (nonlinear) {
      rhot += rho;
      Tet += Te;
      Tit += Ti;
      P += rho * (Te + Ti);

      eta = eta0 * pow(Tet / Te0, -1.5); // Update resistivity based on Te

      tau_e = tau_enorm * pow(Tet, 1.5) / rhot; // Update electron collision rate
    }

    if (flux_method) {
      {
        TRACE("Flux vExB");
        // ExB velocity
        vExB = (cross(B0vec, Grad_perp(phi))) / (B0 * B0);
        vExB.applyBoundary();
      }

      ////////// Density equation ////////////////

      {
        TRACE("Flux Density");

        // Diffusive flux (perpendicular)
        vD = -D_perp * Grad_perp(rho);
        vD.applyBoundary();

        ddt(rho) = -Div(vExB + vD, rhot);
      }

      ////////// Temperature equations ////////////

      {
        TRACE("Flux Te");

        vD = -chi_eperp * Grad_perp(Te) - Grad_par(Te, CELL_YLOW) * chi_epar * B0vec;
        vD.applyBoundary();

        ddt(Te) = -b0xGrad_dot_Grad(phi, Tet) / B0 - (2. / 3.) * Tet * Div(vExB)
                  - Div(vD, Te) / rhot;
      }

      {
        TRACE("Flux Ti");

        vD = -chi_iperp * Grad_perp(Ti) - Grad_par(Ti, CELL_YLOW) * chi_ipar * B0vec;
        vD.applyBoundary();

        ddt(Ti) = -b0xGrad_dot_Grad(phi, Tit) / B0 - (2. / 3.) * Tit * Div(vExB)
                  - Div(vD, Ti) / rhot;
      }
    } else {
      // Use analytic expressions, expand terms

      // Divergence of ExB velocity (neglecting parallel term)
      {
        TRACE("divExB");
        divExB = b0xcv * Grad(phi) / B0 - b0xGrad_dot_Grad(1. / B0, phi);
      }

      {
        TRACE("density");
        ddt(rho) = -bracket(phi, rhot, bm)              // ExB advection
                   - divExB * rhot                      // Divergence of ExB (compression)
                   - Vpar_Grad_par(Vpar, rho)           // Parallel advection
                   - rhot * Div_parP(Vpar, CELL_CENTRE) // Parallel compression
                   + D_perp * Delp2(rho)                // Perpendicular diffusion
            ;

        if (electron_density) {
          // Using electron parallel velocity rather than ion
          ddt(rho) += (Mi / (Charge * sqrt(MU0 * rhonorm))) * Div_parP(Jpar, CELL_CENTRE);
        }

        if (low_pass_z > 0)
          ddt(rho) = lowPass(ddt(rho), low_pass_z);
      }

      {
        TRACE("Te");
        ddt(Te) = -bracket(phi, Tet, bm) - Vpar_Grad_par(Vpar, Tet) // advection
                  - (2. / 3.) * Tet
                        * (divExB + Div_parP(Vpar, CELL_CENTRE)) // Divergence of flow
                  + Div_par_K_Grad_par(chi_epar, Te) / rhot      // Parallel diffusion
                  + chi_eperp * Delp2(Te) / rhot // Perpendicular diffusion
            ;

        if (ohmic_heating)
          ddt(Te) += (2. / 3) * eta * Jpar * Jpar / rhot; // Ohmic heating
      }

      {
        TRACE("Ti");
        ddt(Ti) = -bracket(phi, Tit, bm) - Vpar_Grad_par(Vpar, Tit)
                  - (2. / 3.) * Tit * (divExB + Div_parP(Vpar, CELL_CENTRE))
                  + Div_par_K_Grad_par(chi_ipar, Ti) / rhot
                  + chi_iperp * Delp2(Ti) / rhot;
      }

      if (Wei > 0.0) {
        TRACE("Wei");
        // electron-ion collision term
        // Calculate Wi * (2/3)/rho term. Wei is a scaling factor from options
        Field3D Tei = Wei * 2. * Me_Mi * (Te - Ti) / tau_e;

        ddt(Ti) += Tei;
        ddt(Te) -= Tei;
      }

      if (low_pass_z > 0) {
        ddt(Te) = lowPass(ddt(Te), low_pass_z);
        ddt(Ti) = lowPass(ddt(Ti), low_pass_z);
      }
    }

    ////////// Vorticity equation ////////////

    if (vorticity_momentum) {
      TRACE("vorticity_momentum");
      // Vorticity is b dot curl(rho * v)

      ddt(U) = SQ(B0) * Grad_parP(Jpar / B0, CELL_CENTRE) // b0 dot J
               + 2. * b0xcv * Grad(P)
               - rhot * (divExB + Div_parP(Vpar, CELL_CENTRE)) * Delp2(phi)
                     / B0 // drho/dt term
          ;

      // b dot J0
      if (full_bfield) {
        Vector3D Btilde = Curl(B0vec * Apar / B0);
        ddt(U) += B0 * Btilde * Grad(J0 / B0);
      } else {
        ddt(U) -= SQ(B0) * bracket(Apar, J0 / B0, BRACKET_ARAKAWA);
      }

      if (electron_density) {
        // drho/dt jpar term
        ddt(U) += (Mi / (Charge * sqrt(MU0 * rhonorm))) * rhot
                  * Div_parP(Jpar / rhot, CELL_CENTRE) * Delp2(phi) / B0;
      }

      if (include_profiles) {
        ddt(U) += SQ(B0) * Grad_par(J0 / B0) // b0 dot J0
                  + 2. * b0xcv * Grad(P0);
      }

      if (nonlinear) {
        ddt(U) -= bracket(phi, U, bm);    // Advection
        ddt(U) -= Vpar_Grad_par(Vpar, U); // Parallel advection
      }

      // Viscosity terms
      if (viscos_par > 0.0)
        ddt(U) += viscos_par * Grad2_par2(U); // Parallel viscosity

      if (viscos_perp > 0.0)
        ddt(U) += viscos_perp * rhot * Delp2(U / rhot); // Perpendicular viscosity

    } else {
      TRACE("vorticity");
      // Vorticity is b dot curl(v)
      ddt(U) = (SQ(B0) * Grad_parP(Jpar / B0, CELL_CENTRE)
                + 2. * b0xcv * Grad(P) // curvature term
                )
               / rhot;

      // b dot J0
      if (full_bfield) {
        Vector3D Btilde = Curl(B0vec * Apar / B0);
        ddt(U) += B0 * Btilde * Grad(J0 / B0) / rhot;
      } else {
        ddt(U) -= SQ(B0) * bracket(Apar, J0 / B0, BRACKET_ARAKAWA) / rhot;
      }

      if (include_profiles) {
        ddt(U) += (SQ(B0) * Grad_par(J0 / B0) // b0 dot J0
                   + 2. * b0xcv * Grad(P0))
                  / rhot;
      }

      if (nonlinear) {
        ddt(U) -= bracket(phi, U);        // Advection
        ddt(U) -= Vpar_Grad_par(Vpar, U); // Parallel advection
      }

      // Viscosity terms
      if (viscos_par > 0.0)
        ddt(U) += viscos_par * Grad2_par2(U) / rhot; // Parallel viscosity

      if (viscos_perp > 0.0)
        ddt(U) += viscos_perp * Delp2(U) / rhot; // Perpendicular viscosity

      // Collisional viscosity
      if (viscos_coll > 0.0)
        ddt(U) += viscos_coll / MU0 * eta * Delp2(U) / rhot;
    }

    if (low_pass_z > 0)
      ddt(U) = lowPass(ddt(U), low_pass_z);

    ////////// Parallel velocity equation ////////////

    {
      TRACE("Vpar");

      ddt(Vpar) = -Grad_parP(P + P0, CELL_YLOW);
      if (nonlinear) {
        ddt(Vpar) -= bracket(phi, Vpar);        // Advection
        ddt(Vpar) -= Vpar_Grad_par(Vpar, Vpar); // Parallel advection
      }

      if (low_pass_z > 0)
        ddt(Vpar) = lowPass(ddt(Vpar), low_pass_z);
    }

    ////////// Magnetic potential equation ////////////

    {
      TRACE("Apar");
      ddt(Apar) = -Grad_parP(phi, CELL_YLOW) - eta * Jpar;

      if (hyperresist > 0.0) {
        ddt(Apar) += eta * hyperresist * Delp2(Jpar);
      }
    }

    if (low_pass_z > 0)
      ddt(Apar) = lowPass(ddt(Apar), low_pass_z);

    return 0;
  }
};

BOUTMAIN(Jorek);
