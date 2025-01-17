/****************************************************************
 * DALF3 model
 * 
 * Four-field model for electron pressure, vorticity, A|| and
 * parallel velocity
 *
 * References:
 *
 *   B.Scott, Plasma Phys. Contr. Fusion 39 (1997) 1635
 *
 *   B.Scott, "Drift Wave versus Interchange Turbulence in
 *             Tokamak Geometry: Linear versus Nonlinear
 *             Mode Structure"
 *             arXiv:physics/0207126  Feb 2001
 *
 * NOTE: The normalisation used here is different to in the above
 *       papers. See manual in doc/ subdirectory for details
 *
 ****************************************************************/

#include <bout/physicsmodel.hxx>

#include <bout/interpolation.hxx>
#include <bout/invert/laplacexy.hxx>
#include <bout/invert_laplace.hxx>
#include <bout/utils.hxx>
#include <math.h>
#include <memory>

#include <bout/constants.hxx>

// Constants
const BoutReal MU0 = 4.0e-7 * PI;
const BoutReal Charge = SI::qe;   // electron charge e (C)
const BoutReal Mi = 2.0 * SI::Mp; // Ion mass
const BoutReal Me = SI::Me;       // Electron mass
const BoutReal Me_Mi = Me / Mi;   // Electron mass / Ion mass

class DALF3 : public PhysicsModel {
private:
  // Normalisation factors
  BoutReal Tenorm, Nenorm, Bnorm;
  BoutReal Cs, rho_s, wci;

  // Poisson brackets: b0 x Grad(f) dot Grad(g) / B = [f, g]
  // Method to use: BRACKET_ARAKAWA, BRACKET_STD or BRACKET_SIMPLE
  BRACKET_METHOD bm; // Bracket method for advection terms

  // Evolving quantities
  Field3D Vort, Ajpar, Pe, Vpar;

  Field3D phi, apar, jpar;

  Field2D B0, Pe0, Jpar0;
  Vector2D b0xcv;

  Field2D eta; // Collisional damping (resistivity)
  BoutReal beta_hat, mu_hat;
  BoutReal viscosity_par;

  bool split_n0;
  bool ZeroElMass, estatic;
  bool curv_kappa;
  bool flat_resist;
  BoutReal mul_resist;
  bool parallel_lc;
  bool nonlinear;
  bool jpar_noderiv; // Don't take Delp2(apar) to get jpar

  bool filter_z;

  BoutReal viscosity, hyper_viscosity;

  bool smooth_separatrix;

  FieldGroup comms;

  std::unique_ptr<Laplacian> phiSolver{nullptr};  // Laplacian solver in X-Z
  std::unique_ptr<Laplacian> aparSolver{nullptr}; // Laplacian solver in X-Z for Apar
  std::unique_ptr<LaplaceXY> laplacexy{nullptr};  // Laplacian solver in X-Y (n=0)
  Field2D phi2D; // Axisymmetric potential, used when split_n0=true

protected:
  int init(bool UNUSED(restarting)) override {

    /////////////////////////////////////////////////////
    // Load data from the grid

    GRID_LOAD(Jpar0);
    SAVE_ONCE(Jpar0);

    Field2D Ni0, Te0;
    GRID_LOAD2(Ni0, Te0);
    Ni0 *= 1e20;                   // To m^-3
    Pe0 = 2. * Charge * Ni0 * Te0; // Electron pressure in Pascals
    SAVE_ONCE(Pe0);

    // Load curvature term
    b0xcv.covariant = false;  // Read contravariant components
    mesh->get(b0xcv, "bxcv"); // mixed units x: T y: m^-2 z: m^-2

    // Metric coefficients
    Field2D Rxy, Bpxy, Btxy, hthe;
    Field2D I; // Shear factor

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

    //////////////////////////////////////////////////////////////
    // Options

    auto& globalOptions = Options::root();
    auto& options = globalOptions["dalf3"];

    split_n0 = options["split_n0"].withDefault(false);
    estatic = options["estatic"].withDefault(false);
    ZeroElMass = options["ZeroElMass"].withDefault(false);
    jpar_noderiv = options["jpar_noderiv"].withDefault(true);
    curv_kappa = options["curv_kappa"].withDefault(false);
    flat_resist = options["flat_resist"].withDefault(false);
    mul_resist = options["mul_resist"].withDefault(1.0);
    viscosity = options["viscosity"].withDefault(-1.0);
    hyper_viscosity = options["hyper_viscosity"].withDefault(-1.0);
    viscosity_par = options["viscosity_par"].withDefault(-1.0);
    smooth_separatrix = options["smooth_separatrix"].withDefault(false);

    filter_z = options["filter_z"].withDefault(false);

    parallel_lc = options["parallel_lc"].withDefault(true);
    nonlinear = options["nonlinear"].withDefault(true);

    const int bracket_method = options["bracket_method"].withDefault(0);
    switch (bracket_method) {
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

    Coordinates* coord = mesh->getCoordinates();

    // SHIFTED RADIAL COORDINATES

    // Check type of parallel transform
    std::string ptstr =
        Options::root()["mesh"]["paralleltransform"]["type"].withDefault("identity");

    if (lowercase(ptstr) == "shifted") {
      // Dimits style, using local coordinate system
      b0xcv.z += I * b0xcv.x;
      I = 0.0; // I disappears from metric
    }

    ///////////////////////////////////////////////////
    // Normalisation

    Tenorm = max(Te0, true);
    Nenorm = max(Ni0, true);
    Bnorm = max(B0, true);

    // Sound speed in m/s
    Cs = sqrt(Charge * Tenorm / Mi);

    // drift scale
    rho_s = Cs * Mi / (Charge * Bnorm);

    // Ion cyclotron frequency
    wci = Charge * Bnorm / Mi;

    beta_hat = 4.e-7 * PI * Charge * Tenorm * Nenorm / (Bnorm * Bnorm);

    if (ZeroElMass) {
      mu_hat = 0.;
    } else {
      mu_hat = Me / Mi;
    }

    SAVE_ONCE3(Tenorm, Nenorm, Bnorm);
    SAVE_ONCE3(Cs, rho_s, wci);
    SAVE_ONCE2(beta_hat, mu_hat);

    // Spitzer resistivity
    if (flat_resist) {
      // eta in Ohm-m. NOTE: ln(Lambda) = 20
      eta = 0.51 * 1.03e-4 * 20. * pow(Tenorm, -1.5);
    } else {
      eta = 0.51 * 1.03e-4 * 20. * pow(Te0, -1.5);
    }
    if (mul_resist < 0.0) {
      mul_resist = 0.0;
    }
    eta *= mul_resist;

    // Plasma quantities
    Jpar0 /= Nenorm * Charge * Cs;
    Pe0 /= Nenorm * Charge * Tenorm;

    // Coefficients
    eta *= Charge * Nenorm / Bnorm;

    viscosity /= wci * SQ(rho_s);
    hyper_viscosity /= wci * SQ(SQ(rho_s));
    viscosity_par /= wci * SQ(rho_s);

    b0xcv.x /= Bnorm;
    b0xcv.y *= rho_s * rho_s;
    b0xcv.z *= rho_s * rho_s;

    // Metrics
    Rxy /= rho_s;
    hthe /= rho_s;
    I *= rho_s * rho_s * Bnorm;
    Bpxy /= Bnorm;
    Btxy /= Bnorm;
    B0 /= Bnorm;

    coord->dx /= rho_s * rho_s * Bnorm;

    ///////////////////////////////////////////////////
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

    SOLVE_FOR3(Vort, Pe, Vpar);
    comms.add(Vort, Pe, Vpar);
    if (!(estatic && ZeroElMass)) {
      SOLVE_FOR(Ajpar);
      // Never differentiate Ajpar -> don't communicate
    }
    if (estatic) {
      comms.add(jpar);
    } else {
      // Need to communicate apar first then jpar
      comms.add(apar);
    }

    comms.add(phi);

    phi.setBoundary("phi");
    apar.setBoundary("apar");
    jpar.setBoundary("jpar");

    SAVE_REPEAT3(jpar, apar, phi);

    if (nonlinear) {
      SAVE_REPEAT(eta);
    } else {
      SAVE_ONCE(eta);
    }

    // Create a solver for the Laplacian
    phiSolver = Laplacian::create(&options["phiSolver"]);

    // LaplaceXY for n=0 solve
    if (split_n0) {
      // Create an XY solver for n=0 component
      laplacexy = std::make_unique<LaplaceXY>(mesh);
      phi2D = 0.0; // Starting guess
    }

    // Solver for Apar
    // ajpar = beta_hat*apar + mu_hat*jpar
    aparSolver = Laplacian::create(&options["aparSolver"]);
    aparSolver->setCoefA(beta_hat);
    aparSolver->setCoefD(-mu_hat);

    return 0;
  }

  // Curvature operator
  Field3D Kappa(const Field3D& f) {
    if (curv_kappa) {
      // Use the b0xcv vector from grid file
      return -2. * b0xcv * Grad(f) / B0;
    }

    return 2. * bracket(log(B0), f, bm);
  }

  Field3D Grad_parP(const Field3D& f) {
    if (nonlinear) {
      return ::Grad_parP(apar * beta_hat, f);
    }
    return Grad_par(f);
  }

  int rhs(BoutReal UNUSED(time)) override {

    // Invert vorticity to get electrostatic potential
    if (split_n0) {
      Field2D Vort2D = DC(Vort); // n=0 component
      phi2D = laplacexy->solve(Vort2D, phi2D);

      // Solve non-axisymmetric part using X-Z solver
      phi = phiSolver->solve((Vort - Vort2D) * B0);
      phi += phi2D; // Add axisymmetric part
    } else {
      phi = phiSolver->solve(Vort * B0);
    }
    phi.applyBoundary();

    // Calculate apar and jpar
    if (estatic) {
      // Electrostatic
      apar = 0.;
      if (ZeroElMass) {
        // Not evolving Ajpar
        jpar = Grad_par(Pe - phi) / eta;
        jpar.applyBoundary();
      } else {
        jpar = Ajpar / mu_hat;
      }
      mesh->communicate(comms);
    } else {
      // Electromagnetic
      if (ZeroElMass) {
        // Ignore electron inertia term
        apar = Ajpar / beta_hat;

        mesh->communicate(comms);
        jpar = -Delp2(apar);
        jpar.applyBoundary();
        mesh->communicate(jpar);
      } else {
        // All terms - solve Helmholtz equation
        // ajpar = beta_hat*apar + mu_hat*jpar
        apar = aparSolver->solve(Ajpar);
        apar.applyBoundary();

        mesh->communicate(comms);
        if (jpar_noderiv) {
          // Already applied boundaries on Ajpar and apar
          jpar = (Ajpar - beta_hat * apar) / mu_hat;
        } else {
          jpar = -Delp2(apar);
          jpar.applyBoundary();
          mesh->communicate(jpar);
        }
      }
    }

    Field3D Pet = Pe0;
    if (nonlinear) {
      Pet += Pe;
    }

    // Boundary in jpar
    if (mesh->firstX()) {
      for (int i = 4; i >= 0; i--) {
        for (int j = 0; j < mesh->LocalNy; j++) {
          for (int k = 0; k < mesh->LocalNz; k++) {
            jpar(i, j, k) = 0.5 * jpar(i + 1, j, k);
          }
        }
      }
    }
    if (mesh->lastX()) {
      for (int i = mesh->LocalNx - 5; i < mesh->LocalNx; i++) {
        for (int j = 0; j < mesh->LocalNy; j++) {
          for (int k = 0; k < mesh->LocalNz; k++) {
            jpar(i, j, k) = 0.5 * jpar(i - 1, j, k);
          }
        }
      }
    }

    // Vorticity equation
    ddt(Vort) = B0 * B0 * Grad_parP(jpar / B0) - B0 * Kappa(Pe);

    if (nonlinear) {
      ddt(Vort) -= bracket(phi, Vort, bm); // ExB advection
    }

    if (viscosity > 0.0) {
      ddt(Vort) += viscosity * Delp2(Vort);
    }
    if (hyper_viscosity > 0.0) {
      Field3D delp2_vort = Delp2(Vort);
      delp2_vort.applyBoundary("neumann");
      mesh->communicate(delp2_vort);

      ddt(Vort) += hyper_viscosity * Delp2(delp2_vort);
    }

    if (filter_z) {
      ddt(Vort) = filter(ddt(Vort), 1);
    }

    // Parallel Ohm's law
    if (!(estatic && ZeroElMass)) {
      // beta_hat*apar + mu_hat*jpar
      ddt(Ajpar) = Grad_parP(Pe - phi) - beta_hat * bracket(apar, Pe0, BRACKET_ARAKAWA)
                   - eta * jpar;

      if (nonlinear) {
        ddt(Ajpar) -= mu_hat * bracket(phi, jpar, bm);
      }

      if (filter_z) {
        ddt(Ajpar) = filter(ddt(Ajpar), 1);
      }
    }

    // Parallel velocity
    ddt(Vpar) = -Grad_parP(Pe) + beta_hat * bracket(apar, Pe0);

    if (nonlinear) {
      ddt(Vpar) -= bracket(phi, Vpar, bm);
    }

    if (viscosity_par > 0.) {
      ddt(Vpar) += viscosity_par * Grad2_par2(Vpar);
    }

    if (filter_z) {
      ddt(Vpar) = filter(ddt(Vpar), 1);
    }

    // Electron pressure
    ddt(Pe) = -bracket(phi, Pet, bm)
              + Pet * (Kappa(phi - Pe) + B0 * Grad_parP(jpar - Vpar) / B0);

    if (smooth_separatrix) {
      // Experimental smoothing across separatrix
      ddt(Vort) += mesh->smoothSeparatrix(Vort);
    }

    if (filter_z) {
      ddt(Pe) = filter(ddt(Pe), 1);
    }

    // Boundary in Vpar and vorticity

    if (mesh->firstX()) {
      for (int i = 3; i >= 0; i--) {
        for (int j = 0; j < mesh->LocalNy; j++) {
          for (int k = 0; k < mesh->LocalNz; k++) {
            ddt(Vpar)(i, j, k) = ddt(Vpar)(i + 1, j, k);
            ddt(Vort)(i, j, k) = ddt(Vort)(i + 1, j, k);
          }
        }
      }

      // Subtract DC component
      for (int i = 0; i < 10; i++) {
        for (int j = 0; j < mesh->LocalNy; j++) {
          BoutReal avg = 0.;
          for (int k = 0; k < mesh->LocalNz; k++) {
            avg += ddt(Vort)(i, j, k);
          }
          avg /= (BoutReal)mesh->LocalNz;
          for (int k = 0; k < mesh->LocalNz; k++) {
            ddt(Vort)(i, j, k) -= avg;
          }
        }
      }
    }
    if (mesh->lastX()) {
      for (int i = mesh->LocalNx - 3; i < mesh->LocalNx; i++) {
        for (int j = 0; j < mesh->LocalNy; j++) {
          for (int k = 0; k < mesh->LocalNz; k++) {
            ddt(Vpar)(i, j, k) = ddt(Vpar)(i - 1, j, k);
            ddt(Vort)(i, j, k) = ddt(Vort)(i - 1, j, k);
          }
        }
      }
    }

    return 0;
  }
};

BOUTMAIN(DALF3);
