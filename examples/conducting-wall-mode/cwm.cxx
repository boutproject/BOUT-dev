/*******************************************************************************
 * Written by Brett Friedman: friedman@physics.ucla.edu
 * Linear Conducting Wall Mode Instability
 * Mode discoverd by H.L. Berk et. al. 1993
 * Model version in the code created by M. Umansky and J. Myra.
 *******************************************************************************/
#include <bout/physicsmodel.hxx>

#include <bout/derivs.hxx>
#include <bout/initialprofiles.hxx>
#include <bout/interpolation.hxx>
#include <bout/invert_laplace.hxx>

class CWM : public PhysicsModel {
private:
  // 2D initial profiles
  Field2D Ni0, Te0;

  // 3D evolving fields
  Field3D rho, te;

  // Derived 3D variables
  Field3D phi;

  // e-i Collision frequency
  Field3D nu;

  // Phi boundary conditions
  Field3D dphi_bc_ydown, dphi_bc_yup;

  // Metric coefficients
  Field2D Rxy, Bpxy, Btxy, hthe, Zxy;

  // parameters
  BoutReal Te_x, Ni_x, Vi_x, bmag, rho_s, fmei, AA, ZZ;
  BoutReal lambda_ei;
  BoutReal nu_hat, wci, nueix;

  bool bout_exb; // Use BOUT-06 expression for ExB velocity

  BoutReal zeff, nu_perp;
  BoutReal ShearFactor;

  bool filter_z;
  int filter_z_mode;

  // Coefficients for linear sheath problem
  Field2D LAMBDA1, LAMBDA2;

  // Coordinate system
  Coordinates* coord;

  // Inverts a Laplacian to get potential
  std::unique_ptr<Laplacian> phiSolver{nullptr};

  int init(bool UNUSED(restarting)) override {
    Field2D I; // Shear factor

    /************* LOAD DATA FROM GRID FILE ****************/

    // Load 2D profiles (set to zero if not found)
    GRID_LOAD(Ni0, Te0);

    coord = mesh->getCoordinates();

    // Load metrics
    GRID_LOAD(Rxy, Zxy, Bpxy, Btxy, hthe);
    mesh->get(coord->dx, "dpsi");
    mesh->get(I, "sinty");

    // Load normalisation values
    GRID_LOAD(Te_x, Ni_x, bmag);

    Ni_x *= 1.0e14;
    bmag *= 1.0e4;

    /*************** READ OPTIONS *************************/
    // Read some parameters

    auto& globalOptions = Options::root();
    auto& options = globalOptions["2fluid"];
    AA = options["AA"].withDefault(4.0);
    ZZ = options["ZZ"].withDefault(1.0);

    zeff = options["zeff"].withDefault(1.0);
    nu_perp = options["nu_perp"].withDefault(0.0);
    ShearFactor = options["ShearFactor"].withDefault(1.0);
    bout_exb = options["bout_exb"].withDefault(false);

    // Toroidal filtering
    filter_z = options["filter_z"].withDefault(false); // Filter a single n
    filter_z_mode = options["filter_z_mode"].withDefault(1);

    /************* SHIFTED RADIAL COORDINATES ************/
    // Check type of parallel transform
    std::string ptstr =
        Options::root()["mesh"]["paralleltransform"]["type"].withDefault<std::string>(
            "identity");

    if (lowercase(ptstr) == "shifted") {
      ShearFactor = 0.0; // I disappears from metric
    }

    /************** CALCULATE PARAMETERS *****************/

    rho_s = 1.02 * sqrt(AA * Te_x) / ZZ / bmag;
    fmei = 1. / 1836.2 / AA;

    lambda_ei = 24. - log(sqrt(Ni_x) / Te_x);
    wci = 9.58e3 * ZZ * bmag / AA;
    nueix = 2.91e-6 * Ni_x * lambda_ei / pow(Te_x, 1.5);
    nu_hat = zeff * nueix / wci;

    Vi_x = wci * rho_s;

    output.write("Collisions: nueix = {:e}, nu_hat = {:e}\n", nueix, nu_hat);

    /************** PRINT Z INFORMATION ******************/

    BoutReal hthe0;
    if (mesh->get(hthe0, "hthe0") == 0) {
      output.write(
          "    ****NOTE: input from BOUT, Z length needs to be divided by {:e}\n",
          hthe0 / rho_s);
    }

    /************** NORMALISE QUANTITIES *****************/

    output.write("\tNormalising to rho_s = {:e}\n", rho_s);

    // Normalise profiles
    Ni0 /= Ni_x / 1.0e14;
    Te0 /= Te_x;

    // Normalise geometry
    Rxy /= rho_s;
    hthe /= rho_s;
    I *= rho_s * rho_s * (bmag / 1e4) * ShearFactor;
    coord->dx /= rho_s * rho_s * (bmag / 1e4);

    // Normalise magnetic field
    Bpxy /= (bmag / 1.e4);
    Btxy /= (bmag / 1.e4);
    coord->Bxy /= (bmag / 1.e4);

    // Set nu
    nu = nu_hat * Ni0 / pow(Te0, 1.5);

    /**************** CALCULATE METRICS ******************/

    coord->g11 = SQ(Rxy * Bpxy);
    coord->g22 = 1.0 / SQ(hthe);
    coord->g33 = SQ(I) * coord->g11 + SQ(coord->Bxy) / coord->g11;
    coord->g12 = 0.0;
    coord->g13 = -I * coord->g11;
    coord->g23 = -Btxy / (hthe * Bpxy * Rxy);

    coord->J = hthe / Bpxy;

    coord->g_11 = 1.0 / coord->g11 + SQ(I * Rxy);
    coord->g_22 = SQ(coord->Bxy * hthe / Bpxy);
    coord->g_33 = Rxy * Rxy;
    coord->g_12 = Btxy * hthe * I * Rxy / Bpxy;
    coord->g_13 = I * Rxy * Rxy;
    coord->g_23 = Btxy * hthe * Rxy / Bpxy;

    coord->geometry();

    /**************** SET EVOLVING VARIABLES *************/

    // Tell BOUT++ which variables to evolve
    // add evolving variables to the communication object
    SOLVE_FOR(rho, te);

    SAVE_ONCE(Rxy, Bpxy, Btxy, Zxy, hthe);
    SAVE_ONCE(nu_hat, hthe0);

    // Create a solver for the Laplacian
    phiSolver = Laplacian::create();

    return 0;
  }

  //////////////////////////////////////////////////////////////////

  /// Add variables to the output. This can be used to calculate
  /// diagnostics
  ///
  /// @param[inout] state  A nested dictionary that can be added to
  void outputVars(Options& state) override {
    // Set time-varying quantity (assignRepeat)
    state["phi"].assignRepeat(phi).setAttributes({{"units", "V"},
                                                  {"conversion", Te_x},
                                                  {"standard_name", "potential"},
                                                  {"long_name", "Plasma potential"}});

    // Force updates to non-varying quantities
    state["Ni0"].force(Ni0).setAttributes(
        {{"units", "m^-3"},
         {"conversion", Ni_x},
         {"standard_name", "background ion density"},
         {"long_name", "Background ion number density"}});

    state["Te0"].force(Te0).setAttributes(
        {{"units", "eV"},
         {"conversion", Te_x},
         {"standard_name", "background electron temperature"},
         {"long_name", "Background electron temperature"}});

    state["rho_s"].force(rho_s).setAttributes(
        {{"units", "m"},
         {"conversion", 1},
         {"standard_name", "length normalisation"},
         {"long_name", "Gyro-radius length normalisation"}});

    state["wci"].force(wci).setAttributes(
        {{"units", "s^-1"},
         {"conversion", 1},
         {"standard_name", "frequency normalisation"},
         {"long_name", "Cyclotron frequency normalisation"}});

    state["Te_x"].force(Te_x).setAttributes(
        {{"units", "eV"},
         {"conversion", 1},
         {"standard_name", "temperature normalisation"},
         {"long_name", "Temperature normalisation"}});

    state["Ni_x"].force(Ni_x).setAttributes(
        {{"units", "m^-3"},
         {"conversion", 1},
         {"standard_name", "density temperature"},
         {"long_name", "Number density normalisation"}});

    state["bmag"].force(bmag).setAttributes({{"units", "G"},
                                             {"conversion", 1},
                                             {"standard_name", "B field normalisation"},
                                             {"long_name", "B field normalisation"}});

    state["AA"].force(AA).setAttributes({{"units", ""},
                                         {"conversion", 1},
                                         {"standard_name", "amu"},
                                         {"long_name", "Atomic mass number"}});

    state["ZZ"].force(ZZ).setAttributes({{"units", ""},
                                         {"conversion", 1},
                                         {"standard_name", "ion charge"},
                                         {"long_name", "Ion charge"}});

    state["zeff"].force(zeff).setAttributes({{"units", ""},
                                             {"conversion", 1},
                                             {"standard_name", "Zeff"},
                                             {"long_name", "Effective ion charge"}});
  }

  ///////////////////////////////////////////////////////////////////
  // Function called at each time step
  // Time derivatives calculated here
  int rhs(BoutReal UNUSED(t)) override {
    // Invert vorticity to get phi

    // Solves \nabla^2_\perp x + (1./c)*\nabla_perp c\cdot\nabla_\perp x + a x = b
    phi = phiSolver->solve(rho / Ni0);

    // Communicate variables
    mesh->communicate(phi, te, rho);

    // 'initial guess' for phi boundary values, before applying sheath boundary conditions
    // to set the parallel current.
    phi.applyBoundary();
    phi_sheath_bndryconds();

    // Evolve rho and te
    ddt(rho) = -Ni0 / (fmei * 0.51 * nu) * Grad2_par2(phi);

    ddt(te) = -vE_Grad(Te0, phi);

    // Z filtering
    if (filter_z) {
      // Filter out all except filter_z_mode

      ddt(rho) = filter(ddt(rho), filter_z_mode);
      ddt(te) = filter(ddt(te), filter_z_mode);
    }

    return 0;
  }

  /////////////////////////////////////////////////////////////////

  /****************BOUNDARY FUNCTIONS*****************************/
  // Sheath Boundary Conditions on Phi
  // Linearized
  void phi_sheath_bndryconds() {
    LAMBDA1 = 0.0;
    LAMBDA2 = 1.0;
    // LAMBDA1 = 1.0;
    // LAMBDA2 = log(2.0*sqrt(PI*fmei));

    dphi_bc_ydown = fmei * 0.51 * nu * (LAMBDA1 * phi + LAMBDA2 * te);
    dphi_bc_yup = -fmei * 0.51 * nu * (LAMBDA1 * phi + LAMBDA2 * te);

    bndry_ydown_Grad_par(phi, dphi_bc_ydown);
    bndry_yup_Grad_par(phi, dphi_bc_yup);
  }

  // Boundary gradient to specified Field3D object
  void bndry_yup_Grad_par(Field3D& var, const Field3D& value) {

    RangeIterator xrup = mesh->iterateBndryUpperY();

    for (xrup.first(); !xrup.isDone(); xrup.next()) {
      for (int jy = mesh->yend + 1; jy < mesh->LocalNy; jy++) {
        for (int jz = 0; jz < mesh->LocalNz; jz++) {

          var(xrup.ind, jy, jz) = var(xrup.ind, jy - 1, jz)
                                  + coord->dy(xrup.ind, jy, jz)
                                        * sqrt(coord->g_22(xrup.ind, jy, jz))
                                        * value(xrup.ind, jy, jz);
        }
      }
    }
  }

  void bndry_ydown_Grad_par(Field3D& var, const Field3D& value) {

    RangeIterator xrdn = mesh->iterateBndryLowerY();

    for (xrdn.first(); !xrdn.isDone(); xrdn.next()) {
      for (int jy = mesh->ystart - 1; jy >= 0; jy--) {
        for (int jz = 0; jz < mesh->LocalNz; jz++) {

          var(xrdn.ind, jy, jz) = var(xrdn.ind, jy + 1, jz)
                                  - coord->dy(xrdn.ind, jy, jz)
                                        * sqrt(coord->g_22(xrdn.ind, jy, jz))
                                        * value(xrdn.ind, jy, jz);
        }
      }
    }
  }

  /////////////////////////////////////////////////////////////////
  // ExB terms. These routines allow comparisons with BOUT-06
  // if bout_exb=true is set in BOUT.inp
  /////////////////////////////////////////////////////////////////
  const Field3D vE_Grad(const Field2D& f, const Field3D& p) {
    Field3D result;
    if (bout_exb) {
      // Use a subset of terms for comparison to BOUT-06
      result = VDDX(DDZ(p), f);
    } else {
      // Use full expression with all terms
      result = b0xGrad_dot_Grad(p, f) / coord->Bxy;
    }
    return result;
  }

  const Field3D vE_Grad(const Field3D& f, const Field2D& p) {
    Field3D result;
    if (bout_exb) {
      // Use a subset of terms for comparison to BOUT-06
      result = VDDZ(-DDX(p), f);
    } else {
      // Use full expression with all terms
      result = b0xGrad_dot_Grad(p, f) / coord->Bxy;
    }
    return result;
  }

  const Field3D vE_Grad(const Field3D& f, const Field3D& p) {
    Field3D result;
    if (bout_exb) {
      // Use a subset of terms for comparison to BOUT-06
      result = VDDX(DDZ(p), f) + VDDZ(-DDX(p), f);
    } else {
      // Use full expression with all terms
      result = b0xGrad_dot_Grad(p, f) / coord->Bxy;
    }
    return result;
  }
};

BOUTMAIN(CWM);
