/*******************************************************************************
 * 2-fluid equations for drift-wave tests
 *
 * Settings:
 *  - ZeroElMass
 *  - AparInEpar
 *******************************************************************************/

#include <bout/physicsmodel.hxx>

#include <derivs.hxx>
#include <initialprofiles.hxx>
#include <invert_laplace.hxx>

class EMdrift : public PhysicsModel {
private:
  // 2D initial profiles
  Field2D Ni0, Ti0, Te0;

  // 3D evolving fields
  Field3D rho, Ni, Ajpar;

  // Derived 3D variables
  Field3D phi, Apar, Ve, jpar;

  // Non-linear coefficients
  Field3D nu, mu_i;

  // Metric coefficients
  Field2D Rxy, Bpxy, Btxy, hthe;

  // parameters
  BoutReal Te_x, Ti_x, Ni_x, Vi_x, bmag, rho_s, fmei, AA, ZZ;
  BoutReal lambda_ei, lambda_ii;
  BoutReal nu_hat, mui_hat, wci, nueix, nuiix;
  BoutReal beta_p;

  // settings
  bool estatic, ZeroElMass; // Switch for electrostatic operation (true = no Apar)
  bool AparInEpar;
  BoutReal zeff, nu_perp;
  bool evolve_ajpar;
  BoutReal ShearFactor;
  BoutReal nu_factor;

  int phi_flags, apar_flags; // Inversion flags

  // Communication object
  FieldGroup comms;

  // Inverts a Laplacian to get potential
  Laplacian *phiSolver;
  
  // Solves the electromagnetic potential
  Laplacian *aparSolver;
  Field2D acoef; // Coefficient in the Helmholtz equation
  
  int init(bool UNUSED(restarting)) override {
    Field2D I; // Shear factor

    output.write("Solving 6-variable 2-fluid equations\n");

    /************* LOAD DATA FROM GRID FILE ****************/

    Coordinates* coord = mesh->getCoordinates();

    // Load 2D profiles (set to zero if not found)
    mesh->get(Ni0, "Ni0");
    mesh->get(Ti0, "Ti0");
    mesh->get(Te0, "Te0");

    // Load metrics
    mesh->get(Rxy, "Rxy");
    mesh->get(Bpxy, "Bpxy");
    mesh->get(Btxy, "Btxy");
    mesh->get(hthe, "hthe");
    mesh->get(coord->dx, "dpsi");
    mesh->get(I, "sinty");

    // Load normalisation values
    mesh->get(Te_x, "Te_x");
    mesh->get(Ti_x, "Ti_x");
    mesh->get(Ni_x, "Ni_x");
    mesh->get(bmag, "bmag");

    Ni_x *= 1.0e14;
    bmag *= 1.0e4;

    /*************** READ OPTIONS *************************/

    auto globalOptions = Options::root();
    auto options = globalOptions["2fluid"];

    AA = options["AA"].withDefault(2.0);
    ZZ = options["ZZ"].withDefault(1.0);

    estatic = options["estatic"].withDefault(false);
    ZeroElMass = options["ZeroElMass"].withDefault(false);
    AparInEpar = options["AparInEpar"].withDefault(true);

    zeff = options["zeff"].withDefault(1.0);
    nu_perp = options["nu_perp"].withDefault(0.0);
    ShearFactor = options["ShearFactor"].withDefault(1.0);
    nu_factor = options["nu_factor"].withDefault(1.0);

    phi_flags = options["phi_flags"].withDefault(0);
    apar_flags = options["apar_flags"].withDefault(0);

    evolve_ajpar = globalOptions["Ajpar"]["evolve"].withDefault(true);

    if (ZeroElMass) {
      evolve_ajpar = 0; // Don't need ajpar - calculated from ohm's law
    }
    
    /************* SHIFTED RADIAL COORDINATES ************/

    // Check type of parallel transform
    std::string ptstr =  Options::root()["mesh"]["paralleltransform"].withDefault<std::string>("identity");

    if (lowercase(ptstr) == "shifted") {
      ShearFactor = 0.0; // I disappears from metric
    }

    /************** CALCULATE PARAMETERS *****************/

    rho_s = 1.02 * sqrt(AA * Te_x) / ZZ / bmag;
    fmei = 1. / 1836.2 / AA;

    lambda_ei = 24. - log(sqrt(Ni_x) / Te_x);
    lambda_ii = 23. - log(ZZ * ZZ * ZZ * sqrt(2. * Ni_x) / pow(Ti_x, 1.5));
    wci = 9.58e3 * ZZ * bmag / AA;
    nueix = 2.91e-6 * Ni_x * lambda_ei / pow(Te_x, 1.5);
    nuiix = 4.78e-8 * pow(ZZ, 4.) * Ni_x * lambda_ii / pow(Ti_x, 1.5) / sqrt(AA);
    nu_hat = nu_factor * zeff * nueix / wci;

    if (nu_perp < 1.e-10) {
      mui_hat = (3. / 10.) * nuiix / wci;
    } else
      mui_hat = nu_perp;

    if (estatic) {
      beta_p = 1.e-29;
    } else {
      beta_p = 4.03e-11 * Ni_x * Te_x / bmag / bmag;
    }

    Vi_x = wci * rho_s;

    output.write("Normalisation: rho_s = %e  wci = %e  beta_p = %e\n", rho_s, wci,
                 beta_p);

    /************** PRINT Z INFORMATION ******************/

    BoutReal hthe0;
    if (mesh->get(hthe0, "hthe0") == 0) {
      output.write("    ****NOTE: input from BOUT, Z length needs to be divided by %e\n",
                   hthe0 / rho_s);
    }

    /************** NORMALISE QUANTITIES *****************/

    output.write("\tNormalising to rho_s = %e\n", rho_s);

    // Normalise profiles
    Ni0 /= Ni_x / 1.0e14;
    Ti0 /= Te_x;
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

    /**************** SET EVOLVING VARIABLES *************/

    // Tell BOUT++ which variables to evolve
    // add evolving variables to the communication object
    SOLVE_FOR(rho);
    comms.add(rho);

    SOLVE_FOR(Ni);
    comms.add(Ni);

    if (evolve_ajpar) {
      SOLVE_FOR(Ajpar);
      comms.add(Ajpar);
      output.write("=> Evolving ajpar\n");
    } else {
      output.write("=> Not evolving Apar\n");
      initial_profile("Ajpar", Ajpar);
      if (ZeroElMass) {
        dump.add(Ajpar, "Ajpar", 1); // output calculated Ajpar
      }
    }

    jpar.setBoundary("jpar");

    /************** SETUP COMMUNICATIONS **************/

    // add extra variables to communication
    comms.add(phi, Apar);

    // Add any other variables to be dumped to file
    SAVE_REPEAT(phi, Apar, jpar);
    
    SAVE_ONCE(Ni0, Te0, Ti0);
    SAVE_ONCE(Te_x, Ti_x, Ni_x, rho_s, wci, zeff, AA);
    
    // Create a solver for the Laplacian
    phiSolver = Laplacian::create();
    phiSolver->setFlags(phi_flags);

    if (! (estatic || ZeroElMass)) {
      // Create a solver for the electromagnetic potential
      aparSolver = Laplacian::create();
      aparSolver->setFlags(apar_flags);
      acoef = (-0.5 * beta_p / fmei) * Ni0;
      aparSolver->setCoefA(acoef);
    }
    
    return 0;
  }

  // just define a macro for V_E dot Grad
#define vE_Grad(f, p) (b0xGrad_dot_Grad(p, f) / coord->Bxy)

  int rhs(BoutReal UNUSED(t)) override {

    Coordinates* coord = mesh->getCoordinates();

    // Solve EM fields
    phi = phiSolver->solve(rho / Ni0);

    if (estatic || ZeroElMass) {
      // Electrostatic operation

      Apar = 0.0;
    } else {
      Apar = aparSolver->solve(-acoef * Ajpar);
    }

    // Communicate variables
    mesh->communicate(comms);

    // Update non-linear coefficients on the mesh
    nu = nu_hat * Ni0 / pow(Te0, 1.5);
    mu_i = mui_hat * Ni0 / sqrt(Ti0);

    if (ZeroElMass) {
      // Set jpar,Ve,Ajpar neglecting the electron inertia term

      jpar = ((Te0 * Grad_par(Ni)) - (Ni0 * Grad_par(phi))) / (fmei * 0.51 * nu);

      // Set boundary conditions on jpar
      jpar.applyBoundary();

      // Need to communicate jpar
      mesh->communicate(jpar);

      Ve = -jpar / Ni0;
      Ajpar = Ve;
    } else {
      // Evolving electron parallel velocity

      if (AparInEpar) {
        // Include Apar term in Eparallel
        Ve = Ajpar + Apar;
      } else {
        Ve = Ajpar;
      }
      jpar = -Ni0 * Ve;
    }

    // DENSITY EQUATION

    ddt(Ni) = -vE_Grad(Ni0, phi);

    // VORTICITY

    ddt(rho) = SQ(coord->Bxy) * Div_par(jpar);

    // AJPAR

    ddt(Ajpar) = 0.0;
    if (evolve_ajpar) {
      ddt(Ajpar) += (1. / fmei) * Grad_par(phi);
      ddt(Ajpar) -= (1. / fmei) * (Te0 / Ni0) * Grad_par(Ni);
      ddt(Ajpar) += 0.51 * nu * jpar / Ni0;
    }

    return 0;
  }
};

BOUTMAIN(EMdrift);
