/*******************************************************************************
 * 2-fluid equations
 * Same as Maxim's version of BOUT - simplified 2-fluid for benchmarking
 *******************************************************************************/

#include <bout/physicsmodel.hxx>

#include <initialprofiles.hxx>
#include <invert_laplace.hxx>
#include <derivs.hxx>

using bout::globals::mesh;

class ShearAlfven : public PhysicsModel {
private:
  // 2D initial profiles
  Field2D Ni0, Ti0, Te0, Vi0, phi0, Ve0, rho0, Ajpar0;
  Vector2D b0xcv; // for curvature terms

  // 3D evolving fields
  Field3D rho, Ajpar;

  // Derived 3D variables
  Field3D phi, Apar, jpar;

  // pressures
  Field3D pei, pe;
  Field2D pei0, pe0;

  // Metric coefficients
  Field2D Rxy, Bpxy, Btxy, hthe;

  // parameters
  BoutReal Te_x, Ti_x, Ni_x, Vi_x, bmag, rho_s, fmei, AA, ZZ;
  BoutReal wci;
  BoutReal beta_p;

  // settings
  bool ZeroElMass;
  BoutReal zeff, nu_perp;
  BoutReal ShearFactor;

  int phi_flags, apar_flags; // Inversion flags

  // Group fields together for communication
  FieldGroup comms;

  // Coordinate system
  Coordinates *coord;

  /// Solver for inverting Laplacian
  std::unique_ptr<Laplacian> phiSolver;
  std::unique_ptr<Laplacian> aparSolver;
  
protected:
  int init(bool UNUSED(restarting)) override {
    Field2D I; // Shear factor

    output << "Solving 6-variable 2-fluid equations\n";

    /************* LOAD DATA FROM GRID FILE ****************/

    // Load 2D profiles (set to zero if not found)
    GRID_LOAD(Ni0);
    GRID_LOAD(Ti0);
    GRID_LOAD(Te0);
    GRID_LOAD(Vi0);
    GRID_LOAD(Ve0);
    GRID_LOAD(phi0);
    GRID_LOAD(rho0);
    GRID_LOAD(Ajpar0);

    // Coordinate system
    coord = mesh->getCoordinates();

    // Load magnetic curvature term
    b0xcv.covariant = false;  // Read contravariant components
    mesh->get(b0xcv, "bxcv"); // b0xkappa terms

    b0xcv *= -1.0; // NOTE: THIS IS FOR 'OLD' GRID FILES ONLY

    // Load metrics
    GRID_LOAD(Rxy);
    GRID_LOAD(Bpxy);
    GRID_LOAD(Btxy);
    GRID_LOAD(hthe);
    mesh->get(coord->Bxy, "Bxy");
    mesh->get(coord->dx, "dpsi");
    mesh->get(I, "sinty");
    
    // Load normalisation values
    GRID_LOAD(Te_x);
    GRID_LOAD(Ti_x);
    GRID_LOAD(Ni_x);
    GRID_LOAD(bmag);

    Ni_x *= 1.0e14;
    bmag *= 1.0e4;

    /*************** READ OPTIONS *************************/

    // Read some parameters
    auto globalOptions = Options::root();
    auto options = globalOptions["2fluid"];
    AA = options["AA"].withDefault(2.0);
    ZZ = options["ZZ"].withDefault(1.0);

    ZeroElMass = options["ZeroElMass"].withDefault(false);
    zeff = options["zeff"].withDefault(1.0);
    nu_perp = options["nu_perp"].withDefault(0.0);
    ShearFactor = options["ShearFactor"].withDefault(1.0);

    /************* SHIFTED RADIAL COORDINATES ************/

    // Check type of parallel transform
    std::string ptstr = Options::root()["mesh"]["paralleltransform"]["type"]
                        .withDefault<std::string>("identity");
    
    if (lowercase(ptstr) == "shifted") {
      ShearFactor = 0.0; // I disappears from metric
      b0xcv.z += I * b0xcv.x;
    }

    /************** CALCULATE PARAMETERS *****************/

    rho_s = 1.02 * sqrt(AA * Te_x) / ZZ / bmag;
    fmei = 1. / 1836.2 / AA;

    wci = 9.58e3 * ZZ * bmag / AA;

    beta_p = 4.03e-11 * Ni_x * Te_x / bmag / bmag;

    Vi_x = wci * rho_s;

    /************** PRINT Z INFORMATION ******************/

    BoutReal hthe0;
    if (mesh->get(hthe0, "hthe0") == 0) {
      output.write("    ****NOTE: input from BOUT, Z length needs to be divided by {:e}\n",
                   hthe0 / rho_s);
    }

    /************** SHIFTED GRIDS LOCATION ***************/

    // Velocities defined on cell boundaries
    Ajpar.setLocation(CELL_YLOW);

    // Apar and jpar too
    Apar.setLocation(CELL_YLOW);
    jpar.setLocation(CELL_YLOW);

    /************** NORMALISE QUANTITIES *****************/

    output.write("\tNormalising to rho_s = {:e}\n", rho_s);

    // Normalise profiles
    Ni0 /= Ni_x / 1.0e14;
    Ti0 /= Te_x;
    Te0 /= Te_x;
    phi0 /= Te_x;
    Vi0 /= Vi_x;

    // Normalise curvature term
    b0xcv.x /= (bmag / 1e4);
    b0xcv.y *= rho_s * rho_s;
    b0xcv.z *= rho_s * rho_s;

    // Normalise geometry
    Rxy /= rho_s;
    hthe /= rho_s;
    I *= rho_s * rho_s * (bmag / 1e4) * ShearFactor;
    coord->dx /= rho_s * rho_s * (bmag / 1e4);

    // Normalise magnetic field
    Bpxy /= (bmag / 1.e4);
    Btxy /= (bmag / 1.e4);
    coord->Bxy /= (bmag / 1.e4);

    // calculate pressures
    pei0 = (Ti0 + Te0) * Ni0;
    pe0 = Te0 * Ni0;

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

    SOLVE_FOR(rho, Ajpar);
    comms.add(rho, Ajpar);

    // Set boundary conditions
    jpar.setBoundary("jpar");

    /************** SETUP COMMUNICATIONS **************/

    // add extra variables to communication
    comms.add(phi);
    if (!ZeroElMass) {
      comms.add(Apar);
    }

    // Add any other variables to be dumped to file
    SAVE_REPEAT(phi, Apar, jpar);

    SAVE_ONCE(Ni0, Te0, Ti0);
    SAVE_ONCE(Te_x, Ti_x, Ni_x, rho_s, wci);
    
    // Create a solver for the Laplacian
    phiSolver = Laplacian::create(&options["phiSolver"]);

    aparSolver = Laplacian::create(&options["aparSolver"]);
    aparSolver->setCoefA((-0.5 * beta_p / fmei) * Ni0);

    return 0;
  }

  int rhs(BoutReal UNUSED(t)) override {
    // Solve EM fields

    phi = phiSolver->solve(rho / Ni0);

    if (ZeroElMass) {
      mesh->communicate(comms);
      
      Apar = -Ajpar;
      jpar = -Delp2(Apar);
      
      mesh->communicate(jpar);
    } else {

      Apar = aparSolver->solve((0.5 * beta_p / fmei) * Ni0 * Ajpar);

      // Communicate variables
      mesh->communicate(comms);

      jpar = -Ni0 * (Ajpar + Apar);
    }

    // VORTICITY
    ddt(rho) = SQ(coord->Bxy) * Div_par(jpar, CELL_CENTRE);

    // AJPAR

    ddt(Ajpar) = (1. / fmei) * Grad_par(phi, CELL_YLOW);

    return 0;
  }
};

BOUTMAIN(ShearAlfven);
