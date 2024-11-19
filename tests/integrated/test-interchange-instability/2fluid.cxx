/*******************************************************************************
 * 2-fluid equations
 * Same as Maxim's version of BOUT - simplified 2-fluid for benchmarking
 *******************************************************************************/

#include <bout/tokamak_coordinates_factory.hxx>

#include <bout/derivs.hxx>
#include <bout/initialprofiles.hxx>
#include <bout/invert_laplace.hxx>
#include <bout/unused.hxx>

#include <cmath>
#include <cstdio>
#include <cstdlib>

class Interchange : public PhysicsModel {

  // 2D initial profiles
  Field2D Ni0, Ti0, Te0;

  // 3D evolving fields
  Field3D rho, Ni;

  // Derived 3D variables
  Field3D phi;

  // Parameters
  BoutReal Te_x, Ti_x, Ni_x, bmag, rho_s, AA, ZZ, wci;

  // Laplacian inversion
  std::unique_ptr<Laplacian> phi_solver;

  Coordinates* coord;

  TokamakCoordinatesFactory tokamak_coordinates_factory = TokamakCoordinatesFactory(*mesh);

protected:
  int init(bool UNUSED(restarting)) override {

    output << "Solving 2-variable equations\n";

    /************* LOAD DATA FROM GRID FILE ****************/

    // Load 2D profiles (set to zero if not found)
    GRID_LOAD(Ni0);
    GRID_LOAD(Ti0);
    GRID_LOAD(Te0);

    // Load normalisation values
    GRID_LOAD(Te_x);
    GRID_LOAD(Ti_x);
    GRID_LOAD(Ni_x);
    GRID_LOAD(bmag);

    Ni_x *= 1.0e14;
    bmag *= 1.0e4;

    /*************** READ OPTIONS *************************/

    // Read some parameters
    Options* globalOptions = Options::getRoot();
    Options* options = globalOptions->getSection("2fluid");
    OPTION(options, AA, 2.0);
    OPTION(options, ZZ, 1.0);

    BoutReal ShearFactor;
    OPTION(options, ShearFactor, 1.0);

    /*************** INITIALIZE LAPLACE SOLVER ***********/
    phi_solver = Laplacian::create();

    /************* SHIFTED RADIAL COORDINATES ************/

    bool noshear = false;
    const bool ShiftXderivs = (*globalOptions)["ShiftXderivs"].withDefault(false);
    if (ShiftXderivs) {
      ShearFactor = 0.0; // I disappears from metric
      noshear = true;
    }

    /************** CALCULATE PARAMETERS *****************/

    rho_s = 1.02 * sqrt(AA * Te_x) / ZZ / bmag;
    wci = 9.58e3 * ZZ * bmag / AA;

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
    Ti0 /= Te_x;
    Te0 /= Te_x;

//    b0xcv *= -1.0; // NOTE: THIS IS FOR 'OLD' GRID FILES ONLY  // TODO: Check if needed
    tokamak_coordinates_factory.normalise(rho_s, bmag / 1e4, ShearFactor);
    coord = tokamak_coordinates_factory.make_tokamak_coordinates(noshear, true);

    // Tell BOUT++ which variables to evolve
    SOLVE_FOR2(rho, Ni);

    // Add any other variables to be dumped to file
    SAVE_REPEAT(phi);
    SAVE_ONCE3(Ni0, Te0, Ti0);
    SAVE_ONCE5(Te_x, Ti_x, Ni_x, rho_s, wci);

    // Initialise aux fields
    phi = 0.;

    return (0);
  }

  int rhs(BoutReal UNUSED(t)) override {
    // Solve EM fields
    phi = phi_solver->solve(rho / Ni0, phi);

    // Communicate variables
    mesh->communicate(rho, Ni, phi);
    Field3D pei = (Te0 + Ti0) * Ni;

    // DENSITY EQUATION
    ddt(Ni) = -b0xGrad_dot_Grad(phi, Ni0) / coord->Bxy();

    // VORTICITY
    ddt(rho) = 2.0 * coord->Bxy() * tokamak_coordinates_factory.get_b0xcv() * Grad(pei);

    return (0);
  }
};

BOUTMAIN(Interchange);
