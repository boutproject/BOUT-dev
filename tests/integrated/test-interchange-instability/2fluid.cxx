/*******************************************************************************
 * 2-fluid equations
 * Same as Maxim's version of BOUT - simplified 2-fluid for benchmarking
 *******************************************************************************/

#include <bout.hxx>

#include <derivs.hxx>
#include <initialprofiles.hxx>
#include <invert_laplace.hxx>
#include <unused.hxx>

#include <cmath>
#include <cstdio>
#include <cstdlib>

class Interchange : public PhysicsModel {

  // 2D initial profiles
  Field2D Ni0, Ti0, Te0;
  Vector2D b0xcv; // for curvature terms

  // 3D evolving fields
  Field3D rho, Ni;

  // Derived 3D variables
  Field3D phi;

  // Metric coefficients
  Field2D Rxy, Bpxy, Btxy, hthe;

  // Parameters
  BoutReal Te_x, Ti_x, Ni_x, bmag, rho_s, AA, ZZ, wci;

  // Laplacian inversion
  Laplacian* phi_solver;

  Coordinates *coord;
protected:
  int init(bool UNUSED(restarting)) override {
    Field2D I; // Shear factor

    output << "Solving 2-variable equations\n";

    /************* LOAD DATA FROM GRID FILE ****************/

    // Load 2D profiles (set to zero if not found)
    GRID_LOAD(Ni0);
    GRID_LOAD(Ti0);
    GRID_LOAD(Te0);

    // Load magnetic curvature term
    b0xcv.covariant = false;  // Read contravariant components
    mesh->get(b0xcv, "bxcv"); // b0xkappa terms

    b0xcv *= -1.0; // NOTE: THIS IS FOR 'OLD' GRID FILES ONLY

    // Coordinate system
    coord = mesh->getCoordinates();

    // Load metrics
    GRID_LOAD(Rxy);
    GRID_LOAD(Bpxy);
    GRID_LOAD(Btxy);
    GRID_LOAD(hthe);
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
    Options *globalOptions = Options::getRoot();
    Options *options = globalOptions->getSection("2fluid");
    OPTION(options, AA, 2.0);
    OPTION(options, ZZ, 1.0);

    BoutReal ShearFactor;
    OPTION(options, ShearFactor, 1.0);

    /*************** INITIALIZE LAPLACE SOLVER ***********/
    phi_solver = Laplacian::create();

    /************* SHIFTED RADIAL COORDINATES ************/
    bool ShiftXderivs;
    globalOptions->get("shiftXderivs", ShiftXderivs, false); // Read global flag
    if (ShiftXderivs) {
      ShearFactor = 0.0; // I disappears from metric
      b0xcv.z += I * b0xcv.x;
    }

    /************** CALCULATE PARAMETERS *****************/

    rho_s = 1.02 * sqrt(AA * Te_x) / ZZ / bmag;
    wci       = 9.58e3*ZZ*bmag/AA;
    
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
    ddt(Ni) = - b0xGrad_dot_Grad(phi, Ni0) / coord->Bxy;

    // VORTICITY
    ddt(rho) = 2.0 * coord->Bxy * b0xcv * Grad(pei);

    return (0);
  }
};

BOUTMAIN(Interchange);
