/*******************************************************************************
 * UEDGE benchmark case
 *
 * Solves equations for
 *  density Ni
 *  parallel ion velocity Vi
 *  electron and ion temperatures Te, Ti
 *
 * Intended to be run for NZ=1 (i.e. X and Y only) for comparison with UEDGE
 *
 *******************************************************************************/

#include <bout.hxx>
#include <bout/physicsmodel.hxx>
#include <derivs.hxx>

#include <cmath>

class UedgeBenchmark : public PhysicsModel {
private:
  // 2D initial profiles
  Field2D Ni0, Ti0, Te0, Vi0;
  
  // 3D evolving fields
  Field3D Te, Ni, Vi, Ti;
  
  // Non-linear coefficients
  Field3D kapa_Te, kapa_Ti;
  
  // 3D total values
  Field3D Nit, Tit, Tet, Vit;
  
  // pressures
  Field3D peit, pe;
  Field2D pei0, pe0;

  // Metric coefficients
  Field2D Rxy, Bpxy, Btxy, hthe;
  
  // parameters
  BoutReal Te_x, Ti_x, Ni_x, Vi_x, bmag, rho_s, fmei, AA, ZZ;
  BoutReal lambda_ei, lambda_ii;
  BoutReal nu_hat, mui_hat, wci, nueix, nuiix;
  
  BoutReal chi_perp, D_perp, mu_perp;
  
protected:
  
  int init(bool UNUSED(restarting)) {
    Field2D I; // Shear factor
    
    output.write("Solving transport equations for Ni, Vi, Ti, Te\n");
    
    /////////////// LOAD DATA FROM GRID FILE //////////////
    
    // Load 2D profiles (set to zero if not found)
    GRID_LOAD(Ni0, Ti0, Te0, Vi0);

    // Load metrics
    GRID_LOAD(Rxy);        // Major radius [m]
    GRID_LOAD(Bpxy, Btxy); // Poloidal, Toroidal B field [T]
    GRID_LOAD(hthe);       // Poloidal arc length [m / radian]
    mesh->get(mesh->getCoordinates()->dx, "dpsi");

    // Load normalisation values
    GRID_LOAD(Te_x, Ti_x, Ni_x, bmag);

    Ni_x *= 1.0e14;
    bmag *= 1.0e4;

    /////////////// READ OPTIONS //////////////////////////

    // Read some parameters
    auto& globalOptions = Options::root();
    auto& options = globalOptions["uedge"];
    AA = options["AA"].withDefault(2.0);
    ZZ = options["ZZ"].withDefault(1.0);

    chi_perp = options["chi_perp"].withDefault(0.6); // Read in m^2 / s
    D_perp = options["D_perp"].withDefault(0.6);
    mu_perp = options["mu_perp"].withDefault(0.6);

    ////////////// CALCULATE PARAMETERS ///////////////////

    rho_s = 1.02 * sqrt(AA * Te_x) / ZZ / bmag;
    fmei = 1. / 1836.2 / AA;

    lambda_ei = 24. - log(sqrt(Ni_x) / Te_x);
    lambda_ii = 23. - log(ZZ * ZZ * ZZ * sqrt(2. * Ni_x) / pow(Ti_x, 1.5));
    wci = 9.58e3 * ZZ * bmag / AA;
    nueix = 2.91e-6 * Ni_x * lambda_ei / pow(Te_x, 1.5);
    nuiix = 4.78e-8 * pow(ZZ, 4.) * Ni_x * lambda_ii / pow(Ti_x, 1.5) / sqrt(AA);

    Vi_x = wci * rho_s;

    ///////////// PRINT Z INFORMATION /////////////////////

    BoutReal hthe0;
    if (GRID_LOAD1(hthe0) == 0) {
      output.write("    ****NOTE: input from BOUT, Z length needs to be divided by %e\n",
                   hthe0 / rho_s);
    }

    ///////////// NORMALISE QUANTITIES ////////////////////
    
    output.write("\tNormalising to rho_s = %e\n", rho_s);
    
    auto* coords = mesh->getCoordinates();
    // Normalise profiles
    Ni0 /= Ni_x / 1.0e14;
    Ti0 /= Te_x;
    Te0 /= Te_x;
    Vi0 /= Vi_x;
    
    // Normalise geometry
    Rxy /= rho_s;
    hthe /= rho_s;
    coords->dx /= rho_s * rho_s * (bmag / 1e4);

    // Normalise magnetic field
    Bpxy /= (bmag / 1e4);
    Btxy /= (bmag / 1e4);
    coords->Bxy /= (bmag / 1e4);

    // calculate pressures
    pei0 = (Ti0 + Te0) * Ni0;
    pe0 = Te0 * Ni0;

    // Normalise coefficients
    chi_perp /= rho_s * rho_s * wci;
    D_perp /= rho_s * rho_s * wci;
    mu_perp /= rho_s * rho_s * wci;

    chi_perp = 0.1;
    D_perp = 0.1;
    mu_perp = 0.1;

    output.write("Diffusion coefficients: chi %e D %e Mu %e\n", chi_perp, D_perp, mu_perp);

    /////////////// CALCULATE METRICS /////////////////

    coords->g11 = pow(Rxy * Bpxy, 2.0);
    coords->g22 = 1.0 / pow(hthe, 2.0);
    coords->g33 = pow(coords->Bxy, 2.0) / coords->g11;
    coords->g12 = 0.0;
    coords->g13 = 0.0;
    coords->g23 = -Btxy / (hthe * Bpxy * Rxy);

    coords->J = hthe / Bpxy;

    coords->g_11 = 1.0 / coords->g11;
    coords->g_22 = pow(coords->Bxy * hthe / Bpxy, 2.0);
    coords->g_33 = Rxy * Rxy;
    coords->g_12 = 0.0;
    coords->g_13 = 0.0;
    coords->g_23 = Btxy * hthe * Rxy / Bpxy;

    coords->geometry(); // Calculate other metrics

    //////////////// BOUNDARIES ///////////////////////
    //
    // We want to apply the relaxing boundries to total density,
    // temperature etc.

    Nit.setBoundary("Ni");
    Tet.setBoundary("Te");
    Tit.setBoundary("Ti");
    Vit.setBoundary("Vi");
    
    Ni0.applyBoundary("neumann");
    Te0.applyBoundary("neumann");
    Ti0.applyBoundary("neumann");
    
    ///////////// SET EVOLVING VARIABLES //////////////
    //
    // Tell BOUT++ which variables to evolve

    SOLVE_FOR(Ni, Vi, Te, Ti);

    ///////////// ADD OUTPUT VARIABLES ////////////////
    //
    // Add any other variables to be dumped to file

    SAVE_ONCE(Ni0, Te0, Ti0);                // Background quantities
    SAVE_ONCE(Te_x, Ti_x, Ni_x, rho_s, wci); // Normalisation factors

    return 0;
  }

  // Operator for radial diffusive flux
  /*
    Field3D Div_X_K_Grad_X(const Field3D &difVi, const Field3D &Vi) {
    Field2D sg = 1./sqrt(mesh->g_11);
    return difVi * D2DX2(Vi)/mesh->g_11
    + DDX( difVi * sg ) * DDX(Vi) * sg;
    }
  */
  
  // This version the same as in BOUT-06. Note the R's moved,and hthe added
  Field3D Div_X_K_Grad_X(const Field3D& difFi, const Field3D& Fi) {
    Field3D result;
    
    result = difFi * (pow(Rxy * Bpxy, 2.0)) * D2DX2(Fi)
      + (Bpxy / hthe) * DDX(difFi * Rxy * Rxy * Bpxy * hthe) * DDX(Fi);
    
    return result;
  }
  
  int rhs(BoutReal UNUSED(t)) {
    // Communicate variables
    mesh->communicate(Ni, Vi, Te, Ti);

    // Update profiles
    Nit = Ni0 + DC(Ni);
    Tit = Ti0 + DC(Ti);
    Tet = Te0 + DC(Te);
    Vit = Vi0 + DC(Vi);

    // Apply boundary conditions to total fields
    Nit.applyBoundary();
    Tet.applyBoundary();
    Tit.applyBoundary();
    Vit.applyBoundary();
    
    // Update non-linear coefficients on the mesh
    kapa_Te = 3.2 * (1. / fmei) * (wci / nueix) * pow(Tet, 2.5);
    kapa_Ti = 3.9 * (wci / nuiix) * pow(Tit, 2.5);

    peit = (Tet + Tit) * Nit;

    // DENSITY EQUATION
    ddt(Ni) = -Vpar_Grad_par(Vit, Nit) - Nit * Div_par(Vit)
      + Div_X_K_Grad_X(D_perp * (Nit * 0.0 + 1.0), Nit);

    // ION VELOCITY
    ddt(Vi) = (-Grad_par(peit) + Div_X_K_Grad_X(mu_perp * Nit, Vit)) / Nit
      - Vpar_Grad_par(Vit, Nit * Vit) / Nit - ddt(Ni) * Vit / Nit;

    // ELECTRON TEMPERATURE
    ddt(Te) = (Div_par_K_Grad_par(kapa_Te, Tet) + Div_X_K_Grad_X(chi_perp * Nit, Tet))
      / (1.5 * Nit)
      - ddt(Ni) * Tet / Nit;

    // ION TEMPERATURE
    ddt(Ti) = (Div_par_K_Grad_par(kapa_Ti, Tit) + Div_X_K_Grad_X(chi_perp * Nit, Tit))
      / (1.5 * Nit)
      - ddt(Ni) * Tit / Nit;

    return 0;
  }
};

BOUTMAIN(UedgeBenchmark);
