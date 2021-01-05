/*******************************************************************************
 * High-Beta Flute-Reduced MHD
 * see elm_reduced.pdf
 * Basically the same as Hazeltine-Meiss but different normalisations.
 * Can also include the Vpar compressional term
 *
 * July 2014: Ben Dudson <benjamin.dudson@york.ac.uk>
 *      o Simplified and cleaned up for MMS testing
 *
 *******************************************************************************/

#include <bout.hxx>
#include <initialprofiles.hxx>
#include <invert_laplace.hxx>
#include <invert_parderiv.hxx>
#include <interpolation.hxx>
#include <derivs.hxx>
#include <sourcex.hxx>
#include <bout/physicsmodel.hxx>
#include <bout/constants.hxx>
#include <msg_stack.hxx>
#include <utils.hxx>

#include <field_factory.hxx>

#include <math.h>

class ELMpb : public PhysicsModel {
private:
  // 2D inital profiles
  Field2D J0, P0; // Current and pressure
  Vector2D b0xcv; // Curvature term
  Field2D beta, gradparB;   // Used for Vpar terms
  Field2D phi0;   // When diamagnetic terms used

  // B field vectors
  Vector2D B0vec; // B0 field vector

  // 3D evolving variables
  Field3D U, Psi, P, Vpar;

  // Derived 3D variables
  Field3D Jpar, phi; // Parallel current, electric potential

  Field3D Jpar2; //  Delp2 of Parallel current

  Field3D PsiExact;

  // Parameters
  BoutReal density; // Number density [m^-3]
  BoutReal Bbar, Lbar, Tbar, Va; // Normalisation constants
  BoutReal dnorm; // For diamagnetic terms: 1 / (2. * wci * Tbar)
  BoutReal dia_fact; // Multiply diamagnetic term by this
  BoutReal delta_i; // Normalized ion skin depth
  BoutReal omega_i; // ion gyrofrequency

  BoutReal diffusion_par;  // Parallel pressure diffusion

  BoutReal viscos_perp; // Perpendicular viscosity

  // options
  bool include_curvature, include_jpar0, compress0;

  BoutReal vacuum_pressure;
  BoutReal vacuum_trans; // Transition width
  Field3D vac_mask;

  std::unique_ptr<Laplacian> phi_solver{nullptr};

  bool nonlinear;
  BoutReal g; // Only if compressible
  bool phi_curv;

  // Poisson brackets: b0 x Grad(f) dot Grad(g) / B = [f, g]
  // Method to use: BRACKET_ARAKAWA, BRACKET_STD or BRACKET_SIMPLE
  /*
   * Bracket method
   *
   * BRACKET_STD      - Same as b0xGrad_dot_Grad, methods in BOUT.inp
   * BRACKET_SIMPLE   - Subset of terms, used in BOUT-06
   * BRACKET_ARAKAWA  - Arakawa central differencing (2nd order)
   * BRACKET_CTU      - 1st order upwind method
   *
   */

  // Bracket method for advection terms
  BRACKET_METHOD bm_exb;
  BRACKET_METHOD bm_mag;
  int bm_exb_flag;
  int bm_mag_flag;
  /* BRACKET_METHOD bm_ExB = BRACKET_STD;
     BRACKET_METHOD bm_mflutter = BRACKET_STD; */

  bool diamag;
  bool diamag_grad_t; // Grad_par(Te) term in Psi equation
  bool diamag_phi0;   // Include the diamagnetic equilibrium phi0

  bool eHall;
  BoutReal AA; // ion mass in units of the proton mass; AA=Mi/Mp

  bool filter_z;
  int filter_z_mode;
  int low_pass_z;
  bool zonal_flow;
  bool zonal_field;
  bool zonal_bkgd;

  BoutReal vac_lund, core_lund;       // Lundquist number S = (Tau_R / Tau_A). -ve -> infty
  BoutReal vac_resist,  core_resist;  // The resistivities (just 1 / S)
  Field3D eta;                    // Resistivity profile (1 / S)
  bool spitzer_resist;  // Use Spitzer formula for resistivity
  BoutReal Zeff;            // Z effective for resistivity formula

  BoutReal hyperresist;    // Hyper-resistivity coefficient (in core only)
  BoutReal ehyperviscos;   // electron Hyper-viscosity coefficient

  // Metric coefficients
  Field2D Rxy, Bpxy, Btxy, B0, hthe;
  Field2D I; // Shear factor

  bool mms; // True if testing with Method of Manufactured Solutions

  const BoutReal MU0 = 4.0e-7*PI;
  const BoutReal Mi = 2.0*1.6726e-27; // Ion mass

  // Communication objects
  FieldGroup comms;

  // Parallel gradient along perturbed field-line
  const Field3D Grad_parP(const Field3D &f, CELL_LOC loc = CELL_DEFAULT) {
    Field3D result;

    result = Grad_par(f, loc);

    if(nonlinear) {
      result -= bracket(Psi, f, bm_mag)*B0;
    }

    return result;
  }

public:
  int init(bool restarting) {
    output.write("Solving high-beta flute reduced equations\n");
    output.write("\tFile    : {:s}\n", __FILE__);
    output.write("\tCompiled: {:s} at {:s}\n", __DATE__, __TIME__);

    //////////////////////////////////////////////////////////////
    // Load data from the grid

    // Load 2D profiles
    mesh->get(J0, "Jpar0");    // A / m^2
    mesh->get(P0, "pressure"); // Pascals

    // Load curvature term
    b0xcv.covariant = false; // Read contravariant components
    mesh->get(b0xcv, "bxcv"); // mixed units x: T y: m^-2 z: m^-2

    // Load metrics
    if(mesh->get(Rxy,  "Rxy")) { // m
      output_error.write("Error: Cannot read Rxy from grid\n");
      return 1;
    }
    if(mesh->get(Bpxy, "Bpxy")) { // T
      output_error.write("Error: Cannot read Bpxy from grid\n");
      return 1;
    }
    mesh->get(Btxy, "Btxy"); // T
    mesh->get(B0,   "Bxy");  // T
    mesh->get(hthe, "hthe"); // m
    mesh->get(I,    "sinty");// m^-2 T^-1

    Coordinates *coords = mesh->getCoordinates();

    //////////////////////////////////////////////////////////////
    // Read parameters from the options file
    //
    // Options.get ( NAME,    VARIABLE,    DEFAULT VALUE)
    //
    // or if NAME = "VARIABLE" then just
    //
    // OPTION(VARIABLE, DEFAULT VALUE)
    //
    // Prints out what values are assigned
    /////////////////////////////////////////////////////////////

    Options *globalOptions = Options::getRoot();
    Options *options = globalOptions->getSection("highbeta");

    OPTION(options, density,           1.0e19); // Number density [m^-3]

    // Effects to include/exclude
    OPTION(options, include_curvature, true);
    OPTION(options, include_jpar0,     true);

    OPTION(options, compress0,          false);
    OPTION(options, nonlinear,         false);

    // option for ExB Poisson Bracket
    OPTION(options, bm_exb_flag,         0);
    switch(bm_exb_flag) {
    case 0: {
      bm_exb = BRACKET_STD;
      output << "\tBrackets for ExB: default differencing\n";
      break;
    }
    case 1: {
      bm_exb = BRACKET_SIMPLE;
      output << "\tBrackets for ExB: simplified operator\n";
      break;
    }
    case 2: {
      bm_exb = BRACKET_ARAKAWA;
      output << "\tBrackets for ExB: Arakawa scheme\n";
      break;
    }
    case 3: {
      bm_exb = BRACKET_CTU;
      output << "\tBrackets for ExB: Corner Transport Upwind method\n";
      break;
    }
    default:
      output << "ERROR: Invalid choice of bracket method. Must be 0 - 3\n";
      return 1;
    }

    // option for magnetic flutter Poisson Bracket
    OPTION(options, bm_mag_flag,         0);
    switch(bm_mag_flag) {
    case 0: {
      bm_mag = BRACKET_STD;
      output << "\tBrackets: default differencing\n";
      break;
    }
    case 1: {
      bm_mag = BRACKET_SIMPLE;
      output << "\tBrackets: simplified operator\n";
      break;
    }
    case 2: {
      bm_mag = BRACKET_ARAKAWA;
      output << "\tBrackets: Arakawa scheme\n";
      break;
    }
    case 3: {
      bm_mag = BRACKET_CTU;
      output << "\tBrackets: Corner Transport Upwind method\n";
      break;
    }
    default:
      output << "ERROR: Invalid choice of bracket method. Must be 0 - 3\n";
      return 1;
    }

    OPTION(options, eHall,            false);  // electron Hall or electron parallel pressue gradient effects?
    OPTION(options, AA,               1.0);    // ion mass in units of proton mass

    OPTION(options, diamag,            false);  // Diamagnetic effects?
    OPTION(options, diamag_grad_t,     diamag); // Grad_par(Te) term in Psi equation
    OPTION(options, diamag_phi0,       diamag); // Include equilibrium phi0
    OPTION(options, dia_fact,          1.0);    // Scale diamagnetic effects by this factor

    // Toroidal filtering
    OPTION(options, filter_z,      false);      // Filter a single n
    OPTION(options, filter_z_mode,     1);
    OPTION(options, low_pass_z,       -1);      // Low-pass filter
    OPTION(options, zonal_flow,    false);      // zonal flow filter
    OPTION(options, zonal_field,   false);      // zonal field filter
    OPTION(options, zonal_bkgd,    false);      // zonal background P filter

    // Vacuum region control
    OPTION(options, vacuum_pressure,   0.02);   // Fraction of peak pressure
    OPTION(options, vacuum_trans,      0.005);  // Transition width in pressure

    // Resistivity and hyper-resistivity options
    OPTION(options, vac_lund,          0.0);    // Lundquist number in vacuum region
    OPTION(options, core_lund,         0.0);    // Lundquist number in core region
    OPTION(options, hyperresist,       -1.0);
    OPTION(options, ehyperviscos,      -1.0);
    OPTION(options, spitzer_resist,    false);  // Use Spitzer resistivity
    OPTION(options, Zeff,              2.0);    // Z effective

    // Viscosity and hyper-viscosity
    OPTION(options, viscos_perp,       -1.0);  // Perpendicular viscosity

    // parallel pressure diffusion
    OPTION(options, diffusion_par,        -1.0);  // Parallel pressure diffusion

    // Compressional terms
    OPTION(options, phi_curv,          true);
    options->get("gamma",             g,                 5.0/3.0);

    OPTION(globalOptions->getSection("solver"), mms, false);

    if(!include_curvature)
      b0xcv = 0.0;

    if(!include_jpar0)
      J0 = 0.0;

    //////////////////////////////////////////////////////////////
    // INITIALIZE LAPLACIAN SOLVER

    phi_solver = Laplacian::create();

    //////////////////////////////////////////////////////////////
    // SHIFTED RADIAL COORDINATES

    bool ShiftXderivs;
    globalOptions->get("shiftXderivs", ShiftXderivs, false); // Read global flag
    if(ShiftXderivs) {
      if(mesh->IncIntShear) {
        // BOUT-06 style, using d/dx = d/dpsi + I * d/dz
        coords->IntShiftTorsion = I;

      }else {
        // Dimits style, using local coordinate system
        if(include_curvature)
      b0xcv.z += I*b0xcv.x;
        I = 0.0;  // I disappears from metric
      }
    }

    //////////////////////////////////////////////////////////////
    // NORMALISE QUANTITIES

    if(mesh->get(Bbar, "bmag")) // Typical magnetic field
      Bbar = 1.0;
    if(mesh->get(Lbar, "rmag")) // Typical length scale
      Lbar = 1.0;

    Va = sqrt(Bbar*Bbar / (MU0*density*Mi));

    Tbar = Lbar / Va;

    dnorm = dia_fact * Mi / (2.*1.602e-19*Bbar*Tbar);

    delta_i = AA*60.67*5.31e5/sqrt(density/1e6)/(Lbar*100.0);

    output.write("Normalisations: Bbar = {:e} T   Lbar = {:e} m\n", Bbar, Lbar);
    output.write("                Va = {:e} m/s   Tbar = {:e} s\n", Va, Tbar);
    output.write("                dnorm = {:e}\n", dnorm);
    output.write("    Resistivity\n");

    if(eHall)
      output.write("                delta_i = {:e}   AA = {:e} \n", delta_i, AA);

    if(vac_lund > 0.0) {
      output.write("        Vacuum  Tau_R = {:e} s   eta = {:e} Ohm m\n", vac_lund * Tbar,
           MU0 * Lbar * Lbar / (vac_lund * Tbar));
      vac_resist = 1. / vac_lund;
    }else {
      output.write("        Vacuum  - Zero resistivity -\n");
      vac_resist = 0.0;
    }

    if(core_lund > 0.0) {
      output.write("        Core    Tau_R = {:e} s   eta = {:e} Ohm m\n", core_lund * Tbar,
           MU0 * Lbar * Lbar / (core_lund * Tbar));
      core_resist = 1. / core_lund;
    }else {
      output.write("        Core    - Zero resistivity -\n");
      core_resist = 0.0;
    }

    if(ehyperviscos > 0.0) {
      output.write("    electron Hyper-viscosity coefficient: {:e}\n", ehyperviscos);
    }

    Field2D Te;
    Te = P0 / (2.0*density * 1.602e-19); // Temperature in eV

    J0 = - MU0*Lbar * J0 / B0;
    P0 = 2.0*MU0 * P0 / (Bbar*Bbar);

    b0xcv.x /= Bbar;
    b0xcv.y *= Lbar*Lbar;
    b0xcv.z *= Lbar*Lbar;

    Rxy  /= Lbar;
    Bpxy /= Bbar;
    Btxy /= Bbar;
    B0   /= Bbar;
    hthe /= Lbar;
    coords->dx   /= Lbar*Lbar*Bbar;
    I    *= Lbar*Lbar*Bbar;

    BoutReal pnorm = max(P0, true); // Maximum over all processors

    vacuum_pressure *= pnorm; // Get pressure from fraction
    vacuum_trans *= pnorm;

    // Transitions from 0 in core to 1 in vacuum
    vac_mask = (1.0 - tanh( (P0 - vacuum_pressure) / vacuum_trans )) / 2.0;

    if(spitzer_resist) {
      // Use Spitzer resistivity
      output.write("\tTemperature: {:e} -> {:e} [eV]\n", min(Te), max(Te));
      eta = 0.51*1.03e-4*Zeff*20.*pow(Te, -1.5); // eta in Ohm-m. NOTE: ln(Lambda) = 20
      output.write("\tSpitzer resistivity: {:e} -> {:e} [Ohm m]\n", min(eta), max(eta));
      eta /= MU0 * Va * Lbar;
      output.write("\t -> Lundquist {:e} -> {:e}\n", 1.0/max(eta), 1.0/min(eta));
    }else {
      // transition from 0 for large P0 to resistivity for small P0
      eta = core_resist + (vac_resist - core_resist) * vac_mask;
    }

    dump.add(eta, "eta", 0);

    /**************** CALCULATE METRICS ******************/

    coords->g11 = SQ(Rxy*Bpxy);
    coords->g22 = 1.0 / SQ(hthe);
    coords->g33 = SQ(I)*coords->g11 + SQ(B0)/coords->g11;
    coords->g12 = 0.0;
    coords->g13 = -I*coords->g11;
    coords->g23 = -Btxy/(hthe*Bpxy*Rxy);

    coords->J = hthe / Bpxy;
    coords->Bxy = B0;

    coords->g_11 = 1.0/coords->g11 + (SQ(I*Rxy));
    coords->g_22 = SQ(B0*hthe/Bpxy);
    coords->g_33 = Rxy*Rxy;
    coords->g_12 = Btxy*hthe*I*Rxy/Bpxy;
    coords->g_13 = I*Rxy*Rxy;
    coords->g_23 = Btxy*hthe*Rxy/Bpxy;

    coords->geometry(); // Calculate quantities from metric tensor

    // Set B field vector

    B0vec.covariant = false;
    B0vec.x = 0.;
    B0vec.y = Bpxy / hthe;
    B0vec.z = 0.;

    /**************** SET VARIABLE LOCATIONS *************/

    P.setLocation(CELL_CENTRE);
    U.setLocation(CELL_CENTRE);
    phi.setLocation(CELL_CENTRE);
    Psi.setLocation(CELL_YLOW);
    Jpar.setLocation(CELL_YLOW);
    Vpar.setLocation(CELL_YLOW);

    /**************** SET EVOLVING VARIABLES *************/

    // Tell BOUT which variables to evolve
    SOLVE_FOR3(U, P, Psi);
    dump.add(Jpar, "jpar", 1);

    if(compress0) {
      output.write("Including compression (Vpar) effects\n");

      SOLVE_FOR(Vpar);

      beta = B0*B0 / ( 0.5 + (B0*B0 / (g*P0)));
      gradparB = Grad_par(B0) / B0;

      output.write("Beta in range {:e} -> {:e}\n",
                   min(beta), max(beta));
    }

    // Phi solved in RHS (explicitly)
    dump.add(phi, "phi", 1);

    // Diamagnetic phi0
    if(diamag_phi0) {
      phi0 = -0.5*dnorm*P0/B0;
      SAVE_ONCE(phi0);
    }

    // Add some equilibrium quantities and normalisations
    // everything needed to recover physical units
    SAVE_ONCE2(J0, P0);
    SAVE_ONCE4(density, Lbar, Bbar, Tbar);
    SAVE_ONCE2(Va, B0);

    /////////////// CHECK VACUUM ///////////////////////
    // In vacuum region, initial vorticity should equal zero

    if(!restarting) {
      // Only if not restarting: Check initial perturbation

      // Set U to zero where P0 < vacuum_pressure
      U = where(P0 - vacuum_pressure, U, 0.0);

      // Phi should be consistent with U
      phi = phi_solver->solve(U);

      //if(diamag) {
      //phi -= 0.5*dnorm * P / B0;
      //}
    }

    /************** SETUP COMMUNICATIONS **************/

    comms.add(U, P, Psi);

    phi.setBoundary("phi"); // Set boundary conditions
    Jpar.setBoundary("J");
    Jpar2.setBoundary("J");

    PsiExact.setBoundary("Psi");

    return 0;
  }

  int rhs(BoutReal t) {
    // Perform communications
    mesh->communicate(comms);

    Coordinates *coords = mesh->getCoordinates();

    ////////////////////////////////////////////
    // Transitions from 0 in core to 1 in vacuum
    if(nonlinear) {
      vac_mask = (1.0 - tanh( ((P0 + P) - vacuum_pressure) / vacuum_trans )) / 2.0;

      // Update resistivity
      if(spitzer_resist) {
        // Use Spitzer formula
        Field3D Te;
        Te = (P0+P)*Bbar*Bbar/(4.*MU0) / (density * 1.602e-19); // eV
        eta = 0.51*1.03e-4*Zeff*20.*pow(Te, -1.5); // eta in Ohm-m. ln(Lambda) = 20
        eta /= MU0 * Va * Lbar; // Normalised eta
      }else {
        // Use specified core and vacuum Lundquist numbers
        eta = core_resist + (vac_resist - core_resist) * vac_mask;
      }
    }

    ////////////////////////////////////////////
    // Inversion
    if(mms) {
      // Solve for potential, adding a source term
      Field3D phiS = FieldFactory::get()->create3D("phi:source", Options::getRoot(), mesh, CELL_CENTRE, t);
      phi = phi_solver->solve(U + phiS, phi);
    }else {
      phi = phi_solver->solve(U, phi);
    }

    if(diamag) {
      phi -= 0.5*dnorm * P / B0;
    }

    // Apply a boundary condition on phi for target plates
    //phi.applyBoundary();
    mesh->communicate(phi);

    // Get J from Psi
    //PsiExact = FieldFactory::get()->create3D("psi:solution", Options::getRoot(), mesh, CELL_CENTRE, t);
    //PsiExact.applyBoundary(t);
    //Jpar = Delp2(PsiExact);
    Jpar = Delp2(Psi);

    Jpar.applyBoundary(t);
    mesh->communicate(Jpar);

    //Jpar = FieldFactory::get()->create3D("j:solution", Options::getRoot(), mesh, CELL_CENTRE, t);
    //Jpar.applyBoundary(t);

    //phi = FieldFactory::get()->create3D("phi:solution", Options::getRoot(), mesh, CELL_CENTRE, t);

    // Get Delp2(J) from J
    Jpar2 = Delp2(Jpar);

    Jpar2.applyBoundary(t);
    mesh->communicate(Jpar2);

    ////////////////////////////////////////////////////
    // Parallel electric field
    // Evolving vector potential

    ddt(Psi) = -Grad_parP(phi, CELL_CENTRE) + eta*Jpar;

    if(eHall) {
      ddt(Psi) +=  0.25*delta_i*(Grad_parP(B0*P, CELL_CENTRE) / B0
                                 +b0xGrad_dot_Grad(P0, Psi));   // electron parallel pressure
    }

    if(diamag_phi0)
      ddt(Psi) -= b0xGrad_dot_Grad(phi0, Psi);   // Equilibrium flow

    if(diamag_grad_t) {
      // grad_par(T_e) correction

      ddt(Psi) += 1.71 * dnorm * 0.5 * Grad_parP(P, CELL_YLOW) / B0;
    }

    // Hyper-resistivity
    if(hyperresist > 0.0) {
      ddt(Psi) -= eta*hyperresist * Delp2(Jpar);
    }

    // electron Hyper-viscosity coefficient
    if(ehyperviscos > 0.0) {
      ddt(Psi) -= eta*ehyperviscos * Delp2(Jpar2);
    }

    ////////////////////////////////////////////////////
    // Vorticity equation

    ddt(U) = SQ(B0) * b0xGrad_dot_Grad(Psi, J0, CELL_CENTRE); // Grad j term

    ddt(U) += b0xcv*Grad(P);  // curvature term

    // Parallel current term
    ddt(U) -= SQ(B0)*Grad_parP(Jpar, CELL_CENTRE); // b dot grad j

    if(diamag_phi0)
      ddt(U) -= b0xGrad_dot_Grad(phi0, U);   // Equilibrium flow

    if(nonlinear) {
      ddt(U) -= bracket(phi, U, bm_exb)*B0;    // Advection
    }

    // Viscosity terms

    if(viscos_perp > 0.0)
      ddt(U) += viscos_perp * Delp2(U);     // Perpendicular viscosity

    ddt(U) -= 10*(SQ(SQ(coords->dx))*D4DX4(U) + SQ(SQ(coords->dz))*D4DZ4(U));

    ////////////////////////////////////////////////////
    // Pressure equation

    ddt(P) = -b0xGrad_dot_Grad(phi, P0);


    if(diamag_phi0)
      ddt(P) -= b0xGrad_dot_Grad(phi0, P);   // Equilibrium flow

    if(nonlinear)
      ddt(P) -= bracket(phi, P, bm_exb)*B0;    // Advection

    // Parallel diffusion terms
    if(diffusion_par > 0.0)
      ddt(P) += diffusion_par * Grad2_par2(P); // Parallel diffusion

    ddt(P) -= 10*(SQ(SQ(coords->dx))*D4DX4(P) + SQ(SQ(coords->dz))*D4DZ4(P));

    ////////////////////////////////////////////////////
    // Compressional effects

    if(compress0) {

      //ddt(P) += beta*( - Grad_parP(Vpar, CELL_CENTRE) + Vpar*gradparB );
      ddt(P) -= beta*Div_par_CtoL(Vpar);

      if(phi_curv) {
        ddt(P) -= 2.*beta*b0xcv*Grad(phi);
      }

      // Vpar equation

      //ddt(Vpar) = -0.5*Grad_parP(P + P0, CELL_YLOW);
      ddt(Vpar) = -0.5*Grad_par_LtoC(P + P0);

      if(nonlinear)
        ddt(Vpar) -= bracket(phi, Vpar, bm_exb)*B0; // Advection
    }

    if(filter_z) {
      // Filter out all except filter_z_mode

      ddt(Psi) = filter(ddt(Psi), filter_z_mode);
      ddt(U) = filter(ddt(U), filter_z_mode);
      ddt(P) = filter(ddt(P), filter_z_mode);
    }

    if(low_pass_z > 0) {
      // Low-pass filter, keeping n up to low_pass_z

      ddt(Psi) = lowPass(ddt(Psi), low_pass_z, zonal_field);
      ddt(U) = lowPass(ddt(U), low_pass_z, zonal_flow);
      ddt(P) = lowPass(ddt(P), low_pass_z, zonal_bkgd);
    }

    return 0;
  }
};

BOUTMAIN(ELMpb);
