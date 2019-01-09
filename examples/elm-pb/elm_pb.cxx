/*******************************************************************************
 * High-Beta Flute-Reduced MHD
 * see elm_reduced.pdf
 * Basically the same as Hazeltine-Meiss but different normalisations.
 * Can also include the Vpar compressional term
 *******************************************************************************/

#include <bout/constants.hxx>
#include <bout.hxx>
#include <boutmain.hxx>
#include <derivs.hxx>
#include <initialprofiles.hxx>
#include <interpolation.hxx>
#include <invert_laplace.hxx>
#include <invert_parderiv.hxx>
#include <msg_stack.hxx>
#include <sourcex.hxx>
#include <utils.hxx>

#include <math.h>

// 2D inital profiles
Field2D J0, P0;         // Current and pressure
Vector2D b0xcv;         // Curvature term
Field2D beta, gradparB; // Used for Vpar terms
Field2D phi0;           // When diamagnetic terms used
Field2D U0, Psixy, x;   // 0th vorticity of equilibrium flow,
// radial flux coordinate, normalized radial flux coordinate

bool constn0;
// the total height, average width and center of profile of N0
BoutReal n0_height, n0_ave, n0_width, n0_center, n0_bottom_x, Nbar, Tibar, Tebar;

BoutReal Tconst; // the ampitude of constant temperature

Field2D N0, Ti0, Te0, Ne0; // number density and temperature
Field2D Pi0, Pe0;
Field2D q95;
Field3D ubyn;
BoutReal q95_input;
bool n0_fake_prof, T0_fake_prof;
BoutReal Zi; // charge number of ion

// B field vectors
Vector2D B0vec; // B0 field vector

// V0 field vectors
Vector2D V0net; // net flow

// 3D evolving variables
Field3D U, Psi, P, Vpar;

// Derived 3D variables
Field3D Jpar, phi; // Parallel current, electric potential

Field3D Jpar2; //  Delp2 of Parallel current

Field3D tmpP2; // Grad2_par2new of pressure
Field3D tmpU2; // Grad2_par2new of Parallel vorticity
Field3D tmpA2; // Grad2_par2new of Parallel vector potential

// Constraint
Field3D C_phi;

// Parameters
BoutReal density;              // Number density [m^-3]
BoutReal Bbar, Lbar, Tbar, Va; // Normalisation constants
BoutReal dnorm;                // For diamagnetic terms: 1 / (2. * wci * Tbar)
BoutReal dia_fact;             // Multiply diamagnetic term by this
BoutReal delta_i;              // Normalized ion skin depth
BoutReal omega_i;              // ion gyrofrequency

BoutReal diffusion_p4; // xqx: parallel hyper-viscous diffusion for pressure
BoutReal diffusion_u4; // xqx: parallel hyper-viscous diffusion for vorticity
BoutReal diffusion_a4; // xqx: parallel hyper-viscous diffusion for vector potential

BoutReal diffusion_par; // Parallel pressure diffusion
BoutReal heating_P;     // heating power in pressure
BoutReal hp_width;      // heating profile radial width in pressure
BoutReal hp_length;     // heating radial domain in pressure
BoutReal sink_P;        // sink in pressure
BoutReal sp_width;      // sink profile radial width in pressure
BoutReal sp_length;     // sink radial domain in pressure

BoutReal sink_Ul;    // left edge sink in vorticity
BoutReal su_widthl;  // left edge sink profile radial width in vorticity
BoutReal su_lengthl; // left edge sink radial domain in vorticity

BoutReal sink_Ur;    // right edge sink in vorticity
BoutReal su_widthr;  // right edge sink profile radial width in vorticity
BoutReal su_lengthr; // right edge sink radial domain in vorticity

BoutReal viscos_par;  // Parallel viscosity
BoutReal viscos_perp; // Perpendicular viscosity
BoutReal hyperviscos; // Hyper-viscosity (radial)
Field3D hyper_mu_x;   // Hyper-viscosity coefficient

Field3D Dperp2Phi0, Dperp2Phi, GradPhi02, GradPhi2; // Temporary variables for gyroviscous
Field3D GradparPhi02, GradparPhi2, GradcPhi, GradcparPhi;
Field3D Dperp2Pi0, Dperp2Pi, bracketPhi0P, bracketPhiP0, bracketPhiP;
BoutReal Upara2;

// options
bool include_curvature, include_jpar0, compress0;
bool evolve_pressure, gyroviscous;

BoutReal vacuum_pressure;
BoutReal vacuum_trans; // Transition width
Field3D vac_mask;

int phi_flags, apar_flags;
bool nonlinear;
bool evolve_jpar;
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

// net flow, Er=-R*Bp*Dphi0,Dphi0=-D_min-0.5*D_0*(1.0-tanh(D_s*(x-x0)))
Field2D V0;    // net flow amplitude
Field2D Dphi0; // differential potential to flux
BoutReal D_0;  // potential amplitude
BoutReal D_s;  // shear parameter
BoutReal x0;   // velocity peak location
BoutReal sign; // direction of flow
BoutReal Psiaxis, Psibndry;
bool withflow;
bool K_H_term;  // Kelvin-Holmhotz term
Field2D perp;   // for test
BoutReal D_min; // constant in flow

// for C_mod
bool experiment_Er; // read in total Er from experiment

bool nogradparj;
bool filter_z;
int filter_z_mode;
int low_pass_z;
int zonal_flow;
int zonal_field;
int zonal_bkgd;
bool relax_j_vac;
BoutReal relax_j_tconst; // Time-constant for j relax
Field3D Psitarget;       // The (moving) target to relax to

bool smooth_j_x; // Smooth Jpar in the x direction

int jpar_bndry_width; // Zero jpar in a boundary region

bool sheath_boundaries; // Apply sheath boundaries in Y

bool parallel_lr_diff; // Use left and right shifted stencils for parallel differences

bool phi_constraint; // Solver for phi using a solver constraint

bool include_rmp;     // Include RMP coil perturbation
bool simple_rmp;      // Just use a simple form for the perturbation
int rmp_n, rmp_m;     // toroidal and poloidal mode numbers
BoutReal rmp_polwid;  // Poloidal width (-ve -> full, fraction of 2pi)
BoutReal rmp_polpeak; // Peak poloidal location (fraction of 2pi)
BoutReal rmp_factor;  // Multiply amplitude by this factor
BoutReal rmp_ramp;    // Ramp-up time for RMP [s]. negative -> instant
BoutReal rmp_freq;    // Amplitude oscillation frequency [Hz] (negative -> no oscillation)
BoutReal rmp_rotate;  // Rotation rate [Hz]
bool rmp_vac_mask;    // Should a vacuum mask be applied?
Field3D
    rmp_Psi0; // Parallel vector potential from Resonant Magnetic Perturbation (RMP) coils
Field3D rmp_Psi;   // Value used in calculations
Field3D rmp_dApdt; // Time variation

BoutReal vac_lund, core_lund;     // Lundquist number S = (Tau_R / Tau_A). -ve -> infty
BoutReal vac_resist, core_resist; // The resistivities (just 1 / S)
Field3D eta;                      // Resistivity profile (1 / S)
bool spitzer_resist;              // Use Spitzer formula for resistivity
BoutReal Zeff;                    // Z effective for resistivity formula

BoutReal hyperresist;  // Hyper-resistivity coefficient (in core only)
BoutReal ehyperviscos; // electron Hyper-viscosity coefficient

int damp_width;        // Width of inner damped region
BoutReal damp_t_const; // Timescale of damping

// Metric coefficients
Field2D Rxy, Bpxy, Btxy, B0, hthe;
Field2D I; // Shear factor

const BoutReal MU0 = 4.0e-7 * PI;
const BoutReal Mi = 2.0 * 1.6726e-27; // Ion mass
const BoutReal Me = 9.1094e-31;       // Electron mass
const BoutReal mi_me = Mi / Me;

// Communication objects
FieldGroup comms;

/// Solver for inverting Laplacian
Laplacian* phiSolver;
Laplacian* aparSolver;

int precon(BoutReal t, BoutReal cj, BoutReal delta); // Preconditioner
int jacobian(BoutReal t);                            // Jacobian-vector multiply

int precon_phi(BoutReal t, BoutReal cj,
               BoutReal delta); // Preconditioner with phi constraint

const Field3D Grad2_par2new(const Field3D& f); // for 4th order diffusion

const Field2D N0tanh(BoutReal n0_height, BoutReal n0_ave, BoutReal n0_width,
                     BoutReal n0_center, BoutReal n0_bottom_x) {
  Field2D result;
  result.allocate();

  BoutReal Grid_NX, Grid_NXlimit; // the grid number on x, and the
  BoutReal Jysep;
  mesh->get(Grid_NX, "nx");
  mesh->get(Jysep, "jyseps1_1");
  Grid_NXlimit = n0_bottom_x * Grid_NX;
  output.write("Jysep1_1 = %i   Grid number = %e\n", int(Jysep), Grid_NX);

  if (Jysep > 0.) { // for single null geometry
    BoutReal Jxsep, Jysep2;
    mesh->get(Jxsep, "ixseps1");
    mesh->get(Jysep2, "jyseps2_2");
    // output.write("Jysep2_2 = %i   Ixsep1 = %i\n", int(Jysep2), int(Jxsep));

    for (auto i : result) {
      BoutReal mgx = mesh->GlobalX(i.x());
      BoutReal xgrid_num = (Jxsep + 1.) / Grid_NX;
      // output.write("mgx = %e xgrid_num = %e\n", mgx);

      int globaly = mesh->YGLOBAL(i.y());
      // output.write("local y = %i;   global y: %i\n", i.y, globaly);
      if (mgx > xgrid_num || (globaly <= int(Jysep) - 4) || (globaly > int(Jysep2)))
        mgx = xgrid_num;
      BoutReal rlx = mgx - n0_center;
      BoutReal temp = exp(rlx / n0_width);
      BoutReal dampr = ((temp - 1.0 / temp) / (temp + 1.0 / temp));
      result[i] = 0.5 * (1.0 - dampr) * n0_height + n0_ave;
    }
  } else { // circular geometry
    for (auto i : result) {
      BoutReal mgx = mesh->GlobalX(i.x());
      BoutReal xgrid_num = Grid_NXlimit / Grid_NX;
      if (mgx > xgrid_num)
        mgx = xgrid_num;
      BoutReal rlx = mgx - n0_center;
      BoutReal temp = exp(rlx / n0_width);
      BoutReal dampr = ((temp - 1.0 / temp) / (temp + 1.0 / temp));
      result[i] = 0.5 * (1.0 - dampr) * n0_height + n0_ave;
    }
  }

  mesh->communicate(result);

  return result;
}

/*!
 * This function implements d2/dy2 where y is the poloidal coordinate theta
 */
const Field3D Grad2_par2new(const Field3D& f) {
  TRACE("Grad2_par2new( Field3D )");

  Field3D result = D2DY2(f);

#ifdef TRACK
  result.name = "Grad2_par2new(" + f.name + ")";
#endif

  return result;
}

int physics_init(bool restarting) {
  bool noshear;

  Coordinates* metric = mesh->getCoordinates();

  output.write("Solving high-beta flute reduced equations\n");
  output.write("\tFile    : %s\n", __FILE__);
  output.write("\tCompiled: %s at %s\n", __DATE__, __TIME__);

  //////////////////////////////////////////////////////////////
  // Load data from the grid

  // Load 2D profiles
  mesh->get(J0, "Jpar0");    // A / m^2
  mesh->get(P0, "pressure"); // Pascals

  // Load curvature term
  b0xcv.covariant = false;  // Read contravariant components
  mesh->get(b0xcv, "bxcv"); // mixed units x: T y: m^-2 z: m^-2

  // Load metrics
  if (mesh->get(Rxy, "Rxy")) { // m
    output_error.write("Error: Cannot read Rxy from grid\n");
    return 1;
  }
  if (mesh->get(Bpxy, "Bpxy")) { // T
    output_error.write("Error: Cannot read Bpxy from grid\n");
    return 1;
  }
  mesh->get(Btxy, "Btxy");          // T
  mesh->get(B0, "Bxy");             // T
  mesh->get(hthe, "hthe");          // m
  mesh->get(I, "sinty");            // m^-2 T^-1
  mesh->get(Psixy, "psixy");        // get Psi
  mesh->get(Psiaxis, "psi_axis");   // axis flux
  mesh->get(Psibndry, "psi_bndry"); // edge flux

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

  auto globalOptions = Options::root();
  auto options = globalOptions["highbeta"];

  OPTION(options, constn0, true);
  OPTION(options, n0_fake_prof, false); // use the hyperbolic profile of n0. If both  n0_fake_prof and
                 // T0_fake_prof are false, use the profiles from grid file
  OPTION(options, n0_height, 0.4); // the total height of profile of N0, in percentage of Ni_x
  OPTION(options, n0_ave, 0.01);  // the center or average of N0, in percentage of Ni_x
  OPTION(options, n0_width, 0.1); // the width of the gradient of N0,in percentage of x
  OPTION(options, n0_center, 0.633); // the grid number of the center of N0, in percentage of x
  OPTION(options, n0_bottom_x, 0.81); // the start of flat region of N0 on SOL side, in percentage of x
  OPTION(options, T0_fake_prof, false);
  OPTION(options, Tconst, -1.0); // the amplitude of constant temperature, in percentage

  OPTION(options, density, 1.0e19); // Number density [m^-3]

  OPTION(options, evolve_jpar, false);    // If true, evolve J raher than Psi
  OPTION(options, phi_constraint, false); // Use solver constraint for phi

  // Effects to include/exclude
  OPTION(options, include_curvature, true);
  OPTION(options, include_jpar0, true);
  OPTION(options, evolve_pressure, true);
  OPTION(options, nogradparj, false);

  OPTION(options, compress0, false);
  OPTION(options, gyroviscous, false);
  OPTION(options, nonlinear, false);

  // option for ExB Poisson Bracket
  OPTION(options, bm_exb_flag, 0);
  switch (bm_exb_flag) {
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
  OPTION(options, bm_mag_flag, 0);
  switch (bm_mag_flag) {
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

  OPTION(options, eHall, false); // electron Hall or electron parallel pressue gradient effects?
  OPTION(options, AA, 1.0); // ion mass in units of proton mass

  OPTION(options, diamag, false);         // Diamagnetic effects?
  OPTION(options, diamag_grad_t, diamag); // Grad_par(Te) term in Psi equation
  OPTION(options, diamag_phi0, diamag);   // Include equilibrium phi0
  OPTION(options, dia_fact, 1.0);         // Scale diamagnetic effects by this factor

  OPTION(options, withflow, false); // withflow or not
  OPTION(options, K_H_term, true);  // keep K-H term
  OPTION(options, D_0, 0.0);        // velocity magnitude
  OPTION(options, D_s, 0.0);        // flowshear
  OPTION(options, x0, 0.0);         // flow location
  OPTION(options, sign, 1.0);       // flow direction, -1 means negative electric field
  OPTION(options, D_min, 3000.0);   // a constant

  OPTION(options, experiment_Er, false);

  OPTION(options, noshear, false);

  OPTION(options, relax_j_vac, false); // Relax vacuum current to zero
  OPTION(options, relax_j_tconst, 0.1);

  // Toroidal filtering
  OPTION(options, filter_z, false); // Filter a single n
  OPTION(options, filter_z_mode, 1);
  OPTION(options, low_pass_z, -1);  // Low-pass filter
  OPTION(options, zonal_flow, -1);  // zonal flow filter
  OPTION(options, zonal_field, -1); // zonal field filter
  OPTION(options, zonal_bkgd, -1);  // zonal background P filter

  // Radial smoothing
  OPTION(options, smooth_j_x, false); // Smooth Jpar in x

  // Jpar boundary region
  OPTION(options, jpar_bndry_width, -1);

  OPTION(options, sheath_boundaries, false);

  // Parallel differencing
  OPTION(options, parallel_lr_diff, false);

  // RMP-related options
  OPTION(options, include_rmp, false); // Read RMP data from grid

  OPTION(options, simple_rmp, false); // Include a simple RMP model
  OPTION(options, rmp_factor, 1.0);
  OPTION(options, rmp_ramp, -1.0);
  OPTION(options, rmp_freq, -1.0);
  OPTION(options, rmp_rotate, 0.0);

  // Vacuum region control
  OPTION(options, vacuum_pressure, 0.02); // Fraction of peak pressure
  OPTION(options, vacuum_trans, 0.005);   // Transition width in pressure

  // Resistivity and hyper-resistivity options
  OPTION(options, vac_lund, 0.0);  // Lundquist number in vacuum region
  OPTION(options, core_lund, 0.0); // Lundquist number in core region
  OPTION(options, hyperresist, -1.0);
  OPTION(options, ehyperviscos, -1.0);
  OPTION(options, spitzer_resist, false); // Use Spitzer resistivity
  OPTION(options, Zeff, 2.0);             // Z effective

  // Inner boundary damping
  OPTION(options, damp_width, 0);
  OPTION(options, damp_t_const, 0.1);

  // Viscosity and hyper-viscosity
  OPTION(options, viscos_par, -1.0);  // Parallel viscosity
  OPTION(options, viscos_perp, -1.0); // Perpendicular viscosity
  OPTION(options, hyperviscos, -1.0); // Radial hyperviscosity

  // parallel pressure diffusion
  OPTION(options, diffusion_par, -1.0); // Parallel pressure diffusion
  OPTION(options, diffusion_p4, -1.0); // xqx: parallel hyper-viscous diffusion for pressure
  OPTION(options, diffusion_u4, -1.0); // xqx: parallel hyper-viscous diffusion for vorticity
  OPTION(options, diffusion_a4, -1.0); // xqx: parallel hyper-viscous diffusion for vector potential

  // heating factor in pressure
  OPTION(options, heating_P, -1.0); //  heating power in pressure
  OPTION(options, hp_width, 0.1);   //  the percentage of radial grid points for heating
                                    //  profile radial width in pressure
  OPTION(options, hp_length, 0.04); //  the percentage of radial grid points for heating
                                    //  profile radial domain in pressure

  // sink factor in pressure
  OPTION(options, sink_P, -1.0);    //  sink in pressure
  OPTION(options, sp_width, 0.05);  //  the percentage of radial grid points for sink
                                    //  profile radial width in pressure
  OPTION(options, sp_length, 0.04); //  the percentage of radial grid points for sink
                                    //  profile radial domain in pressure

  // left edge sink factor in vorticity
  OPTION(options, sink_Ul, -1.0);   //  left edge sink in vorticity
  OPTION(options, su_widthl, 0.06); //  the percentage of left edge radial grid points for
                                    //  sink profile radial width in vorticity
  OPTION(options, su_lengthl, 0.15); //  the percentage of left edge radial grid points
                                     //  for sink profile radial domain in vorticity

  // right edge sink factor in vorticity
  OPTION(options, sink_Ur, -1.0);    //  right edge sink in vorticity
  OPTION(options, su_widthr, 0.06);  //  the percentage of right edge radial grid points
                                     //  for sink profile radial width in vorticity
  OPTION(options, su_lengthr, 0.15); //  the percentage of right edge radial grid points
                                     //  for sink profile radial domain in vorticity

  // Compressional terms
  OPTION(options, phi_curv, true);
  g = options["gamma"].withDefault(5.0 / 3.0);

  // Field inversion flags
  OPTION(options, phi_flags, 0);
  OPTION(options, apar_flags, 0);

  x = (Psixy - Psiaxis) / (Psibndry - Psiaxis);

  if (experiment_Er) { // get er from experiment
    mesh->get(Dphi0, "Epsi");
    diamag_phi0 = false;
    K_H_term = false;
  } else {
    Dphi0 = -D_min - 0.5 * D_0 * (1.0 - tanh(D_s * (x - x0)));
  }

  if (sign < 0) // change flow direction
    Dphi0 *= -1;

  V0 = -Rxy * Bpxy * Dphi0 / B0;

  if (simple_rmp)
    include_rmp = true;

  if (include_rmp) {
    // Including external field coils.
    if (simple_rmp) {
      // Use a fairly simple form for the perturbation

      Field2D pol_angle;
      if (mesh->get(pol_angle, "pol_angle")) {
        output_warn.write("     ***WARNING: need poloidal angle for simple RMP\n");
        include_rmp = false;
      } else {
        OPTION(options, rmp_n, 3);
        OPTION(options, rmp_m, 9);
        OPTION(options, rmp_polwid, -1.0);
        OPTION(options, rmp_polpeak, 0.5);
        OPTION(options, rmp_vac_mask, true);
        // Divide n by the size of the domain
        int zperiod = globalOptions["zperiod"].withDefault(1);
        if ((rmp_n % zperiod) != 0)
          output_warn.write("     ***WARNING: rmp_n (%d) not a multiple of zperiod (%d)\n", rmp_n,
                            zperiod);

        output.write("\tMagnetic perturbation: n = %d, m = %d, magnitude %e Tm\n", rmp_n,
                     rmp_m, rmp_factor);

        rmp_Psi0 = 0.0;
        if (mesh->lastX()) {
          // Set the outer boundary
          for (int jx = mesh->LocalNx - 4; jx < mesh->LocalNx; jx++)
            for (int jy = 0; jy < mesh->LocalNy; jy++)
              for (int jz = 0; jz < mesh->LocalNz; jz++) {

                BoutReal angle = rmp_m * pol_angle(jx, jy)
                                 + rmp_n * ((BoutReal)jz) * mesh->getCoordinates()->dz;
                rmp_Psi0(jx, jy, jz) =
                    (((BoutReal)(jx - 4)) / ((BoutReal)(mesh->LocalNx - 5))) * rmp_factor
                    * cos(angle);
                if (rmp_polwid > 0.0) {
                  // Multiply by a Gaussian in poloidal angle
                  BoutReal gx =
                      ((pol_angle(jx, jy) / (2. * PI)) - rmp_polpeak) / rmp_polwid;
                  rmp_Psi0(jx, jy, jz) *= exp(-gx * gx);
                }
              }
        }

        // Now have a simple model for Psi due to coils at the outer boundary
        // Need to calculate Psi inside the domain, enforcing j = 0

        Jpar = 0.0;
        Laplacian* psiLap = Laplacian::create();
        psiLap->setInnerBoundaryFlags(INVERT_AC_GRAD); // Zero gradient inner BC
        psiLap->setOuterBoundaryFlags(INVERT_SET); // Set to rmp_Psi0 on outer boundary
        rmp_Psi0 = psiLap->solve(Jpar, rmp_Psi0);
        mesh->communicate(rmp_Psi0);
      }
    } else {
      // Load perturbation from grid file.
      include_rmp = !mesh->get(rmp_Psi0, "rmp_A"); // Only include if found
      if (!include_rmp) {
        output_warn.write("WARNING: Couldn't read 'rmp_A' from grid file\n");
      }
      // Multiply by factor
      rmp_Psi0 *= rmp_factor;
    }
  }

  if (!include_curvature)
    b0xcv = 0.0;

  if (!include_jpar0)
    J0 = 0.0;

  if (noshear) {
    if (include_curvature)
      b0xcv.z += I * b0xcv.x;
    I = 0.0;
  }

  //////////////////////////////////////////////////////////////
  // SHIFTED RADIAL COORDINATES

  if (mesh->IncIntShear) {
    // BOUT-06 style, using d/dx = d/dpsi + I * d/dz
    metric->IntShiftTorsion = I;

  } else {
    // Dimits style, using local coordinate system
    if (include_curvature)
      b0xcv.z += I * b0xcv.x;
    I = 0.0; // I disappears from metric
  }

  //////////////////////////////////////////////////////////////
  // NORMALISE QUANTITIES

  if (mesh->get(Bbar, "bmag")) // Typical magnetic field
    Bbar = 1.0;
  if (mesh->get(Lbar, "rmag")) // Typical length scale
    Lbar = 1.0;

  Va = sqrt(Bbar * Bbar / (MU0 * density * Mi));

  Tbar = Lbar / Va;

  dnorm = dia_fact * Mi / (2. * 1.602e-19 * Bbar * Tbar);

  delta_i = AA * 60.67 * 5.31e5 / sqrt(density / 1e6) / (Lbar * 100.0);

  output.write("Normalisations: Bbar = %e T   Lbar = %e m\n", Bbar, Lbar);
  output.write("                Va = %e m/s   Tbar = %e s\n", Va, Tbar);
  output.write("                dnorm = %e\n", dnorm);
  output.write("    Resistivity\n");

  if (gyroviscous) {
    omega_i = 9.58e7 * Zeff * Bbar;
    Upara2 = 0.5 / (Tbar * omega_i);
    // Upara3 = 1.0;
    output.write("Upara2 = %e     Omega_i = %e\n", Upara2, omega_i);
  }

  if (eHall)
    output.write("                delta_i = %e   AA = %e \n", delta_i, AA);

  if (vac_lund > 0.0) {
    output.write("        Vacuum  Tau_R = %e s   eta = %e Ohm m\n", vac_lund * Tbar,
                 MU0 * Lbar * Lbar / (vac_lund * Tbar));
    vac_resist = 1. / vac_lund;
  } else {
    output.write("        Vacuum  - Zero resistivity -\n");
    vac_resist = 0.0;
  }
  if (core_lund > 0.0) {
    output.write("        Core    Tau_R = %e s   eta = %e Ohm m\n", core_lund * Tbar,
                 MU0 * Lbar * Lbar / (core_lund * Tbar));
    core_resist = 1. / core_lund;
  } else {
    output.write("        Core    - Zero resistivity -\n");
    core_resist = 0.0;
  }

  if (hyperresist > 0.0) {
    output.write("    Hyper-resistivity coefficient: %e\n", hyperresist);
  }

  if (ehyperviscos > 0.0) {
    output.write("    electron Hyper-viscosity coefficient: %e\n", ehyperviscos);
  }

  if (hyperviscos > 0.0) {
    output.write("    Hyper-viscosity coefficient: %e\n", hyperviscos);
    dump.add(hyper_mu_x, "hyper_mu_x", 1);
  }

  if (diffusion_par > 0.0) {
    output.write("    diffusion_par: %e\n", diffusion_par);
    dump.add(diffusion_par, "diffusion_par", 1);
  }

  // xqx: parallel hyper-viscous diffusion for pressure
  if (diffusion_p4 > 0.0) {
    output.write("    diffusion_p4: %e\n", diffusion_p4);
    dump.add(diffusion_p4, "diffusion_p4", 1);
  }

  // xqx: parallel hyper-viscous diffusion for vorticity
  if (diffusion_u4 > 0.0) {
    output.write("    diffusion_u4: %e\n", diffusion_u4);
    dump.add(diffusion_u4, "diffusion_u4", 1);
  }

  // xqx: parallel hyper-viscous diffusion for vector potential
  if (diffusion_a4 > 0.0) {
    output.write("    diffusion_a4: %e\n", diffusion_a4);
    dump.add(diffusion_a4, "diffusion_a4", 1);
  }

  if (heating_P > 0.0) {
    output.write("    heating_P(watts): %e\n", heating_P);
    dump.add(heating_P, "heating_P", 1);

    output.write("    hp_width(%): %e\n", hp_width);
    dump.add(hp_width, "hp_width", 1);

    output.write("    hp_length(%): %e\n", hp_length);
    dump.add(hp_length, "hp_length", 1);
  }

  if (sink_P > 0.0) {
    output.write("    sink_P(rate): %e\n", sink_P);
    dump.add(sink_P, "sink_P", 1);

    output.write("    sp_width(%): %e\n", sp_width);
    dump.add(sp_width, "sp_width", 1);

    output.write("    sp_length(%): %e\n", sp_length);
    dump.add(sp_length, "sp_length", 1);
  }

  if (K_H_term)
    output.write("    keep K-H term\n");
  else
    output.write("   drop K-H term\n");

  Field2D Te;
  Te = P0 / (2.0 * density * 1.602e-19); // Temperature in eV

  J0 = -MU0 * Lbar * J0 / B0;
  P0 = 2.0 * MU0 * P0 / (Bbar * Bbar);
  V0 = V0 / Va;
  Dphi0 *= Tbar;

  b0xcv.x /= Bbar;
  b0xcv.y *= Lbar * Lbar;
  b0xcv.z *= Lbar * Lbar;

  Rxy /= Lbar;
  Bpxy /= Bbar;
  Btxy /= Bbar;
  B0 /= Bbar;
  hthe /= Lbar;
  metric->dx /= Lbar * Lbar * Bbar;
  I *= Lbar * Lbar * Bbar;

  if (constn0) {
    T0_fake_prof = false;
    n0_fake_prof = false;
  } else {
    Nbar = 1.0;
    Tibar = 1000.0;
    Tebar = 1000.0;

    if ((!T0_fake_prof) && n0_fake_prof) {
      N0 = N0tanh(n0_height * Nbar, n0_ave * Nbar, n0_width, n0_center, n0_bottom_x);

      Ti0 = P0 / N0 / 2.0;
      Te0 = Ti0;
    } else if (T0_fake_prof) {
      Ti0 = Tconst;
      Te0 = Ti0;
      N0 = P0 / (Ti0 + Te0);
    } else {
      if (mesh->get(N0, "Niexp")) { // N_i0
        output_error.write("Error: Cannot read Ni0 from grid\n");
        return 1;
      }

      if (mesh->get(Ti0, "Tiexp")) { // T_i0
        output_error.write("Error: Cannot read Ti0 from grid\n");
        return 1;
      }

      if (mesh->get(Te0, "Teexp")) { // T_e0
        output_error.write("Error: Cannot read Te0 from grid\n");
        return 1;
      }
      N0 /= Nbar;
      Ti0 /= Tibar;
      Te0 /= Tebar;
    }
  }

  if (gyroviscous) {
    Dperp2Phi0.setLocation(CELL_CENTRE);
    Dperp2Phi0.setBoundary("phi");
    Dperp2Phi.setLocation(CELL_CENTRE);
    Dperp2Phi.setBoundary("phi");
    GradPhi02.setLocation(CELL_CENTRE);
    GradPhi02.setBoundary("phi");
    GradcPhi.setLocation(CELL_CENTRE);
    GradcPhi.setBoundary("phi");
    Dperp2Pi0.setLocation(CELL_CENTRE);
    Dperp2Pi0.setBoundary("P");
    Dperp2Pi.setLocation(CELL_CENTRE);
    Dperp2Pi.setBoundary("P");
    bracketPhi0P.setLocation(CELL_CENTRE);
    bracketPhi0P.setBoundary("P");
    bracketPhiP0.setLocation(CELL_CENTRE);
    bracketPhiP0.setBoundary("P");
    if (nonlinear) {
      GradPhi2.setLocation(CELL_CENTRE);
      GradPhi2.setBoundary("phi");
      bracketPhiP.setLocation(CELL_CENTRE);
      bracketPhiP.setBoundary("P");
    }
  }

  BoutReal pnorm = max(P0, true); // Maximum over all processors

  vacuum_pressure *= pnorm; // Get pressure from fraction
  vacuum_trans *= pnorm;

  // Transitions from 0 in core to 1 in vacuum
  vac_mask = (1.0 - tanh((P0 - vacuum_pressure) / vacuum_trans)) / 2.0;

  if (spitzer_resist) {
    // Use Spitzer resistivity
    output.write("\tTemperature: %e -> %e [eV]\n", min(Te), max(Te));
    eta = 0.51 * 1.03e-4 * Zeff * 20.
          * pow(Te, -1.5); // eta in Ohm-m. NOTE: ln(Lambda) = 20
    output.write("\tSpitzer resistivity: %e -> %e [Ohm m]\n", min(eta), max(eta));
    eta /= MU0 * Va * Lbar;
    output.write("\t -> Lundquist %e -> %e\n", 1.0 / max(eta), 1.0 / min(eta));
  } else {
    // transition from 0 for large P0 to resistivity for small P0
    eta = core_resist + (vac_resist - core_resist) * vac_mask;
  }

  SAVE_ONCE(eta);

  if (include_rmp) {
    // Normalise RMP quantities

    rmp_Psi0 /= Bbar * Lbar;

    rmp_ramp /= Tbar;
    rmp_freq *= Tbar;
    rmp_rotate *= Tbar;

    rmp_Psi = rmp_Psi0;
    rmp_dApdt = 0.0;

    bool apar_changing = false;

    output.write("Including magnetic perturbation\n");
    if (rmp_ramp > 0.0) {
      output.write("\tRamping up over period t = %e (%e ms)\n", rmp_ramp,
                   rmp_ramp * Tbar * 1000.);
      apar_changing = true;
    }
    if (rmp_freq > 0.0) {
      output.write("\tOscillating with frequency f = %e (%e kHz)\n", rmp_freq,
                   rmp_freq / Tbar / 1000.);
      apar_changing = true;
    }
    if (rmp_rotate != 0.0) {
      output.write("\tRotating with a frequency f = %e (%e kHz)\n", rmp_rotate,
                   rmp_rotate / Tbar / 1000.);
      apar_changing = true;
    }

    if (apar_changing) {
      SAVE_REPEAT(rmp_Psi, rmp_dApdt);
    } else {
      SAVE_ONCE(rmp_Psi);
    }
  } else
    rmp_Psi = 0.0;

  /**************** CALCULATE METRICS ******************/

  metric->g11 = SQ(Rxy * Bpxy);
  metric->g22 = 1.0 / SQ(hthe);
  metric->g33 = SQ(I) * metric->g11 + SQ(B0) / metric->g11;
  metric->g12 = 0.0;
  metric->g13 = -I * metric->g11;
  metric->g23 = -Btxy / (hthe * Bpxy * Rxy);

  metric->J = hthe / Bpxy;
  metric->Bxy = B0;

  metric->g_11 = 1.0 / metric->g11 + SQ(I * Rxy);
  metric->g_22 = SQ(B0 * hthe / Bpxy);
  metric->g_33 = Rxy * Rxy;
  metric->g_12 = Btxy * hthe * I * Rxy / Bpxy;
  metric->g_13 = I * Rxy * Rxy;
  metric->g_23 = Btxy * hthe * Rxy / Bpxy;

  metric->geometry(); // Calculate quantities from metric tensor

  // Set B field vector

  B0vec.covariant = false;
  B0vec.x = 0.;
  B0vec.y = Bpxy / hthe;
  B0vec.z = 0.;

  V0net.covariant = false; // presentation for net flow
  V0net.x = 0.;
  V0net.y = Rxy * Btxy * Bpxy / (hthe * B0 * B0) * Dphi0;
  V0net.z = -Dphi0;

  U0 = B0vec * Curl(V0net) / B0; // get 0th vorticity for Kelvin-Holmholtz term

  /**************** SET EVOLVING VARIABLES *************/

  // Tell BOUT which variables to evolve
  SOLVE_FOR(U, P);

  if (evolve_jpar) {
    output.write("Solving for jpar: Inverting to get Psi\n");
    SOLVE_FOR(Jpar);
    dump.add(Psi, "Psi", 1);
  } else {
    output.write("Solving for Psi, Differentiating to get jpar\n");
    SOLVE_FOR(Psi);
    dump.add(Jpar, "jpar", 1);
  }

  if (compress0) {
    output.write("Including compression (Vpar) effects\n");

    SOLVE_FOR(Vpar);
    comms.add(Vpar);

    beta = B0 * B0 / (0.5 + (B0 * B0 / (g * P0)));
    gradparB = Grad_par(B0) / B0;

    output.write("Beta in range %e -> %e\n", min(beta), max(beta));
  } else {
    Vpar = 0.0;
  }

  if (phi_constraint) {
    // Implicit Phi solve using IDA

    if (!bout_constrain(phi, C_phi, "phi")) {
      output_error.write("ERROR: Cannot constrain. Run again with phi_constraint=false\n");
      throw BoutException("Aborting.\n");
    }

    // Set preconditioner
    solver->setPrecon(precon_phi);

  } else {
    // Phi solved in RHS (explicitly)
    dump.add(phi, "phi", 1);

    // Set preconditioner
    solver->setPrecon(precon);

    // Set Jacobian
    solver->setJacobian(jacobian);
  }

  // Diamagnetic phi0
  if (diamag_phi0) {
    if (constn0)
      phi0 = -0.5 * dnorm * P0 / B0;
    else
      // Stationary equilibrium plasma. ExB velocity balances diamagnetic drift
      phi0 = -0.5 * dnorm * P0 / B0 / N0;
    SAVE_ONCE(phi0);
  }

  // Add some equilibrium quantities and normalisations
  // everything needed to recover physical units
  SAVE_ONCE(J0, P0);
  SAVE_ONCE(density, Lbar, Bbar, Tbar);
  SAVE_ONCE(Va, B0);
  SAVE_ONCE(Dphi0, U0);
  SAVE_ONCE(V0);
  if (!constn0)
    SAVE_ONCE(Ti0, Te0, N0);

  // Create a solver for the Laplacian
  phiSolver = Laplacian::create();
  phiSolver->setFlags(phi_flags);

  aparSolver = Laplacian::create();
  aparSolver->setFlags(apar_flags);
  aparSolver->setCoefA(-delta_e_inv * N0 * N0);
  
  /////////////// CHECK VACUUM ///////////////////////
  // In vacuum region, initial vorticity should equal zero

  // ubyn.setLocation(CELL_YLOW);
  // ubyn.setBoundary("U");

  if (!restarting) {
    // Only if not restarting: Check initial perturbation

    // Set U to zero where P0 < vacuum_pressure
    U = where(P0 - vacuum_pressure, U, 0.0);

    if (constn0) {
      ubyn = U;
      // Phi should be consistent with U
      phi = phiSolver->solve(ubyn);
    } else {
      ubyn = U / N0;
      phiSolver->setCoefC(N0);
      phi = phiSolver->solve(ubyn);
    }

    // if(diamag) {
    // phi -= 0.5*dnorm * P / B0;
    //}
  }

  /************** SETUP COMMUNICATIONS **************/

  comms.add(U, P);

  phi.setBoundary("phi"); // Set boundary conditions
  tmpU2.setBoundary("U");
  tmpP2.setBoundary("P");
  tmpA2.setBoundary("J");

  if (evolve_jpar) {
    comms.add(Jpar);
  } else {
    comms.add(Psi);
    // otherwise Need to communicate Jpar separately
    Jpar.setBoundary("J");
  }
  Jpar2.setBoundary("J");

  return 0;
}

// Parallel gradient along perturbed field-line
const Field3D Grad_parP(const Field3D& f, CELL_LOC loc = CELL_DEFAULT) {
  Field3D result;

  if (parallel_lr_diff) {
    // Use left/right biased stencils. NOTE: First order only!
    if (loc == CELL_YLOW) {
      result = Grad_par_CtoL(f);
    } else
      result = Grad_par_LtoC(f);
  } else
    result = Grad_par(f, loc);

  if (nonlinear) {
    result -= bracket(Psi, f, bm_mag) * B0;

    if (include_rmp) {
      result -= bracket(rmp_Psi, f, bm_mag) * B0;
    }
  }

  return result;
}

bool first_run = true; // For printing out some diagnostics first time around

int physics_run(BoutReal t) {
  // Perform communications
  mesh->communicate(comms);

  Coordinates* metric = mesh->getCoordinates();

  ////////////////////////////////////////////
  // Transitions from 0 in core to 1 in vacuum
  if (nonlinear) {
    vac_mask = (1.0 - tanh(((P0 + P) - vacuum_pressure) / vacuum_trans)) / 2.0;

    // Update resistivity
    if (spitzer_resist) {
      // Use Spitzer formula
      Field3D Te;
      Te = (P0 + P) * Bbar * Bbar / (4. * MU0) / (density * 1.602e-19); // eV
      eta = 0.51 * 1.03e-4 * Zeff * 20. * pow(Te, -1.5); // eta in Ohm-m. ln(Lambda) = 20
      eta /= MU0 * Va * Lbar;                            // Normalised eta
    } else {
      // Use specified core and vacuum Lundquist numbers
      eta = core_resist + (vac_resist - core_resist) * vac_mask;
    }
  }

  ////////////////////////////////////////////
  // Resonant Magnetic Perturbation code

  if (include_rmp) {

    if ((rmp_ramp > 0.0) || (rmp_freq > 0.0) || (rmp_rotate != 0.0)) {
      // Need to update the RMP terms

      if ((rmp_ramp > 0.0) && (t < rmp_ramp)) {
        // Still in ramp phase

        rmp_Psi = (t / rmp_ramp) * rmp_Psi0; // Linear ramp

        rmp_dApdt = rmp_Psi0 / rmp_ramp;
      } else {
        rmp_Psi = rmp_Psi0;
        rmp_dApdt = 0.0;
      }

      if (rmp_freq > 0.0) {
        // Oscillating the amplitude

        rmp_dApdt = rmp_dApdt * sin(2. * PI * rmp_freq * t)
                    + rmp_Psi * (2. * PI * rmp_freq) * cos(2. * PI * rmp_freq * t);

        rmp_Psi *= sin(2. * PI * rmp_freq * t);
      }

      if (rmp_rotate != 0.0) {
        // Rotate toroidally at given frequency

        shiftZ(rmp_Psi, 2 * PI * rmp_rotate * t);
        shiftZ(rmp_dApdt, 2 * PI * rmp_rotate * t);

        // Add toroidal rotation term. CHECK SIGN

        rmp_dApdt += DDZ(rmp_Psi) * 2 * PI * rmp_rotate;
      }

      // Set to zero in the core
      if (rmp_vac_mask)
        rmp_Psi *= vac_mask;
    } else {
      // Set to zero in the core region
      if (rmp_vac_mask)
        rmp_Psi =
            rmp_Psi0 * vac_mask; // Only in vacuum -> skin current -> diffuses inwards
    }

    mesh->communicate(rmp_Psi);
  }

  ////////////////////////////////////////////
  // Inversion

  if (evolve_jpar) {
    // Invert laplacian for Psi
    Psi = aparSolver->solve(Jpar);
    mesh->communicate(Psi);
  }

  if (phi_constraint) {
    // Phi being solved as a constraint

    Field3D Ctmp = phi;
    Ctmp.setBoundary("phi"); // Look up boundary conditions for phi
    Ctmp.applyBoundary();
    Ctmp -= phi; // Now contains error in the boundary

    C_phi = Delp2(phi) - U; // Error in the bulk
    C_phi.setBoundaryTo(Ctmp);

  } else {

    if (constn0) {
      phi = phiSolver->solve(U);

      if (diamag) {
        phi -= 0.5 * dnorm * P / B0;
      }
    } else {
      ubyn = U / N0;
      if (diamag) {
        ubyn -= 0.5 * dnorm / (N0 * B0) * Delp2(P);
        mesh->communicate(ubyn);
      }
      // Invert laplacian for phi
      phiSolver->setCoefC(N0);
      phi = phiSolver->solve(ubyn);
    }
    // Apply a boundary condition on phi for target plates
    phi.applyBoundary();
    mesh->communicate(phi);
  }

  if (!evolve_jpar) {
    // Get J from Psi
    Jpar = Delp2(Psi);
    if (include_rmp)
      Jpar += Delp2(rmp_Psi);

    Jpar.applyBoundary();
    mesh->communicate(Jpar);

    if (jpar_bndry_width > 0) {
      // Zero j in boundary regions. Prevents vorticity drive
      // at the boundary

      for (int i = 0; i < jpar_bndry_width; i++)
        for (int j = 0; j < mesh->LocalNy; j++)
          for (int k = 0; k < mesh->LocalNz; k++) {
            if (mesh->firstX())
              Jpar(i, j, k) = 0.0;
            if (mesh->lastX())
              Jpar(mesh->LocalNx - 1 - i, j, k) = 0.0;
          }
    }

    // Smooth j in x
    if (smooth_j_x) {
      Jpar = smooth_x(Jpar);
      Jpar.applyBoundary();

      // Recommunicate now smoothed
      mesh->communicate(Jpar);
    }
    
    // Get Delp2(J) from J
    Jpar2 = Delp2(Jpar);

    Jpar2.applyBoundary();
    mesh->communicate(Jpar2);

    if (jpar_bndry_width > 0) {
      // Zero jpar2 in boundary regions. Prevents vorticity drive
      // at the boundary

      for (int i = 0; i < jpar_bndry_width; i++)
        for (int j = 0; j < mesh->LocalNy; j++)
          for (int k = 0; k < mesh->LocalNz; k++) {
            if (mesh->firstX())
              Jpar2(i, j, k) = 0.0;
            if (mesh->lastX())
              Jpar2(mesh->LocalNx - 1 - i, j, k) = 0.0;
          }
    }
  }

  ////////////////////////////////////////////////////
  // Sheath boundary conditions
  // Normalised and linearised, since here we have only pressure
  // rather than density and temperature. Applying a boundary
  // to Jpar so that Jpar = sqrt(mi/me)/(2*pi) * phi
  //

  if (sheath_boundaries) {

    // At y = ystart (lower boundary)

    for (RangeIterator r = mesh->iterateBndryLowerY(); !r.isDone(); r++) {
      for (int jz = 0; jz < mesh->LocalNz; jz++) {

        // Zero-gradient potential
        BoutReal phisheath = phi(r.ind, mesh->ystart, jz);

        BoutReal jsheath = -(sqrt(mi_me) / (2. * sqrt(PI))) * phisheath;

        // Apply boundary condition half-way between cells
        for (int jy = mesh->ystart - 1; jy >= 0; jy--) {
          // Neumann conditions
          P(r.ind, jy, jz) = P(r.ind, mesh->ystart, jz);
          phi(r.ind, jy, jz) = phisheath;
          // Dirichlet condition on Jpar
          Jpar(r.ind, jy, jz) = 2. * jsheath - Jpar(r.ind, mesh->ystart, jz);
        }
      }
    }

    // At y = yend (upper boundary)

    for (RangeIterator r = mesh->iterateBndryUpperY(); !r.isDone(); r++) {
      for (int jz = 0; jz < mesh->LocalNz; jz++) {

        // Zero-gradient potential
        BoutReal phisheath = phi(r.ind, mesh->yend, jz);

        BoutReal jsheath = (sqrt(mi_me) / (2. * sqrt(PI))) * phisheath;

        // Apply boundary condition half-way between cells
        for (int jy = mesh->yend - 1; jy >= 0; jy--) {
          // Neumann conditions
          P(r.ind, jy, jz) = P(r.ind, mesh->yend, jz);
          phi(r.ind, jy, jz) = phisheath;
          // Dirichlet condition on Jpar
          Jpar(r.ind, jy, jz) = 2. * jsheath - Jpar(r.ind, mesh->yend, jz);
        }
      }
    }
  }

  ////////////////////////////////////////////////////
  // Parallel electric field

  if (evolve_jpar) {
    // Jpar
    Field3D B0U = B0 * U;
    mesh->communicate(B0U);
    ddt(Jpar) = -Grad_parP(B0U, CELL_YLOW) / B0 + eta * Delp2(Jpar);

    if (relax_j_vac) {
      // Make ddt(Jpar) relax to zero.

      ddt(Jpar) -= vac_mask * Jpar / relax_j_tconst;
    }
  } else {
    // Vector potential
    ddt(Psi) = -Grad_parP(phi, CELL_CENTRE) + eta * Jpar;

    if (eHall) {
      ddt(Psi) += 0.25 * delta_i
                  * (Grad_parP(P, CELL_CENTRE)
                     + bracket(P0, Psi, bm_mag)); // electron parallel pressure
    }

    if (diamag_phi0)
      ddt(Psi) -= bracket(phi0, Psi, bm_exb); // Equilibrium flow

    if (withflow) // net flow
      ddt(Psi) -= V_dot_Grad(V0net, Psi);

    if (diamag_grad_t) {
      // grad_par(T_e) correction

      ddt(Psi) += 1.71 * dnorm * 0.5 * Grad_parP(P, CELL_YLOW) / B0;
    }

    // Hyper-resistivity
    if (hyperresist > 0.0) {
      ddt(Psi) -= eta * hyperresist * Delp2(Jpar);
    }

    // electron Hyper-viscosity coefficient
    if (ehyperviscos > 0.0) {
      ddt(Psi) -= eta * ehyperviscos * Delp2(Jpar2);
    }

    // xqx: parallel hyper-viscous diffusion for vector potential
    if (diffusion_a4 > 0.0) {
      tmpA2 = Grad2_par2new(Psi);
      mesh->communicate(tmpA2);
      tmpA2.applyBoundary();
      ddt(Psi) -= diffusion_a4 * Grad2_par2new(tmpA2);
    }

    // Vacuum solution
    if (relax_j_vac) {
      // Calculate the J and Psi profile we're aiming for
      Field3D Jtarget = Jpar * (1.0 - vac_mask); // Zero in vacuum

      // Invert laplacian for Psi
      Psitarget = aparSolver->solve(Jtarget);

      // Add a relaxation term in the vacuum
      ddt(Psi) = ddt(Psi) * (1. - vac_mask) - (Psi - Psitarget) * vac_mask / relax_j_tconst;
    }
  }

  ////////////////////////////////////////////////////
  // Vorticity equation

  // Grad j term
  ddt(U) = SQ(B0) * b0xGrad_dot_Grad(Psi, J0, CELL_CENTRE);
  if (include_rmp) {
    ddt(U) += SQ(B0) * b0xGrad_dot_Grad(rmp_Psi, J0, CELL_CENTRE);
  }

  ddt(U) += b0xcv * Grad(P); // curvature term

  if (!nogradparj) {
    // Parallel current term
    ddt(U) -= SQ(B0) * Grad_parP(Jpar, CELL_CENTRE); // b dot grad j
  }

  if (withflow && K_H_term) // K_H_term
    ddt(U) -= b0xGrad_dot_Grad(phi, U0);

  if (diamag_phi0)
    ddt(U) -= b0xGrad_dot_Grad(phi0, U); // Equilibrium flow

  if (withflow) // net flow
    ddt(U) -= V_dot_Grad(V0net, U);

  if (nonlinear) {
    ddt(U) -= bracket(phi, U, bm_exb) * B0; // Advection
  }

  // Viscosity terms
  if (viscos_par > 0.0)
    ddt(U) += viscos_par * Grad2_par2(U); // Parallel viscosity

  // xqx: parallel hyper-viscous diffusion for vorticity
  if (diffusion_u4 > 0.0) {
    tmpU2 = Grad2_par2new(U);
    mesh->communicate(tmpU2);
    tmpU2.applyBoundary();
    //    tmpU2.applyBoundary("neumann");
    ddt(U) -= diffusion_u4 * Grad2_par2new(tmpU2);
  }

  if (viscos_perp > 0.0)
    ddt(U) += viscos_perp * Delp2(U); // Perpendicular viscosity

  // Hyper-viscosity
  if (hyperviscos > 0.0) {
    // Calculate coefficient.

    hyper_mu_x = hyperviscos * metric->g_11 * SQ(metric->dx) * abs(metric->g11 * D2DX2(U))
                 / (abs(U) + 1e-3);
    hyper_mu_x.applyBoundary("dirichlet"); // Set to zero on all boundaries

    ddt(U) += hyper_mu_x * metric->g11 * D2DX2(U);

    if (first_run) { // Print out maximum values of viscosity used on this processor
      output.write("   Hyper-viscosity values:\n");
      output.write("      Max mu_x = %e, Max_DC mu_x = %e\n", max(hyper_mu_x),
                   max(DC(hyper_mu_x)));
    }
  }

  if (gyroviscous) {

    Field3D Pi;
    Field2D Pi0;
    Pi = 0.5 * P;
    Pi0 = 0.5 * P0;

    Dperp2Phi0 = Field3D(Delp2(B0 * phi0));
    Dperp2Phi0.applyBoundary();
    mesh->communicate(Dperp2Phi0);

    Dperp2Phi = Delp2(B0 * phi);
    Dperp2Phi.applyBoundary();
    mesh->communicate(Dperp2Phi);

    Dperp2Pi0 = Field3D(Delp2(Pi0));
    Dperp2Pi0.applyBoundary();
    mesh->communicate(Dperp2Pi0);

    Dperp2Pi = Delp2(Pi);
    Dperp2Pi.applyBoundary();
    mesh->communicate(Dperp2Pi);

    bracketPhi0P = bracket(B0 * phi0, Pi, bm_exb);
    bracketPhi0P.applyBoundary();
    mesh->communicate(bracketPhi0P);

    bracketPhiP0 = bracket(B0 * phi, Pi0, bm_exb);
    bracketPhiP0.applyBoundary();
    mesh->communicate(bracketPhiP0);

    ddt(U) -= 0.5 * Upara2 * bracket(Pi, Dperp2Phi0, bm_exb) / B0;
    ddt(U) -= 0.5 * Upara2 * bracket(Pi0, Dperp2Phi, bm_exb) / B0;
    Field3D B0phi = B0 * phi;
    mesh->communicate(B0phi);
    Field3D B0phi0 = B0 * phi0;
    mesh->communicate(B0phi0);
    ddt(U) += 0.5 * Upara2 * bracket(B0phi, Dperp2Pi0, bm_exb) / B0;
    ddt(U) += 0.5 * Upara2 * bracket(B0phi0, Dperp2Pi, bm_exb) / B0;
    ddt(U) -= 0.5 * Upara2 * Delp2(bracketPhi0P) / B0;
    ddt(U) -= 0.5 * Upara2 * Delp2(bracketPhiP0) / B0;

    if (nonlinear) {
      Field3D B0phi = B0 * phi;
      mesh->communicate(B0phi);
      bracketPhiP = bracket(B0phi, Pi, bm_exb);
      bracketPhiP.applyBoundary();
      mesh->communicate(bracketPhiP);

      ddt(U) -= 0.5 * Upara2 * bracket(Pi, Dperp2Phi, bm_exb) / B0;
      ddt(U) += 0.5 * Upara2 * bracket(B0phi, Dperp2Pi, bm_exb) / B0;
      ddt(U) -= 0.5 * Upara2 * Delp2(bracketPhiP) / B0;
    }
  }

  // left edge sink terms
  if (sink_Ul > 0.0) {
    ddt(U) -= sink_Ul * sink_tanhxl(P0, U, su_widthl, su_lengthl); // core sink
  }

  // right edge sink terms
  if (sink_Ur > 0.0) {
    ddt(U) -= sink_Ur * sink_tanhxr(P0, U, su_widthr, su_lengthr); //  sol sink
  }

  ////////////////////////////////////////////////////
  // Pressure equation

  ddt(P) = 0.0;
  if (evolve_pressure) {
    ddt(P) -= b0xGrad_dot_Grad(phi, P0);

    if (diamag_phi0)
      ddt(P) -= b0xGrad_dot_Grad(phi0, P); // Equilibrium flow

    if (withflow) // net flow
      ddt(P) -= V_dot_Grad(V0net, P);

    if (nonlinear)
      ddt(P) -= bracket(phi, P, bm_exb) * B0; // Advection
  }

  // Parallel diffusion terms
  if (diffusion_par > 0.0)
    ddt(P) += diffusion_par * Grad2_par2(P); // Parallel diffusion

  // xqx: parallel hyper-viscous diffusion for pressure
  if (diffusion_p4 > 0.0) {
    tmpP2 = Grad2_par2new(P);
    mesh->communicate(tmpP2);
    tmpP2.applyBoundary();
    ddt(P) = diffusion_p4 * Grad2_par2new(tmpP2);
  }

  // heating source terms
  if (heating_P > 0.0) {
    BoutReal pnorm = P0(0, 0);
    ddt(P) += heating_P * source_expx2(P0, 2. * hp_width, 0.5 * hp_length)
              * (Tbar / pnorm); // heat source
    ddt(P) += (100. * source_tanhx(P0, hp_width, hp_length) + 0.01) * metric->g11
              * D2DX2(P) * (Tbar / Lbar / Lbar); // radial diffusion
  }

  // sink terms
  if (sink_P > 0.0) {
    ddt(P) -= sink_P * sink_tanhxr(P0, P, sp_width, sp_length) * Tbar; // sink
  }

  ////////////////////////////////////////////////////
  // Compressional effects

  if (compress0) {

    // ddt(P) += beta*( - Grad_parP(Vpar, CELL_CENTRE) + Vpar*gradparB );
    ddt(P) -= beta * Div_par_CtoL(Vpar);

    if (phi_curv) {
      ddt(P) -= 2. * beta * b0xcv * Grad(phi);
    }

    // Vpar equation

    // ddt(Vpar) = -0.5*Grad_parP(P + P0, CELL_YLOW);
    ddt(Vpar) = -0.5 * (Grad_par_LtoC(P) + Grad_par_LtoC(P0));

    if (nonlinear)
      ddt(Vpar) -= bracket(phi, Vpar, bm_exb) * B0; // Advection
  }

  if (filter_z) {
    // Filter out all except filter_z_mode

    if (evolve_jpar) {
      ddt(Jpar) = filter(ddt(Jpar), filter_z_mode);
    } else
      ddt(Psi) = filter(ddt(Psi), filter_z_mode);

    ddt(U) = filter(ddt(U), filter_z_mode);
    ddt(P) = filter(ddt(P), filter_z_mode);
  }

  if (low_pass_z > 0) {
    // Low-pass filter, keeping n up to low_pass_z
    if (evolve_jpar) {
      ddt(Jpar) = lowPass(ddt(Jpar), low_pass_z, zonal_field);
    } else
      ddt(Psi) = lowPass(ddt(Psi), low_pass_z, zonal_field);

    ddt(U) = lowPass(ddt(U), low_pass_z, zonal_flow);
    ddt(P) = lowPass(ddt(P), low_pass_z, zonal_bkgd);
  }

  if (damp_width > 0) {
    for (int i = 0; i < damp_width; i++) {
      for (int j = 0; j < mesh->LocalNy; j++)
        for (int k = 0; k < mesh->LocalNz; k++) {
          if (mesh->firstX())
            ddt(U)(i, j, k) -= U(i, j, k) / damp_t_const;
          if (mesh->lastX())
            ddt(U)(mesh->LocalNx - 1 - i, j, k) -=
                U(mesh->LocalNx - 1 - i, j, k) / damp_t_const;
        }
    }
  }

  first_run = false;

  return 0;
}

/*******************************************************************************
 * Preconditioner
 *
 * o System state in variables (as in rhs function)
 * o Values to be inverted in time derivatives
 *
 * o Return values should be in time derivatives
 *
 * enable by setting solver / use_precon = true in BOUT.inp
 *******************************************************************************/

int precon(BoutReal UNUSED(t), BoutReal gamma, BoutReal UNUSED(delta)) {
  // First matrix, applying L
  mesh->communicate(ddt(Psi));
  Field3D Jrhs = Delp2(ddt(Psi));
  Jrhs.applyBoundary("neumann");

  if (jpar_bndry_width > 0) {
    // Boundary in jpar
    if (mesh->firstX()) {
      for (int i = jpar_bndry_width; i >= 0; i--)
        for (int j = 0; j < mesh->LocalNy; j++)
          for (int k = 0; k < mesh->LocalNz; k++) {
            Jrhs(i, j, k) = 0.5 * Jrhs(i + 1, j, k);
          }
    }
    if (mesh->lastX()) {
      for (int i = mesh->LocalNx - jpar_bndry_width - 1; i < mesh->LocalNx; i++)
        for (int j = 0; j < mesh->LocalNy; j++)
          for (int k = 0; k < mesh->LocalNz; k++) {
            Jrhs(i, j, k) = 0.5 * Jrhs(i - 1, j, k);
          }
    }
  }

  mesh->communicate(Jrhs, ddt(P));

  Field3D U1 = ddt(U);
  U1 += (gamma * B0 * B0) * Grad_par(Jrhs, CELL_CENTRE) + (gamma * b0xcv) * Grad(P);

  // Second matrix, solving Alfven wave dynamics
  static InvertPar* invU = 0;
  if (!invU)
    invU = InvertPar::Create();

  invU->setCoefA(1.);
  invU->setCoefB(-SQ(gamma) * B0 * B0);
  ddt(U) = invU->solve(U1);
  ddt(U).applyBoundary();

  // Third matrix, applying U
  Field3D phi3 = phiSolver->solve(ddt(U));
  mesh->communicate(phi3);
  phi3.applyBoundary("neumann");
  Field3D B0phi3 = B0 * phi3;
  mesh->communicate(B0phi3);
  ddt(Psi) = ddt(Psi) - gamma * Grad_par(B0phi3) / B0;
  ddt(Psi).applyBoundary();

  return 0;
}

/*******************************************************************************
 * Jacobian-vector multiply
 *
 * Input
 *   System state is in (P, Psi, U)
 *   Vector v is in (F_P, F_Psi, F_U)
 * Output
 *   Jacobian-vector multiplied Jv should be in (P, Psi, U)
 *
 * NOTE: EXPERIMENTAL
 * enable by setting solver / use_jacobian = true in BOUT.inp
 *******************************************************************************/

int jacobian(BoutReal UNUSED(t)) {
  // NOTE: LINEAR ONLY!

  // Communicate
  mesh->communicate(ddt(P), ddt(Psi), ddt(U));

  phi = phiSolver->solve(ddt(U));

  Jpar = Delp2(ddt(Psi));

  mesh->communicate(phi, Jpar);

  Field3D JP = -b0xGrad_dot_Grad(phi, P0);
  JP.setBoundary("P");
  JP.applyBoundary();
  Field3D B0phi = B0 * phi;
  mesh->communicate(B0phi);
  Field3D JPsi = -Grad_par(B0phi, CELL_YLOW) / B0;
  JPsi.setBoundary("Psi");
  JPsi.applyBoundary();

  Field3D JU = b0xcv * Grad(ddt(P)) - SQ(B0) * Grad_par(Jpar, CELL_CENTRE)
               + SQ(B0) * b0xGrad_dot_Grad(ddt(Psi), J0, CELL_CENTRE);
  JU.setBoundary("U");
  JU.applyBoundary();

  // Put result into time-derivatives

  ddt(P) = JP;
  ddt(Psi) = JPsi;
  ddt(U) = JU;

  return 0;
}

/*******************************************************************************
 * Preconditioner for when phi solved as a constraint
 * Currently only possible with the IDA solver
 *
 * o System state in variables (as in rhs function)
 * o Values to be inverted in F_vars
 *
 * o Return values should be in vars (overwriting system state)
 *******************************************************************************/

int precon_phi(BoutReal UNUSED(t), BoutReal UNUSED(cj), BoutReal UNUSED(delta)) {
  ddt(phi) = phiSolver->solve(C_phi - ddt(U));
  return 0;
}
