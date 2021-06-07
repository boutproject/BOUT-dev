/*******************************************************************************
 * High-Beta Flute-Reduced MHD with 6-field of (N_i, T_e, T_i, U, Psi, Vipar)
 * This is a simple model for BOUT++ workshop 2013 hand-on exercise
 * Basically the same as Hazeltine-Meiss but different normalisations.
 * diffusion_par can turn on the parallel thermal conductivity
 * T. Xia
 *******************************************************************************/

#include "bout.hxx"
#include "derivs.hxx"
#include "initialprofiles.hxx"
#include "interpolation_xz.hxx"
#include "invert_laplace.hxx"
#include "invert_parderiv.hxx"
#include "msg_stack.hxx"
#include "sourcex.hxx"
#include "bout/constants.hxx"
#include "bout/physicsmodel.hxx"

#include <cmath>

constexpr BoutReal eV_K = 11605.0; // 1eV = 11605K

class Elm_6f : public PhysicsModel {

  /// The total height, average width and center of profile of N0
  BoutReal n0_height, n0_ave, n0_width, n0_center, n0_bottom_x;
  /// The ampitude of congstant temperature
  BoutReal Tconst;

  /// Test the effect of first order term of invert Laplace function
  BoutReal laplace_alpha;
  /// The ratio of Ti0/Te0
  BoutReal Tau_ie;

  // 2D inital profiles
  /// Current and pressure
  Field2D J0, P0;
  /// Curvature term
  Vector2D b0xcv;
  /// When diamagnetic terms used
  Field2D phi0;

  /// Number density and temperature for ions and electrons
  Field2D N0, Ti0, Te0, Ne0;
  Field2D Pi0, Pe0;
  Field2D q95;
  BoutReal q95_input;
  bool local_q;
  BoutReal q_alpha;
  bool n0_fake_prof, T0_fake_prof;
  /// Charge number of ion
  BoutReal Zi;

  /// B0 field vector
  Vector2D B0vec;

  // V0 field vectors
  /// V0 field vector in convection
  Vector2D V0vec;
  /// Effective V0 field vector in Ohm's law
  Vector2D V0eff;

  // 3D evolving variables
  Field3D U, Psi, P, Pi, Pe;
  Field3D Ni, Te, Ti, Ne;
  Field3D Vipar, Vepar;

  // Derived 3D variables
  Field3D Jpar, phi; // Parallel current, electric potential

  Field3D Ajpar; // Parallel current, electric potential
  bool emass;
  BoutReal emass_inv; // inverse of electron mass
  BoutReal coef_jpar;
  BoutReal delta_e;     // Normalized electron skin depth
  BoutReal delta_e_inv; // inverse normalized electron skin depth
  BoutReal gyroAlv;     // Normalized ion current coef
  Field3D ubyn;

  Field3D Jpar2;                         //  Delp2 of Parallel current
  Field3D tmpA2;                         // Grad2_par2new of Parallel vector potential
  Field3D tmpN2, tmpTi2, tmpTe2, tmpVp2; // Grad2_par2new of Parallel density

  // Constraint
  Field3D C_phi;

  // Parameters
  BoutReal density;              // Number density [m^-3]
  BoutReal Bbar, Lbar, Tbar, Va; // Normalisation constants
  BoutReal Nbar, Tibar, Tebar;
  BoutReal dia_fact; // Multiply diamagnetic term by this

  BoutReal diffusion_par;  // Parallel thermal conductivity
  BoutReal diffusion_perp; // Perpendicular thermal conductivity (>0 open)
  BoutReal diffusion_n4, diffusion_ti4,
      diffusion_te4; // M: 4th Parallel density diffusion
  BoutReal diffusion_v4;
  BoutReal diffusion_u4; // xqx: parallel hyper-viscous diffusion for vorticity

  BoutReal heating_P; // heating power in pressure
  BoutReal hp_width;  // heating profile radial width in pressure
  BoutReal hp_length; // heating radial domain in pressure
  BoutReal sink_vp;   // sink in pressure
  BoutReal sp_width;  // sink profile radial width in pressure
  BoutReal sp_length; // sink radial domain in pressure

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

  Field3D Dperp2Phi0, Dperp2Phi, GradPhi02,
      GradPhi2; // Temporary variables for gyroviscous
  Field3D GradparPhi02, GradparPhi2, GradcPhi, GradcparPhi;
  Field3D Dperp2Pi0, Dperp2Pi, bracketPhi0P, bracketPhiP0, bracketPhiP;

  BoutReal Psipara1, Upara0,
      Upara1; // Temporary normalization constants for all the equations
  BoutReal Upara2, Upara3, Nipara1;
  BoutReal Tipara1, Tipara2;
  BoutReal Tepara1, Tepara2, Tepara3, Tepara4;
  BoutReal Vepara, Vipara;
  BoutReal Low_limit; // To limit the negative value of total density and temperatures

  Field3D Te_tmp, Ti_tmp, N_tmp;   // to avoid the negative value of total value
  BoutReal gamma_i_BC, gamma_e_BC; // sheath energy transmission factors
  int Sheath_width;
  Field3D c_se, Jpar_sh, q_se, q_si, vth_et,
      c_set; // variables for sheath boundary conditions
  Field2D vth_e0, c_se0, Jpar_sh0;
  BoutReal const_cse;

  // options
  bool include_curvature, include_jpar0, compress0;
  bool evolve_pressure, continuity, gyroviscous;
  Field3D diff_radial, ddx_ni, ddx_n0;
  BoutReal diffusion_coef_Hmode0, diffusion_coef_Hmode1;
  Field3D eta_i0, pi_ci;

  BoutReal vacuum_pressure;
  BoutReal vacuum_trans; // Transition width
  Field3D vac_mask;

  bool nonlinear;
  bool evolve_jpar;
  BoutReal g; // Only if compressible
  bool phi_curv;

  // Poisson brackets: b0 x Grad(f) dot Grad(g) / B = [f, g]
  // Method to use: BRACKET_ARAKAWA, BRACKET_STD or BRACKET_SIMPLE
  BRACKET_METHOD bm_exb, bm_mag; // Bracket method for advection terms
  int bracket_method_exb, bracket_method_mag;

  bool diamag;
  bool energy_flux, energy_exch; // energy flux term
  bool diamag_phi0;              // Include the diamagnetic equilibrium phi0
  bool thermal_force;            // Include the thermal flux term in Ohm's law
  bool eHall;
  BoutReal AA; // ion mass in units of the proton mass; AA=Mi/Mp

  BoutReal Vt0; // equilibrium toroidal flow normalized to Alfven velocity
  BoutReal Vp0; // equilibrium poloidal flow normalized to Alfven velocity

  bool experiment_Er; // read in phi_0 from experiment
  Field2D V0, Dphi0;  // net flow amplitude, differential potential to flux
  Vector2D V0net;     // net flow

  bool nogradparj;
  bool filter_z;
  int filter_z_mode;
  int low_pass_z;
  bool zonal_flow;
  bool zonal_field;
  bool zonal_bkgd;
  bool relax_j_vac;
  BoutReal relax_j_tconst; // Time-constant for j relax
  Field3D Psitarget;       // The (moving) target to relax to

  bool smooth_j_x; // Smooth Jpar in the x direction
  BoutReal filter_nl;

  int jpar_bndry_width; // Zero jpar in a boundary region

  bool parallel_lagrange; // Use (semi-) Lagrangian method for parallel derivatives
  bool parallel_project;  // Use Apar to project field-lines

  //********************

  Field3D Xip_x, Xip_z; // Displacement of y+1 (in cell index space)

  Field3D Xim_x, Xim_z; // Displacement of y-1 (in cell index space)

  bool phi_constraint; // Solver for phi using a solver constraint

  BoutReal vac_lund, core_lund;     // Lundquist number S = (Tau_R / Tau_A). -ve -> infty
  BoutReal vac_resist, core_resist; // The resistivities (just 1 / S)
  Field3D eta;                      // Resistivity profile (1 / S)
  bool spitzer_resist;              // Use Spitzer formula for resistivity

  Field3D eta_spitzer;        // Resistivity profile (kg*m^3 / S / C^2)
  Field3D nu_i;               // Ion collision frequency profile (1 / S)
  Field3D nu_e;               // Electron collision frequency profile (1 / S)
  Field3D vth_i;              // Ion Thermal Velocity profile (M / S)
  Field3D vth_e;              // Electron Thermal Velocity profile (M / S)
  Field3D kappa_par_i;        // Ion Thermal Conductivity profile (kg&M / S^2)
  Field3D kappa_par_e;        // Electron Thermal Conductivity profile (kg*M / S^2)
  Field2D omega_ci, omega_ce; // cyclotron frequency
  Field3D kappa_perp_i; // Ion perpendicular Thermal Conductivity profile (kg&M / S^2)
  // Electron perpendicular Thermal Conductivity profile (kg*M / S^2)
  Field3D kappa_perp_e;

  bool output_transfer; // output the results of energy transfer
  bool output_ohm;      // output the results of the terms in Ohm's law
  bool output_flux_par; // output the results of parallel particle and heat flux
  // Maxwell stress, Reynolds stress, ion diamagbetic and curvature term
  Field3D T_M, T_R, T_ID, T_C, T_G;
  Field3D ohm_phi, ohm_hall, ohm_thermal;
  // particle flux, ion and elelctron heat flux
  Field3D gamma_par_i, heatf_par_i, heatf_par_e;

  BoutReal hyperresist;  // Hyper-resistivity coefficient (in core only)
  BoutReal ehyperviscos; // electron Hyper-viscosity coefficient
  Field3D hyper_eta_x;   // Radial resistivity profile
  Field3D hyper_eta_z;   // Toroidal resistivity profile

  int damp_width;        // Width of inner damped region
  BoutReal damp_t_const; // Timescale of damping

  // Metric coefficients
  Field2D Rxy, Bpxy, Btxy, B0, hthe;
  Field2D I;         // Shear factor
  BoutReal LnLambda; // ln(Lambda)

  /// Ion mass
  BoutReal Mi = SI::amu;

  /// Communication objects
  FieldGroup comms;

  /// Solver for inverting Laplacian
  std::unique_ptr<Laplacian> phiSolver{nullptr};
  std::unique_ptr<Laplacian> aparSolver{nullptr};

  /// For printing out some diagnostics first time around
  bool first_run = true;

  Field3D field_larger(const Field3D& f, const BoutReal limit) {
    Field3D result;
    result.allocate();

    for (auto i : result) {
      if (f[i] >= limit) {
        result[i] = f[i];
      } else {
        result[i] = limit;
      }
    }
    mesh->communicate(result);
    return result;
  }

  Field3D Grad2_par2new(const Field3D& f) {
    /*
     * This function implements d2/dy2 where y is the poloidal coordinate theta
     */

    TRACE("Grad2_par2new( Field3D )");

    Field3D result = D2DY2(f);

#if BOUT_USE_TRACK
    result.name = "Grad2_par2new(" + f.name + ")";
#endif

    return result;
  }

  Field2D N0tanh(BoutReal n0_height, BoutReal n0_ave, BoutReal n0_width,
                 BoutReal n0_center, BoutReal n0_bottom_x) {
    Field2D result;
    result.allocate();

    BoutReal Grid_NX, Grid_NXlimit; // the grid number on x, and the
    BoutReal Jysep;
    mesh->get(Grid_NX, "nx");
    mesh->get(Jysep, "jyseps1_1");
    Grid_NXlimit = n0_bottom_x * Grid_NX;
    output.write("Jysep1_1 = {:d}   Grid number = {:e}\n", int(Jysep), Grid_NX);

    if (Jysep > 0.) { // for single null geometry

      BoutReal Jxsep, Jysep2;
      mesh->get(Jxsep, "ixseps1");
      mesh->get(Jysep2, "jyseps2_2");

      for (int jx = 0; jx < mesh->LocalNx; jx++) {
        BoutReal mgx = mesh->GlobalX(jx);
        BoutReal xgrid_num = (Jxsep + 1.) / Grid_NX;
        // output.write("mgx = {:e} xgrid_num = {:e}\n", mgx);
        for (int jy = 0; jy < mesh->LocalNy; jy++) {
          int globaly = mesh->getGlobalYIndex(jy);
          // output.write("local y = {:d};   global y: {:d}\n", jy, globaly);
          if (mgx > xgrid_num || (globaly <= int(Jysep) - 2)
              || (globaly > int(Jysep2) + 2)) {
            mgx = xgrid_num;
          }
          BoutReal rlx = mgx - n0_center;
          BoutReal temp = exp(rlx / n0_width);
          BoutReal dampr = ((temp - 1.0 / temp) / (temp + 1.0 / temp));
          result(jx, jy) = 0.5 * (1.0 - dampr) * n0_height + n0_ave;
        }
      }
    } else { // circular geometry
      for (int jx = 0; jx < mesh->LocalNx; jx++) {
        BoutReal mgx = mesh->GlobalX(jx);
        BoutReal xgrid_num = Grid_NXlimit / Grid_NX;
        if (mgx > xgrid_num) {
          mgx = xgrid_num;
        }
        BoutReal rlx = mgx - n0_center;
        BoutReal temp = exp(rlx / n0_width);
        BoutReal dampr = ((temp - 1.0 / temp) / (temp + 1.0 / temp));
        for (int jy = 0; jy < mesh->LocalNy; jy++) {
          result(jx, jy) = 0.5 * (1.0 - dampr) * n0_height + n0_ave;
        }
      }
    }

    mesh->communicate(result);

    return result;
  }

  // Parallel gradient along perturbed field-line
  Field3D Grad_parP(const Field3D& f, CELL_LOC loc = CELL_DEFAULT) {
    TRACE("Grad_parP");

    Field3D result;

    if (parallel_lagrange || parallel_project) {
      // Moving stencil locations

      ASSERT2((not mesh->StaggerGrids) or loc == CELL_DEFAULT or loc == f.getLocation());

      Field3D fp, fm; // Interpolated on + and - y locations

      fp = interpolate(f, Xip_x, Xip_z);
      fm = interpolate(f, Xim_x, Xim_z);

      Coordinates* coord = mesh->getCoordinates();

      result.allocate();
      for (auto i : result) {
        result[i] =
            (fp[i.yp()] - fm[i.ym()]) / (2. * coord->dy[i] * sqrt(coord->g_22[i]));
      }
    } else {
      result = Grad_par(f, loc);

      if (nonlinear) {
        result -= bracket(Psi, f, bm_mag) * B0;
      }
    }

    return result;
  }

protected:
  int init(bool restarting) override {
    bool noshear;

    // Get the metric tensor
    Coordinates* coord = mesh->getCoordinates();

    output.write("Solving high-beta flute reduced equations\n");
    output.write("\tFile    : {:s}\n", __FILE__);
    output.write("\tCompiled: {:s} at {:s}\n", __DATE__, __TIME__);

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
    mesh->get(Btxy, "Btxy"); // T
    mesh->get(B0, "Bxy");    // T
    mesh->get(hthe, "hthe"); // m
    mesh->get(I, "sinty");   // m^-2 T^-1

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

    auto& globalOptions = Options::root();
    auto& options = globalOptions["highbeta"];

    n0_fake_prof = options["n0_fake_prof"]
                       .doc("use the hyperbolic profile of n0. If both  n0_fake_prof and "
                            "T0_fake_prof are "
                            "false, use the profiles from grid file")
                       .withDefault(false);
    n0_height = options["n0_height"]
                    .doc("the total height of profile of N0, in percentage of Ni_x")
                    .withDefault(0.4);
    n0_ave = options["n0_ave"]
                 .doc("the center or average of N0, in percentage of Ni_x")
                 .withDefault(0.01);
    n0_width = options["n0_width"]
                   .doc("the width of the gradient of N0,in percentage of x")
                   .withDefault(0.1);
    n0_center = options["n0_center"]
                    .doc("the grid number of the center of N0, in percentage of x")
                    .withDefault(0.633);
    n0_bottom_x =
        options["n0_bottom_x"]
            .doc("the start of flat region of N0 on SOL side, in percentage of x")
            .withDefault(0.81);
    T0_fake_prof = options["T0_fake_prof"].doc("").withDefault(false);
    Tconst = options["Tconst"]
                 .doc("the amplitude of constant temperature, in percentage")
                 .withDefault(-1.0);

    experiment_Er = options["experiment_Er"].withDefault(false);

    laplace_alpha = options["laplace_alpha"]
                        .doc("test parameter for the cross term of invert Lapalace")
                        .withDefault(1.0);
    Low_limit = options["Low_limit"]
                    .doc("limit the negative value of total quantities")
                    .withDefault(1.0e-10);
    q95_input = options["q95_input"]
                    .doc("input q95 as a constant, if <0 use profile from grid")
                    .withDefault(5.0);
    local_q = options["local_q"]
                  .doc("Using magnetic field to calculate q profile?")
                  .withDefault(false);
    q_alpha = options["q_alpha"]
                  .doc("flux-limiting coefficient, typical value is [0.03, 3]")
                  .withDefault(1.0);

    gamma_i_BC = options["gamma_i_BC"]
                     .doc("sheath energy transmission factor for ion")
                     .withDefault(-1.0);
    gamma_e_BC = options["gamma_e_BC"]
                     .doc("sheath energy transmission factor for electron")
                     .withDefault(-1.0);
    Sheath_width = options["Sheath_width"]
                       .doc("Sheath boundary width in grid number")
                       .withDefault(1);

    density = options["density"].doc("Number density [m^-3]").withDefault(1.0e19);
    Zi = options["Zi"].doc("ion charge number").withDefault(1);
    continuity = options["continuity"].doc("use continuity equation").withDefault(false);

    evolve_jpar =
        options["evolve_jpar"].doc("If true, evolve J raher than Psi").withDefault(false);
    phi_constraint =
        options["phi_constraint"].doc("Use solver constraint for phi").withDefault(false);

    // Effects to include/exclude
    include_curvature = options["include_curvature"].withDefault(true);
    include_jpar0 = options["include_jpar0"].withDefault(true);
    evolve_pressure = options["evolve_pressure"].withDefault(true);

    compress0 = options["compress0"].withDefault(false);
    nonlinear = options["nonlinear"].withDefault(false);

    //  int bracket_method;
    bracket_method_exb = options["bracket_method_exb"].withDefault(0);
    switch (bracket_method_exb) {
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

    //  int bracket_method;
    bracket_method_mag = options["bracket_method_mag"].withDefault(2);
    switch (bracket_method_mag) {
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

    AA = options["AA"].doc("ion mass in units of proton mass").withDefault(1.0);
    Mi *= AA;

    emass = options["emass"]
                .doc("including electron inertial, electron mass")
                .withDefault(false);
    emass_inv = options["emass_inv"].doc("inverse of electron mass").withDefault(1.0);

    diamag = options["diamag"].doc("Diamagnetic effects?").withDefault(false);
    diamag_phi0 =
        options["diamag_phi0"].doc("Include equilibrium phi0").withDefault(diamag);
    dia_fact = options["dia_fact"]
                   .doc("Scale diamagnetic effects by this factor")
                   .withDefault(1.0);

    noshear = options["noshear"].withDefault(false);

    relax_j_vac =
        options["relax_j_vac"].doc("Relax vacuum current to zero").withDefault(false);
    relax_j_tconst = options["relax_j_tconst"].withDefault(0.1);

    // Toroidal filtering
    filter_z = options["filter_z"].doc("Filter a single n").withDefault(false);
    filter_z_mode = options["filter_z_mode"].withDefault(1);
    low_pass_z = options["low_pass_z"].doc("Low pass filter. < 0 -> off").withDefault(-1);
    zonal_flow = options["zonal_flow"].withDefault(false);   // zonal flow filter
    zonal_field = options["zonal_field"].withDefault(false); // zonal field filter
    zonal_bkgd = options["zonal_bkgd"].withDefault(false);   // zonal background P filter

    filter_nl = options["filter_nl"].doc("zonal background P filter").withDefault(-1);

    // Radial smoothing
    smooth_j_x = options["smooth_j_x"].doc("Smooth Jpar in x").withDefault(false);

    // Jpar boundary region
    jpar_bndry_width = options["jpar_bndry_width"].withDefault(-1);

    // Parallel differencing
    parallel_lagrange = options["parallel_lagrange"]
                            .doc("Use a (semi-) Lagrangian method for Grad_parP")
                            .withDefault(false);
    parallel_project = options["parallel_project"].withDefault(false);

    // Vacuum region control
    vacuum_pressure =
        options["vacuum_pressure"].doc("Fraction of peak pressure").withDefault(0.02);
    vacuum_trans =
        options["vacuum_trans"].doc("Transition width in pressure").withDefault(0.005);

    // Resistivity and hyper-resistivity options
    vac_lund =
        options["vac_lund"].doc("Lundquist number in vacuum region").withDefault(0.0);
    core_lund =
        options["core_lund"].doc("Lundquist number in core region").withDefault(0.0);
    hyperresist = options["hyperresist"].withDefault(-1.0);
    ehyperviscos = options["ehyperviscos"].withDefault(-1.0);
    spitzer_resist =
        options["spitzer_resist"].doc("Use Spitzer resistivity?").withDefault(false);

    // Inner boundary damping
    damp_width = options["damp_width"].withDefault(0);
    damp_t_const = options["damp_t_const"].withDefault(0.1);

    // Viscosity and hyper-viscosity
    viscos_par = options["viscos_par"].doc("Parallel viscosity").withDefault(-1.0);
    viscos_perp = options["viscos_perp"].doc("Perpendicular viscosity").withDefault(-1.0);
    hyperviscos = options["hyperviscos"].doc("Radial hyperviscosity").withDefault(-1.0);

    diffusion_par =
        options["diffusion_par"].doc("Parallel temperature diffusion").withDefault(-1.0);
    diffusion_n4 =
        options["diffusion_n4"].doc("4th Parallel density diffusion").withDefault(-1.0);
    diffusion_ti4 = options["diffusion_ti4"]
                        .doc("4th Parallel ion temperature diffusion")
                        .withDefault(-1.0);
    diffusion_te4 = options["diffusion_te4"]
                        .doc("4th Parallel electron temperature diffusion")
                        .withDefault(-1.0);
    diffusion_u4 = options["diffusion_u4"]
                       .doc("parallel hyper-viscous diffusion for vorticity")
                       .withDefault(-1.0);
    diffusion_v4 = options["diffusion_v4"]
                       .doc("4th order Parallel ion velocity diffusion (< 0 = none)")
                       .withDefault(-1.0);

    // heating factor in pressure
    heating_P = options["heating_P"].doc("heating power in pressure").withDefault(-1.0);
    //
    hp_width = options["hp_width"]
                   .doc("the percentage of radial grid points for heating profile radial "
                        "width in pressure")
                   .withDefault(0.1);
    hp_length = options["hp_length"]
                    .doc("the percentage of radial grid points for heating profile "
                         "radial domain in pressure")
                    .withDefault(0.04);

    // sink factor in pressure
    sink_vp = options["sink_vp"].doc("sink in pressure").withDefault(-1.0);
    sp_width = options["sp_width"]
                   .doc("the percentage of radial grid points for sink profile radial "
                        "width in pressure")
                   .withDefault(0.05);
    sp_length = options["sp_length"]
                    .doc("the percentage of radial grid points for sink profile radial "
                         "domain in pressure")
                    .withDefault(0.04);

    // left edge sink factor in vorticity
    sink_Ul = options["sink_Ul"].doc("left edge sink in vorticity").withDefault(-1.0);
    su_widthl = options["su_widthl"]
                    .doc("the percentage of left edge radial grid points for sink "
                         "profile radial width in vorticity")
                    .withDefault(0.06);
    su_lengthl = options["su_lengthl"]
                     .doc("the percentage of left edge radial grid points for sink "
                          "profile radial domain in vorticity")
                     .withDefault(0.15);

    // right edge sink factor in vorticity
    // right edge sink in vorticity
    sink_Ur = options["sink_Ur"].doc("").withDefault(-1.0);
    // the percentage of right edge radial grid points for sink profile
    // radial width in vorticity
    su_widthr = options["su_widthr"].doc("").withDefault(0.06);
    // the percentage of right edge radial grid points for sink profile
    // radial domain in vorticity
    su_lengthr = options["su_lengthr"].doc("").withDefault(0.15);

    // Compressional terms
    phi_curv = options["phi_curv"].doc("Compressional ExB terms").withDefault(true);
    g = options["gamma"].doc("Ratio of specific heats").withDefault(5.0 / 3.0);

    if (!include_curvature) {
      b0xcv = 0.0;
    }

    if (!include_jpar0) {
      J0 = 0.0;
    }

    if (noshear) {
      if (include_curvature) {
        b0xcv.z += I * b0xcv.x;
      }
      I = 0.0;
    }

    //////////////////////////////////////////////////////////////
    // SHIFTED RADIAL COORDINATES

    if (mesh->IncIntShear) {
      // BOUT-06 style, using d/dx = d/dpsi + I * d/dz
      coord->IntShiftTorsion = I;

    } else {
      // Dimits style, using local coordinate system
      if (include_curvature) {
        b0xcv.z += I * b0xcv.x;
      }
      I = 0.0; // I disappears from metric
    }

    //////////////////////////////////////////////////////////////
    // NORMALISE QUANTITIES

    if (mesh->get(Bbar, "bmag")) { // Typical magnetic field
      Bbar = 1.0;
    }
    if (mesh->get(Lbar, "rmag")) { // Typical length scale
      Lbar = 1.0;
    }

    if (mesh->get(Tibar, "Ti_x")) { // Typical ion temperature scale
      Tibar = 1.0;
    }

    if (mesh->get(Tebar, "Te_x")) { // Typical electron temperature scale
      Tebar = 1.0;
    }

    if (mesh->get(Nbar, "Nixexp")) { // Typical ion density scale
      Nbar = 1.0;
    }
    Nbar *= 1.e20 / density;

    Tau_ie = Tibar / Tebar;

    Va = sqrt(Bbar * Bbar / (SI::mu0 * Mi * Nbar * density));

    Tbar = Lbar / Va;

    output.write("Normalisations: Bbar = {:e} T   Lbar = {:e} m\n", Bbar, Lbar);
    output.write("                Va = {:e} m/s   Tbar = {:e} s\n", Va, Tbar);
    output.write("                Nbar = {:e} * {:e} m^-3\n", Nbar, density);
    output.write("Tibar = {:e} eV   Tebar = {:e} eV    Ti/Te = {:e}\n", Tibar, Tebar,
                 Tau_ie);
    output.write("    Resistivity\n");

    Upara0 = SI::kb * Tebar * eV_K / (Zi * SI::qe * Bbar * Va * Lbar);
    Upara1 = SI::kb * Tebar * eV_K / Mi / Va / Va;
    output.write("vorticity cinstant: Upara0 = {:e}     Upara1 = {:e}\n", Upara0, Upara1);

    if (diamag) {
      Nipara1 = SI::kb * Tibar * eV_K / (Zi * SI::qe * Bbar * Lbar * Va);
      Tipara2 = Nipara1;
      Tepara2 = SI::kb * Tebar * eV_K / (SI::qe * Bbar * Lbar * Va);
      Tepara3 = Bbar / (SI::qe * SI::mu0 * Nbar * density * Lbar * Va);
      output.write("Nipara1 = {:e}     Tipara2 = {:e}\n", Nipara1, Tipara2);
      output.write("Tepara2 = {:e}     Tepara3 = {:e}\n", Tepara2, Tepara3);
    }

    if (compress0) {
      output.write("Including compression (Vipar) effects\n");
      Vipara = SI::mu0 * SI::kb * Nbar * density * Tebar * eV_K / (Bbar * Bbar);
      Vepara = Bbar / (SI::mu0 * Zi * SI::qe * Nbar * density * Lbar * Va);
      output.write("Normalized constant for Vipar :   Vipara = {:e}\n", Vipara);
      output.write("Normalized constant for Vepar :   Vepara = {:e}\n", Vepara);
    }

    if (diffusion_par > 0.0) {
      Tipara1 = 2.0 / 3.0 / (Lbar * Va);
      Tepara1 = Tipara1 / Zi;
    }

    if (vac_lund > 0.0) {
      output.write("        Vacuum  Tau_R = {:e} s   eta = {:e} Ohm m\n", vac_lund * Tbar,
                   SI::mu0 * Lbar * Lbar / (vac_lund * Tbar));
      vac_resist = 1. / vac_lund;
    } else {
      output.write("        Vacuum  - Zero resistivity -\n");
      vac_resist = 0.0;
    }
    if (core_lund > 0.0) {
      output.write("        Core    Tau_R = {:e} s   eta = {:e} Ohm m\n",
                   core_lund * Tbar, SI::mu0 * Lbar * Lbar / (core_lund * Tbar));
      core_resist = 1. / core_lund;
    } else {
      output.write("        Core    - Zero resistivity -\n");
      core_resist = 0.0;
    }

    if (hyperresist > 0.0) {
      output.write("    Hyper-resistivity coefficient: {:e}\n", hyperresist);
      dump.add(hyper_eta_x, "hyper_eta_x", 1);
      dump.add(hyper_eta_z, "hyper_eta_z", 1);
    }

    if (ehyperviscos > 0.0) {
      output.write("    electron Hyper-viscosity coefficient: {:e}\n", ehyperviscos);
    }

    if (hyperviscos > 0.0) {
      output.write("    Hyper-viscosity coefficient: {:e}\n", hyperviscos);
      dump.add(hyper_mu_x, "hyper_mu_x", 1);
    }

    if (diffusion_par > 0.0) {
      output.write("    diffusion_par: {:e}\n", diffusion_par);
      dump.add(diffusion_par, "diffusion_par", 0);
    }

    // M: 4th order diffusion of p
    if (diffusion_n4 > 0.0) {
      output.write("    diffusion_n4: {:e}\n", diffusion_n4);
      dump.add(diffusion_n4, "diffusion_n4", 0);
    }

    // M: 4th order diffusion of Ti
    if (diffusion_ti4 > 0.0) {
      output.write("    diffusion_ti4: {:e}\n", diffusion_ti4);
      dump.add(diffusion_ti4, "diffusion_ti4", 0);
    }

    // M: 4th order diffusion of Te
    if (diffusion_te4 > 0.0) {
      output.write("    diffusion_te4: {:e}\n", diffusion_te4);
      dump.add(diffusion_te4, "diffusion_te4", 0);
    }

    // M: 4th order diffusion of Vipar
    if (diffusion_v4 > 0.0) {
      output.write("    diffusion_v4: {:e}\n", diffusion_v4);
      dump.add(diffusion_v4, "diffusion_v4", 0);
    }

    // xqx: parallel hyper-viscous diffusion for vorticity
    if (diffusion_u4 > 0.0) {
      output.write("    diffusion_u4: {:e}\n", diffusion_u4);
      dump.add(diffusion_u4, "diffusion_u4", 0);
    }

    if (sink_vp > 0.0) {
      output.write("    sink_vp(rate): {:e}\n", sink_vp);
      dump.add(sink_vp, "sink_vp", 1);

      output.write("    sp_width(%%): {:e}\n", sp_width);
      dump.add(sp_width, "sp_width", 1);

      output.write("    sp_length(%%): {:e}\n", sp_length);
      dump.add(sp_length, "sp_length", 1);
    }

    J0 = SI::mu0 * Lbar * J0 / B0;
    P0 = P0 / (SI::kb * (Tibar + Tebar) * eV_K / 2. * Nbar * density);

    b0xcv.x /= Bbar;
    b0xcv.y *= Lbar * Lbar;
    b0xcv.z *= Lbar * Lbar;

    Rxy /= Lbar;
    Bpxy /= Bbar;
    Btxy /= Bbar;
    B0 /= Bbar;
    hthe /= Lbar;
    coord->dx /= Lbar * Lbar * Bbar;
    I *= Lbar * Lbar * Bbar;

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

    Ne0 = Zi * N0; // quasi-neutral condition
    Pi0 = N0 * Ti0;
    Pe0 = Ne0 * Te0;

    nu_e.setBoundary("kappa");
    if (spitzer_resist) {
      eta_spitzer.setBoundary("kappa");
    }
    if (diffusion_par > 0.0) {
      nu_i.setBoundary("kappa");
      vth_i.setBoundary("kappa");
      vth_e.setBoundary("kappa");
      kappa_par_i.setBoundary("kappa");
      kappa_par_e.setBoundary("kappa");
      kappa_perp_i.setBoundary("kappa");
      kappa_perp_e.setBoundary("kappa");
    }

    if (compress0) {
      eta_i0.setBoundary("Ti");
      pi_ci.setBoundary("Ti");
    }

    BoutReal pnorm = max(P0, true); // Maximum over all processors

    vacuum_pressure *= pnorm; // Get pressure from fraction
    vacuum_trans *= pnorm;

    // Transitions from 0 in core to 1 in vacuum
    vac_mask = (1.0 - tanh((P0 - vacuum_pressure) / vacuum_trans)) / 2.0;

    if (diffusion_par > 0.0) {
      if (q95_input > 0) {
        q95 = q95_input; // use a constant for test
      } else {
        if (local_q) {
          q95 = abs(hthe * Btxy / (Bpxy)) * q_alpha;
        } else {
          output.write("\tUsing q profile from grid.\n");
          if (mesh->get(q95, "q")) {
            output.write(
                "Cannot get q profile from grid!\nPlease run addqprofile.pro first\n");
            return 1;
          }
        }
      }
      output.write("\tlocal max q: {:e}\n", max(q95));
      output.write("\tlocal min q: {:e}\n", min(q95));
    }

    LnLambda =
        24.0
        - log(pow(Zi * Nbar * density / 1.e6, 0.5) * pow(Tebar, -1.0)); // xia: ln Lambda
    output.write("\tlog Lambda: {:e}\n", LnLambda);

    nu_e = 2.91e-6 * LnLambda * ((N0)*Nbar * density / 1.e6)
           * pow(Te0 * Tebar, -1.5); // nu_e in 1/S.
    output.write("\telectron collision rate: {:e} -> {:e} [1/s]\n", min(nu_e), max(nu_e));
    // nu_e.applyBoundary();
    // mesh->communicate(nu_e);

    if (diffusion_par > 0.0) {

      output.write("\tion thermal noramlized constant: Tipara1 = {:e}\n", Tipara1);
      output.write("\telectron normalized thermal constant: Tepara1 = {:e}\n", Tepara1);
      // xqx addition, begin
      // Use Spitzer thermal conductivities
      nu_i = 4.80e-8 * (Zi * Zi * Zi * Zi / sqrt(AA)) * LnLambda
             * ((N0)*Nbar * density / 1.e6) * pow(Ti0 * Tibar, -1.5); // nu_i in 1/S.
      // output.write("\tCoulomb Logarithm: {:e} \n", max(LnLambda));
      output.write("\tion collision rate: {:e} -> {:e} [1/s]\n", min(nu_i), max(nu_i));

      // nu_i.applyBoundary();
      // mesh->communicate(nu_i);

      vth_i = 9.79e3 * sqrt((Ti0)*Tibar / AA); // vth_i in m/S.
      output.write("\tion thermal velocity: {:e} -> {:e} [m/s]\n", min(vth_i),
                   max(vth_i));
      // vth_i.applyBoundary();
      // mesh->communicate(vth_i);
      vth_e = 4.19e5 * sqrt((Te0)*Tebar); // vth_e in m/S.
      output.write("\telectron thermal velocity: {:e} -> {:e} [m/s]\n", min(vth_e),
                   max(vth_e));
      // vth_e.applyBoundary();
      // mesh->communicate(vth_e);
    }

    if (compress0) {
      eta_i0 = 0.96 * Pi0 * Tau_ie * nu_i * Tbar;
      output.write("\tCoefficients of parallel viscocity: {:e} -> {:e} [kg/(m s)]\n",
                   min(eta_i0), max(eta_i0));
    }

    if (diffusion_par > 0.0) {
      kappa_par_i = 3.9 * vth_i * vth_i / nu_i; // * 1.e4;
      kappa_par_e = 3.2 * vth_e * vth_e / nu_e; // * 1.e4;

      output.write("\tion thermal conductivity: {:e} -> {:e} [m^2/s]\n", min(kappa_par_i),
                   max(kappa_par_i));
      output.write("\telectron thermal conductivity: {:e} -> {:e} [m^2/s]\n",
                   min(kappa_par_e), max(kappa_par_e));

      output.write("\tnormalized ion thermal conductivity: {:e} -> {:e} \n",
                   min(kappa_par_i * Tipara1), max(kappa_par_i * Tipara1));
      output.write("\tnormalized electron thermal conductivity: {:e} -> {:e} \n",
                   min(kappa_par_e * Tepara1), max(kappa_par_e * Tepara1));

      Field3D kappa_par_i_fl, kappa_par_e_fl;

      kappa_par_i_fl = vth_i * (q95 * Lbar); // * 1.e2;
      kappa_par_e_fl = vth_e * (q95 * Lbar); // * 1.e2;

      kappa_par_i *= kappa_par_i_fl / (kappa_par_i + kappa_par_i_fl);
      kappa_par_i *= Tipara1 * N0;
      output.write("\tUsed normalized ion thermal conductivity: {:e} -> {:e} \n",
                   min(kappa_par_i), max(kappa_par_i));
      // kappa_par_i.applyBoundary();
      // mesh->communicate(kappa_par_i);
      kappa_par_e *= kappa_par_e_fl / (kappa_par_e + kappa_par_e_fl);
      kappa_par_e *= Tepara1 * N0 / Zi;
      output.write("\tUsed normalized electron thermal conductivity: {:e} -> {:e} \n",
                   min(kappa_par_e), max(kappa_par_e));
      // kappa_par_e.applyBoundary();
      // mesh->communicate(kappa_par_e);

      dump.add(kappa_par_i, "kappa_par_i", 1);
      dump.add(kappa_par_e, "kappa_par_e", 1);
    }

    if (spitzer_resist) {
      // Use Spitzer resistivity
      output.write("\n\tSpizter parameters");
      // output.write("\tTemperature: {:e} -> {:e} [eV]\n", min(Te), max(Te));
      eta_spitzer = 0.51 * 1.03e-4 * Zi * LnLambda
                    * pow(Te0 * Tebar, -1.5); // eta in Ohm-m. NOTE: ln(Lambda) = 20
      output.write("\tSpitzer resistivity: {:e} -> {:e} [Ohm m]\n", min(eta_spitzer),
                   max(eta_spitzer));
      eta_spitzer /= SI::mu0 * Va * Lbar;
      // eta_spitzer.applyBoundary();
      // mesh->communicate(eta_spitzer);
      output.write("\t -> Lundquist {:e} -> {:e}\n", 1.0 / max(eta_spitzer),
                   1.0 / min(eta_spitzer));
      dump.add(eta_spitzer, "eta_spitzer", 1);
    } else {
      // transition from 0 for large P0 to resistivity for small P0
      eta = core_resist + (vac_resist - core_resist) * vac_mask;
      eta_spitzer = 0.;
      dump.add(eta, "eta", 0);
    }

    /**************** CALCULATE METRICS ******************/

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

    // Set V0vec field vector

    V0vec.covariant = false;
    V0vec.x = 0.;
    V0vec.y = Vp0 / hthe;
    V0vec.z = Vt0 / Rxy;

    // Set V0eff field vector

    V0eff.covariant = false;
    V0eff.x = 0.;
    V0eff.y = -(Btxy / (B0 * B0)) * (Vp0 * Btxy - Vt0 * Bpxy) / hthe;
    V0eff.z = (Bpxy / (B0 * B0)) * (Vp0 * Btxy - Vt0 * Bpxy) / Rxy;

    Pe.setBoundary("P");
    Pi.setBoundary("P");

    /**************** SET EVOLVING VARIABLES *************/

    // Tell BOUT which variables to evolve
    SOLVE_FOR(U);
    SOLVE_FOR(Ni);
    SOLVE_FOR(Ti);
    SOLVE_FOR(Te);
    SOLVE_FOR(Psi);

    SAVE_REPEAT(Jpar, P, Vepar);

    if (parallel_lagrange) {
      // Evolving the distortion of the flux surfaces (Ideal-MHD only!)

      solver->add(Xip_x, "Xip_x");
      solver->add(Xip_z, "Xip_z");

      solver->add(Xim_x, "Xim_x");
      solver->add(Xim_z, "Xim_z");
    }

    if (parallel_project) {
      // Add Xi to the dump file
      SAVE_REPEAT(Xip_x, Xip_z, Xim_x, Xim_z);
    }

    if (compress0) {
      SOLVE_FOR(Vipar);
      if (!restarting) {
        Vipar = 0.0;
      }
    }

    if (phi_constraint) {
      // Implicit Phi solve using IDA

      solver->constraint(phi, C_phi, "phi");

    } else {
      // Phi solved in RHS (explicitly)
      SAVE_REPEAT(phi);
    }

    // Diamagnetic phi0
    if (diamag && diamag_phi0) {
      if (experiment_Er) { // get phi0 from grid file
        mesh->get(phi0, "Phi_0");
        phi0 /= B0 * Lbar * Va;
      } else {
        // Stationary equilibrium plasma. ExB velocity balances diamagnetic drift
        phi0 = -Upara0 * Pi0 / B0 / N0;
      }
      SAVE_ONCE(phi0);
    }

    // Add some equilibrium quantities and normalisations
    // everything needed to recover physical units
    SAVE_ONCE(J0, P0);
    SAVE_ONCE(density, Lbar, Bbar, Tbar);
    SAVE_ONCE(Tibar, Tebar, Nbar);
    SAVE_ONCE(Va, B0);
    SAVE_ONCE(Ti0, Te0, N0);

    // Create a solver for the Laplacian
    phiSolver = Laplacian::create(&options["phiSolver"]);

    aparSolver = Laplacian::create(&options["aparSolver"]);

    /////////////// CHECK VACUUM ///////////////////////
    // In vacuum region, initial vorticity should equal zero

    ubyn.setBoundary("U");

    if (!restarting) {
      // Only if not restarting: Check initial perturbation

      // Set U to zero where P0 < vacuum_pressure
      U = where(P0 - vacuum_pressure, U, 0.0);

      //    Field2D lap_temp = 0.0;
      Field2D logn0 = laplace_alpha * N0;
      Field3D Ntemp;
      Ntemp = N0;
      ubyn = U * B0 / Ntemp;
      // Phi should be consistent with U
      if (laplace_alpha <= 0.0) {
        phi = phiSolver->solve(ubyn) / B0;
      } else {
        phiSolver->setCoefC(logn0);
        phi = phiSolver->solve(ubyn) / B0;
      }
    }

    /************** SETUP COMMUNICATIONS **************/

    comms.add(U);
    comms.add(Ni);
    comms.add(Ti);
    comms.add(Te);
    if (!emass) {
      comms.add(Psi);
    } else {
      comms.add(Ajpar);
    }

    if (compress0) {
      comms.add(Vipar);
      Vepar.setBoundary("Vipar");
    }

    if (diffusion_u4 > 0.0) {
      tmpA2.setBoundary("J");
    }

    if (diffusion_n4 > 0.0) {
      tmpN2.setBoundary("Ni");
    }

    if (diffusion_ti4 > 0.0) {
      tmpTi2.setBoundary("Ti");
    }

    if (diffusion_te4 > 0.0) {
      tmpTe2.setBoundary("Te");
    }

    if (diffusion_v4 > 0.0) {
      tmpVp2.setBoundary("Vipar");
    }

    phi.setBoundary("phi"); // Set boundary conditions

    P.setBoundary("P");
    Jpar.setBoundary("J");
    Jpar2.setBoundary("J");

    return 0;
  }
  int rhs(BoutReal UNUSED(t)) override {

    Coordinates* coord = mesh->getCoordinates();

    // Perform communications
    mesh->communicate(comms);

    // Inversion
    Pi = Ni * Ti0 + N0 * Ti;
    if (nonlinear) {
      Pi += Ni * Ti;
    }
    mesh->communicate(Pi);

    Pe = Zi * (Ni * Te0 + N0 * Te);
    if (nonlinear) {
      Pe += Zi * Ni * Te;
    }
    mesh->communicate(Pe);

    P = Tau_ie * Pi + Pe;
    mesh->communicate(P);

    //  Field2D lap_temp=0.0;
    Field2D logn0 = laplace_alpha * N0;
    ubyn = U * B0 / N0;
    if (diamag) {
      ubyn -= Upara0 / N0 * Delp2(Pi) / B0;
      mesh->communicate(ubyn);
      ubyn.applyBoundary();
    }
    // Invert laplacian for phi
    if (laplace_alpha > 0.0) {
      phiSolver->setCoefC(logn0);
    }
    phi = phiSolver->solve(ubyn) / B0;

    mesh->communicate(phi);

    if (emass) {
      Field2D acoeff = -delta_e_inv * N0 * N0;
      if (compress0) {
        Psi = aparSolver->solve(acoeff * Ajpar - gyroAlv * Vipar);
      } else {
        Psi = aparSolver->solve(acoeff * Ajpar);
      }
      mesh->communicate(Psi);
    }

    BoutReal N_tmp1;
    N_tmp1 = Low_limit;
    N_tmp = field_larger(N0 + Ni, N_tmp1);

    BoutReal Te_tmp1, Ti_tmp1;
    Te_tmp1 = Low_limit;
    Ti_tmp1 = Low_limit;

    Ti_tmp = field_larger(Ti0 + Ti, Ti_tmp1);
    Te_tmp = field_larger(Te0 + Te, Te_tmp1);

    // vac_mask transitions from 0 in core to 1 in vacuum
    if (nonlinear) {
      vac_mask = (1.0 - tanh(((P0 + P) - vacuum_pressure) / vacuum_trans)) / 2.0;
      // Update resistivity
      if (spitzer_resist) {
        // Use Spitzer formula
        eta_spitzer = 0.51 * 1.03e-4 * Zi * LnLambda
                      * pow(Te_tmp * Tebar, -1.5); // eta in Ohm-m. ln(Lambda) = 20
        eta_spitzer /= SI::mu0 * Va * Lbar;
      } else {
        eta = core_resist + (vac_resist - core_resist) * vac_mask;
      }

      nu_e = 2.91e-6 * LnLambda * (N_tmp * Nbar * density / 1.e6)
             * pow(Te_tmp * Tebar, -1.5); // nu_e in 1/S.

      if (diffusion_par > 0.0) {
        // Use Spitzer thermal conductivities

        nu_i = 4.80e-8 * (Zi * Zi * Zi * Zi / sqrt(AA)) * LnLambda
               * (N_tmp * Nbar * density / 1.e6)
               * pow(Ti_tmp * Tibar, -1.5);         // nu_i in 1/S.
        vth_i = 9.79e3 * sqrt(Ti_tmp * Tibar / AA); // vth_i in m/S.
        vth_e = 4.19e5 * sqrt(Te_tmp * Tebar);      // vth_e in m/S.
      }

      if (diffusion_par > 0.0) {
        kappa_par_i = 3.9 * vth_i * vth_i / nu_i; // * 1.e4;
        kappa_par_e = 3.2 * vth_e * vth_e / nu_e; // * 1.e4;

        Field3D kappa_par_i_fl, kappa_par_e_fl;

        kappa_par_i_fl = vth_i * (q95 * Lbar); // * 1.e2;
        kappa_par_e_fl = vth_e * (q95 * Lbar); // * 1.e2;

        kappa_par_i *= kappa_par_i_fl / (kappa_par_i + kappa_par_i_fl);
        kappa_par_i *= Tipara1 * N_tmp;
        kappa_par_e *= kappa_par_e_fl / (kappa_par_e + kappa_par_e_fl);
        kappa_par_e *= Tepara1 * N_tmp * Zi;
      }
    }

    Jpar = -Delp2(Psi);
    Jpar.applyBoundary();
    mesh->communicate(Jpar);

    if (jpar_bndry_width > 0) {
      // Zero j in boundary regions. Prevents vorticity drive
      // at the boundary

      for (int i = 0; i < jpar_bndry_width; i++) {
        for (int j = 0; j < mesh->LocalNy; j++) {
          for (int k = 0; k < mesh->LocalNz; k++) {
            if (mesh->firstX()) {
              Jpar(i, j, k) = 0.0;
            }
            if (mesh->lastX()) {
              Jpar(mesh->LocalNx - 1 - i, j, k) = 0.0;
            }
          }
        }
      }
    }

    // Smooth j in x
    if (smooth_j_x) {
      Jpar = smooth_x(Jpar);
    }

    if (compress0) {
      if (nonlinear) {
        Vepar = Vipar - B0 * (Jpar) / N_tmp * Vepara;
      } else {
        Vepar = Vipar - B0 * (Jpar) / N0 * Vepara;
        Vepar.applyBoundary();
        mesh->communicate(Vepar);
      }
    }

    // Get Delp2(J) from J
    Jpar2 = -Delp2(Jpar);

    Jpar2.applyBoundary();
    mesh->communicate(Jpar2);

    if (jpar_bndry_width > 0) {
      // Zero jpar2 in boundary regions. Prevents vorticity drive
      // at the boundary

      for (int i = 0; i < jpar_bndry_width; i++) {
        for (int j = 0; j < mesh->LocalNy; j++) {
          for (int k = 0; k < mesh->LocalNz; k++) {
            if (mesh->firstX()) {
              Jpar2(i, j, k) = 0.0;
            }
            if (mesh->lastX()) {
              Jpar2(mesh->LocalNx - 1 - i, j, k) = 0.0;
            }
          }
        }
      }
    }

    ////////////////////////////////////////////////////
    // Parallel electric field
    {
      TRACE("ddt(Psi)");

      ddt(Psi) = 0.0;

      if (spitzer_resist) {
        ddt(Psi) = -Grad_parP(B0 * phi) / B0 - eta_spitzer * Jpar;
      } else {
        ddt(Psi) = -Grad_parP(B0 * phi) / B0 - eta * Jpar;
      }

      if (diamag) {
        ddt(Psi) -= bracket(B0 * phi0, Psi, bm_exb); // Equilibrium flow
      }

      // Hyper-resistivity
      if (hyperresist > 0.0) {
        ddt(Psi) += hyperresist * Delp2(Jpar);
      }
    }

    ////////////////////////////////////////////////////
    // Vorticity equation

    {
      TRACE("ddt(U)");

      ddt(U) = 0.0;

      ddt(U) = -SQ(B0) * bracket(Psi, J0, bm_mag) * B0; // Grad j term

      ddt(U) += 2.0 * Upara1 * b0xcv * Grad(P); // curvature term

      ddt(U) += SQ(B0) * Grad_parP(Jpar); // b dot grad j

      if (diamag) {
        ddt(U) -= bracket(B0 * phi0, U, bm_exb); // Equilibrium flow
      }

      if (nonlinear) {
        ddt(U) -= bracket(B0 * phi, U, bm_exb); // Advection
      }

      // parallel hyper-viscous diffusion for vector potential
      if (diffusion_u4 > 0.0) {
        tmpA2 = Grad2_par2new(Psi);
        mesh->communicate(tmpA2);
        tmpA2.applyBoundary();
        ddt(U) -= diffusion_u4 * Grad2_par2new(tmpA2);
      }

      // Viscosity terms
      if (viscos_par > 0.0) {
        ddt(U) += viscos_par * Grad2_par2(U); // Parallel viscosity
      }

      if (hyperviscos > 0.0) {
        // Calculate coefficient.

        hyper_mu_x = hyperviscos * coord->g_11 * SQ(coord->dx)
                     * abs(coord->g11 * D2DX2(U)) / (abs(U) + 1e-3);
        hyper_mu_x.applyBoundary("dirichlet"); // Set to zero on all boundaries

        ddt(U) += hyper_mu_x * coord->g11 * D2DX2(U);

        if (first_run) {
          // Print out maximum values of viscosity used on this processor
          output.write("   Hyper-viscosity values:\n");
          output.write("      Max mu_x = {:e}, Max_DC mu_x = {:e}\n", max(hyper_mu_x),
                       max(DC(hyper_mu_x)));
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
    }

    ///////////////////////////////////////////////
    // number density equation

    {
      TRACE("ddt(Ni)");

      ddt(Ni) = 0.0;

      ddt(Ni) -= bracket(B0 * phi, N0, bm_exb);

      if (diamag) {
        ddt(Ni) -= bracket(B0 * phi0, Ni, bm_exb); // Equilibrium flow
      }

      if (nonlinear) {
        ddt(Ni) -= bracket(B0 * phi, Ni, bm_exb); // Advection
      }

      if (compress0) {
        ddt(Ni) -= N0 * B0 * Grad_parP(Vipar / B0);
      }

      // 4th order Parallel diffusion terms
      if (diffusion_n4 > 0.0) {
        tmpN2 = Grad2_par2new(Ni);
        mesh->communicate(tmpN2);
        tmpN2.applyBoundary();
        ddt(Ni) -= diffusion_n4 * Grad2_par2new(tmpN2);
      }
    }

    ///////////////////////////////////////////////
    // ion temperature equation
    {
      TRACE("ddt(Ti)");

      ddt(Ti) = 0.0;

      ddt(Ti) -= bracket(B0 * phi, Ti0, bm_exb);

      if (diamag) {
        ddt(Ti) -= bracket(phi0 * B0, Ti, bm_exb); // Equilibrium flow
      }

      if (nonlinear) {
        ddt(Ti) -= bracket(phi * B0, Ti, bm_exb); // Advection
      }

      if (compress0) {
        ddt(Ti) -= 2.0 / 3.0 * Ti0 * B0 * Grad_parP(Vipar / B0);
      }

      if (diffusion_par > 0.0) {
        ddt(Ti) += kappa_par_i * Grad2_par2(Ti) / N0; // Parallel diffusion
        ddt(Ti) += Grad_par(kappa_par_i) * Grad_par(Ti) / N0;
      }

      // 4th order Parallel diffusion terms
      if (diffusion_ti4 > 0.0) {
        tmpTi2 = Grad2_par2new(Ti);
        mesh->communicate(tmpTi2);
        tmpTi2.applyBoundary();
        ddt(Ti) -= diffusion_ti4 * Grad2_par2new(tmpTi2);
      }
    }

    ///////////////////////////////////////////////
    // electron temperature equation

    {
      TRACE("ddt(Te)");

      ddt(Te) = 0.0;

      ddt(Te) -= bracket(B0 * phi, Te0, bm_exb);

      if (diamag) {
        ddt(Te) -= bracket(B0 * phi0, Te, bm_exb); // Equilibrium flow
      }

      if (nonlinear) {
        ddt(Te) -= bracket(B0 * phi, Te, bm_exb); // Advection
      }

      if (compress0) {
        ddt(Te) -= 2.0 / 3.0 * Te0 * B0 * Grad_parP(Vepar / B0);
      }

      if (diffusion_par > 0.0) {
        ddt(Te) += kappa_par_e * Grad2_par2(Te) / N0; // Parallel diffusion
        ddt(Te) += Grad_par(kappa_par_e) * Grad_par(Te) / N0;
      }

      if (diffusion_te4 > 0.0) {
        tmpTe2 = Grad2_par2new(Te);
        mesh->communicate(tmpTe2);
        tmpTe2.applyBoundary();
        ddt(Te) -= diffusion_te4 * Grad2_par2new(tmpTe2);
      }
    }

    //////////////////////////////////////////////////////////////////////
    if (compress0) { // parallel velocity equation
      TRACE("ddt(Vipar)");

      ddt(Vipar) = 0.0;

      ddt(Vipar) -= Vipara * Grad_parP(P) / N0;
      ddt(Vipar) += Vipara * bracket(Psi, P0, bm_mag) * B0 / N0;

      if (diamag) {
        ddt(Vipar) -= bracket(B0 * phi0, Vipar, bm_exb);
      }

      if (nonlinear) {
        ddt(Vipar) -= bracket(B0 * phi, Vipar, bm_exb);
      }

      // xqx: parallel hyper-viscous diffusion for vector potential
      if (diffusion_v4 > 0.0) {
        tmpVp2 = Grad2_par2new(Vipar);
        mesh->communicate(tmpVp2);
        tmpVp2.applyBoundary();
        ddt(Vipar) -= diffusion_v4 * Grad2_par2new(tmpVp2);
      }

      if (sink_vp > 0.0) {
        Field2D V0tmp = 0.;
        ddt(Vipar) -= sink_vp * sink_tanhxl(V0tmp, Vipar, sp_width, sp_length); // sink
      }
    }

    ///////////////////////////////////////////////////////////////////////

    if (filter_z) {
      // Filter out all except filter_z_mode
      TRACE("filter_z");

      ddt(Psi) = filter(ddt(Psi), filter_z_mode);

      ddt(U) = filter(ddt(U), filter_z_mode);

      ddt(Ni) = filter(ddt(Ni), filter_z_mode);

      ddt(Ti) = filter(ddt(Ti), filter_z_mode);

      ddt(Te) = filter(ddt(Te), filter_z_mode);

      if (compress0) {
        ddt(Vipar) = filter(ddt(Vipar), filter_z_mode);
      }
    }

    ///////////////////////////////////////////////////////////////////////

    if (low_pass_z > 0) {
      // Low-pass filter, keeping n up to low_pass_z
      TRACE("low_pass_z");

      if (!emass) {
        ddt(Psi) = lowPass(ddt(Psi), low_pass_z, zonal_field);
      } else {
        ddt(Ajpar) = lowPass(ddt(Ajpar), low_pass_z, zonal_field);
      }

      ddt(U) = lowPass(ddt(U), low_pass_z, zonal_flow);

      ddt(Ti) = lowPass(ddt(Ti), low_pass_z, zonal_bkgd);
      ddt(Te) = lowPass(ddt(Te), low_pass_z, zonal_bkgd);
      ddt(Ni) = lowPass(ddt(Ni), low_pass_z, zonal_bkgd);

      if (compress0) {
        ddt(Vipar) = lowPass(ddt(Vipar), low_pass_z, zonal_bkgd);
      }
    }

    if (damp_width > 0) {
      for (int i = 0; i < damp_width; i++) {
        for (int j = 0; j < mesh->LocalNy; j++) {
          for (int k = 0; k < mesh->LocalNz; k++) {
            if (mesh->firstX()) {
              ddt(U)(i, j, k) -= U(i, j, k) / damp_t_const;
            }
            if (mesh->lastX()) {
              ddt(U)(mesh->LocalNx - 1 - i, j, k) -=
                  U(mesh->LocalNx - 1 - i, j, k) / damp_t_const;
            }
          }
        }
      }
    }

    if (filter_nl > 0) {
      TRACE("filter_nl");
      ddt(Ni) = nl_filter(ddt(Ni), filter_nl);
    }

    first_run = false;

    return 0;
  }
};

BOUTMAIN(Elm_6f)
