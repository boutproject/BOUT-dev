/*******************************************************************************
 * High-Beta Flute-Reduced MHD with 6-field of (N_i, T_e, T_i, U, Psi, Vipar)
 * Basically the same as Hazeltine-Meiss but different normalisations.
 * diffusion_par can turn on the parallel thermal conductivity
 * T.Y. Xia, Ben Zhu, Nami Li, Malamas Tsagkaridis
 *******************************************************************************/

/* merging 6f module from bout-v3 to bout_v5 by NamiLi Dec 2023
 * 1. add the neutral model with gas puffing and fixed-fraction impurity, test needed
 * 2. flux driven source
 *      fixed amplitude one, PI controller TBD
 * 3. for zonal field solver, hypre needed
 * */

// debug mode to output flags
#define DEBUG_6F 0

#include "bout/bout.hxx"
#include "bout/derivs.hxx"
#include "bout/initialprofiles.hxx"
#include "bout/interpolation_xz.hxx"
#include "bout/invert_laplace.hxx"
#include "bout/invert_parderiv.hxx"
#include "bout/msg_stack.hxx"
#include "bout/constants.hxx"
#include "bout/physicsmodel.hxx"
#include <bout/fft.hxx>
#include <bout/invert/laplacexy.hxx>
#include "bout/interpolation.hxx"
#include "bout/sourcex.hxx"
#include "bout/neutral.hxx"
#include <bout/utils.hxx>
//#include "bout/integrops.hxx"

#include <math.h>

#if BOUT_HAS_HYPRE
#include <bout/invert/laplacexy2_hypre.hxx>
#endif

#include <bout/field_factory.hxx>

CELL_LOC loc = CELL_CENTRE;

/// Set default options
/// This sets sensible defaults for when it's not set in the input
BOUT_OVERRIDE_DEFAULT_OPTION("phi:bndry_target", "neumann");
BOUT_OVERRIDE_DEFAULT_OPTION("phi:bndry_xin", "none");
BOUT_OVERRIDE_DEFAULT_OPTION("phi:bndry_xout", "none");

class Elm_6f : public PhysicsModel {
private:


  /********** Physical constants ************************/
  const BoutReal MU0 = 4.0e-7 * PI;              // [m kg s-2 A-2] permeability of free space
  const BoutReal Mp = 1.6726e-27;                // [kg] proton mass  
  const BoutReal KB = 1.38065e-23;               // Boltamann constant
  const BoutReal ee = 1.602e-19;                 
  const BoutReal ratio_pe = 1836.2;              // proton/electron mass ratio
  const BoutReal eV_K = 11605.0;                 // 1eV = 11605K

  /********** Magnetic configuration *******************************************/
  int mag_config;                                // magnetic geometry: 1-circular, 2-circular with limiter, 3-single null, 4-double null
  Field2D Rxy, Bpxy, Btxy, B0, hthe, I;          // I: shear factor
  BoutReal ixsep, ixsep2;
  BoutReal jysep1, jysep2, jysep1_2, jysep2_1;   // index for x-point on y direction
  Vector2D B0vec;                                // B0 field vector

  /********** Primary variables ************************************************/
  BoutReal AA, Zi, Mi;                           // main ion info (atomic mass, charge and  mass)

  Field2D J0, P0;                                // Current and pressure
  BoutReal J0_factor, P0_factor;
  Vector2D b0xcv;                                // Curvature term
  Field2D phi0;                                  // When diamagnetic terms used

  Field2D N0, Ti0, Te0, Ne0, N_imp0, T_imp0;     // number density and temperature
  Field2D Pi0, Pe0, P_imp0, density_tmp;
  Field2D q95;
  BoutReal q95_input;
  bool local_q;
  bool n0_fake_prof, n0_p0_0p3, T0_fake_prof, Nimp_lowlimit, quasi_neutral_Ni;
  Field2D LnLambda;                              // ln(Lambda)

  // V0 field vectors
  Vector2D Ve0;                                  // equilibrium ExB velocity: Ve0 = Ve0_net + Ve0_dia
  Vector2D Ve0_net;                              // net flow
  Vector2D Ve0_dia;                              // diamagnetic drift velocity.
  Vector2D V0vec;                                // V0 field vector in convection
  Vector2D V0eff;                                // effective V0 field vector in Ohm's law
  Vector3D Vexb;                                 // total ExB velocity
  Vector3D Vbtilde;                              // perturbed B vec: Vbtilde = -b0vec cross Grad Psi

  // Er0 field vectors
  Vector2D Er0;                                  // total equilibrium Er0
  BoutReal Er0_factor;                           // Er0 *= Er0_factor
  Vector2D Er0_dia;
  Vector2D Er0_net;                              // Er0 = Er0_net + Er0_dia
  
  // Equilibrium/background profiles for Vipar, Vepar, and U.
  // NOTE(malamast): Only the perturbed part of Vipar, Vepar, and U was considered in the 6-field model. 
  //                 See (Xia T.Y. et al, 2013)
  Field2D Vipar0, Vepar0, U0; 

  // 3D evolving variables
  Field3D Ni, Te, Ti, P, Pi, Pe;
  Field3D U, Vipar, Vepar, Apar, Psi;
  // derived variables
  Field3D Jpar, phi;                             // Parallel current, electric potential
  Field2D phiDC,VortDC;

  Field3D ubyn;

  Field3D Jpar2;                                 // Delp2 of Parallel current
  Field3D tmpU2, tmpA2;                          // Grad2_par2new of Parallel vector potential
  Field3D tmpN2, tmpTi2, tmpTe2, tmpVp2;         // Grad2_par2new of Parallel density

  Field3D nu_e, nu_i;                            // Electron/ion collision frequency profile (1 / S)
  Field3D vth_i, vth_e;                          // Electron/ion Thermal Velocity profile (M / S)
  Field3D kappa_par_i;                           // Ion Thermal Conductivity profile (kg*M / S^2)
  Field3D kappa_par_e;                           // Electron Thermal Conductivity profile (kg*M / S^2)
  Field2D omega_ci, omega_ce;                    // cyclotron frequency
  Field3D kappa_perp_i;                          // Ion perpendicular Thermal Conductivity profile (kg*M / S^2)
  Field3D kappa_perp_e;                          // Electron perpendicular Thermal Conductivity profile (kg*M / S^2)
  Field3D kappa_par_i_fl, kappa_par_e_fl;
  Field3D kappa_perp_i_fl, kappa_perp_e_fl;
  Field3D q_par_i, q_par_e;
  Field3D q_par_fl, q_par_landau;

  /********** Model options and additional variables ***************************/
  bool evolve_psi;
  bool emass;
  Field3D Ajpar;                                 // Parallel current, electric potential
  BoutReal emass_inv;                            // inverse of electron mass
  BoutReal coef_jpar;
  BoutReal delta_e;                              // Normalized electron skin depth
  BoutReal delta_e_inv;                          // inverse normalized electron skin depth
  BoutReal gyroAlv;                              // Normalized ion current coef

  bool diamag;
  BoutReal dia_fact;                             // Multiply diamagnetic term by this

  bool energy_flux, energy_exch;                 // energy flux term
  bool diamag_phi0;                              // Include the diamagnetic equilibrium phi0
  bool thermal_force;                            // Include the thermal flux term in Ohm's law
  bool eHall;
  bool diff_par_flutter;

  bool experiment_Er, KH_term;                   // read in phi_0 from experiment
  Field2D phi0_net, U0_net;                      // calculate Kelvin-Helmholtz term
  Field2D V0, Dphi0;                             // net flow amplitude, differential potential to flux
  bool diamag_er;                                // switch phi0 to Er

  bool nonlinear;
  bool evolve_jpar;
  BoutReal g;                                    // Only if compressible
  bool phi_curv;

  bool include_curvature, include_jpar0, compress0;
  bool include_vipar, include_vpar0, include_U0;
  bool evolve_pressure, continuity;
  bool parallel_viscous;
  
  Field3D eta_i0, pi_ci;
  
  bool gyroviscous;
  /// Temporary variables for gyroviscous
  Field3D Dperp2Phi0, Dperp2Phi, GradPhi02, GradPhi2; 
  Field3D GradparPhi02, GradparPhi2, GradcPhi, GradcparPhi;
  Field3D Dperp2Pi0, Dperp2Pi, bracketPhi0P, bracketPhiP0, bracketPhiP;
  
  /// Test the effect of first order term of invert Laplace function
  BoutReal laplace_alpha;

  /// Bootsctrap current
  bool BScurrent;
  bool radial_diffusion;
  Field3D Jpar_BS0, nu_estar, nu_istar;
  BoutReal Aratio;
  Field3D diff_radial, ddx_ni, ddx_n0;
  BoutReal diffusion_coef_Hmode0, diffusion_coef_Hmode1;

  /// neoclassical effects
  bool neoclassic_i, neoclassic_e;
  Field3D xii_neo, xie_neo, Dri_neo, rho_i, rho_e, tmpddx2;
  Field3D heatf_neo_i, heatf_neo_e, partf_neo_i;
  BoutReal epsilon;

  /// resistivity
  bool spitzer_resist;                           // Use Spitzer formula for resistivity
  BoutReal vac_lund, core_lund;                  // Lundquist number S = (Tau_R / Tau_A). -ve -> infty
  BoutReal vac_resist, core_resist;              // The resistivities (just 1 / S)
  Field3D eta;                                   // Resistivity profile (1 / S)
  BoutReal FZ;                                   // correction coefficient in Spitzer model
  
  BoutReal vacuum_pressure;
  BoutReal vacuum_trans;                         // Transition width
  Field3D vac_mask;

  /// parallel heat flux
  // case1: flux limited expression
  bool fluxlimit;                                // flux limited condition
  BoutReal q_alpha;                              // flux limiting coefficient

  // // case2: Landau damping closure
  // bool Landau;           // Use Gyro Landau Fluid closure instead of thermal conductivity
  // bool Landau_coll;      // collisional Landau Damping
  // BoutReal Landau_coeff; // Coefficient for Landau Damping
  // int nLorentzian;       // number of Lorentzians, collisional: [3, 7, 12], collisionless: >=7
  // BoutReal kappa_0;      // collisionless Landau damping coefficient
  // Field3D kappa_i;       // ion collisional Landau damping coefficient (0.5*nu_i/vth_i)
  // Field3D kappa_e;       // electron collisional Landau damping coefficient (0.5*nu_e/vthe_e)  

  Field3D SBC_value_i, SBC_value_e;
  bool full_sbc;                                 // full-f version of sheath BC for ion parallel velocity

  /// impurity
  bool impurity_prof, load_impurity, impurity_gyro;
  BoutReal Z_imp, A_imp;                         // impurity ion info
  Field3D Dperp2Pimp0, bracketPhiPimp0;
  Field2D Upara_imp;

  /// source
  bool source;
  Field3D NiSource, TeSource, TiSource;          // Axisymmetric 2D/3D sources
  BoutReal NiAmp, TeAmp, TiAmp;                  // Amplitude of the Gaussian shape sources
  int NiLoc, TeLoc, TiLoc, NiSig, TeSig, TiSig;  // Center locations and standard deviation of the sources

  /// neutral
  bool neutral, Solving_Eq_Nn, Solving_Eq_Vn;
  bool initial_Nn;
  bool full_sbc_Vn;                              // full-f version of sheath BC for neutral parallel velocity
  bool read_collision_rate;
  bool constent_Dn, Gradperp_Dn, Gradpar_Dn, Gradpar_etan, Gradperp_etan, external_source;
  bool with_vipar, with_viperp;                  // Vi~Vn
  Field3D Nn, Vn, Pn;
  Field3D Dn, Dn_fl, Dn1, etan, etan_perp, Sn, Sn_ext, Sv, S_tmp;
  bool fl_Dn;
  BoutReal Lnn_min;
  Field3D nu_iz, nu_cx, nu_rc, sigma_cx;
  // used for initialize neutral profile
  BoutReal NnAmp, NnLoc, NnSig, NnLoc_y, NnSig_y, fac_A, Diff_n, fac_Dn, fac_etan;  
  BoutReal SnAmp, SnLoc, SnSig, SnLoc_y, SnSig_y;// used for gas puffing profile
  Field3D Gamma_nn;
  BoutReal Rcyc_Nn, Rcyc_Vn;                     // recycle coefficients
  bool Nn_recyc_BC, Vn_recyc_BC;

  // Fixed-fraction model for impurity radiation
  bool fix_fraction_imp, Limp_carbon, Limp_carbon_adas, Limp_nitro, Limp_nitro_adas, Limp_Neon, Limp_Neon_adas, Limp_Argon, Limp_Argon_adas;
  Field3D Limp, Srad, Wrad;
  Field2D N_tmp0, Ne_tmp0, Te_tmp0, Limp0, Srad0, Wrad0; // Used for the linearization of the fixed-fraction radiation model

  BoutReal Wiz, Wrc;                              // ionization and recombination energy in eV
  BoutReal frac_imp;
  
  /// 3D const fueling flux
  bool with_fueling;
  bool initial_Nm, gas_puffing;
  BoutReal CF_BC_x0, CF_BC_y0, CF_BC_y1, CF_BC_z0, CF_BC_z1;
  Field3D Nm, Vmx, pm, Nm_tmp, S_diss, nu_diss, Sgas;
  Vector3D Vm;
  BoutReal Nm0, Vm0, Tm_x, Mm;

  /// parallel and perpendicular hyperdiffusion
  // M: 4th Parallel density diffusion
  BoutReal hyperdiff_par_n4, hyperdiff_par_ti4, hyperdiff_par_te4;
  // aprallel hyper-viscous diffusion for vorticity
  BoutReal hyperdiff_par_v4, hyperdiff_par_apar4, hyperdiff_par_u4; 
  // M: 4th Perpendicular density diffusion
  BoutReal hyperdiff_perp_n4, hyperdiff_perp_ti4, hyperdiff_perp_te4;
  BoutReal hyperdiff_perp_v4, hyperdiff_perp_apar4, hyperdiff_perp_u4;
  BoutReal viscos_par;                           // Parallel viscosity
  BoutReal viscos_perp;                          // Perpendicular viscosity
  BoutReal hyperviscos;                          // Hyper-viscosity (radial)
  Field3D hyper_mu_x;                            // Hyper-viscosity coefficient
  
  /// position filter
  bool pos_filter, pos_filter2, keep_zonalPF;
  BoutReal filter_position_ni, filter_position_ti, filter_position_te;
  int position_tmpi, position_tmpe, position_tmp;
  bool pos_filter_zf;
  BoutReal pos_sink_zf, Grid_NX, Grid_NY;
  BoutReal pos_filter_width, pos_filter_length, sink_pos_zf;
  Field2D pos_sink_ti, pos_sink_te, pos_sink_ni;

  /// filter low/high-n mode
  bool filter_z, filter_z_nonlinear;
  int filter_z_mode;
  int low_pass_z;
  bool zonal_flow;  //NOTE(malamast): This used to be int instead of bool. CHECK
  bool zonal_field; //NOTE(malamast): This used to be int instead of bool. CHECK
  bool zonal_bkgd;  //NOTE(malamast): This used to be int instead of bool. CHECK

  // Load background/equilibrium profiles from 2D transport similation data
  bool load_2d_bkgd;

  /********** Normalization coefficients ***************************************/
  BoutReal density;                              // density normalization factor [m^-3]
  BoutReal density_unit;                         // density unit for grid [m^-3]
  BoutReal Bbar, Lbar, Tbar, Va;                 // Normalization constants
  BoutReal Nbar, Tibar, Tebar, Tau_ie;
  // coefficients in the equations
  BoutReal Psipara1, Upara0, Upara1;
  BoutReal Upara2, Upara3, Nipara1;
  BoutReal Tipara1, Tipara2, Tipara3;
  BoutReal Tepara1, Tepara2, Tepara3, Tepara4;
  BoutReal Vepara, Vipara;
  
  BoutReal Vt0;                                  // equilibrium toroidal flow normalized to Alfven velocity
  BoutReal Vp0;                                  // equilibrium poloidal flow normalized to Alfven velocity

  /*********** Numerical/solver realted variables ************************/
  int jx, jy, jz, ncz;                           // index varriable
  int xind, indx, indy;
  BoutReal dindx;
  Field3D fp, fm;                                // Interpolated on + and - y locations
  Field2D F2D_tmp;

  /// Poisson brackets: b0 x Grad(f) dot Grad(g) / B = [f, g]
  // Method to use: BRACKET_ARAKAWA, BRACKET_STD or BRACKET_SIMPLE
  BRACKET_METHOD bm_exb, bm_mag;                 // Bracket method for advection terms
  int bracket_method_exb, bracket_method_mag;
  
  /// For preconditioner
  Field3D ni_tmp, ti_tmp, te_tmp, vi_tmp, psi_tmp, u_tmp, p_tmp, jpar1, phi_tmp;
  Field3D u_tmp1, ti_tmp2, te_tmp2;
  Field2D kappa_par_i_lin, kappa_par_e_lin;      // for preconditioner
  BoutReal Te_tmp1, Ti_tmp1, N_tmp1;

  bool limit_jacobi;                             // limit the infinity value of jacobi at x-point
  BoutReal bpxy_constraint, const_bp;
  BoutReal hthe_constraint;
  bool PF_limit;                                 // filter the instability in PF region
  BoutReal PF_limit_range;
  BoutReal PF_sink, PFs_width, PFs_length;       // sink at inner boundary of PF

  bool relax_j_vac;
  BoutReal relax_j_tconst;                       // Time-constant for j relax
  Field3D Psitarget;                             // The (moving) target to relax to

  bool smooth_j_x;                               // Smooth Jpar in the x direction
  bool mask_j_x, mask_phi_x;                     // Mask Jpar at the inner boundary of x
  Field3D mask_jx1d, mask_px1d;                  // the variable of mask function, normalized to 1.
  int mask_flag_j, mask_flag_phi;
  BoutReal mask_width, mask_length;
  BoutReal filter_nl;
  
  int jpar_bndry_width;                          // Zero jpar in a boundary region

  bool parallel_lr_diff;                         // Use left and right shifted stencils for parallel differences
  bool parallel_lagrange;                        // Use (semi-) Lagrangian method for parallel derivatives
  bool parallel_project;                         // Use Apar to project field-lines

  /// for debug purpose
  Field3D term1, term2, term3, term4, term5;

  /****************************************************************************/
  /****************************************************************************/

  /// The total height, average width and center of profile of N0
  BoutReal n0_height, n0_ave, n0_width, n0_center, n0_bottom_x;
  /// The ampitude of congstant temperature
  BoutReal Tconst;

  /// Constraint
  Field3D C_phi;

  /// Parameters
  BoutReal diffusion_par;                        // Parallel thermal conductivity
  BoutReal diffusion_perp;                       // Perpendicular thermal conductivity (>0 open)

  BoutReal heating_P;                            // heating power in pressure
  BoutReal hp_width;                             // heating profile radial width in pressure
  BoutReal hp_length;                            // heating radial domain in pressure
  BoutReal sink_vp;                              // sink in pressure
  BoutReal sp_width;                             // sink profile radial width in pressure
  BoutReal sp_length;                            // sink radial domain in pressure

  BoutReal sink_Ul;                              // left edge sink in vorticity
  BoutReal su_widthl;                            // left edge sink profile radial width in vorticity
  BoutReal su_lengthl;                           // left edge sink radial domain in vorticity

  BoutReal sink_Ur;                              // right edge sink in vorticity
  BoutReal su_widthr;                            // right edge sink profile radial width in vorticity
  BoutReal su_lengthr;                           // right edge sink radial domain in vorticity

  BoutReal sink_Ter;                             // right edge sink in Psi
  BoutReal ste_widthr;                           // right edge sink profile radial width in Psi
  BoutReal ste_lengthr;                          // right edge sink radial domain in Psi

  BoutReal Low_limit;                            // To limit the negative value of total density and temperatures

  Field3D Te_tmp, Ti_tmp, N_tmp, Ne_tmp;         // to avoid the negative value of total value
  BoutReal gamma_i_BC, gamma_e_BC;               // sheath energy transmission factors
  int Sheath_width;
  bool SBC_phi;
  // variables for sheath boundary conditions
  Field3D c_se, Jpar_sh, q_se, q_si, vth_et, c_set, phi_sh; 
  Field2D vth_e0, c_se0, Jpar_sh0, phi_sh0;
  BoutReal const_cse;
  
  /*****************************************************************************/

  Field3D Xip_x, Xip_z;                          // Displacement of y+1 (in cell index space)
  Field3D Xim_x, Xim_z;                          // Displacement of y-1 (in cell index space)

  bool phi_constraint;                           // Solver for phi using a solver constraint

  bool output_transfer;                          // output the results of energy transfer
  bool output_ohm;                               // output the results of the terms in Ohm's law
  bool output_flux_par;                          // output the results of parallel particle and heat flux
  bool output_vradial;                           // output the results of radial velocity, induced by ExB and magnetic flutter
  bool output_Teterms, output_Titerms, output_Tevegradte, output_qparcompare;

  /// Maxwell stress, Reynolds stress, ion diamagbetic and curvature term
  Field3D T_M, T_R, T_ID, T_C, T_G;
  Field3D ohm_phi, ohm_hall, ohm_thermal;
  // particle flux, ion and elelctron heat flux
  Field3D gamma_par_i, heatf_par_i, heatf_par_e;
  Field3D heatf_par_flutter_i, heatf_par_flutter_e;
  // temp variable for perturbed parallel thermal conduction
  Field3D bracket1i, gradpar_ti, bracket1e, gradpar_te; 
  Field3D Vbti_par, Vbte_par;

  BoutReal hyperresist;                          // Hyper-resistivity coefficient (in core only)
  BoutReal ehyperviscos;                         // electron Hyper-viscosity coefficient

  int damp_width;                                // Width of inner damped region
  BoutReal damp_t_const;                         // Timescale of damping

  /// Communication objects
  FieldGroup comms;

#if BOUT_HAS_HYPRE
  std::unique_ptr<LaplaceXY2Hypre> laplacexy{nullptr}; // Laplacian solver in X-Y (n=0)
#else
  std::unique_ptr<LaplaceXY> laplacexy{nullptr};       // Laplacian solver in X-Y (n=0)
#endif

  /// Solver for inverting Laplacian
  std::unique_ptr<Laplacian> phiSolver{nullptr};
  std::unique_ptr<Laplacian> aparSolver{nullptr};

  const Field3D ret_const_flux_BC(const Field3D &var, const BoutReal value) {
    Field3D result;
    result.allocate();
    for (int jx = 0; jx < mesh->LocalNx; jx++) {
      int x_glb = mesh->getGlobalXIndex(jx);
      for (int jy = 0; jy < mesh->LocalNy; jy++) {
	      int y_glb = mesh->getGlobalYIndex(jy);
        for (int jz = 0; jz < mesh->LocalNz; jz++) {
          if(x_glb >= CF_BC_x0 && y_glb >= CF_BC_y0 && y_glb <= CF_BC_y1) {
	          result(jx,jy,jz) = value;
	        } else {
	          result(jx,jy,jz) = var(jx,jy,jz);
	        }
	      }	  
      }
    }
    mesh->communicate(result);
    return(result);
  }

  const Field2D tanhxl_core(const int filter_index) {
    Field2D result;
    result.allocate();
    result = Field2D(0.);
    BoutReal xpos, width, length, tanh_tmp;
    int indy;
    xpos = filter_index;
    width = xpos * pos_filter_width;   // the width of the tanh filter function
    length = xpos * pos_filter_length; // the middle point of the tanh filter function
    
    for (jx = 0; jx < mesh->LocalNx; jx++) {
      indx = mesh->getGlobalXIndex(jx);
      tanh_tmp = (1. - tanh((indx - length) / width)) / 2.;
      if (tanh_tmp < 0.005)
        tanh_tmp = 0.;
      for (jy = 0; jy < mesh->LocalNy; jy++) {
        indy = mesh->getGlobalYIndex(jy);
        if (((indy > int(jysep1) - 2) && (indy <= int(jysep2_1))) || ((indy > int(jysep1_2)) && (indy <= int(jysep2) + 2))) {
          result(jx,jy) = tanh_tmp;
        } else
          result(jx,jy) = 0.;
      }
    }
  #ifdef CHECK
    msg_stack.pop();
  #endif
    return result;
  }

  const Field3D lowPass_pos2(const Field3D &var, const Field3D &prof) {
    Field3D tmp_result, result;
    result.allocate();
    tmp_result.allocate();
    Field2D prof2d = DC(prof);
    BoutReal y_ind;
    bool zonal = false;
    
    result = var;
    tmp_result = lowPass(var, low_pass_z, zonal);

    for (int jx = 0; jx < mesh->LocalNx; jx++) {
      for (int jy = 0; jy < mesh->LocalNy; jy++) {
        y_ind = mesh->getGlobalYIndex(jy);
        if (keep_zonalPF && (((y_ind > int(jysep1) - 2) && (y_ind <= int(jysep2_1))) || ((y_ind > int(jysep1_2)) && (y_ind <= int(jysep2) + 2)))) {
          if (prof2d(jx,jy) < 0.) {
            for (int jz = 0; jz < mesh->LocalNz; jz++) {
              result(jx,jy,jz) = tmp_result(jx,jy,jz);
	    }
	  }
        }
      }
    }
    return result;
  }

  const Field3D lowPass_pos(const Field3D &var, int filter_index) {
    Field3D tmp_result, result;
    result.allocate();
    tmp_result.allocate();
    int y_ind;
    bool zonal = false;

    result = var;
    tmp_result = lowPass(var, low_pass_z, zonal);

    for (int jx = 0; jx < mesh->LocalNx; jx++) {
      if (mesh->getGlobalXIndex(jx) <= filter_index) {
        for (int jy = 0; jy < mesh->LocalNy; jy++) {
          y_ind = mesh->getGlobalYIndex(jy);
          if (keep_zonalPF && (((y_ind > int(jysep1) - 2) && (y_ind <= int(jysep2_1))) || ((y_ind > int(jysep1_2)) && (y_ind <= int(jysep2) + 2)))) {
            for (int jz = 0; jz < mesh->LocalNz; jz++) {
              result(jx,jy,jz) = tmp_result(jx,jy,jz);
            }
          }
        }
      }
    }
    return result;
  }

  const Field3D filter_z_non(const Field3D &var, int N0, int N1) {
    // ASSERT1(var.isAllocated());
    static dcomplex *f = (dcomplex *)NULL;
    ncz = mesh->LocalNz - 1;
    if (f == (dcomplex *)NULL) {
      // Allocate memory
      f = new dcomplex[ncz / 2 + 1];
    }

    Field3D result;
    result.allocate();
    for (jx = 0; jx < mesh->LocalNx; jx++) {
      for (jy = 0; jy < mesh->LocalNy; jy++) {

        rfft(var(jx,jy), ncz, f); // Forward FFT

        for (jz = 0; jz <= ncz / 2; jz++) {

          if ((jz != N0) && (jz != N1)) {
            // Zero this component
            f[jz] = 0.0;
          }
        }

        irfft(f, ncz, result(jx,jy)); // Reverse FFT

        result(jx,jy,ncz) = result(jx,jy,0);
      }
    }
  #ifdef TRACK
    result.name = "filter(" + var.name + ")";
  #endif
    // result.location = var.location;
    return result;
  }

  BoutReal TanH(BoutReal a)
  {
    BoutReal temp = exp(a);
    return ((temp - 1.0 / temp) / (temp + 1.0 / temp));
  }

  const Field3D mask_x_1d(bool BoutRealspace, int mask_flag, BoutReal mask_width, BoutReal mask_length) {
    Field3D result;
    result.allocate();

    BoutReal Grid_NX; // the grid number on x, and the
    mesh->get(Grid_NX, "nx");

    // create a radial buffer zone to set jpar zero near radial boundary
    BoutReal min_tmp = (TanH((4. / Grid_NX - mask_length) / mask_width) + 1.) / 2.;
    BoutReal max_tmp = (TanH(((1. - (Grid_NX - 5.) / Grid_NX) - mask_length) / mask_width) + 1.) / 2.;
    for (jx = 0; jx < mesh->LocalNx; jx++)
      for (jy = 0; jy < mesh->LocalNy; jy++)
        for (jz = 0; jz < mesh->LocalNz; jz++) {
          BoutReal lx = mesh->GlobalX(jx);
          BoutReal dampl = (TanH((lx - mask_length) / mask_width) + 1.) / 2. - min_tmp;
          BoutReal dampr = (TanH(((1. - lx) - mask_length) / mask_width) + 1.) / 2. - max_tmp;
          if (mask_flag == 0) // left mask
            result(jx,jy,jz) = dampl;
          else if (mask_flag == 1) // right mask
            result(jx,jy,jz) = dampr;
          else // mask on both boundary
            result(jx,jy,jz) = dampl * dampr;
          if (result(jx,jy,jz) < 0)
            result(jx,jy,jz) = 0.;
        }
    result /= max(result, true);
    if (BoutRealspace)
      //result = result.shiftZ(false); // Shift back
    
    // Need to communicate boundaries
    mesh->communicate(result);
 
    return result;
  }
  
  const Field2D field_larger(const Field2D &f, const BoutReal limit) {
    Field2D result;
    result.allocate();

    // pragma omp parallel for
    for (jx = 0; jx < mesh->LocalNx; jx++)
      for (jy = 0; jy < mesh->LocalNy; jy++) {
        if (f(jx,jy) >= limit)
          result(jx,jy) = f(jx,jy);
        else {
          result(jx,jy) = 0.9 * limit + 0.1 * f(jx,jy);
        }
      }
    for (jx = 1; jx < mesh->LocalNx - 1; jx++)
      for (jy = 1; jy < mesh->LocalNy - 1; jy++) {
        {
          if (f(jx,jy) <= 1.2 * limit) {
            result(jx,jy) = 0.5 * result(jx,jy) + 0.25 * (result(jx - 1,jy) + result(jx + 1,jy));
            result(jx,jy) = 0.5 * result(jx,jy) + 0.25 * (result(jx - 1,jy) + result(jx + 1,jy));
          }
        // Constraint
        Field3D C_phi;
     }
      }
    mesh->communicate(result);
    return result;
  }

  const Field3D field_larger(const Field3D& f, const BoutReal limit) {
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

  const Field3D PF_filter(const Field3D &input, const BoutReal PF_limit_range) {
    Field3D result;
    result.allocate();
    for (int jx = 0; jx < mesh->LocalNx; jx++) {
      indx = mesh->getGlobalXIndex(jx);
      dindx = indx / ixsep;
      for (int jy = 0; jy < mesh->LocalNy; jy++) {
        indy = mesh->getGlobalYIndex(jy);
        if ((dindx < PF_limit_range) && ((indy <= int(jysep1) - 2) || ((indy > int(jysep2_1)) && (indy <= int(jysep1_2))) || (indy > int(jysep2) + 2))) {
          for (int jz = 0; jz < mesh->LocalNz; jz++)
            result(jx,jy,jz) = 0.;
        } else {
          for (int jz = 0; jz < mesh->LocalNz; jz++)
            result(jx,jy,jz) = input(jx,jy,jz);
        }
      }
    }

    mesh->communicate(result);
    return result;
  }

  const Field3D PF_filter_2(const Field3D &input, const BoutReal PF_limit_range) {
    Field3D result;
    result.allocate();
    Field2D prof2d = DC(input);

    for (int jx = 0; jx < mesh->LocalNx; jx++) {
      indx = mesh->getGlobalXIndex(jx);
      dindx = indx / ixsep;
      for (int jy = 0; jy < mesh->LocalNy; jy++) {
        indy = mesh->getGlobalYIndex(jy);
        if ((dindx < PF_limit_range) && ((indy <= int(jysep1) - 2) || ((indy > int(jysep2_1)) && (indy <= int(jysep1_2))) || (indy > int(jysep2) + 2))) {
          BoutReal value_dc = prof2d(jx,jy);
          for (int jz = 0; jz < mesh->LocalNz; jz++)
            result(jx,jy,jz) = value_dc;
        } else {
          for (int jz = 0; jz < mesh->LocalNz; jz++)
            result(jx,jy,jz) = input(jx,jy,jz);
        }
      }
    }

    mesh->communicate(result);
    return result;
  }

  const Field3D sink_zonal_core(const Field3D &var, int filter_index) {
    Field3D result;
    result.allocate();
    static dcomplex *f = NULL, *f2 = NULL;
    int indx;

  #ifdef CHECK
    msg_stack.push("sink_zonal_core(Field3D, int)", filter_index);
  #endif

    BoutReal xpos, width, length, tanh_tmp;
    xpos = filter_index;
    width = xpos * pos_filter_width;   // the width of the tanh filter function
    length = xpos * pos_filter_length; // the middle point of the tanh filter function

    // output.write("Error zonal_part~\n");
    if (!var.isAllocated()) {
      return DC(var);
    }
    int ncz = mesh->LocalNz - 1;

    if (f == NULL)
      f = new dcomplex[ncz / 2 + 1];

    if (f2 == NULL)
      f2 = new dcomplex[ncz / 2 + 1];

    for (jx = 0; jx < mesh->LocalNx; jx++) {
      indx = mesh->getGlobalXIndex(jx);
      tanh_tmp = (1. + TanH((indx - length) / width)) / 2.;
      if (tanh_tmp < 0.005)
        tanh_tmp = 0.;
      for (jy = 0; jy < mesh->LocalNy; jy++) {
        indy = mesh->getGlobalYIndex(jy);
        //if ( ((indy > int(jysep1)) && (indy <= int(jysep2_1))) || ((indy > int(jysep1_2)) && (indy <= int(jysep2))) )
        {
          // Take FFT in the Z direction
          rfft(var(jx,jy), ncz, f);
          // Filter the zonal component based on the filter_index
          // f[0] *= (1.+TanH( (indx-length)/width ))/2.;
          f[0] *= tanh_tmp;
        }

        irfft(f, ncz, result(jx,jy)); // Reverse FFT
        result(jx,jy,ncz) = result(jx,jy,0);
      }
    }

  #ifdef CHECK
    msg_stack.pop();
  #endif
    // mesh->communicate(result);
    return result;
  }

  const Field3D sink_PF(const Field2D &f0, const Field3D &f, const BoutReal width, const BoutReal length) {
    Field3D result;
    result.allocate();

    result = sink_tanhxl(f0, f, width, length);
    for (jy = 0; jy < mesh->LocalNy; jy++) {
      indy = mesh->getGlobalYIndex(jy);
      if (((indy > int(jysep1) - 2) && (indy <= int(jysep2_1))) || ((indy > int(jysep1_2)) && (indy <= int(jysep2) + 2))) {
        for (jx = 0; jx < mesh->LocalNx; jx++)
          for (jz = 0; jz < mesh->LocalNz; jz++)
            result(jx,jy,jz) = 0.;
      } else {
        for (jx = 0; jx < mesh->LocalNx; jx++) {
          if (mesh->getGlobalXIndex(jx) >= ixsep)
            for (jz = 0; jz < mesh->LocalNz; jz++)
              result(jx,jy,jz) = 0.;
        }
      }
    }
    mesh->communicate(result);
    return result;
  }

  const Field3D Grad2_par2new(const Field3D& f) {
    // This function implements d2/dy2 where y is the poloidal coordinate theta

    TRACE("Grad2_par2new( Field3D )");

    Field3D result = D2DY2(f);

  #if BOUT_USE_TRACK
    result.name = "Grad2_par2new(" + f.name + ")";
  #endif

    return result;
  }

  const Field2D N0tanh(BoutReal n0_height, BoutReal n0_ave, BoutReal n0_width, BoutReal n0_center, BoutReal n0_bottom_x) {
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
        for (int jy = 0; jy < mesh->LocalNy; jy++) {
          int globaly = mesh->getGlobalYIndex(jy);
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
  const Field3D Grad_parP(const Field3D& f, CELL_LOC loc = CELL_DEFAULT) {
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
        if (evolve_psi)
	  result -= bracket(Psi, f, bm_mag) * B0;
	else
	  result -= bracket(Apar, f, bm_mag);
      }
    }

    return result;
  }

  const Field3D BS_ft(const int index) {
    Field3D result, result1;
    result.allocate();
    result1.allocate();
    result1 = 0.;

    BoutReal xlam, dxlam;
    dxlam = 1. / max(B0) / index;
    xlam = 0.;

    for (int i = 0; i < index; i++) {
      result1 += xlam * dxlam / sqrt(1. - xlam * B0);
      xlam += dxlam;
    }
    result = 1. - 0.75 * B0 * B0 * result1;

    return result;
  }

  const Field3D F31(const Field3D input) {
    Field3D result;
    result.allocate();

    result = (1 + 1.4 / (Zi + 1.)) * input;
    result -= 1.9 / (Zi + 1.) * input * input;
    result += 0.3 / (Zi + 1.) * input * input * input;
    result += 0.2 / (Zi + 1.) * input * input * input * input;

    return result;
  }

  const Field3D F32ee(const Field3D input) {
    Field3D result;
    result.allocate();

    result = (0.05 + 0.62 * Zi) / (Zi * (1 + 0.44 * Zi)) * (input - input * input * input * input);
    result += 1. / (1. + 0.22 * Zi) * (input * input - input * input * input * input - 1.2 * (input * input * input - input * input * input * input));
    result += 1.2 / (1. + 0.5 * Zi) * input * input * input * input;

    return result;
  }

  const Field3D F32ei(const Field3D input) {
    Field3D result;
    result.allocate();

    result = -(0.56 + 1.93 * Zi) / (Zi * (1 + 0.44 * Zi)) * (input - input * input * input * input);
    result += 4.95 / (1. + 2.48 * Zi) * (input * input - input * input * input * input - 0.55 * (input * input * input - input * input * input * input));
    result -= 1.2 / (1. + 0.5 * Zi) * input * input * input * input;

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

    // Load data from the grid
    mesh->get(ixsep, "ixseps1");
    mesh->get(ixsep2, "ixseps2");
    mesh->get(jysep1, "jyseps1_1");
    mesh->get(jysep2, "jyseps2_2");
    mesh->get(jysep1_2, "jyseps1_2");
    mesh->get(jysep2_1, "jyseps2_1");

    mesh->get(Grid_NX, "nx");
    mesh->get(Grid_NY, "ny");

    if (ixsep==Grid_NX) {
      output.write("Cicular geometry without limiter!\n");
      mag_config = 1;
    } else if (ixsep2<Grid_NX) {
      output.write("Double null geometry!\n");
      mag_config = 4;
    } else if (jysep1<0) {
      output.write("Circular geometry with limiter!\n");
      mag_config = 2;
    } else if (jysep1>0) {
      output.write("Single null geometry!\n");
      mag_config = 3;
    } else {
      output.write("CAUTION: magnetic configuration cannot be determined!\n");
    }

    // Load 2D profiles
    mesh->get(J0, "Jpar0");                      // A / m^2

    if (mesh->get(P0, "pressure_s")) {           // Pascals
      mesh->get(P0, "pressure");
      output.write("Using pressure as P0.\n");
    } else {
      output.write("Using pressure_s as P0.\n");
    }

    // Load curvature term
    b0xcv.covariant = false;                     // Read contravariant components
    mesh->get(b0xcv, "bxcv");                    // mixed units x: T y: m^-2 z: m^-2

    // Load metrics
    if (mesh->get(Rxy, "Rxy")) { // m
      output_error.write("Error: Cannot read Rxy from grid\n");
      return 1;
    }
    if (mesh->get(Bpxy, "Bpxy")) { // T
      output_error.write("Error: Cannot read Bpxy from grid\n");
      return 1;
    }
    mesh->get(Btxy, "Btxy");                     // T
    mesh->get(B0, "Bxy");                        // T
    mesh->get(hthe, "hthe");                     // m
    mesh->get(I, "sinty");                       // m^-2 T^-1

    /////////////////////////////////////////////////////////////////
    // Read parameters from the options file
    //
    // Options.get ( NAME,    VARIABLE,    DEFAULT VALUE)
    //
    // or if NAME = "VARIABLE" then just
    //
    // OPTION(VARIABLE, DEFAULT VALUE)
    //
    // Prints out what values are assigned
    /////////////////////////////////////////////////////////////////

    auto& globalOptions = Options::root();
    auto& options = globalOptions["highbeta"];

    n0_fake_prof = options["n0_fake_prof"]
                   .doc("use the hyperbolic profile of n0. If both  n0_fake_prof and "
                        "T0_fake_prof are "
                        "false, use the profiles from grid file")
                   .withDefault(false);
    n0_p0_0p3 = options["n0_p0_0p3"].doc("use n0 ~ P0^0.3").withDefault(false);
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
    n0_bottom_x = options["n0_bottom_x"]
                   .doc("the start of flat region of N0 on SOL side, in percentage of x")
                   .withDefault(0.81);
    T0_fake_prof = options["T0_fake_prof"].withDefault(false);
    Tconst = options["Tconst"]
                   .doc("the amplitude of constant temperature, in percentage")
                   .withDefault(-1.0);

                   
    load_2d_bkgd = options["load_2d_bkgd"]
	           .doc("Load background/equilibrium profiles from 2D transport similation data.")
		   .withDefault(false);

    impurity_prof = options["impurity_prof"]
	           .doc("Include the profile of impurity in the equilibrium")
		   .withDefault(false);
    impurity_gyro = options["impurity_gyro"]
	           .doc("Include the gyro viscous of impurity")
		   .withDefault(false);
    load_impurity = options["load_impurity"]
	           .doc("if load impurity from grid")
		   .withDefault(false);
    Z_imp = options["Z_imp"].doc("The charge number of impurity").withDefault(6.0);
    A_imp = options["A_imp"].doc("The mass number of impurity").withDefault(12.0);
    Nimp_lowlimit = options["Nimp_lowlimit"]
	           .doc("The switch to limit the lowest value of impurity density")
		   .withDefault(false);
    quasi_neutral_Ni = options["quasi_neutral_Ni"]
	           .doc("The switch to use quasi neutral condition to calculate N0 instead of reading from grid")
		   .withDefault(false);

    experiment_Er = options["experiment_Er"].withDefault(false);
    Er0_factor = options["Er0_factor"].withDefault(1.0);     // change Er0 *= Er0_factor
    KH_term = options["KH_term"].withDefault(false);         // switch to Kelvin-Helmholtz term
    J0_factor = options["J0_factor"].withDefault(1.0);
    P0_factor = options["P0_factor"].withDefault(1.0);

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
    SBC_phi = options["SBC_phi"]
                  .doc("use sheath boundary on phi instead of Jpar")
                  .withDefault(false);

    density = options["density"].doc("number density normalization factor [m^-3]").withDefault(1.0e20);
    density_unit = options["density_unit"].doc("Number density unit for grid [m^-3]").withDefault(1.0e20);
    Zi = options["Zi"].doc("ion charge number").withDefault(1);

    evolve_jpar =
        options["evolve_jpar"].doc("If true, evolve J raher than Psi").withDefault(false);
    phi_constraint =
        options["phi_constraint"].doc("Use solver constraint for phi").withDefault(false);

    // Effects to include/exclude
    nonlinear = options["nonlinear"].withDefault(false);
    include_curvature = options["include_curvature"].withDefault(true);
    include_jpar0 = options["include_jpar0"].withDefault(true);
    evolve_pressure = options["evolve_pressure"].withDefault(true);
    evolve_psi = options["evolve_psi"].withDefault(true);

    continuity = options["continuity"].doc("use continuity equation").withDefault(false);
    compress0 = options["compress0"].withDefault(false);
    gyroviscous = options["gyroviscous"].withDefault(false);
    parallel_viscous = options["parallel_viscous"].withDefault(false);

    include_vipar = options["include_vipar"].withDefault(false);
    include_vpar0 = options["include_vpar0"].withDefault(false);    // include Vpar0 terms. Vpar0  is the ecuilibrium 2D backgraound from 2D trasnport simulations
    include_U0 = options["include_U0"].withDefault(false);          // include U0 terms. U0 is the ecuilibrium 2D vorticity backgraound from 2D trasnport simulations

    BScurrent = options["BScurrent"].withDefault(false);
    Aratio = options["Aratio"].withDefault(0.35);

    radial_diffusion = options["radial_diffusion"].withDefault(false);
    diffusion_coef_Hmode0 = options["diffusion_coef_Hmode0"].withDefault(1.0);    // default value of radial diffusion coefficient
    diffusion_coef_Hmode1 = options["diffusion_coef_Hmode1"].withDefault(10.0);   // upper limit of radial diffusion coefficient

    //path = options["path"].withDefault("./");                                   // The path of the original Vexb data
    pos_filter = options["pos_filter"].withDefault(false);                        // switch to turn on the filter of the negative value of zonal background
    pos_filter2 = options["pos_filter2"].withDefault(false);                      // switch to turn on the filter inside certain position
    pos_filter_zf = options["pos_filter_zf"].withDefault(false);                  // switch to turn on the filter of the dc profiles inside certain postion with tanh function
    keep_zonalPF = options["keep_zonalPF"].withDefault(false);                    // keep the zonal component in PF region when zonal filter is turned on
    filter_position_ni = options["filter_position_ni"].withDefault(100);          // radial index of the filter. Zonal component of Ni in the x range of 0 - filter_position will be filtered.
    filter_position_ti = options["filter_position_ti"].withDefault(100);
    filter_position_te = options["filter_position_te"].withDefault(100);
    pos_filter_width = options["pos_filter_width"].withDefault(0.02);             // width of po_filter_zf with the normalization of filter_position
    pos_filter_length = options["pos_filter_length"].withDefault(0.95);           // length of po_filter_zf with the normalization of filter_position
    pos_sink_zf = options["pos_sink_zf"].withDefault(-10);                        // switch to turn on the sink of the dc profiles inside certain postion with tanh function

    if (pos_filter2 || pos_filter_zf || (pos_sink_zf > 0.)) {
      if (filter_position_ni <= 0) {
        filter_position_ni = ixsep;
        output << "\tWarning: filter_position_ni has a negtive value!\n";
        output << "\tfilter_position_ni is forced to be isxeps1=" << ixsep << ".\n";
      }
      if (filter_position_ti <= 0) {
        filter_position_ti = ixsep;
	      output << "\tWarning: filter_position_ti has a negtive value!\n";
	      output << "\tfilter_position_ti is forced to be isxeps1=" << ixsep << ".\n";	
      }
      if (filter_position_te <= 0) {
        filter_position_te = ixsep;
	      output << "\tWarning: filter_position_te has a negtive value!\n";
	      output << "\tfilter_position_te is forced to be isxeps1=" << ixsep << ".\n";
      }
      position_tmpi = std::min(filter_position_ni, filter_position_ti);
      position_tmpe = std::min(filter_position_ni, filter_position_te);
      position_tmp = std::min(position_tmpi, position_tmpe);
    }     

    //source 
    NiAmp = options["NiAmp"].withDefault(-1.0);    // Amplitude of the explicit particle sourcing
    TeAmp = options["TeAmp"].withDefault(-1.0);    // Amplitude of the explicit electron energy sourcing
    TiAmp = options["TiAmp"].withDefault(-1.0);    // Amplitude of the explicit ion energy sourcing
    NiLoc = options["NiLoc"].withDefault(floor(ixsep/4.));
    TeLoc = options["TeLoc"].withDefault(floor(ixsep/4.));
    TiLoc = options["TiLoc"].withDefault(floor(ixsep/4.));
    NiSig = options["NiSig"].withDefault(floor(ixsep/12.));
    TeSig = options["TeSig"].withDefault(floor(ixsep/12.));
    TiSig = options["TiSig"].withDefault(floor(ixsep/12.));

    neutral = options["neutral"].withDefault(false);
    Solving_Eq_Nn = options["Solving_Eq_Nn"].withDefault(false);
    Solving_Eq_Vn = options["Solving_Eq_Vn"].withDefault(false);
    with_vipar = options["with_vipar"].withDefault(false);
    initial_Nn = options["initial_Nn"].withDefault(false);
    fl_Dn = options["fl_Dn"].withDefault(false);
    constent_Dn = options["constent_Dn"].withDefault(false);
    Gradperp_Dn = options["Gradperp_Dn"].withDefault(false);
    Gradpar_Dn = options["Gradpar_Dn"].withDefault(false);
    Gradpar_etan = options["Gradpar_etan"].withDefault(false);
    Gradperp_etan = options["Gradperp_etan"].withDefault(false);
    external_source = options["external_source"].withDefault(false);
    read_collision_rate = options["read_collision_rate"].withDefault(false);
    NnAmp = options["NnAmp"].withDefault(0.2);
    NnLoc = options["NnLoc"].withDefault(ixsep);
    NnSig = options["NnSig"].withDefault(floor(Grid_NX/4.));
    NnLoc_y = options["NnLoc_y"].withDefault(31.0);
    NnSig_y = options["NnSig_y"].withDefault(4.0);
    Nn_recyc_BC = options["Nn_recyc_BC"].withDefault(false);
    Vn_recyc_BC = options["Vn_recyc_BC"].withDefault(false);
    full_sbc_Vn = options["full_sbc_Vn"].withDefault(false);
    Rcyc_Nn = options["Rcyc_Nn"].withDefault(1.0);
    Rcyc_Vn = options["Rcyc_Vn"].withDefault(1.0);
    fac_A = options["fac_A"].withDefault(0.4);
    fac_Dn = options["fac_Dn"].withDefault(1.0);
    fac_etan = options["fac_etan"].withDefault(1.0);
    Lnn_min = options["Lnn_min"].withDefault(1.0e-3);
    Diff_n = options["Diff_n"].withDefault(5.0);
    SnAmp = options["SnAmp"].withDefault(0.1);
    SnLoc = options["SnLoc"].withDefault(ixsep);
    SnSig = options["SnSig"].withDefault(floor(Grid_NX/4.));
    SnLoc_y = options["SnLoc_y"].withDefault(31);
    SnSig_y = options["SnSig_y"].withDefault(4);
    Wiz = options["Wiz"].withDefault(13.6);    // ionization energy in eV   
    Wrc = options["Wrc"].withDefault(4.5);     //recombination energy in eV

    fix_fraction_imp = options["fix_fraction_imp"].withDefault(false);
    Limp_carbon = options["Limp_carbon"].withDefault(false);
    Limp_nitro = options["Limp_nitro"].withDefault(false);
    Limp_Neon = options["Limp_Neon"].withDefault(false);
    Limp_Argon = options["Limp_Argon"].withDefault(false);
    Limp_carbon_adas = options["Limp_carbon_adas"].withDefault(false);
    Limp_nitro_adas = options["Limp_nitro_adas"].withDefault(false);
    Limp_Neon_adas = options["Limp_Neon_adas"].withDefault(false);
    Limp_Argon_adas = options["Limp_Argon_adas"].withDefault(false);
    frac_imp = options["frac_imp"].withDefault(0.0);

    with_fueling = options["with_fueling"].withDefault(false);
    initial_Nm = options["initial_Nm"].withDefault(false);
    gas_puffing = options["gas_puffing"].withDefault(false);
    CF_BC_x0 = options["CF_BC_x0"].withDefault(1.01);
    CF_BC_y0 = options["CF_BC_y0"].withDefault(0.47);
    CF_BC_y1 = options["CF_BC_y1"].withDefault(0.53);
    CF_BC_z0 = options["CF_BC_z0"].withDefault(0.0);
    CF_BC_z1 = options["CF_BC_z1"].withDefault(2.0);  //default smbi in whole Z
    Vm0 = options["Vm0"].withDefault(-500.0);         // Read in m/s
    Nm0 = options["Nm0"].withDefault(1.0e7);          // Read in in 1.e^20 m^-3
    Tm_x = options["Tm_x"].withDefault(0.0258);       // eV
    Mm = options["Mm"].withDefault(2.0);              // in Mi

    // int bracket_method;
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

    // int bracket_method;
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

    eHall = options["eHall"]
	    .doc("electron Hall or electron parallel pressue gradient effects?")
	    .withDefault(false);
    thermal_force = options["thermal_force"]
	    .doc("Thermal flux in Ohm's Law")
	    .withDefault(false);
    AA = options["AA"]
	    .doc("ion mass in units of proton mass")
	    .withDefault(2.0);
    Mi = Mp * AA;

    emass = options["emass"]
            .doc("including electron inertial, electron mass")
            .withDefault(false);
    emass_inv = options["emass_inv"].doc("inverse of electron mass").withDefault(1.0);

    diamag = options["diamag"].doc("Diamagnetic effects?").withDefault(false);
    diamag_phi0 = options["diamag_phi0"].doc("Include equilibrium phi0").withDefault(diamag);
    dia_fact = options["dia_fact"].doc("Scale diamagnetic effects by this factor").withDefault(1.0);
    diamag_er = options["diamag_er"].doc("switch from phi0 to Er0").withDefault(false);
    energy_flux = options["energy_flux"].doc("energy flux").withDefault(false);
    energy_exch = options["energy_exch"].doc("energy exchange").withDefault(false);

    noshear = options["noshear"].withDefault(false);

    relax_j_vac =
        options["relax_j_vac"].doc("Relax vacuum current to zero").withDefault(false);
    relax_j_tconst = options["relax_j_tconst"].withDefault(0.1);

    // Toroidal filtering
    filter_z = options["filter_z"].doc("Filter a single n").withDefault(false);
    filter_z_mode = options["filter_z_mode"].withDefault(1);
    filter_z_nonlinear = options["filter_z_nonlinear"]
	                .doc("Filter a single n and zonal").withDefault(false);
    low_pass_z = options["low_pass_z"].doc("Low pass filter. < 0 -> off").withDefault(-1);
    zonal_flow = options["zonal_flow"].withDefault(false);           // zonal flow filter
    zonal_field = options["zonal_field"].withDefault(false);         // zonal field filter
    zonal_bkgd = options["zonal_bkgd"].withDefault(false);           // zonal background P filter
    if (zonal_flow) {
      // Create an XY solver for n=0 component
#if BOUT_HAS_HYPRE
      laplacexy = bout::utils::make_unique<LaplaceXY2Hypre>(mesh);
#else
      laplacexy = bout::utils::make_unique<LaplaceXY>(mesh);
#endif
      // Set coefficients for Boussinesq solve
      laplacexy->setCoefs(1.0, 0.0);
      phiDC = 0.0; // Starting guess
      phiDC.setBoundary("phi");
    }

    filter_nl = options["filter_nl"].doc("zonal background P filter").withDefault(-1);

    limit_jacobi = options["limit_jacobi"].withDefault(false);       // limit the value of jacobi at x-point region
    bpxy_constraint = options["bpxy_constraint"].withDefault(0.04);
    hthe_constraint = options["hthe_constraint"].withDefault(0.04);
    const_bp = options["const_bp"].withDefault(-1.);
    PF_limit = options["PF_limit"].withDefault(false);               // filter the instability in PF
    PF_limit_range = options["PF_limit_range"].withDefault(0.1);     // range of filter in PF
    PF_sink = options["PF_sink"].withDefault(-1.);                   // the coefficents of PF sink
    PFs_width = options["PFs_width"].withDefault(0.2);               // the percentage of radial grid points for sink profile radial width in PF
    PFs_length = options["PFs_length"].withDefault(0.4);             // the percentage of radial grid points for sink profile radial domain in PF

    // Radial smoothing
    smooth_j_x = options["smooth_j_x"]
	          .doc("Smooth Jpar in x")
		  .withDefault(false);
    mask_j_x = options["mask_j_x"]
	          .doc("mask Jpar in x at boundary with tanh function")
	          .withDefault(false);
    mask_flag_j = options["mask_flag_j"]
	          .doc("mask flag, 0: mask on left boundary; 1: right; others: both on left and right")
		  .withDefault(1);
    mask_phi_x = options["mask_phi_x"]
	          .doc("mask phi in x at boundary with tanh function")
		  .withDefault(false);
    mask_flag_phi = options["mask_flag_phi"]
	          .doc("mask flag, 0: mask on left boundary; 1: right; others: both on left and right")
		  .withDefault(0);
    mask_length = options["mask_length"]
	          .doc("the center of tanh function")
		  .withDefault(0.1);
    mask_width = options["mask_width"]
	          .doc("the width of tanh function")
		  .withDefault(0.03);

    // Jpar boundary region
    jpar_bndry_width = options["jpar_bndry_width"].withDefault(-1);

    // Parallel differencing
    parallel_lr_diff = options["parallel_lr_diff"].withDefault(false);
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
    spitzer_resist = options["spitzer_resist"].doc("Use Spitzer resistivity?").withDefault(false);

    // Inner boundary damping
    damp_width = options["damp_width"].withDefault(0);
    damp_t_const = options["damp_t_const"].withDefault(0.1);

    // Viscosity and hyper-viscosity
    viscos_par = options["viscos_par"].doc("Parallel viscosity").withDefault(-1.0);
    viscos_perp = options["viscos_perp"].doc("Perpendicular viscosity").withDefault(-1.0);
    hyperviscos = options["hyperviscos"].doc("Radial hyperviscosity").withDefault(-1.0);

    diffusion_par = options["diffusion_par"].doc("Parallel temperature diffusion").withDefault(-1.0);
    diffusion_perp = options["diffusion_perp"].doc("Perpendicular temperature diffusion").withDefault(-1.0);
    diff_par_flutter =  options["diff_par_flutter"].doc("add magnetic flutter terms").withDefault(false);
    full_sbc =  options["full_sbc"].withDefault(false);
    fluxlimit =  options["fluxlimit"].withDefault(false);
    hyperdiff_par_n4 =  options["hyperdiff_par_n4"]
	                .doc("4th Parallel density diffusion")
			.withDefault(-1.0);
    hyperdiff_par_ti4 = options["hyperdiff_par_ti4"]
                        .doc("4th Parallel ion temperature diffusion")
                        .withDefault(-1.0);
    hyperdiff_par_te4 = options["hyperdiff_par_te4"]
                        .doc("4th Parallel electron temperature diffusion")
                        .withDefault(-1.0);
    hyperdiff_par_u4 = options["hyperdiff_par_u4"]
                       .doc("parallel hyper-viscous diffusion for vorticity")
                       .withDefault(-1.0);
    hyperdiff_par_v4 = options["hyperdiff_par_v4"]
                       .doc("4th order Parallel ion velocity diffusion (< 0 = none)")
                       .withDefault(-1.0);
    hyperdiff_par_apar4 = options["hyperdiff_par_apar4"].withDefault(-1.0);
    hyperdiff_perp_n4 = options["hyperdiff_perp_n4"]                  
	               .doc("4th order Perpendicular density diffusion")
	               .withDefault(-1.0);
    hyperdiff_perp_ti4 = options["hyperdiff_perp_ti4"]                  
	               .doc("4th order Perpendicular ion temperature diffusion")
	               .withDefault(-1.0);
    hyperdiff_perp_te4 = options["hyperdiff_perp_te4"]                  
	               .doc("4th order Perpendicular electron temperature diffusion")
	               .withDefault(-1.0);
    hyperdiff_perp_v4 = options["hyperdiff_perp_v4"]                  
	               .doc("4th order Perpendicular ion parallel velocity diffusion")
	               .withDefault(-1.0);
    hyperdiff_perp_u4 = options["hyperdiff_perp_u4"]                  
	               .doc("4th order Perpendicular vorticity diffusion")
	               .withDefault(-1.0);
    hyperdiff_perp_apar4 = options["hyperdiff_perp_apar4"]                  
	               .withDefault(-1.0);
    neoclassic_i = options["neoclassic_i"].doc("switch for ion neoclassical transport").withDefault(false);
    neoclassic_e = options["neoclassic_e"].doc("switch for electron neoclassical transport").withDefault(false);
    epsilon  = options["epsilon"].withDefault(0.2);    // the value of reverse aspect ratio

    // output terms
    output_Teterms  = options["output_Teterms"].withDefault(false);
    output_Titerms  = options["output_Titerms"].withDefault(false);
    output_Tevegradte  = options["output_Tevegradte"].withDefault(false);
    output_transfer  = options["output_transfer"].withDefault(false);
    output_ohm  = options["output_ohm"].withDefault(false);
    output_flux_par  = options["output_flux_par"].withDefault(false);
    output_vradial  = options["output_vradial"].withDefault(false);
    output_Tevegradte  = options["output_Tevegradte"].withDefault(false);
    output_qparcompare  = options["output_qparcompare"].withDefault(false);

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

    // right edge sink factor in Te
    sink_Ter = options["sink_Ter"].withDefault(-1.0);
    ste_widthr = options["ste_widthr"].withDefault(0.06);
    ste_lengthr = options["ste_lengthr"].withDefault(0.15);

    // Compressional terms
    phi_curv = options["phi_curv"].doc("Compressional ExB terms").withDefault(true);
    g = options["gamma"].doc("Ratio of specific heats").withDefault(5.0 / 3.0);

    if (diffusion_par < 0. && output_flux_par) {
      output_flux_par = false;
      output.write("No parallel thermal conduction. Set 'output_flux_par' to be false.\n");
      if (diff_par_flutter) {
        diff_par_flutter = false;
        output.write("No parallel thermal conduction. Set 'diff_par_flutter' to be false.\n");
      }
    }

    if (!nonlinear) {
      if (output_transfer) {
        output_transfer = false;
        output.write("Linear simulation! Set output_transfer to false.\n");
      }
      if (output_ohm) {
        output_ohm = false;
        output.write("Linear simulation! Set output_ohm to false.\n");
      }
      if (output_flux_par) {
        output_flux_par = false;
        output.write("Linear simulation! Set output_flux_par to false.\n");
      }
      if (output_vradial) {
        output_vradial = false;
        output.write("Linear simulation! Set output_vradial to false.\n");
      }
    }

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

    /////////////////////////////////////////////////////////////////
    // Initialize variables
    U0 = 0.0, Vipar0 = 0.0, Vepar0 = 0.0;

    U = 0.0, Ni = 0.0, Ti = 0.0, Te = 0.0;
    if (emass) {
      Ajpar = 0.0;
    } else {
      if (evolve_psi) {
        Psi = 0.0;
      } else {
        Apar = 0.0;
      }
    }

    /////////////////////////////////////////////////////////////////
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

    /////////////////////////////////////////////////////////////////
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
    Nbar *= density_unit / density;

    Tau_ie = Tibar / Tebar;

    Va = sqrt(Bbar * Bbar / (MU0 * Mi * Nbar * density));

    Tbar = Lbar / Va;

    output.write("Normalisations: Bbar = {:e} T   Lbar = {:e} m\n", Bbar, Lbar);
    output.write("\t Va = {:e} m/s   Tbar = {:e} s\n", Va, Tbar);
    output.write("\t Nbar = {:e} * {:e} m^-3\n", Nbar, density);
    output.write("\t Pibar = {:e} Pa   Pebar = {:e} Pa\n", ee * Tibar * Nbar * density, ee * Tebar * Nbar * density);
    output.write("\t Tibar = {:e} eV   Tebar = {:e} eV    Ti/Te = {:e}\n", Tibar, Tebar, Tau_ie);
    output.write("Resistivity\n");
    output.write("\t etabar = {:e} [Ohm m]\n", MU0 * Va * Lbar);
    output.write("\t mi = {:e} \n", Mi);

    if (emass) {
      delta_e = 5.31e5 / sqrt(Nbar * density / 1e6) / (Lbar * 100.0) * emass_inv;
      delta_e_inv = 1.e0 / delta_e / delta_e;
      gyroAlv = 1.602e-19 * Bbar * Tbar / Mi;
      output.write("delta_e = {:e}    wci*T_A = {:e}\n", delta_e, gyroAlv);
    }

    if (thermal_force || eHall) {
      Psipara1 = KB * Tebar * eV_K / ee / Bbar / Lbar / Va;
      output.write("Psipara1 = {:e}   AA = {:e}\n", Psipara1, AA);
    }

    Upara0 = KB * Tebar * eV_K / (Zi * ee * Bbar * Va * Lbar);
    Upara1 = KB * Tebar * eV_K / Mi / Va / Va;
    output.write("vorticity cinstant: Upara0 = {:e}     Upara1 = {:e}\n", Upara0, Upara1);

    if (gyroviscous) {
      Upara2 = KB * Tibar * eV_K / (Zi * ee * Bbar * Lbar * Va);
      Upara3 = 1.0;
      output.write("Upara2 = {:e}     Upara3 = {:e}\n", Upara2, Upara3);
    }

    if ((diamag && continuity) || energy_flux) {
      Nipara1 = KB * Tibar * eV_K / (Zi * ee * Bbar * Lbar * Va);
      Tipara2 = Nipara1;
      Tipara3 = Mi / Zi / ee / Bbar / Tbar;
      Tepara2 = KB * Tebar * eV_K / (ee * Bbar * Lbar * Va);
      Tepara3 = Bbar / (ee * MU0 * Nbar * density * Lbar * Va);
      output.write("Nipara1 = {:e}     Tipara2 = {:e}\n", Nipara1, Tipara2);
      output.write("Tepara2 = {:e}     Tepara3 = {:e}\n", Tepara2, Tepara3);
    }

    if (energy_exch) {
      Tepara4 = Bbar * Bbar / (MU0 * KB * Nbar * density * Tebar * eV_K);
      output.write("energy exchange constant:   Tepara4 = {:e}\n", Tepara4);
    }

    if (compress0) {
      output.write("Including compression (Vipar) effects\n");
      Vipara = MU0 * KB * Nbar * density * Tebar * eV_K / (Bbar * Bbar);
      // Vepara = Bbar / (MU0 * Zi * ee * Nbar * density * Lbar * Va);
      Vepara = Bbar / (MU0 * ee * Nbar * density * Lbar * Va);
      output.write("Normalized constant for Vipar :   Vipara = {:e}\n", Vipara);
      output.write("Normalized constant for Vepar :   Vepara = {:e}\n", Vepara);
    } else {
      include_vipar = false;
    }

    if (diffusion_par > 0.0 || diffusion_perp > 0.0) {
      Tipara1 = 2.0 / 3.0 / (Lbar * Va);
      Tepara1 = Tipara1;
    }

    if (vac_lund > 0.0) {
      output.write(" Vacuum  Tau_R = {:e} s   eta = {:e} Ohm m\n", vac_lund * Tbar, MU0 * Lbar * Lbar / (vac_lund * Tbar));
      vac_resist = 1. / vac_lund;
    } else {
      output.write("Vacuum  - Zero resistivity -\n");
      vac_resist = 0.0;
    }
    if (core_lund > 0.0) {
      output.write("Core    Tau_R = {:e} s   eta = {:e} Ohm m\n", core_lund * Tbar, MU0 * Lbar * Lbar / (core_lund * Tbar));
      core_resist = 1. / core_lund;
    } else {
      output.write("Core    - Zero resistivity -\n");
      core_resist = 0.0;
    }

    if (hyperresist > 0.0) {
      output.write("Hyper-resistivity coefficient: {:e}\n", hyperresist);
    }

    if (ehyperviscos > 0.0) {
      output.write("electron Hyper-viscosity coefficient: {:e}\n", ehyperviscos);
    }

    if (hyperviscos > 0.0) {
      output.write("Hyper-viscosity coefficient: {:e}\n", hyperviscos);
      dump.add(hyper_mu_x, "hyper_mu_x", 1);
    }

    if (diffusion_par > 0.0) {
      output.write("diffusion_par: {:e}\n", diffusion_par);
      dump.add(diffusion_par, "diffusion_par", 0);
    }
     
    if (diffusion_perp > 0.0) {
      output.write("diffusion_perp: {:e}\n", diffusion_perp);
      dump.add(diffusion_perp, "diffusion_perp", 0);
    }

    // 4th order diffusion of p
    if (hyperdiff_par_n4 > 0.0) {
      output.write("hyperdiff_par_n4: {:e}\n", hyperdiff_par_n4);
      dump.add(hyperdiff_par_n4, "hyperdiff_par_n4", 0);
    }

    // 4th order diffusion of Ti
    if (hyperdiff_par_ti4 > 0.0) {
      output.write("hyperdiff_par_ti4: {:e}\n", hyperdiff_par_ti4);
      dump.add(hyperdiff_par_ti4, "hyperdiff_par_ti4", 0);
    }

    // 4th order diffusion of Te
    if (hyperdiff_par_te4 > 0.0) {
      output.write("hyperdiff_par_te4: {:e}\n", hyperdiff_par_te4);
      dump.add(hyperdiff_par_te4, "hyperdiff_par_te4", 0);
    }

    // 4th order diffusion of Vipar
    if (hyperdiff_par_v4 > 0.0) {
      output.write("hyperdiff_par_v4: {:e}\n", hyperdiff_par_v4);
      dump.add(hyperdiff_par_v4, "hyperdiff_par_v4", 0);
    }
    
    if (hyperdiff_par_apar4 > 0.0) {
      output.write("parallel hyperdiffusion_apar4: {:e}\n", hyperdiff_par_apar4);
      dump.add(hyperdiff_par_apar4, "hyperdiff_par_apar4", 0);
    }

    // parallel hyper-viscous diffusion for vorticity
    if (hyperdiff_par_u4 > 0.0) {
      output.write("hyperdiff_par_u4: {:e}\n", hyperdiff_par_u4);
      dump.add(hyperdiff_par_u4, "hyperdiff_par_u4", 0);
    }
  
    if (hyperdiff_perp_n4 > 0.0) {
      output.write("perpendicular hyperdiffusion_n4: {:e}\n", hyperdiff_perp_n4);
      dump.add(hyperdiff_perp_n4, "hyperdiff_perp_n4", 0);
    }

    if (hyperdiff_perp_ti4 > 0.0) {
      output.write("perpendicular hyperdiffusion_ti4: {:e}\n", hyperdiff_perp_ti4);
      dump.add(hyperdiff_perp_ti4, "hyperdiff_perp_ti4", 0);
    }
    
    if (hyperdiff_perp_te4 > 0.0) {
      output.write("perpendicular hyperdiffusion_te4: {:e}\n", hyperdiff_perp_te4);
      dump.add(hyperdiff_perp_te4, "hyperdiff_perp_te4", 0);
    }

    if (hyperdiff_perp_v4 > 0.0) {
      output.write("perpendicular hyperdiffusion_v4: {:e}\n", hyperdiff_perp_v4);
      dump.add(hyperdiff_perp_v4, "hyperdiff_perp_v4", 0);
    }

    if (hyperdiff_perp_apar4 > 0.0) {
      output.write("perpendicular hyperdiffusion_apar4: {:e}\n", hyperdiff_perp_apar4);
      dump.add(hyperdiff_perp_apar4, "hyperdiff_perp_apar4", 0);
    }

    if (hyperdiff_perp_u4 > 0.0) {
      output.write("perpendicular hyperdiffusion_u4: {:e}\n", hyperdiff_perp_u4);
      dump.add(hyperdiff_perp_u4, "hyperdiff_perp_u4", 0);
    }

    if (sink_vp > 0.0) {
      output.write("sink_vp(rate): {:e}\n", sink_vp);
      dump.add(sink_vp, "sink_vp", 1);

      output.write("sp_width(%%): {:e}\n", sp_width);
      dump.add(sp_width, "sp_width", 1);

      output.write("sp_length(%%): {:e}\n", sp_length);
      dump.add(sp_length, "sp_length", 1);
    }
    
    if (limit_jacobi) {
      Bpxy = field_larger(Bpxy, bpxy_constraint);
      // Bpxy = smooth_xy(Bpxy);
      dump.add(Bpxy, "Bpxy", 0);
      if (const_bp > 0.)
        Bpxy = const_bp;
      if (hthe_constraint > 0.) {
        hthe = field_larger(hthe, hthe_constraint);
        dump.add(hthe, "hthe", 0);
      }
    }

    J0 = MU0 * Lbar * J0 / Bbar;
    // P0 = P0 / (KB * (Tibar + Tebar) * eV_K / 2. * Nbar * density);
    P0 = P0 / (KB * (Tebar) * eV_K * Nbar * density);

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
      if (n0_p0_0p3) {
        output.write("N0 ~ P0^0.3 used!\n");
        // n0 = n0_height*(P0/P00)^0.3*density
        // this relation is used in the bootstrap current calculation
        // when generating the cbm18_dens_ne* serial grids
        // P00: P0 at magnetic axis, for cbm18_dens6, it's 23072.3
        // densn: P00 = 23072.3*n/6.
        BoutReal P00; // P0 at magnetic axis
        P00 = options["P00"].withDefault(23072.3);
        output.write("n0=n0_height*(P0/P00)^0.3*density is used, P00 is required!\n");
        P00 = P00 / (ee * Tebar * Nbar * density);
        N0 = n0_height * (pow(P0 / P00, 0.3));
      } else {
	      N0 = N0tanh(n0_height * Nbar, n0_ave * Nbar, n0_width, n0_center, n0_bottom_x);
      }
      Ti0 = P0 / N0 / (1.0 + Zi);
      Te0 = Ti0;
      N_imp0 = 0;
      Ne0 = Zi * N0;
    } else if (T0_fake_prof) {
      Ti0 = Tconst;
      Te0 = Ti0;
      N0 = P0 / (Ti0 * Tau_ie + Te0);
      N_imp0 = 0;
      Ne0 = Zi * N0;

  } else if (load_2d_bkgd) {

    if (mesh->get(P0, "P_2D")) { // [Pa]
      output.write("Error: Cannot read P_2D from grid\n");
      return 1;
    }

    if (mesh->get(N0, "Nd+_2D")) { // [1e20 #/m^3]
      output.write("Error: Cannot read Nd+_2D from grid\n");
      return 1;
    }

    if (mesh->get(Ti0, "Td+_2D")) { // [eV]]
      output.write("Error: Cannot read Td+_2D from grid\n");
      return 1;
    }

    if (mesh->get(Ne0, "Ne_2D")) { // [1e20 #/m^3]
      output.write("Error: Cannot read Ne_2D from grid\n");
      return 1;
    }

    if (mesh->get(Te0, "Te_2D")) { // [eV]
      output.write("Error: Cannot read Te_2D from grid\n");
      return 1;
    }

    if (mesh->get(Vipar0, "Vd+_2D")) { // [m/s]
      output.write("Error: Cannot read Vd+_2D from grid\n");
      return 1;
    }

    if (mesh->get(Vepar0, "Ve_2D")) { // [m/s]
      output.write("Error: Cannot read Ve_2D from grid\n");
      return 1;
    }

    // if (mesh->get(phi0, "phi_2D")) { // [V]
    //   output.write("Error: Cannot read phi_2D from grid\n");
    //   return 1;
    // }

    // if (mesh->get(U0, "Vort_2D")) { // [C m^-3] //NOTE(malamast): The units of vorticity used in Hermes-3 and 6-field are different. Be carefull!
    //   output.write("Error: Cannot read Vort_2D from grid\n");
    //   return 1;
    // }

    // if (load_impurity) {

    //   if (mesh->get(N_imp0, "Nd_2D")) { // [1e20 #/m^3]
    //     output.write("Error: Cannot read Nd_2D from grid\n");
    //     return 1;
    //   }

    //   if (mesh->get(T_imp0, "Td_2D")) { // [1e20 #/m^3]
    //     output.write("Error: Cannot read Nd_2D from grid\n");
    //     return 1;
    //   }

    // }

    MPI_Barrier(BoutComm::get());

    P0 = P0 / (KB * (Tebar)*eV_K * Nbar * density);

    N0 /= Nbar;
    Ti0 /= Tibar;
    Ne0 /= Nbar;
    Te0 /= Tebar;

    Vipar0 /= Va;
    Vepar0 /= Va;

    // Ne0 = Zi * N0; // quasi-neutral condition
    N_imp0 = 0;
    // N_imp0 /= Nbar;
    // T_imp0 /= Tibar;
    // P_imp0 = N_imp0 * T_imp0;      

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

      if (impurity_prof) {
        if (mesh->get(Ne0, "Neexp")) { // N_e0
          output.write("Error: Cannot read Ne0 from grid\n");
          return 1;
        }
        Ne0 /= Nbar;
        if (load_impurity) {
          if (mesh->get(N_imp0, "N_imp")) {
            output.write("Error: Cannot read N_imp from grid\n");
            return 1;
          }
          N_imp0 /= Nbar;
        } else{
          N_imp0 = (Ne0 - Zi * N0) / Z_imp;
        }
        if (min(N_imp0, true) <= 0.) {
          output.write("Error: Impurity density has negative value.\n");
          if (Nimp_lowlimit) {
            BoutReal min_Nimp = max(N_imp0, true) * 0.01;
            output.write("Use Nimp_lowlimit of {:d} to modify the smallest value.\n", min_Nimp);
            for (jx = 0; jx < mesh->LocalNx; jx++)
              for (jy = 0; jy < mesh->LocalNy; jy++)
                if (N_imp0(jx,jy) <= min_Nimp) 
                  N_imp0(jx,jy) = min_Nimp;
          } else {
            return 1;
          }
        }
        T_imp0 = Ti0;
        P_imp0 = N_imp0 * T_imp0;
        // P_imp0 = P0 - N0*Ti0*Tau_ie - Ne0*Te0;
              // T_imp0 = P_imp0/N_imp0;
        if (min(T_imp0, true) <= 0.) {
          output.write("Error: Impurity temperature has negative value.\n");
          if (Nimp_lowlimit) {
            BoutReal min_Timp = max(T_imp0, true) * 0.01;
            output.write("Use Timp_lowlimit of {:d} to modify the smallest value.\n", min_Timp);
            for (jx = 0; jx < mesh->LocalNx; jx++)
              for (jy = 0; jy < mesh->LocalNy; jy++)
                if (T_imp0(jx,jy) <= min_Timp)
                  T_imp0(jx,jy) = min_Timp;
            P_imp0 = N_imp0 * T_imp0;
          } else
            return 1;
        }
      } else if (quasi_neutral_Ni) {
        N0 = P0 / (Tau_ie * Ti0 + Te0);
        Ne0 = Zi * N0;
        N_imp0 = 0;
      } else {
        Ne0 = Zi * N0; // quasi-neutral condition
        N_imp0 = 0;
      }
    }

    J0 *= J0_factor;
    P0 *= P0_factor;

    Pi0 = N0 * Ti0;
    Pe0 = Ne0 * Te0;

    jpar1.setBoundary("J");
    u_tmp1.setBoundary("U");
    ti_tmp2.setBoundary("Ti");
    te_tmp2.setBoundary("Te");
    phi_tmp.setBoundary("phi");

    nu_e.setBoundary("kappa");
    if (spitzer_resist) {
      eta.setBoundary("kappa");
    }
    if (diffusion_par > 0.0 || diffusion_perp > 0.0) {
      nu_i.setBoundary("kappa");
      vth_i.setBoundary("kappa");
      vth_e.setBoundary("kappa");
      kappa_par_i.setBoundary("kappa");
      kappa_par_e.setBoundary("kappa");
      kappa_perp_i.setBoundary("kappa");
      kappa_perp_e.setBoundary("kappa");
      if (diff_par_flutter) {
        bracket1i.setBoundary("Ti");
        gradpar_ti.setBoundary("Ti");
        bracket1e.setBoundary("Te");
        gradpar_te.setBoundary("Te");
      }

      // if (Landau) {
      //   q_par_e.setBoundary("q_par_e");
      //   q_par_i.setBoundary("q_par_i");
      //   q_par_e = 0.;
      //   q_par_i = 0.;
      //   dump.add(q_par_e, "q_par_e", 1);
      //   dump.add(q_par_i, "q_par_i", 1);
      //   if (output_qparcompare) {
      //     q_par_fl.setBoundary("q_par_e");
      //     dump.add(q_par_fl, "q_par_fl", 1);
      //     q_par_landau.setBoundary("q_par_e");
      //     dump.add(q_par_landau, "q_par_landau", 1);
      //   }
      //   if (Landau_coll) {
      //     kappa_i.setBoundary("kappa");
      //     kappa_e.setBoundary("kappa");
      //   }
      //   SBC_value_i.setBoundary("kappa");
      //   SBC_value_e.setBoundary("kappa");
      // }      

    }

    if (neoclassic_i || neoclassic_e) {
      xii_neo.setBoundary("kappa");
      xie_neo.setBoundary("kappa");
      rho_i.setBoundary("kappa");
      rho_e.setBoundary("kappa");
      tmpddx2.setBoundary("kappa");
      partf_neo_i.setBoundary("Ni");
      heatf_neo_i.setBoundary("Ti");
      heatf_neo_e.setBoundary("Te");
    }

    if (parallel_viscous && compress0) {
      eta_i0.setBoundary("Ti");
      pi_ci.setBoundary("Ti");
    }

    if (gyroviscous) {
      Dperp2Phi0.setLocation(CELL_CENTRE);
      Dperp2Phi0.setBoundary("phi");
      dump.add(Dperp2Phi0, "Dperp2Phi0", 1);
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
      if (impurity_prof && impurity_gyro) {
        Upara_imp = (A_imp * N_imp0) / (AA * N0);
        output.write("Max Upara_imp = {:e}, Min Upara_imp = {:e}\n", max(Upara_imp, true), min(Upara_imp, true));
        
        Dperp2Pimp0.setLocation(CELL_CENTRE);
        Dperp2Pimp0.setBoundary("P");
        bracketPhiPimp0.setLocation(CELL_CENTRE);
        bracketPhiPimp0.setBoundary("P");
      }
      if (nonlinear) {
        GradPhi2.setLocation(CELL_CENTRE);
        GradPhi2.setBoundary("phi");
	      bracketPhiP.setLocation(CELL_CENTRE);
        bracketPhiP.setBoundary("P");
      }
    }

    if (output_transfer) {
      T_R.setLocation(CELL_YLOW);
      T_R.setBoundary("phi");
      T_M.setLocation(CELL_YLOW);
      T_M.setBoundary("phi");
      T_ID.setLocation(CELL_YLOW);
      T_ID.setBoundary("phi");
      T_C.setLocation(CELL_YLOW);
      T_C.setBoundary("P");
      T_G.setLocation(CELL_YLOW);
      T_G.setBoundary("P");

      dump.add(T_R, "T_R", 1);
      dump.add(T_M, "T_M", 1);
      dump.add(T_ID, "T_ID", 1);
      dump.add(T_C, "T_C", 1);
      dump.add(T_G, "T_G", 1);
    }

    if (output_ohm) {
      ohm_phi.setBoundary("Psi");
      ohm_hall.setBoundary("Psi");
      ohm_thermal.setBoundary("Psi");

      dump.add(ohm_phi, "ohm_phi", 1);
      dump.add(ohm_hall, "ohm_hall", 1);
      dump.add(ohm_thermal, "ohm_thermal", 1);
    }

    if (output_flux_par && diffusion_par > 0.) {
      heatf_par_i.setBoundary("Ti");
      heatf_par_e.setBoundary("Te");

      // dump.add(gamma_par_i, "gamma_i", 1);
      dump.add(heatf_par_i, "heatflux_par_i", 1);
      dump.add(heatf_par_e, "heatflux_par_e", 1);
      if (diff_par_flutter) {
        heatf_par_flutter_i.setBoundary("Ti");
        heatf_par_flutter_e.setBoundary("Te");

        // dump.add(gamma_par_i, "gamma_i", 1);
        dump.add(heatf_par_flutter_i, "heatflux_par_flutter_i", 1);
        dump.add(heatf_par_flutter_e, "heatflux_par_flutter_e", 1);
      }
    }
    
    if (output_vradial) {
      Vexb.covariant = false;
      Vexb.setBoundary("Vipar");
      Vbtilde.covariant = false;
      Vbtilde.setBoundary("Vipar");
      // Vbti_par.setBoundary("Vipar");
      // Vbte_par.setBoundary("Vipar");

      dump.add(Vexb, "Vexb", 1);
      dump.add(Vbtilde, "Vbtild", 1);
    }

    if (neutral) {
      dump.add(nu_iz, "nu_iz", 1);
      dump.add(nu_rc, "nu_rc", 1);
      dump.add(nu_cx, "nu_cx", 1);
      dump.add(Dn, "Dn", 1);
      dump.add(Dn_fl, "Dn_fl", 1);
      dump.add(Dn1, "Dn1", 1);
      dump.add(Gamma_nn, "Gamma_nn", 1);
      dump.add(etan, "etan", 1);
      dump.add(etan_perp, "etan_perp", 1);
      dump.add(Pn, "Pn", 1);
      dump.add(Sn, "Sn", 1);
      dump.add(Sn_ext, "Sn_ext", 1);
      dump.add(Sv, "Sv", 1);
      dump.add(S_diss, "S_diss", 1);
      dump.add(nu_diss, "nu_diss", 1);
      dump.add(Sgas, "Sgas", 1);
      dump.add(sigma_cx, "sigma_cx", 1);
      if (with_vipar) {
        dump.add(term1, "term1", 1);
      }
      dump.add(term2, "term2", 1);
      dump.add(term3, "term3", 1);
    }

    if (fix_fraction_imp) {
      dump.add(Limp, "Limp", 1);
      dump.add(Srad, "Srad", 1);
      dump.add(Wrad, "Wrad", 1);
    }
    
    if (mask_j_x) {
      mask_jx1d = mask_x_1d(false, mask_flag_j, mask_width, mask_length);
      // dump.add(mask_jx1d, "mask_jx1d", 0);
    }
    if (mask_phi_x) {
      mask_px1d = mask_x_1d(false, mask_flag_phi, mask_width, mask_length);
      // dump.add(mask_px1d, "mask_px1d", 0);
    }

    BoutReal pnorm = max(P0, true); // Maximum over all processors

    vacuum_pressure *= pnorm; // Get pressure from fraction
    vacuum_trans *= pnorm;

    // Transitions from 0 in core to 1 in vacuum
    vac_mask = (1.0 - tanh((P0 - vacuum_pressure) / vacuum_trans)) / 2.0;
    LnLambda = 24.0 - log(pow(Zi * Nbar * density * Ne0 / 1.e6, 0.5) / (Tebar * Te0));  
    output.write("\tlog Lambda: {:e} -> {:e} \n", min(LnLambda), max(LnLambda));

    if (spitzer_resist) {
      // Use Spitzer resistivity
      output.write("\n\tSpizter parameters");
      FZ = (1. + 1.198 * Zi + 0.222 * Zi * Zi) / (1. + 2.996 * Zi + 0.753 * Zi * Zi);
      eta = FZ * 1.03e-4 * Zi * LnLambda * (pow(Te0 * Tebar, -1.5)); // eta in Ohm-m. NOTE: ln(Lambda) = 20
      output.write("\tSpitzer resistivity: {:e} -> {:e} [Ohm m]\n", min(eta),
                   max(eta));
      eta /= MU0 * Va * Lbar;
      output.write("\t -> Lundquist {:e} -> {:e}\n", 1.0 / max(eta),
                   1.0 / min(eta));
      dump.add(eta, "eta", 1);
    } else {
      // transition from 0 for large P0 to resistivity for small P0
      if (vac_lund > 0.0) {
        output.write("Vacuum  Tau_R = {:e} s   eta = {:e} Ohm m\n", vac_lund * Tbar, MU0 * Lbar * Lbar / (vac_lund * Tbar));
	      vac_resist = 1. / vac_lund;
      } else {
        output.write("Vacuum  - Zero resistivity -\n");
	      vac_resist = 0.0;
      }
      if (core_lund > 0.0) {
        output.write("Core    Tau_R = {:e} s   eta = {:e} Ohm m\n", core_lund * Tbar, MU0 * Lbar * Lbar / (core_lund * Tbar));
        core_resist = 1. / core_lund;
      } else {
        output.write("Core    - Zero resistivity -\n");
        core_resist = 0.0;
      }
      eta = core_resist + (vac_resist - core_resist) * vac_mask;
      dump.add(eta, "eta", 0);
    }

    if (diffusion_par > 0.0 || diffusion_perp > 0.0 || neoclassic_i || neoclassic_e) {
      if (q95_input > 0) {
        q95 = q95_input; // use a constant for test
      } else {
        if (local_q) {
          q95 = abs(hthe * Btxy / (Bpxy));
        } else {
          output.write("\tUsing q profile from grid.\n");
          if (mesh->get(q95, "q")) {
            output.write("Cannot get q profile from grid!\nPlease run addqprofile.pro first\n");
            return 1;
          }
        }
      }
      output.write("\tlocal max q: {:e}\n", max(q95));
      output.write("\tlocal min q: {:e}\n", min(q95));
    }

    nu_e = 2.91e-6 * LnLambda * ((N0) * Nbar * density / 1.e6) * pow(Te0 * Tebar, -1.5); // nu_e in 1/S.
    output.write("\telectron collision rate: {:e} -> {:e} [1/s]\n", min(nu_e), max(nu_e));
    // nu_e.applyBoundary();
    // mesh->communicate(nu_e);

    if (diffusion_par > 0.0 || diffusion_perp > 0.0 || parallel_viscous || neoclassic_i || neoclassic_e) {

      output.write("\tion thermal noramlized constant: Tipara1 = {:e}\n", Tipara1);
      output.write("\telectron normalized thermal constant: Tepara1 = {:e}\n", Tepara1);
      // Use Spitzer thermal conductivities
      nu_i = 4.80e-8 * (Zi * Zi * Zi * Zi / sqrt(AA)) * LnLambda * ((N0) * Nbar * density / 1.e6) * pow(Ti0 * Tibar, -1.5); // nu_i in 1/S.
      output.write("\tion collision rate: {:e} -> {:e} [1/s]\n", min(nu_i), max(nu_i));

      // nu_i.applyBoundary();
      // mesh->communicate(nu_i);

      vth_i = 9.79e3 * sqrt((Ti0)*Tibar / AA); // vth_i in m/S.
      output.write("\tion thermal velocity: {:e} -> {:e} [m/s]\n", min(vth_i), max(vth_i));
      // vth_i.applyBoundary();
      // mesh->communicate(vth_i);
      vth_e = 4.19e5 * sqrt((Te0)*Tebar); // vth_e in m/S.
      output.write("\telectron thermal velocity: {:e} -> {:e} [m/s]\n", min(vth_e), max(vth_e));
      // vth_e.applyBoundary();
      // mesh->communicate(vth_e);
    }

    if (parallel_viscous && compress0) {
      eta_i0 = 0.96 * Pi0 * Tau_ie / nu_i / Tbar;
      output.write("\tCoefficients of parallel viscocity: {:e} -> {:e} [kg/(m s)]\n", min(eta_i0), max(eta_i0));
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

      kappa_par_i_fl = vth_i * (q95 * Lbar) * q_alpha; // * 1.e2;
      kappa_par_e_fl = vth_e * (q95 * Lbar) * q_alpha; // * 1.e2;
      
      if (fluxlimit) {
        kappa_par_i *= kappa_par_i_fl / (kappa_par_i + kappa_par_i_fl);
        kappa_par_e *= kappa_par_e_fl / (kappa_par_e + kappa_par_e_fl);
      }

      kappa_par_i *= Tipara1 * N0;
      output.write("\tUsed normalized ion thermal conductivity: {:e} -> {:e} \n",
                   min(kappa_par_i), max(kappa_par_i));
      // kappa_par_i.applyBoundary();
      // mesh->communicate(kappa_par_i);
      //kappa_par_e *= Tepara1 * N0 / Zi;
      kappa_par_e *= Tepara1 * Ne0;
      output.write("\tUsed normalized electron thermal conductivity: {:e} -> {:e} \n",
                   min(kappa_par_e), max(kappa_par_e));
      // kappa_par_e.applyBoundary();
      // mesh->communicate(kappa_par_e);

      dump.add(kappa_par_i, "kappa_par_i", 1);
      dump.add(kappa_par_e, "kappa_par_e", 1);

      kappa_par_i_lin = 1.0;
      kappa_par_e_lin = 1.0;
      for (jx = 0; jx < mesh->LocalNx; jx++)
        for (jy = 0; jy < mesh->LocalNy; jy++) {
          kappa_par_i_lin(jx,jy) = kappa_par_i(jx,jy,0);
          kappa_par_e_lin(jx,jy) = kappa_par_e(jx,jy,0);
        }

      // if (Landau) {
      //   if (Landau_coll) {
      //     output.write("Collisional Landau damping closure used!\n");
      //     kappa_i = 0.5 * nu_i / vth_i;
      //     kappa_i *= Lbar; // normalization
      //     kappa_e = 0.5 * nu_e / vth_e;
      //     kappa_e *= Lbar; // normalization
      //   } else {
      //     output.write("Collisionless Landau damping closure used!\n");
      //     kappa_0 = 2. / (20. * (2 * PI * q95_input));
      //     output.write("\tkappa_0 = %e\n", kappa_0);
      //   }
      // }
    }

    if (diffusion_perp > 0.0) {
      omega_ci = Zi * ee * Bbar * B0 / Mi;
      omega_ce = ratio_pe * AA * ee * Bbar * B0 / Mi;
      kappa_perp_i = 2.0 * vth_i * vth_i * nu_i / (omega_ci * omega_ci); // * 1.e4;
      kappa_perp_e = 4.7 * vth_e * vth_e * nu_e / (omega_ce * omega_ce); // * 1.e4;

      output.write("\tion perp thermal conductivity: {:e} -> {:e} [m^2/s]\n", min(kappa_perp_i), max(kappa_perp_i));
      output.write("\telectron perp thermal conductivity: {:e} -> {:e} [m^2/s]\n", min(kappa_perp_e), max(kappa_perp_e));
      output.write("\tnormalized perp ion thermal conductivity: {:e} -> {:e}\n", min(kappa_perp_i * Tipara1), max(kappa_perp_i * Tipara1));
      output.write("\tnormalized perp electron thermal conductivity: {:e} -> {:e}\n", min(kappa_perp_e * Tepara1), max(kappa_perp_e * Tepara1));

      kappa_perp_i_fl = q_alpha * vth_i * q95 * Lbar; // * 1.e4;
      kappa_perp_e_fl = q_alpha * vth_e * q95 * Lbar; // * 1.e4;

      if (fluxlimit) {
        kappa_perp_i *= kappa_perp_i_fl / (kappa_perp_i + kappa_perp_i_fl);
        kappa_perp_e *= kappa_perp_e_fl / (kappa_perp_e + kappa_perp_e_fl);
      }
      kappa_perp_i *= Tipara1 * N0;
      output.write("\tUsed normalized ion perp thermal conductivity: {:e} -> {:e}\n", min(kappa_perp_i), max(kappa_perp_i));
      // kappa_perp_i.applyBoundary();
      // mesh->communicate(kappa_perp_i);
      kappa_perp_e *= Tepara1 * Ne0;
      output.write("\tUsed normalized electron perp thermal conductivity: {:e} -> {:e}\n", min(kappa_perp_e), max(kappa_perp_e));
      // kappa_perp_e.applyBoundary();
      // mesh->communicate(kappa_perp_e);

      dump.add(kappa_perp_i, "kappa_perp_i", 1);
      dump.add(kappa_perp_e, "kappa_perp_e", 1);
    }

    if (neoclassic_i) {
      rho_i = 1.02e-4 * sqrt(AA * Ti0 * Tibar) / B0 / Bbar / Zi;
      // Dri_neo = (1. + 1.6 * q95) * (1. + Tau_ie) * nu_i * rho_i * rho_i;
      xii_neo = (q95 * q95) / epsilon / sqrt(epsilon) * nu_i * (pow(rho_i, 2.));
      output.write("\tion neoclassic heat transport: {:e} -> {:e} [m^2/s]\n", min(xii_neo), max(xii_neo));

      // Dri_neo *= 3./2.* Tipara1;
      xii_neo *= Tipara1;
      // Dri_neo.applyBoundary();
      // xii_neo.applyBoundary();
      output.write("\tNormalized ion neoclassic heat transport: {:e} -> {:e} [m^2/s]\n", min(xii_neo), max(xii_neo));

      // dump.add(Dri_neo, "Dri_neo", 1);
      dump.add(xii_neo, "xii_neo", 1);
      dump.add(heatf_neo_i, "heatf_neo_i", 1);
    }

    if (neoclassic_e) {
      rho_e = 2.38e-6 * sqrt(Te0 * Tebar) / B0 / Bbar;
      xie_neo = (q95 * q95) / epsilon / sqrt(epsilon) * nu_e * (pow(rho_e, 2.));
      output.write("\telectron neoclassic heat transport: {:e} -> {:e} [m^2/s]\n", min(xie_neo), max(xie_neo));

      xie_neo *= Tepara1;
      // xie_neo.applyBoundary();
      output.write("\tNormalized electron neoclassic heat transport: {:e} -> {:e} [m^2/s]\n", min(xie_neo), max(xie_neo));

      dump.add(xie_neo, "xie_neo", 1);
      dump.add(heatf_neo_e, "heatf_neo_e", 1);
    }

    /**************** CALCULATE METRICS ***************************************/

    // BoutReal sbp = 1.0; // Sign of Bp
    // if (min(Bpxy, true) < 0.0) {
    //   sbp = -1.0;
    // }    

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

    // Calculate quantities from metric tensor
    coord->geometry(); 

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
    
    //initialize sourcing profile				
    NiSource = 0.0;
    if (NiAmp > 0.) {
      NiSource.setBoundary("Ni");
      for (jz = 0; jz < mesh->LocalNz; jz++) {
        for (jy = 0; jy < mesh->LocalNy; jy++) {
          indy = mesh->getGlobalYIndex(jy);
          for (jx = 0; jx < mesh->LocalNx; jx++) {
	          indx = mesh->getGlobalXIndex(jx);
            if (mag_config == 1 || mag_config == 2)
              NiSource(jx,jy,jz) = NiAmp * exp( - ((indx - NiLoc) * (indx - NiLoc) / (2. * NiSig * NiSig)));
            if (mag_config == 3) {
              if ( (indy > jysep1 - 2) && (indy <= jysep2 + 2) )
                NiSource(jx,jy,jz) = NiAmp * exp( - ((indx - NiLoc) * (indx - NiLoc) / (2. * NiSig * NiSig)));
            }
            if (mag_config == 4) {
               if ( ((indy > jysep1 - 2) && (indy <= jysep2_1)) || ((indy > jysep1_2) && (indy <= jysep2 + 2)) )
                 NiSource(jx,jy,jz) = NiAmp * exp( - ((indx - NiLoc) * (indx - NiLoc) / (2. * NiSig * NiSig)));
	          }
	        }
	      }
      }
    }
    NiSource.applyBoundary();
    mesh->communicate(NiSource);
    SAVE_ONCE(NiSource);

  #if DEBUG_6F>0
    output.write("NiSource initialization finished\n");
  #endif

    TeSource = 0.0;
    if (TeAmp > 0.) {
      TeSource.setBoundary("Te");
      for (jz = 0; jz < mesh->LocalNz; jz++) {
        for (jy = 0; jy < mesh->LocalNy; jy++) {
          indy = mesh->getGlobalYIndex(jy);
          for (jx = 0; jx < mesh->LocalNx; jx++) {
	          indx = mesh->getGlobalXIndex(jx);
            if (mag_config == 1 || mag_config == 2)
              TeSource(jx,jy,jz) = TeAmp * exp( - ((indx - TeLoc) * (indx - TeLoc) / (2. * TeSig * TeSig)));
            if (mag_config == 3) {
              if ( (indy > jysep1 - 2) && (indy <= jysep2 + 2) )
                TeSource(jx,jy,jz) = TeAmp * exp( - ((indx - TeLoc) * (indx - TeLoc) / (2. * TeSig * TeSig)));
	          }
            if (mag_config == 4) {
               if ( ((indy > jysep1 - 2) && (indy <= jysep2_1)) || ((indy > jysep1_2) && (indy <= jysep2 + 2)) )
                 TeSource(jx,jy,jz) = TeAmp * exp( - ((indx - TeLoc) * (indx - TeLoc) / (2. * TeSig * TeSig)));
	          }
	        }
	      }
      }
    }
    TeSource.applyBoundary();
    mesh->communicate(TeSource);
    SAVE_ONCE(TeSource);
  #if DEBUG_6F>0
    output.write("TeSource initialization finished\n");
  #endif

    TiSource = 0.0;
    if (TiAmp > 0.) {
      TiSource.setBoundary("Ti");
      for (jz = 0;jz < mesh->LocalNz; jz++) {
        for (jy = 0;jy < mesh->LocalNy; jy++) {
          indy = mesh->getGlobalYIndex(jy);
          for (jx = 0; jx < mesh->LocalNx; jx++) {
	          indx = mesh->getGlobalXIndex(jx);
            if (mag_config == 1 || mag_config == 2)
              TiSource(jx,jy,jz) = TiAmp * exp( - ((indx - TiLoc) * (indx - TiLoc) / (2. * TiSig * TiSig)));
            if (mag_config == 3) {
              if ( (indy > jysep1 - 2) && (indy <= jysep2 + 2) )
                TiSource(jx,jy,jz) = TiAmp * exp( - ((indx - TiLoc) * (indx - TiLoc) / (2. * TiSig * TiSig)));
            }
            if (mag_config == 4) {
               if ( ((indy > jysep1 - 2) && (indy <= jysep2_1)) || ((indy > jysep1_2) && (indy <= jysep2 + 2)))
                 TiSource(jx,jy,jz) = TiAmp * exp( - ((indx - TiLoc) * (indx - TiLoc)/(2. * TiSig * TiSig)));
            }
	        }
	      }
      }
    }
    TiSource.applyBoundary();
    mesh->communicate(TiSource);
    SAVE_ONCE(TiSource);
  #if DEBUG_6F>0
    output.write("TiSource initialization finished\n");
  #endif

    //initialize neutral profile				
    if (neutral) {
      output.write("Solving for Nn and Vn\n");
      SOLVE_FOR(Nn);
      SOLVE_FOR(Vn);
      
      if (!restarting) {
        Nn = 1.0e-10;
        Vn = 0.0;
      }

      if (with_fueling) {
        SOLVE_FOR(Nm);
        SOLVE_FOR(Vm);
        if (!restarting) {
          Nm = 1.0e-10;
          Vm = 0.0;
	      }
      }
      
      if (with_fueling) {
        Vmx = Vm.x;
        Vm0 /= Lbar/Tbar;
        //Nm=ret_const_flux_BC(Nm, Nm0);
        //Vmx=ret_const_flux_BC(Vmx, Vm0);
        if (initial_Nm) {
          for (jz = 0; jz < mesh->LocalNz; jz++) {
            for (jy = 0; jy< mesh->LocalNy; jy++) {
              indy = mesh->getGlobalYIndex(jy);
              for (jx = 0;jx < mesh->LocalNx; jx++) {
	              indx = mesh->getGlobalXIndex(jx);
                if (mag_config == 1 || mag_config == 2) {
                  Vmx(jx,jy,jz) = Vm0*(1 - fac_A * (tanh( - (indx - SnLoc) / SnSig) + 1.)) * exp( - (indy - SnLoc_y) * (indy - SnLoc_y) / (2 * SnSig_y * SnSig_y));
                  Nm(jx,jy,jz) = Nm0*(1 - fac_A * (tanh( - (indx - SnLoc) / SnSig) + 1.)) * exp( - (indy - SnLoc_y) * (indy - SnLoc_y) / (2 * SnSig_y * SnSig_y));
	              }
                if (mag_config == 3) {
                  Vmx(jx,jy,jz) = Vm0*(1 - fac_A * (tanh( - (indx - SnLoc) / SnSig) + 1.)) * exp( - (indy - SnLoc_y) * (indy - SnLoc_y) / (2 * SnSig_y * SnSig_y));
                  Nm(jx,jy,jz) = Nm0*(1 - fac_A * (tanh( - (indx - SnLoc) / SnSig) + 1.)) * exp( - (indy - SnLoc_y) * (indy - SnLoc_y) / (2 * SnSig_y * SnSig_y));
	              }
                if (mag_config == 4) {
                  Vmx(jx,jy,jz) = Vm0*(1 - fac_A * (tanh( - (indx - SnLoc) / SnSig) + 1.)) * exp( - (indy - SnLoc_y) * (indy - SnLoc_y) / (2 * SnSig_y * SnSig_y));
                  Nm(jx,jy,jz) = Nm0*(1 - fac_A * (tanh( - (indx - SnLoc) / SnSig) + 1.)) * exp( - (indy - SnLoc_y) * (indy - SnLoc_y) / (2 * SnSig_y * SnSig_y));
	              }
	            }
	          }
	        }
	      } else {
          Vmx = 0.0;
          Nm = 1.0e-10;
	      }
        Vm.x = Vmx;
        Vm.y = 0.0;
        Vm.z = 0.0;
        mesh->communicate(Nm);
        mesh->communicate(Vmx);
        mesh->communicate(Vm);
        Nm.setBoundary("Nm");
        Vm.setBoundary("Vm");
        Nm.applyBoundary();
        Vm.applyBoundary();
        
	      Sgas = 0.;
        if (gas_puffing) {
          for (jz = 0; jz < mesh->LocalNz; jz++) {
            for (jy = 0; jy < mesh->LocalNy; jy++) {
              indy = mesh->getGlobalYIndex(jy);
              for (jx = 0; jx < mesh->LocalNx; jx++) {
	              indx = mesh->getGlobalXIndex(jx);
                if (mag_config == 1 || mag_config == 2) {
                  Sgas(jx,jy,jz) = Nm0 * (1 - fac_A * (tanh( - (indx - SnLoc) / SnSig) + 1.)) * exp( - (indy - SnLoc_y) * (indy - SnLoc_y) / (2 * SnSig_y * SnSig_y));
		            }
		            if (mag_config == 3) {
                  Sgas(jx,jy,jz) = Nm0 * (1 - fac_A * (tanh( - (indx - SnLoc) / SnSig) + 1.)) * exp( - (indy - SnLoc_y) * (indy - SnLoc_y) / (2 * SnSig_y * SnSig_y));
                }
                if (mag_config == 4) {
                  Sgas(jx,jy,jz) = Nm0 * (1 - fac_A * (tanh( - (indx - SnLoc) / SnSig) + 1.)) * exp( - (indy - SnLoc_y) * (indy - SnLoc_y) / (2 * SnSig_y * SnSig_y));
                }
	            }
	          }
	        }
	      }
        Sgas.setBoundary("Nm");
        Sgas.applyBoundary();
        mesh->communicate(Sgas);
      }
      
      // ideally Nn and Vn shall have the neumann bc at radial boundaries
      Nn.setBoundary("Nn");
      //lNn.setBoundary("Ni");
      Vn.setBoundary("Vn");
      // initial neutral profile with tanh function, could be changed easily
      if (initial_Nn) {
       for (jz = 0; jz < mesh->LocalNz; jz++) {
        for (jy = 0; jy < mesh->LocalNy; jy++) {
          indy = mesh->getGlobalYIndex(jy);
          for (jx = 0; jx < mesh->LocalNx; jx++) {
	          indx = mesh->getGlobalXIndex(jx);
            if (mag_config == 1 || mag_config == 2) {
              Nn(jx,jy,jz) = NnAmp * (1 - fac_A * (tanh( - (indx - NnLoc) / NnSig) + 1.)) * (exp((indy - NnLoc_y) / NnSig_y) + exp( - (indy - NnLoc_y) / NnSig_y)) / 2. + 1.0e-10;
	          }
	          if (mag_config == 3) {
              Nn(jx,jy,jz) = NnAmp * (1 - fac_A * (tanh( - (indx - NnLoc) / NnSig) + 1.)) * (exp((indy - NnLoc_y) / NnSig_y) + exp( - (indy - NnLoc_y) / NnSig_y)) / 2. + 1.0e-10;
            }
            if (mag_config == 4) {
              Nn(jx,jy,jz) = NnAmp * (1 - fac_A * (tanh( - (indx - NnLoc) / NnSig) + 1.)) * (exp((indy - NnLoc_y) / NnSig_y) + exp( - (indy - NnLoc_y) / NnSig_y)) / 2. + 1.0e-10;
            }
	        }
        }
       }
      } else {
        Nn = 1.0e-10;
      }
      //lNn.applyBoundary();
      //mesh->communicate(lNn);

      //Nn = NnAmp*(1.-lNn);
      //Nn = NnAmp;//*(lNn+1.0);
      Nn.applyBoundary();
      mesh->communicate(Nn);
      output.write("Max and min of normalized Nn = {:e}, {:e}.\n", max(Nn), min(Nn));
      //SAVE_ONCE(Nn);
      //lNn = log(Nn);
      //output.write("Max and min of normalized ln(Nn) = {:e}, {:e}.\n", max(lNn), min(lNn));
      Sn_ext = 0.;
      if (external_source) {
        for (jz = 0; jz < mesh->LocalNz; jz++) {
          for (jy = 0; jy < mesh->LocalNy; jy++) {
            indy = mesh->getGlobalYIndex(jy);
            for (jx = 0; jx < mesh->LocalNx; jx++) {
	            indx = mesh->getGlobalXIndex(jx);
              if (mag_config == 1 || mag_config == 2)
                Sn_ext(jx,jy,jz) = SnAmp * (1 - fac_A * (tanh( - (indx - SnLoc) / SnSig) + 1.)) * exp( - (indy - SnLoc_y) * (indy - SnLoc_y) / (2 * SnSig_y * SnSig_y));
              if (mag_config == 3) {
                if ( (indy > jysep1 - 2) && (indy <= jysep2 + 2) )
                  Sn_ext(jx,jy,jz) = SnAmp * (1 - fac_A * (tanh( - (indx - SnLoc) / SnSig) + 1.)) * exp( - (indy - SnLoc_y) * (indy - SnLoc_y) / (2 * SnSig_y * SnSig_y));
              }
              if (mag_config == 4) {
                if ( ((indy > jysep1 - 2) && (indy <= jysep2_1)) || ((indy > jysep1_2) && (indy <= jysep2 + 2)))
                  Sn_ext(jx,jy,jz) = SnAmp * (1 - fac_A * (tanh( - (indx - SnLoc) / SnSig) + 1.)) * exp( - (indy - SnLoc_y) * (indy - SnLoc_y) / (2 * SnSig_y * SnSig_y));
              }
            }
          }
        }
        Sn_ext.setBoundary("Nn");
        Sn_ext.applyBoundary();
        mesh->communicate(Sn_ext);
      }
  #if DEBUG_6F>0
    output.write("Neutral initialization finished\n");
  #endif
    }

    

    // Background Impurity radiation source term. Calculated based on the background/equilibrium profiles.
    // To be substracted from the total source term  in the equation for the perturbed Te. 
    if (fix_fraction_imp) {
      N_tmp0 = field_larger(N0, Low_limit);
      Te_tmp0 = field_larger(Te0, Low_limit);
      if (impurity_prof)
        Ne_tmp0 = field_larger(Ne0, Low_limit);
      else
        Ne_tmp0 = Zi * N_tmp0;
      
      Limp0 = 0.0;

      if (Limp_carbon) {
        /// Carbon in coronal equilibrium
        /// From I.H.Hutchinson Nucl. Fusion 34 (10) 1337 - 1348 (1994)
	      Limp0 = 2.0e-31 * (pow(Te_tmp0 * Tebar /10., 3.0)) / ((pow(Te_tmp0 * Tebar / 10., 4.5)) + 1.); // in unit Wm^3 
      }

      if (Limp_carbon_adas) { //NOTE(malamast): I changed the order of iteration to optimize it. C++ data are stored is row-major order!!! 
        Field2D Te_tmp_real = Te_tmp0*Tebar;   //eV 
	      for (jx = 0; jx < mesh->LocalNx; jx++) {
	        for (jy = 0; jy < mesh->LocalNy; jy++) {
            BoutReal Te_tmp_loc  = Te_tmp_real(jx,jy);
            BoutReal logT = log(Te_tmp_loc);
            BoutReal log_out =0;
            if (Te_tmp_loc >= 1.0 && Te_tmp_loc <= 500.) {
              log_out = - 7.87837896e+01 * pow(logT, 0)
                        + 1.55326376e+00 * pow(logT, 1)
                        + 1.65898194e+01 * pow(logT, 2)
                        - 3.23804546e+01 * pow(logT, 3)
                        + 3.12784663e+01 * pow(logT, 4)
                        - 1.74826039e+01 * pow(logT, 5)
                        + 5.91393245e+00 * pow(logT, 6)
                        - 1.22974105e+00 * pow(logT, 7)                 
                        + 1.54004499e-01 * pow(logT, 8)
                        - 1.06797106e-02 * pow(logT, 9)                                                   
                        + 3.15657594e-04 * pow(logT, 10);
              Limp0(jx,jy) = exp(log_out);
            } else if (Te_tmp_loc < 1.0) {
              Limp0(jx,jy) = 6.00623928e-35;
            } else {
              Limp0(jx,jy) = 4.53057707e-33;
            }
          }
        }
      }

      if (Limp_nitro) {
        Field2D Te_tmp_real = Te_tmp0*Tebar;   //eV 
        /// Nitrogen based cooling curve used in Lipschultz 2016
	      for (jx = 0; jx < mesh->LocalNx; jx++) { //NOTE(malamast): I changed the order of iteration to optimize it. C++ data are stored is row-major order!!! 
	        for (jy = 0; jy < mesh->LocalNy; jy++) {
            BoutReal Te_tmp_loc  = Te_tmp_real(jx,jy);
              if (Te_tmp_loc > 1.0 && Te_tmp_loc < 80.) {
	              Limp0(jx,jy) = 5.9e-34 * (sqrt(Te_tmp_loc-1.)) * (80. - Te_tmp_loc) / (3.1e-3 * ((Te_tmp_loc - 1.) * (Te_tmp_loc - 1.)) + 1.0);     //Wm^3
              } // else {
              //  Limp0(jx,jy) = 0.0;
              // }
          }
        }
      }

      if (Limp_nitro_adas) {
        Field2D Te_tmp_real = Te_tmp0*Tebar;   //eV 
	      for (jx = 0;jx < mesh->LocalNx; jx++) {
	        for (jy = 0;jy < mesh->LocalNy; jy++) {
            BoutReal log_out =0;
	          BoutReal Te_tmp_loc  = Te_tmp_real(jx,jy,jz);
            BoutReal logT = log(Te_tmp_loc);
            if (Te_tmp_loc >= 2.0 && Te_tmp_loc <= 500.) {
              log_out = - 5.01649969e+01 * pow(logT, 0)
                        - 1.35749724e+02 * pow(logT, 1)
                        + 2.73509608e+02 * pow(logT, 2)
                        - 2.92109992e+02 * pow(logT, 3)
                        + 1.90120639e+02 * pow(logT, 4)
                        - 7.95164871e+01 * pow(logT, 5)
                        + 2.17762218e+01 * pow(logT, 6)
                        - 3.88334992e+00 * pow(logT, 7)                 
                        + 4.34730098e-01 * pow(logT, 8)
                        - 2.77683605e-02 * pow(logT, 9)                                                   
                        + 7.72720422e-04 * pow(logT, 10);
              Limp0(jx,jy) = exp(log_out);
            } else if (Te_tmp_loc < 2.0) {
              Limp0(jx,jy) = 4.34835380e-34;
            } else {
              Limp0(jx,jy) = 8.11096182e-33;
            }
          }
        }
      }

      if (Limp_Neon) {
        /// Neon based cooling curve produced by Matlab polynominal curve
        /// fitting "polyval" (Ryoko 2020 Nov)
        Field2D Te_tmp_real = Te_tmp0*Tebar;   //eV 
	      for (jx = 0; jx < mesh->LocalNx; jx++) {
	        for (jy = 0; jy < mesh->LocalNy; jy++) {
            BoutReal Te_tmp_loc  = Te_tmp_real(jx,jy);
            if (Te_tmp_loc >= 3 and Te_tmp_loc <= 100) {
              Limp0(jx,jy) =  - 2.0385e-40 * pow(Te_tmp_loc, 5)
                              + 5.4824e-38 * pow(Te_tmp_loc, 4)
                              - 5.1190e-36 * pow(Te_tmp_loc, 3)
                              + 1.7347e-34 * SQ(Te_tmp_loc)
                              - 3.4151e-34 * Te_tmp_loc
                              - 3.2798e-34;
            } else if (Te_tmp_loc >=2 and Te_tmp_loc < 3) {
                Limp0(jx,jy) = 7e-35 * (Te_tmp_loc - 2.) + 1e-35;
            } else if (Te_tmp_loc >=1 and Te_tmp_loc < 2) {
                Limp0(jx,jy) = 1e-35 * (Te_tmp_loc - 1.);
            } else {
                Limp0(jx,jy) = 0.0;
            }
          }
        }
      }

      if (Limp_Neon_adas) {
        Field2D Te_tmp_real = Te_tmp0*Tebar;   //eV 
	      for (jx = 0;jx < mesh->LocalNx; jx++) {
	        for (jy = 0;jy < mesh->LocalNy; jy++) {
            BoutReal log_out =0;
            BoutReal Te_tmp_loc  = Te_tmp_real(jx,jy);
            BoutReal logT = log(Te_tmp_loc);
            if (Te_tmp_loc >= 2.0 && Te_tmp_loc <= 1000.) {
              log_out = - 8.21475117e+01 * pow(logT, 0)
                        + 1.28929854e+01 * pow(logT, 1)
                        - 4.74266289e+01 * pow(logT, 2)
                        + 7.45222324e+01 * pow(logT, 3)
                        - 5.75710722e+01 * pow(logT, 4)
                        + 2.57375965e+01 * pow(logT, 5)
                        - 7.12758563e+00 * pow(logT, 6)
                        + 1.24287546e+00 * pow(logT, 7)
                        - 1.32943407e-01 * pow(logT, 8)
                        + 7.97368445e-03 * pow(logT, 9)
                        - 2.05487897e-04 * pow(logT, 10); 
              Limp0(jx,jy) = exp(log_out);
            } else if (Te_tmp_loc < 2.0) {
              Limp0(jx,jy) = 6.35304113e-36;
            } else {
              Limp0(jx,jy) = 1.17894628e-32;
            }
          }
        }
      }

      if (Limp_Argon) {
        /// Argon based cooling curve produced by Matlab polynominal curve
        /// fitting "polyval" (Ryoko 2020 Nov)
        Field2D Te_tmp_real = Te_tmp0*Tebar;   //eV 
        for (jx=0;jx<mesh->LocalNx;jx++) {
          for (jy=0;jy<mesh->LocalNy;jy++) {
            BoutReal Te_tmp_loc  = Te_tmp_real(jx,jy);
            if (Te_tmp_loc >= 1.5 and Te_tmp_loc <= 100) {
              Limp0(jx,jy) = - 4.9692e-48 * pow(Te_tmp_loc, 10)
                             + 2.8025e-45 * pow(Te_tmp_loc, 9)
                             - 6.7148e-43 * pow(Te_tmp_loc, 8)
                             + 8.8636e-41 * pow(Te_tmp_loc, 7)
                             - 6.9642e-39 * pow(Te_tmp_loc, 6)
                             + 3.2559e-37 * pow(Te_tmp_loc, 5)
                             - 8.3410e-36 * pow(Te_tmp_loc, 4)
                             + 8.6011e-35 * pow(Te_tmp_loc, 3)
                             + 1.9958e-34 * pow(Te_tmp_loc, 2)
                             + 4.9864e-34 * Te_tmp_loc
                             - 9.9412e-34;
            } else if (Te_tmp_loc >= 1.0 and Te_tmp_loc < 1.5) {
              Limp0(jx,jy) = 5e-35 * (Te_tmp_loc - 1.0);
            } else {
              Limp0(jx,jy) = 0.0;
            }
          }
        }
      } 

      if (Limp_Argon_adas) {
        Field2D Te_tmp_real = Te_tmp0*Tebar;   //eV 
        for (jx=0;jx<mesh->LocalNx;jx++) {
          for (jy=0;jy<mesh->LocalNy;jy++) {
            BoutReal log_out =0;
            BoutReal Te_tmp_loc  = Te_tmp_real(jx,jy);
            BoutReal logT = log(Te_tmp_loc);
            if (Te_tmp_loc >= 1.5 && Te_tmp_loc <= 1500.) {
              log_out = - 8.45410692e+01 * pow(logT, 0)
                        + 1.57727040e+01 * pow(logT, 1)
                        - 1.54264860e+01 * pow(logT, 2)
                        + 1.49409902e+01 * pow(logT, 3)
                        - 1.04815113e+01 * pow(logT, 4)
                        + 5.00924595e+00 * pow(logT, 5)
                        - 1.60029106e+00 * pow(logT, 6)
                        + 3.29455609e-01 * pow(logT, 7)
                        - 4.14036827e-02 * pow(logT, 8)
                        + 2.87063206e-03 * pow(logT, 9)
                        - 8.38888002e-05 * pow(logT, 10);
              Limp0(jx,jy) = exp(log_out);
            } else if (Te_tmp_loc < 1.5) {
              Limp0(jx,jy)= 1.95353412e-35;
            } else {
              Limp0(jx,jy) = 1.22649600e-32;
            }
          }
        }
      }

      Srad0 = -Limp0 * frac_imp * Ne_tmp0 * Ne_tmp0 * Nbar * density / (1.602e-19 * Tebar / Tbar);  //Normalized //NOTE(malamast): Should we multiply with Nbar * density?
      mesh->communicate(Srad0);
      Srad0.applyBoundary("neumann");
      Wrad0 = Srad0 / Ne_tmp0;
      mesh->communicate(Wrad0);
      Wrad0.applyBoundary("neumann");
      dump.add(Wrad0, "Wrad0", 0);
    }


    /**************** SET EVOLVING VARIABLES ***********************************/

    // Tell BOUT which variables to evolve
    SOLVE_FOR(U, Ni, Ti, Te);
    SAVE_REPEAT(Jpar, P);

    if (emass) {
      output.write("Solving for Psi, Differentiating to get jpar\n");
      SOLVE_FOR(Ajpar);
    } else {
      if (evolve_psi) {
        output.write("Solving for Psi, Differentiating to get jpar\n");
	      SOLVE_FOR(Psi);
      } else {
        output.write("Solving for Apar, Differentiating to get jpar\n");
	      SOLVE_FOR(Apar);
	      Psi = Apar / B0;
	      SAVE_REPEAT(Psi);
      }
    }

    if (parallel_lagrange) {
      // Evolving the distortion of the flux surfaces (Ideal-MHD only!)
      SOLVE_FOR(Xip_x, Xip_z, Xim_x, Xim_z);
    }

    if (parallel_project) {
      // Add Xi to the dump file
      SAVE_REPEAT(Xip_x, Xip_z, Xim_x, Xim_z);
    }

    if (compress0) {
      SOLVE_FOR(Vipar);
      SAVE_REPEAT(Vepar);
      if (!restarting) {
        Vipar = 0.0;
      }
    }

    if (phi_constraint) {
      // Implicit Phi solve using IDA
      solver->constraint(phi, C_phi, "phi");
      
      // Set preconditioner
      setPrecon(&Elm_6f::precon_phi);
    } else {
      // Phi solved in RHS (explicitly)
      SAVE_REPEAT(phi);

      // Set preconditioner
      setPrecon(&Elm_6f::precon);
    }

    // Diamagnetic phi0
    if (diamag && diamag_phi0) {
      if (n0_p0_0p3) {
        // n0 ~ p0^0.3
	      phi0 = -1 / 0.7 * Upara0 * Pi0 / N0;
      } else {
	      // const density
        phi0 = -Upara0 * Pi0 / N0;
      }
      mesh->communicate(phi0);

      if (experiment_Er) { 
	      if (diamag_er) {
	        // get phi0 from grid file
	        Field2D Er_tmp;
	        mesh->get(Er_tmp, "E_r");
          Er0.x = Er0_factor * Er_tmp / (Bbar * Va * Rxy * Bpxy);
          Er0.y = 0.;
          Er0.z = 0.;

          mesh->communicate(Er0);
          Er0.x.applyBoundary();
          Er0.y.applyBoundary();
          Er0.z.applyBoundary();

          Er0_dia = Upara0 * Grad(Pi0) / N0;
          mesh->communicate(Er0_dia);
          Er0.x.applyBoundary();
          Er0.y.applyBoundary();
          Er0.z.applyBoundary();
 
          Er0_net = Er0 - Er0_dia;
	      } else {
          // WARNING: deprecated
          // get phi0 from grid file
          // TODO: U0_net in the end??
          mesh->get(phi0, "Phi_0");
          phi0 /= Bbar * Lbar * Va;
          // BoutReal N0tmp = max(N0,true);
          // phi0_net =  phi0 + Upara0*Pi0/N0tmp;
          // U0_net = Delp2(phi0_net) * N0tmp/B0;
          U0_net = N0 / B0 * (Delp2(phi0) + Upara0 * Delp2(Pi0) / N0); // NOTE(malamast): WARNING: deprecated. We overwrite U0_net below.
          mesh->communicate(U0_net);
          U0_net.applyBoundary("dirichlet");
        }
      } else {
        if (diamag_er) {
	        // phi0 = -Upara0*Pi0/max(N0,true);
          Er0_dia = Upara0 * Grad(Pi0) / N0;
          mesh->communicate(Er0_dia);
          Er0_dia.x.applyBoundary();
          Er0_dia.y.applyBoundary();
          Er0_dia.z.applyBoundary();

	        Er0 = Er0_dia;
	        mesh->communicate(Er0);
          Er0.x.applyBoundary();
          Er0.y.applyBoundary();
          Er0.z.applyBoundary();

          Er0_net = 0;
	      } else {
	        Er0_dia = -Grad(phi0);
	        mesh->communicate(Er0_dia);
          Er0_dia.x.applyBoundary();
          Er0_dia.y.applyBoundary();
          Er0_dia.z.applyBoundary();

          Er0 = Er0_dia;
          mesh->communicate(Er0);
          (Er0.x).applyBoundary();
          (Er0.y).applyBoundary();
          (Er0.z).applyBoundary();

          Er0_net = 0;
	      }
      }
    } else {
      phi0 = 0.;
      Er0 = 0;
      Er0_net = 0;
      Er0_dia = 0;
      Ve0 = 0;
      Ve0_net = 0;
      Ve0_dia = 0;
    } 

    // Er0_net, Ve0, Ve0_dia, Ve0_net, U0_net
    // applyBoundary(), communicate
    mesh->communicate(Er0_net);
    Er0_net.x.applyBoundary();
    Er0_net.y.applyBoundary();
    Er0_net.z.applyBoundary();

    Ve0 = cross(Er0, B0vec) / (B0 * B0);
    Ve0.setBoundary("Vipar");
    mesh->communicate(Ve0);
    Ve0.x.applyBoundary();
    Ve0.y.applyBoundary();
    Ve0.z.applyBoundary();

    Ve0_dia = cross(Er0_dia, B0vec) / (B0 * B0);
    Ve0_dia.setBoundary("Vipar");
    mesh->communicate(Ve0_dia);
    Ve0_dia.x.applyBoundary();
    Ve0_dia.y.applyBoundary();
    Ve0_dia.z.applyBoundary();

    Ve0_net = Ve0 - Ve0_dia;
    Ve0_net.setBoundary("Vipar");
    mesh->communicate(Ve0_net);
    Ve0_net.x.applyBoundary();
    Ve0_net.y.applyBoundary();
    Ve0_net.z.applyBoundary();

    U0_net = B0vec * Curl(N0 * Ve0_net) / B0;
    mesh->communicate(U0_net);
    U0_net.applyBoundary("dirichlet");

    // Add some equilibrium quantities and normalisations
    // everything needed to recover physical units
    SAVE_ONCE(J0, P0);
    SAVE_ONCE(density, Lbar, Bbar, Tbar);
    SAVE_ONCE(Tibar, Tebar, Nbar);
    SAVE_ONCE(Va, B0);
    SAVE_ONCE(Ti0, Te0, N0, Ne0);
    if (impurity_prof) {
      //SAVE_ONCE(N_imp0, T_imp0, P_imp0);
      SAVE_ONCE(N_imp0);
    }
    if (diamag) {
      SAVE_ONCE(phi0, Er0_dia, Ve0_dia);
      if (experiment_Er) {
        SAVE_ONCE(Er0, Ve0, U0_net);
      }
    }

    if (include_vpar0) {
      SAVE_ONCE2(Vipar0, Vepar0); 
    }

    if (include_U0) {
      SAVE_ONCE(U0); 
    }

    // Create a solver for the Laplacian
    phiSolver = Laplacian::create(&globalOptions["phiSolver"]);

    aparSolver = Laplacian::create(&globalOptions["aparSolver"], loc);

    /////////////// CHECK VACUUM ////////////////////////////////////
    // In vacuum region, initial vorticity should equal zero

    ubyn.setBoundary("U");

    density_tmp = N0 + A_imp / AA * N_imp0;

    if (!restarting) {
      // Only if not restarting: Check initial perturbation

      // Set U to zero where P0 < vacuum_pressure
      U = where(P0 - vacuum_pressure, U, 0.0);

      // Field2D lap_temp = 0.0;
      // TODO: diamag in U?
      if (impurity_prof) {
        Field2D logn0 = laplace_alpha * density_tmp;
        ubyn = U * B0 / density_tmp;
        // Phi should be consistent with U
        if (laplace_alpha <= 0.0) {
          phi = phiSolver->solve(ubyn);
        } else {
          phiSolver->setCoefC(logn0);
          phi = phiSolver->solve(ubyn);
        }
      } else {
        Field2D logn0 = laplace_alpha * N0;
	      ubyn = U * B0 / N0;
        // Phi should be consistent with U
        if (laplace_alpha <= 0.0) {
          phi = phiSolver->solve(ubyn);
        } else {
          phiSolver->setCoefC(logn0);
          phi = phiSolver->solve(ubyn);
        }
      }
    }

    if ((gamma_i_BC > 0.0) && (gamma_e_BC > 0.0)) {
      output.write("Sheath Boundary conditions applied.\n");
      dump.add(c_se, "c_se", 1);
      dump.add(q_si, "q_si", 1);
      dump.add(q_se, "q_se", 1);

      const_cse = sqrt(KB * Tebar * eV_K / Mi);

      c_se0 = sqrt(abs(Tau_ie * Ti0 + Te0));
      c_se0 *= const_cse;
      vth_e0 = 4.19e5 * sqrt(Te0 * Tebar);

      c_se0 /= Va; // NOTE(malamast): Just a test. 
      dump.add(c_se0, "c_se0", 0);
      c_se0 *=Va;

      // output << max(c_se0, true) << "\t" << max(vth_e0, true) << endl;
      dump.add(Jpar_sh, "Jpar_sh", 1);
      Jpar_sh0 = Ne0 * Nbar * density * ee;
      Jpar_sh0 *= c_se0 - vth_e0 / (2.0 * sqrt(PI)) * exp(- ee * (phi0 * Va * Lbar * Bbar) / (KB * Te0 * Tebar * eV_K));
      dump.add(phi_sh, "phi_sh", 1);
      phi_sh.setBoundary("phi");
      // phi_sh0 = -Te0 * Tebar;
      // phi_sh0 *= log(2. * sqrt(PI) * (c_se0-J0 * B0 * Bbar / (MU0 * Lbar) / (Ne0 * Nbar * density * ee)) / vth_e0);
    }

    if (BScurrent) {
      // TODO: check
      Field3D L31, L32, L34;
      Field3D f31, f32ee, f32ei, f34, ft;
      Field3D BSal0, BSal;
      Jpar_BS0.setBoundary("J");

      nu_estar = 100. * nu_e * q95 * Lbar / (vth_e) / pow(Aratio, 1.5);
      nu_istar = 100. * nu_i * q95 * Lbar / (vth_i) / pow(Aratio, 1.5);
      // nu_estar = 0.012 * N0*Nbar*density/1.e20*Zi*Zi*q95*Lbar/(Te0*Tebar/1000. * pow(Aratio, 1.5));
      // nu_istar = 0.012 * N0*Nbar*density/1.e20*Zi*q95*Lbar/(Ti0*Tibar/1000. * pow(Aratio, 1.5));
      output.write("Bootstrap current is included:\n");
      output.write("Normalized electron collisionality: nu_e* = {:e}\n", max(nu_estar));
      output.write("Normalized ion collisionality: nu_i* = {:e}\n", max(nu_istar));
      ft = BS_ft(100);
      output.write("modified collisional trapped particle fraction: ft = {:e}\n", max(ft));
      f31 = ft / (1. + (1. - 0.1 * ft) * sqrt(nu_estar) + 0.5 * (1. - ft) * nu_estar / Zi);
      f32ee = ft / (1. + 0.26 * (1. - ft) * sqrt(nu_estar) + 0.18 * (1. - 0.37 * ft) * nu_estar / sqrt(Zi));
      f32ei = ft / (1. + (1. + 0.6 * ft) * sqrt(nu_estar) + 0.85 * (1. - 0.37 * ft) * nu_estar * (1. + Zi));
      f34 = ft / (1. + (1. - 0.1 * ft) * sqrt(nu_estar) + 0.5 * (1. - 0.5 * ft) * nu_estar / Zi);

      L31 = F31(f31);
      L32 = F32ee(f32ee) + F32ei(f32ei);
      L34 = F31(f34);
      BSal0 = -(1.17 * (1. - ft)) / (1. - 0.22 * ft - 0.19 * ft * ft);
      BSal = (BSal0 + 0.25 * (1 - ft * ft) * sqrt(nu_istar)) / (1. + 0.5 * sqrt(nu_istar)) + 0.31 * nu_istar * nu_istar * ft * ft * ft * ft * ft * ft;
      BSal *= 1. / (1. + 0.15 * nu_istar * nu_istar * ft * ft * ft * ft * ft * ft);

      Jpar_BS0 = L31 * DDX(P0) / Pe0 + L32 * DDX(Te0) / Te0 + L34 * DDX(Ti0) / (Zi * Te0) * BSal;
      Jpar_BS0 *= Field3D(-Rxy * Btxy * Pe0 / (B0 * B0) * (MU0 * KB * Nbar * density * Tebar * eV_K) / (Bbar * Bbar));

      mesh->communicate(Jpar_BS0);
      Jpar_BS0.applyBoundary();
      dump.add(Jpar_BS0, "jpar_BS0", 0);
      dump.add(nu_estar, "nu_estar", 0);
      dump.add(nu_istar, "nu_istar", 0);
    }

    if (radial_diffusion) {
      diffusion_coef_Hmode0 /= Lbar * Lbar / Tbar;
      diffusion_coef_Hmode1 /= Lbar * Lbar / Tbar;

      ddx_ni.setBoundary("Ni");
      diff_radial.setBoundary("Ni");
      ddx_n0 = DDX(N0);
      ddx_n0 = -field_larger(-ddx_n0, Low_limit);
      diff_radial = Field3D(diffusion_coef_Hmode0);

      dump.add(diff_radial, "diff_radial", 1);
    }

    /************** SETUP COMMUNICATIONS ***************************************/

    comms.add(U);
    comms.add(Ni);
    comms.add(Ti);
    comms.add(Te);
    if (!emass) {
      if (evolve_psi) {
        comms.add(Psi);
      } else {
        comms.add(Apar);
        Psi.setBoundary("Psi");
      }
    } else {
      comms.add(Ajpar);
    }

    if (compress0) {
      comms.add(Vipar);
      Vepar.setBoundary("Vipar");
    }

    if (neutral) {
      comms.add(Nn);
      comms.add(Vn);
    }

    if (hyperdiff_par_u4 > 0.0 || hyperdiff_perp_u4 > 0.0) {
      tmpU2.setBoundary("U");
    }

    if (hyperdiff_par_apar4 > 0.0 || hyperdiff_perp_apar4 > 0.0) {
      tmpA2.setBoundary("J");
    }

    if (hyperdiff_par_n4 > 0.0 || hyperdiff_perp_n4 > 0.0) {
      tmpN2.setBoundary("Ni");
    }

    if (hyperdiff_par_ti4 > 0.0 || hyperdiff_perp_ti4 > 0.0) {
      tmpTi2.setBoundary("Ti");
    }

    if (hyperdiff_par_te4 > 0.0 || hyperdiff_perp_te4 > 0.0) {
      tmpTe2.setBoundary("Te");
    }

    if (hyperdiff_par_v4 > 0.0 || hyperdiff_perp_v4 > 0.0) {
      tmpVp2.setBoundary("Vipar");
    }

    phi.setBoundary("phi"); // Set boundary conditions

    P.setBoundary("P");
    Jpar.setBoundary("J");
    Jpar2.setBoundary("J");

    F2D_tmp = 0.;
    
    return 0;
  }

  /// For printing out some diagnostics first time around
  bool first_run = true;

  /*
  --------RHS Function
  */

  int rhs(BoutReal UNUSED(t)) override {

    Coordinates* coord = mesh->getCoordinates();

    // Perform communications
    mesh->communicate(comms);

    if (pos_filter) {
      Ti = lowPass_pos2(Ti, Ti);
      Te = lowPass_pos2(Te, Te);
      Ni = lowPass_pos2(Ni, Ni);
    }

    // Inversion
    Pi = Ni * Ti0 + N0 * Ti;
    if (nonlinear) {
      Pi += Ni * Ti;
    }
    mesh->communicate(Pi);

    Pe = Zi * Ni * Te0 + Ne0 * Te;
    if (nonlinear) {
      Pe += Zi * Ni * Te;
    }
    mesh->communicate(Pe);

    P = Tau_ie * Pi + Pe;
    mesh->communicate(P);

    if (nonlinear && pos_filter) {
      Pi = lowPass_pos2(Pi, Pi);
      Pe = lowPass_pos2(Pe, Pe);
      P = lowPass_pos2(P, P);
    }

    if (nonlinear && pos_filter2) {
      Pi = lowPass_pos(Pi, position_tmpi);
      Pe = lowPass_pos(Pe, position_tmpe);
      P = lowPass_pos(P, position_tmp);
    }

    if (nonlinear && (pos_filter_zf || (pos_sink_zf > 0.))) {
      Pi = sink_zonal_core(Pi, position_tmpi);
      Pe = sink_zonal_core(Pe, position_tmpe);
      P = sink_zonal_core(P, position_tmp);
    }

    if (low_pass_z > 0.) {
      Pi = lowPass(Pi, low_pass_z, zonal_bkgd);
      Pe = lowPass(Pe, low_pass_z, zonal_bkgd);
      P = lowPass(P, low_pass_z, zonal_bkgd);
    }
    
    if ((gamma_i_BC > 0.0) && (gamma_e_BC > 0.0)) {
      SBC_Gradpar(U, 0.0, PF_limit, PF_limit_range);
    }

    // Field2D lap_temp=0.0;
    // Field2D logn0 = laplace_alpha * N0;
    // ubyn = U * B0 / max(N0, true);
    if (impurity_prof) {
      ubyn = U * B0 / density_tmp;
      if (diamag) {
        ubyn -= Upara0 / density_tmp * Delp2(Pi);
        mesh->communicate(ubyn);
        ubyn.applyBoundary();
      }

      // Invert laplacian for phi
      if (laplace_alpha <= 0.0) {
        phi = phiSolver->solve(ubyn);
      } else {
        phiSolver->setCoefC(density_tmp);
        phi = phiSolver->solve(ubyn);
      }
    } else { 
      ubyn = U * B0 / N0;
      if (diamag) {
        ubyn -= Upara0 / N0 * Delp2(Pi);
        mesh->communicate(ubyn);
        ubyn.applyBoundary();
      }
      
      // Invert laplacian for phi
      if (laplace_alpha <= 0.0) {
        phi = phiSolver->solve(ubyn);
      } else {
        phiSolver->setCoefC(N0);
        phi = phiSolver->solve(ubyn); // NOTE(malamast): a) Should this be N0 + Ni  instead of just N0
                                      //                 b) Should we subtract the term 1/N0 Grad_perp(phi0) * Grad_perp(Ni) like we did above with  Upara0 / N0 * Delp2(Pi)?
      }
    }

    if (mask_phi_x) {
      phi *= mask_px1d;
    }

    mesh->communicate(phi);

    if (low_pass_z > 0.) {
      phi = lowPass(phi, low_pass_z, false); // now phi has no DC component
    }
    
    if (zonal_flow) {
      VortDC = DC(ubyn);  // n=0 component	
      if (laplace_alpha <= 0.0) {
        laplacexy->setCoefs(1.0, 0.0);
      } else {
	      if (impurity_prof)
          laplacexy->setCoefs(density_tmp, 0.0);
        else
	        laplacexy->setCoefs(N0, 0.0);
      }
      // Solve axisymmetric (n=0) part
      phiDC = laplacexy->solve(VortDC, phiDC);
      phi += phiDC;
    }

    // Apply a boundary condition on phi for target plates
    // phi.applyBoundary();

    if (emass) {
      Field2D acoeff = -delta_e_inv * N0 * N0;
      if (compress0) {
        aparSolver->setCoefC(acoeff);
        Psi = aparSolver->solve(acoeff * Ajpar - gyroAlv * Vipar);
      } else {
        aparSolver->setCoefC(acoeff);
        Psi = aparSolver->solve(acoeff * Ajpar);
      }
      mesh->communicate(Psi);
    }

    N_tmp1 = Low_limit;
    if (nonlinear) {
      N_tmp = field_larger(N0 + Ni, N_tmp1);
      if (impurity_prof)
        Ne_tmp = field_larger(Ne0 + Zi * Ni, N_tmp1);
      else
        Ne_tmp = Zi * N_tmp;
    }

    // BoutReal Te_tmp1, Ti_tmp1;
    Te_tmp1 = Low_limit;
    Ti_tmp1 = Low_limit;
    if (nonlinear) {
      Ti_tmp = field_larger(Ti0 + Ti, Ti_tmp1);
      Te_tmp = field_larger(Te0 + Te, Te_tmp1);
    }

    if (!nonlinear && (parallel_viscous && compress0)) {
      pi_ci = -eta_i0 * 2. * (pow(B0, -0.5)) * Grad_par((pow(B0, 0.5)) * Vipar);
      pi_ci -= eta_i0 * b0xcv * (Er0_net + Grad(phi)) / B0;
      mesh->communicate(pi_ci);
      pi_ci.applyBoundary();
    }

    // vac_mask transitions from 0 in core to 1 in vacuum
    if (nonlinear) {
      vac_mask = (1.0 - tanh(((P0 + P) - vacuum_pressure) / vacuum_trans)) / 2.0;
      // Update resistivity
      if (spitzer_resist) {
        // Use Spitzer formula
        eta = FZ * 1.03e-4 * Zi * LnLambda * (pow(Te_tmp * Tebar, -1.5)); // eta in Ohm-m. ln(Lambda) = 20
	      eta /= MU0 * Va * Lbar;
      } else {
        eta = core_resist + (vac_resist - core_resist) * vac_mask;
      }

      if (impurity_prof)
        nu_e = 2.91e-6 * LnLambda * (N_tmp * Nbar * density / 1.e6) * (pow(Te_tmp * Tebar, -1.5)); // nu_e in 1/S.
      else 
	      nu_e = 2.91e-6 * LnLambda * (N_tmp * Nbar * density / 1.e6) * (pow(Te_tmp * Tebar, -1.5)); // nu_e in 1/S.
      
      if (diffusion_par > 0.0 || diffusion_perp > 0.0 || parallel_viscous || neoclassic_i || neoclassic_e) {
        // Use Spitzer thermal conductivities

        nu_i = 4.80e-8 * (Zi * Zi * Zi * Zi / sqrt(AA)) * LnLambda * (N_tmp * Nbar * density / 1.e6) * pow(Ti_tmp * Tibar, -1.5);  // nu_i in 1/S.
        vth_i = 9.79e3 * sqrt(Ti_tmp * Tibar / AA);   // vth_i in m/S.
        vth_e = 4.19e5 * sqrt(Te_tmp * Tebar);        // vth_e in m/S.
      }

      if (parallel_viscous && compress0) {
        eta_i0 = 0.96 * (Pi0 + Pi) * Tau_ie * nu_i * Tbar;
        pi_ci = -eta_i0 * 2. / sqrt(B0) * Grad_parP((sqrt(B0) * Vipar));
        pi_ci -= eta_i0 * b0xcv * (Er0_net + Grad(phi)) / B0;
        mesh->communicate(pi_ci);
        pi_ci.applyBoundary();
      }

      if (diffusion_par > 0.0) {
        kappa_par_i = 3.9 * vth_i * vth_i / nu_i; // * 1.e4;
        kappa_par_e = 3.2 * vth_e * vth_e / nu_e; // * 1.e4;

        kappa_par_i_fl = q_alpha * vth_i * (q95 * Lbar); // * 1.e2;
        kappa_par_e_fl = q_alpha * vth_e * (q95 * Lbar); // * 1.e2;

	      if (fluxlimit) {
          kappa_par_i *= kappa_par_i_fl / (kappa_par_i + kappa_par_i_fl);
	        kappa_par_e *= kappa_par_e_fl / (kappa_par_e + kappa_par_e_fl);
	      }
	      kappa_par_i *= Tipara1 * N_tmp;
        kappa_par_e *= Tepara1 * Ne_tmp;
      }

      if (diffusion_perp > 0.0) {
        kappa_perp_i = 2.0 * vth_i * vth_i * nu_i / (omega_ci * omega_ci); // * 1.e4;
        kappa_perp_e = 4.7 * vth_e * vth_e * nu_e / (omega_ce * omega_ce); // * 1.e4;
        kappa_perp_i_fl = q_alpha * vth_i * q95 * Lbar; // * 1.e4;
        kappa_perp_e_fl = q_alpha * vth_e * q95 * Lbar; // * 1.e4;

	      if (fluxlimit) {
          kappa_perp_i *= kappa_perp_i_fl / (kappa_perp_i + kappa_perp_i_fl);
          kappa_perp_e *= kappa_perp_e_fl / (kappa_perp_e + kappa_perp_e_fl);
	      }
	      kappa_perp_i *= Tipara1 * N_tmp;
	      kappa_perp_e *= Tepara1 * Ne_tmp;
      }

      if (neoclassic_i) {
        rho_i = 1.02e-4 * sqrt(AA * Ti_tmp * Tibar) / B0 / Bbar / Zi;
        // Dri_neo = (1.+1.6 * q95) * (1. + Tau_ie) * nu_i * rho_i * rho_i;
        xii_neo = (q95 * q95) / epsilon / sqrt(epsilon) * nu_i * (pow(rho_i, 2.));
        // Dri_neo *= 3./2.* Tipara1;
        xii_neo *= Tipara1;
      	// Dri_neo.applyBoundary();
        // xii_neo.applyBoundary();
      }
      if (neoclassic_e) {
        rho_e = 2.38e-6 * sqrt(Te_tmp * Tebar) / B0 / Bbar;
        xie_neo = (q95 * q95) / epsilon / sqrt(epsilon) * nu_e * rho_e * rho_e;
        xie_neo *= Tepara1;
      }
    }

    // update Landau parallel heat flux
    /*  if (diffusion_par && Landau) {
          SBC_value_i = gamma_i_BC * Pi * c_set / (1.6 * N0  * vth_i);
          SBC_value_e = gamma_e_BC * Pe * c_set / (1.6 * Ne0 * vth_e);
          if (Landau_coll) {
            q_par_i = -Landau_coeff * (2. / 3.) * N0 * 2. * sqrt(2. / PI) * vth_i / Va *eiSign_kpar_wcoll(Ti, kappa_i, 2, nLorentzian, SBC_value_i);
            q_par_e = -Landau_coeff * (2. / 3.) * Ne0 * 2. * sqrt(2. / PI) * vth_e / Va * iSign_kpar_wcoll(Te, kappa_e, 2, nLorentzian, SBC_value_e);
          } else {
            q_par_i = -Landau_coeff * (2. / 3.) * N0 * 2. * sqrt(2. / PI) * vth_i / Va * iSign_kpar(Ti, kappa_0, 2, nLorentzian, SBC_value_i);
            q_par_e = -Landau_coeff * (2. / 3.) * Ne0 * 2. * sqrt(2. / PI) * vth_e / Va * iSign_kpar(Te, kappa_0, 2, nLorentzian, SBC_value_e);
          }
          mesh->communicate(q_par_i);
          mesh->communicate(q_par_e);
          q_par_i.applyBoundary();
          q_par_e.applyBoundary();
        } 
    */

    if (radial_diffusion && nonlinear) {
      ddx_ni = DDX(N_tmp);
      mesh->communicate(ddx_ni);
      ddx_ni.applyBoundary();

      for (jx = 0; jx < mesh->LocalNx; jx++) {
        for (jy = 0; jy < mesh->LocalNy; jy++) {
	        for (jz = 0; jz < mesh->LocalNz; jz++) {
	          if (ddx_ni(jx,jy,jz) > -Low_limit && ddx_ni(jx,jy,jz) < 0.) {
	            ddx_ni(jx,jy,jz) = -Low_limit;
	          } else if (ddx_ni(jx,jy,jz) < Low_limit && ddx_ni(jx,jy,jz) >= 0.) {
	            ddx_ni(jx,jy,jz) = Low_limit;
	          }
	        }
	      }
      }

      diff_radial = diffusion_coef_Hmode0 * ddx_n0 / ddx_ni;

      for (jx = 0; jx < mesh->LocalNx; jx++) {
        for (jy = 0; jy < mesh->LocalNy; jy++) {
	        for (jz = 0; jz < mesh->LocalNz; jz++) {
	          if (diff_radial(jx,jy,jz) > diffusion_coef_Hmode1) {
	            diff_radial(jx,jy,jz) = diffusion_coef_Hmode1;
	          }
	        }
	      }
      }

      diff_radial = nl_filter(diff_radial, 1);
    }

    if (evolve_psi) {
      Jpar = -B0 * Delp2(Psi);
    } else {
      Jpar = -Delp2(Apar);
    }
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
    if (smooth_j_x && !((gamma_i_BC > 0.0) && (gamma_e_BC > 0.0))) {
      Jpar = smooth_x(Jpar);
      /*Jpar = smooth_x(Jpar);
      Jpar = smooth_x(Jpar);
      Jpar = smooth_x(Jpar);
      Jpar = smooth_x(Jpar);
      Jpar = smooth_y(Jpar);
      Jpar = smooth_y(Jpar);
      Jpar = smooth_y(Jpar);
      Jpar = smooth_y(Jpar);
      Jpar = smooth_y(Jpar);
      */
      mesh->communicate(Jpar);
      Jpar.applyBoundary();
    }

    if (mask_j_x) {
      Jpar *= mask_jx1d;
    }

    if ((gamma_i_BC > 0.0) && (gamma_e_BC > 0.0)) {
      if (nonlinear) {
        c_set = sqrt(abs(Tau_ie * Ti_tmp + Te_tmp));
        c_se = c_set - c_se0 / const_cse;
        c_se *= const_cse / Va; // normalized
        c_set *= const_cse;     // not normalized, with unit
        vth_et = 4.19e5 * sqrt(Te_tmp * Tebar);
      } else {
        c_set = sqrt(abs(Tau_ie * Ti0 + Te0));
        c_se = 0.;
        c_se *= const_cse / Va; // normalized
        c_set *= const_cse;     // not normalized, with unit
        vth_et = 4.19e5 * sqrt(Te0 * Tebar);
      }

      if (evolve_psi) {
        phi_sh = -eta * Jpar / B0;

        if (eHall) {
          phi_sh += Psipara1 * Grad_parP(Pe) / B0 / Ne0;
          phi_sh -= Psipara1 * bracket(Psi, Pe0, bm_mag) / Ne0;
	      }
        if (thermal_force) {
          phi_sh += 0.71 * Psipara1 * Grad_parP(Te) / B0;
          phi_sh -= 0.71 * Psipara1 * bracket(Psi, Te0, bm_mag);
        }
      } else {
        phi_sh = -eta * Jpar;

        if (eHall) {
          phi_sh += Psipara1 * Grad_parP(Pe) / Ne0;
          phi_sh -= Psipara1 * bracket(Apar, Pe0, bm_mag) / Ne0;  // NOTE(malamast): Why do we need that?
        }
        if (thermal_force) {
          phi_sh += 0.71 * Psipara1 * Grad_parP(Te);
          phi_sh -= 0.71 * Psipara1 * bracket(Apar, Te0, bm_mag);
        }
      }
      mesh->communicate(phi_sh);
      phi_sh.applyBoundary();
      SBC_Gradpar(phi, phi_sh, PF_limit, PF_limit_range);

      if (nonlinear) {
        Jpar_sh = Ne_tmp * Nbar * density * ee;
        Jpar_sh *= c_set - vth_et / (2.0 * sqrt(PI)) * exp(- ee * ((phi + phi0) * Va * Lbar * Bbar) / (KB * Te_tmp * Tebar * eV_K));
        Jpar_sh -= Jpar_sh0;
        Jpar_sh *= MU0 * Lbar / Bbar;
        SBC_Dirichlet(Jpar, Jpar_sh, PF_limit, PF_limit_range);

        /* phi_sh = -Te_tmp*Tebar;
              if (impurity_prof)
                phi_sh *= log( 2.0*sqrt(PI)*(c_set-(J0+Jpar)*B0*Bbar/(MU0*Lbar)/(Ne_tmp*Nbar*density*ee))/vth_et );
              else
                phi_sh *= log( 2.0*sqrt(PI)*(c_set-(J0+Jpar)*B0*Bbar/(MU0*Lbar)/(N_tmp*Nbar*density*ee*Zi))/vth_et );
              phi_sh -= phi_sh0;
              phi_sh /= Bbar*Va*Lbar;*/
      } else {
        Jpar_sh = Ne0 * Nbar * density * ee;
        Jpar_sh *= vth_et / (2.0 * sqrt(PI)) * (ee * ((phi)*Va * Lbar * Bbar) / (KB * Te0 * Tebar * eV_K));
        Jpar_sh *= MU0 * Lbar / (Bbar);
        SBC_Dirichlet(Jpar, Jpar_sh, PF_limit, PF_limit_range);
      }

      if (diffusion_par > 0.0) {
        if (full_sbc) {
          q_se = -2. / 3. * gamma_e_BC * (Pe + Pe0) * c_set / Va / kappa_par_e; // * (Nbar*density*KB);
          q_si = -2. / 3. * gamma_i_BC * (Pi + Pi0) * c_set / Va / kappa_par_i; // * (Nbar*density*KB);
	      } else {
          q_se = -2. / 3. * gamma_e_BC * Pe * c_se / kappa_par_e; // * (Nbar*density*KB);
          q_si = -2. / 3. * gamma_i_BC * Pi * c_se / kappa_par_i; // * (Nbar*density*KB);
	      }
      } else {
        q_se = 0.;
        q_si = 0.;
      }

      if (compress0) {
        if (full_sbc)
          SBC_Dirichlet(Vipar, c_set/Va, PF_limit, PF_limit_range);
        else
          SBC_Dirichlet(Vipar, c_se, PF_limit, PF_limit_range);
      }

      // SBC_Gradpar(U, 0.0, PF_limit, PF_limit_range);
      SBC_Gradpar(Ni, 0.0, PF_limit, PF_limit_range);
      // if (!Landau) {
        SBC_Gradpar(Ti, q_si, PF_limit, PF_limit_range);
        SBC_Gradpar(Te, q_se, PF_limit, PF_limit_range);
      // } else {
      //   SBC_Gradpar(Ti, 0., PF_limit, PF_limit_range);
      //   SBC_Gradpar(Te, 0., PF_limit, PF_limit_range);
      // }     
    }

  //   // update Landau parallel heat flux
  //   if (diffusion_par && Landau) {
  //     if (full_sbc) {
  //       SBC_value_i = gamma_i_BC * (Pi + Pi0) * c_set / (1.6 * N0  * vth_i);
  //       SBC_value_e = gamma_e_BC * (Pe + Pe0) * c_set / (1.6 * Ne0 * vth_e);
  //     } else {
  //       SBC_value_i = gamma_i_BC * Pi * c_set / (1.6 * N0  * vth_i);
  //       SBC_value_e = gamma_e_BC * Pe * c_set / (1.6 * Ne0 * vth_e);
  //     }
  //     if (Landau_coll) {
  //       q_par_i = -Landau_coeff * (2. / 3.) * N0 * 2. * sqrt(2. / PI) * vth_i / Va * iSign_kpar_wcoll(Ti, kappa_i, 2, nLorentzian, SBC_value_i);
  //       q_par_e = -Landau_coeff * (2. / 3.) * Ne0 * 2. * sqrt(2. / PI) * vth_e / Va * iSign_kpar_wcoll(Te, kappa_e, 2, nLorentzian, SBC_value_e);
  //     } else {
  //       q_par_i = -Landau_coeff * (2. / 3.) * N0 * 2. * sqrt(2. / PI) * vth_i / Va * iSign_kpar(Ti, kappa_0, 2, nLorentzian, SBC_value_i);
  //       q_par_e = -Landau_coeff * (2. / 3.) * Ne0 * 2. * sqrt(2. / PI) * vth_e / Va * iSign_kpar(Te, kappa_0, 2, nLorentzian, SBC_value_e);
  //     }
  //     mesh->communicate(q_par_i);
  //     mesh->communicate(q_par_e);
  //     q_par_i.applyBoundary();
  //     q_par_e.applyBoundary();

  // /*    SBC_value_i = gamma_i_BC * (Pi + Pi0) * c_set / Va;
  //     SBC_value_e = gamma_e_BC * (Pe + Pe0) * c_set / Va;
  //     SBC_Dirichlet(q_par_i, SBC_value_i, PF_limit, PF_limit_range);
  //     SBC_Dirichlet(q_par_e, SBC_value_e, PF_limit, PF_limit_range);
  // */
  //     SBC_Gradpar(q_par_i, 0., PF_limit, PF_limit_range);
  //     SBC_Gradpar(q_par_e, 0., PF_limit, PF_limit_range);
  //   }

    if (neutral) {
      N_tmp = field_larger(N0 + Ni, Low_limit);
      Ti_tmp = field_larger(Ti0 + Ti, Low_limit);
      Te_tmp = field_larger(Te0 + Te, Low_limit);
      if (impurity_prof)
        Ne_tmp = field_larger(Ne0 + Zi * Ni, Low_limit);
      else
        Ne_tmp = Zi * N_tmp;
      
      // calculate from src code in neutral.cxx
      if (read_collision_rate) {
        nu_iz = iz_rate(Te_tmp * Tebar, 1) * Nbar * density * Lbar / Va;
        nu_cx = cx_rate(Ti_tmp * Tibar, 1) * Nbar * density * Lbar / Va;
        nu_rc = rc_rate(Te_tmp * Tebar, N_tmp * Nbar * density, 1) * Nbar * density * Lbar / Va;
        sigma_cx = cx_sect(Ti_tmp * Tibar, 1) * Nbar * density * Lbar;
      } else {
        nu_iz = 3.e4 * Tbar * Nbar * Te_tmp * Te_tmp * Tebar * Tebar / (3. + 0.01 * Te_tmp * Te_tmp * Tebar * Tebar);
        nu_cx = 1.e6 * Tbar * Nbar * (1.7 + 1.9 * ((pow(1.5 * Ti_tmp * Tibar, 0.333)) - 2.4662) / ((pow(150.0 * Ti_tmp * Tibar, 0.333)) - 2.4662));
        Field3D lambda_rec=1.5789e5 / (Te_tmp * Tebar * 1.1604e4);  //Te trasfer unit from eV to Kelvin
        nu_rc = Tbar * Nbar * 5.197 * Zi * sqrt(lambda_rec) * (0.4288 + 0.5 * log(lambda_rec) + 0.469 / pow(lambda_rec, 0.333));
        sigma_cx = cx_sect(Ti_tmp * Tibar, 1) * Nbar * density * Lbar;
        mesh->communicate(nu_iz);
        mesh->communicate(nu_cx);
        mesh->communicate(nu_rc);
	      mesh->communicate(sigma_cx);
      }
      vth_i = 9.79e3 * sqrt((Ti_tmp) * Tibar / AA); // vth_i in m/S.
      // diffusion and viscosity coefficients
      if (constent_Dn) {
        Dn = Diff_n /(Va * Va * Tbar);
      } else {
        Dn = vth_i * vth_i / (Va * Va * (N_tmp * nu_cx))/fac_Dn;
      }
      Dn1 = vth_i * vth_i / (Va * Va * (N_tmp * nu_cx));
      Dn_fl = vth_i /Va * Lnn_min / Lbar;
      if (fl_Dn) {
        Dn *=Dn_fl / (Dn + Dn_fl);
      }
      mesh->communicate(Dn);
      Dn.applyBoundary("neumann");

      etan = 1.0 / (N_tmp + Nn) * vth_i / (Va * sigma_cx);
      etan_perp = etan;
      if (fac_etan > 0.0) {
        etan_perp /= fac_etan;
      }
      mesh->communicate(etan);
      etan.applyBoundary("neumann");

      c_set = sqrt(abs(Tau_ie * Ti_tmp + Te_tmp));
      c_se0 = sqrt(abs(Tau_ie * Ti0 + Te0));
      const_cse = sqrt(KB * Tebar * eV_K / Mi);
      c_se0 *= const_cse;
      c_se = c_set - c_se0 / const_cse;
      c_se *= const_cse / Va; // normalized
      c_set *= const_cse;     // not normalized, with unit
      mesh->communicate(c_set);
      c_set.applyBoundary("neumann");
      mesh->communicate(c_se);
      c_se.applyBoundary("neumann");

      if (full_sbc_Vn) {
        Gamma_nn = Rcyc_Nn * N_tmp * c_set / Va / Dn;
      } else {
        Gamma_nn = Rcyc_Nn * N_tmp * c_se / Dn;
      }
      mesh->communicate(Gamma_nn);
      Gamma_nn.applyBoundary("neumann");
      if (Nn_recyc_BC) {
      SBC_Gradpar(Nn, Gamma_nn, PF_limit, PF_limit_range);
      } else {
        SBC_yup_eq(Nn, Rcyc_Nn * N_tmp, PF_limit, PF_limit_range);
        SBC_ydown_eq(Nn, Rcyc_Nn * N_tmp, PF_limit, PF_limit_range);
      }
      if (Vn_recyc_BC) {
        if (full_sbc_Vn) {
          SBC_Dirichlet(Vn, -Rcyc_Vn * c_set / Va, PF_limit, PF_limit_range);
        } else {
        SBC_Dirichlet(Vn, -Rcyc_Vn * c_se, PF_limit, PF_limit_range);
        }
      } else {
        SBC_Gradpar(Vn, 0.0, PF_limit, PF_limit_range);
      }
    }
    
    if (fix_fraction_imp) {
      N_tmp = field_larger(N0 + Ni, Low_limit);
      Te_tmp = field_larger(Te0 + Te, Low_limit);
      if (impurity_prof)
        Ne_tmp = field_larger(Ne0 + Zi * Ni, Low_limit);
      else
        Ne_tmp = Zi * N_tmp;
      Limp = 0.0;
      if (Limp_carbon) {
        /// Carbon in coronal equilibrium
	      // From I.H.Hutchinson Nucl. Fusion 34 (10) 1337 - 1348 (1994)
	      Limp = 2.0e-31 * (pow(Te_tmp * Tebar /10., 3.0)) / ((pow(Te_tmp * Tebar / 10., 4.5)) + 1.);
      }

      //NOTE(malamast): I changed the order of iteration to optimize it. C++ data are stored is row-major order!!! 
      if (Limp_carbon_adas) {
	      for (jx = 0; jx < mesh->LocalNx; jx++) {
	        for (jy = 0; jy < mesh->LocalNy; jy++) {
            for (jz = 0; jz < mesh->LocalNz; jz++) {
              Field3D Te_tmp_real = Te_tmp * Tebar;   //eV
              BoutReal Te_tmp_loc  = Te_tmp_real(jx,jy,jz);
              BoutReal logT = log(Te_tmp_loc);
              BoutReal log_out =0;
              if (Te_tmp_loc >= 1.0 && Te_tmp_loc <= 500.) {
                log_out = - 7.87837896e+01 * pow(logT, 0)
                          + 1.55326376e+00 * pow(logT, 1)
                          + 1.65898194e+01 * pow(logT, 2)
                          - 3.23804546e+01 * pow(logT, 3)
                          + 3.12784663e+01 * pow(logT, 4)
                          - 1.74826039e+01 * pow(logT, 5)
                          + 5.91393245e+00 * pow(logT, 6)
                          - 1.22974105e+00 * pow(logT, 7)        
                          + 1.54004499e-01 * pow(logT, 8)
                          - 1.06797106e-02 * pow(logT, 9)
                          + 3.15657594e-04 * pow(logT, 10);
                Limp(jx,jy,jz) = exp(log_out);
              } else if (Te_tmp_loc < 1.0) {
                Limp(jx,jy,jz) = 6.00623928e-35;
              } else {
                Limp(jx,jy,jz) = 4.53057707e-33;
              }
	          }
	        }
	      }
      }

      if (Limp_nitro) {
        /// Nitrogen based cooling curve used in Lipschultz 2016
        for (jx = 0; jx < mesh->LocalNx; jx++) {
          for (jy = 0; jy < mesh->LocalNy; jy++) {
            for (jz = 0; jz < mesh->LocalNz; jz++) {
              Field3D Te_tmp_real = Te_tmp * Tebar;
              BoutReal Te_tmp_loc  = Te_tmp_real(jx,jy,jz);
              if (Te_tmp_loc > 1.0 && Te_tmp_loc < 80.) {
                Limp(jx,jy,jz) = 5.9e-34 * (sqrt(Te_tmp_loc-1.)) * (80. - Te_tmp_loc) / (3.1e-3 * ((Te_tmp_loc - 1.) * (Te_tmp_loc - 1.)) + 1.0);     //Wm^3
              } else {
                Limp(jx,jy,jz) = 0.0;
              }
            }
          }
        }
      }

      if (Limp_nitro_adas) {
	      for (jx = 0; jx < mesh->LocalNx; jx++) {
	        for (jy = 0; jy < mesh->LocalNy; jy++) {
            for (jz = 0; jz < mesh->LocalNz; jz++) {
              Field3D Te_tmp_real = Te_tmp * Tebar;
              BoutReal log_out =0;
              BoutReal Te_tmp_loc  = Te_tmp_real(jx,jy,jz);
              BoutReal logT = log(Te_tmp_loc);
              if (Te_tmp_loc >= 2.0 && Te_tmp_loc <= 500.) {
                log_out = - 5.01649969e+01 * pow(logT, 0)
                          - 1.35749724e+02 * pow(logT, 1)
                          + 2.73509608e+02 * pow(logT, 2)
                          - 2.92109992e+02 * pow(logT, 3)
                          + 1.90120639e+02 * pow(logT, 4)
                          - 7.95164871e+01 * pow(logT, 5)
                          + 2.17762218e+01 * pow(logT, 6)
                          - 3.88334992e+00 * pow(logT, 7)        
                          + 4.34730098e-01 * pow(logT, 8)
                          - 2.77683605e-02 * pow(logT, 9)                                                   
                          + 7.72720422e-04 * pow(logT, 10);
                Limp(jx,jy,jz) = exp(log_out);
              } else if (Te_tmp_loc < 2.0) {
                Limp(jx,jy,jz) = 4.34835380e-34;
              } else {
                Limp(jx,jy,jz) = 8.11096182e-33;
              }
            }
          }
        }
      }

      if (Limp_Neon) {
        /// Neon based cooling curve produced by Matlab polynominal curve
	      // fitting "polyval" (Ryoko 2020 Nov)
	      for (jx = 0; jx < mesh->LocalNx; jx++) {
	        for (jy = 0; jy < mesh->LocalNy; jy++) {
            for (jz = 0; jz < mesh->LocalNz; jz++) {
	            Field3D Te_tmp_real = Te_tmp * Tebar;   //eV
              BoutReal Te_tmp_loc  = Te_tmp_real(jx,jy,jz);
	            if (Te_tmp_loc >= 3 and Te_tmp_loc <= 100) {
	              Limp(jx,jy,jz) = - 2.0385e-40 * pow(Te_tmp_loc, 5)
                                 + 5.4824e-38 * pow(Te_tmp_loc, 4)
                                 - 5.1190e-36 * pow(Te_tmp_loc, 3)
                                 + 1.7347e-34 * SQ(Te_tmp_loc)
                                 - 3.4151e-34 * Te_tmp_loc
                                 - 3.2798e-34;
              } else if (Te_tmp_loc >=2 and Te_tmp_loc < 3) {
                Limp(jx,jy,jz) = 7e-35 * (Te_tmp_loc - 2.) + 1e-35;
              } else if (Te_tmp_loc >=1 and Te_tmp_loc < 2) {
                Limp(jx,jy,jz) = 1e-35 * (Te_tmp_loc - 1.);
              } else {
                Limp(jx,jy,jz) = 0.0;
              }
	          }
	        }
	      }
      }

      if (Limp_Neon_adas) {
	      for (jx = 0; jx < mesh->LocalNx; jx++) {
	        for (jy = 0; jy < mesh->LocalNy; jy++) {
            for (jz = 0; jz < mesh->LocalNz; jz++) {
              Field3D Te_tmp_real = Te_tmp * Tebar;   //eV 
              BoutReal log_out = 0;
              BoutReal Te_tmp_loc  = Te_tmp_real(jx,jy,jz);
              BoutReal logT = log(Te_tmp_loc);
              if (Te_tmp_loc >= 2.0 && Te_tmp_loc <= 1000.) {
                log_out = - 8.21475117e+01 * pow(logT, 0)
                          + 1.28929854e+01 * pow(logT, 1)
                          - 4.74266289e+01 * pow(logT, 2)
                          + 7.45222324e+01 * pow(logT, 3)
                          - 5.75710722e+01 * pow(logT, 4)
                          + 2.57375965e+01 * pow(logT, 5)
                          - 7.12758563e+00 * pow(logT, 6)
                          + 1.24287546e+00 * pow(logT, 7)
                          - 1.32943407e-01 * pow(logT, 8)
                          + 7.97368445e-03 * pow(logT, 9)
                          - 2.05487897e-04 * pow(logT, 10);
                Limp(jx,jy,jz) = exp(log_out);
              } else if (Te_tmp_loc < 2.0) {
                Limp(jx,jy,jz) = 6.35304113e-36;
              } else {
                Limp(jx,jy,jz) = 1.17894628e-32;
              }
            }
          }
        }
      }

      if (Limp_Argon) {
        /// Argon based cooling curve produced by Matlab polynominal curve
	      // fitting "polyval" (Ryoko 2020 Nov)
	      for (jx = 0; jx < mesh->LocalNx; jx++) {
	        for (jy = 0; jy < mesh->LocalNy; jy++) {
            for (jz = 0; jz < mesh->LocalNz; jz++) {
              Field3D Te_tmp_real = Te_tmp * Tebar;
              BoutReal Te_tmp_loc  = Te_tmp_real(jx,jy,jz);
              if (Te_tmp_loc >= 1.5 and Te_tmp_loc <= 100) {
                Limp(jx,jy,jz) = - 4.9692e-48 * pow(Te_tmp_loc, 10)
                                      + 2.8025e-45 * pow(Te_tmp_loc, 9)
                                      - 6.7148e-43 * pow(Te_tmp_loc, 8)
                                      + 8.8636e-41 * pow(Te_tmp_loc, 7)
                                      - 6.9642e-39 * pow(Te_tmp_loc, 6)
                                      + 3.2559e-37 * pow(Te_tmp_loc, 5)
                                      - 8.3410e-36 * pow(Te_tmp_loc, 4)
                                      + 8.6011e-35 * pow(Te_tmp_loc, 3)
                                      + 1.9958e-34 * pow(Te_tmp_loc, 2)
                                      + 4.9864e-34 * Te_tmp_loc
                                      - 9.9412e-34;
              } else if (Te_tmp_loc >= 1.0 and Te_tmp_loc < 1.5) {
                Limp(jx,jy,jz) = 5e-35 * (Te_tmp_loc - 1.0);
              } else {
                Limp(jx,jy,jz) = 0.0;
              }
            }
          }
        }
      }

      if (Limp_Argon_adas) {
	      for (jx = 0; jx < mesh->LocalNx; jx++) {
	        for (jy = 0; jy < mesh->LocalNy; jy++) {
            for (jz = 0; jz < mesh->LocalNz; jz++) {
              Field3D Te_tmp_real = Te_tmp * Tebar;   //eV
              BoutReal log_out = 0;
              BoutReal Te_tmp_loc  = Te_tmp_real(jx,jy,jz);
              BoutReal logT = log(Te_tmp_loc);
              if (Te_tmp_loc >= 1.5 && Te_tmp_loc <= 1500.) {
                log_out = - 8.45410692e+01 * pow(logT, 0)
                          + 1.57727040e+01 * pow(logT, 1)
                          - 1.54264860e+01 * pow(logT, 2)
                          + 1.49409902e+01 * pow(logT, 3)
                          - 1.04815113e+01 * pow(logT, 4)
                          + 5.00924595e+00 * pow(logT, 5)
                          - 1.60029106e+00 * pow(logT, 6)
                          + 3.29455609e-01 * pow(logT, 7)
                          - 4.14036827e-02 * pow(logT, 8)
                          + 2.87063206e-03 * pow(logT, 9)
                          - 8.38888002e-05 * pow(logT, 10);
                Limp(jx,jy,jz) = exp(log_out);
              } else if (Te_tmp_loc < 1.5) {
                      Limp(jx,jy,jz) = 1.95353412e-35;
              } else {
                      Limp(jx,jy,jz) = 1.22649600e-32;
              }
            }
          }
        }
      }
      
      Srad = -Limp * frac_imp * Ne_tmp * Ne_tmp * Nbar * density / (1.602e-19 * Tebar / Tbar);        //Normalized //NOTE(malamast): Should we multiply with Nbar * density?
      mesh->communicate(Srad);
      Srad.applyBoundary("neumann");
      Wrad = Srad / Ne_tmp;
      mesh->communicate(Wrad);
      Wrad.applyBoundary("neumann");
    }

    if (jpar_bndry_width > 0) {
      // Zero j in boundary regions. Prevents vorticity drive
      // at the boundary
      for (jx = 0; jx < jpar_bndry_width; jx++)
        for (jy = 0; jy < mesh->LocalNy; jy++)
          for (jz = 0; jz < mesh->LocalNz - 1; jz++) {
            if (mesh->firstX())
              Jpar(jx,jy,jz) = 0.0;
            if (mesh->lastX())
              Jpar(mesh->LocalNx - 1 - jx,jy,jz) = 0.0;
	        }
    }

    if (smooth_j_x && (gamma_i_BC > 0.0) && (gamma_e_BC > 0.0)) {
      Jpar = smooth_x(Jpar);
      /*Jpar = smooth_x(Jpar);
      Jpar = smooth_x(Jpar);
      Jpar = smooth_x(Jpar);
      Jpar = smooth_x(Jpar);
      Jpar = smooth_y(Jpar);
      Jpar = smooth_y(Jpar);
      Jpar = smooth_y(Jpar);
      Jpar = smooth_y(Jpar);
      Jpar = smooth_y(Jpar);
      */
      mesh->communicate(Jpar);
      Jpar.applyBoundary();
    }

    if (mask_j_x) {
      Jpar *= mask_jx1d;
    }
    
    if (compress0) {
      if (nonlinear) {
        if (include_vpar0) {
          Vepar = - Jpar / Ne_tmp * Vepara + N_tmp / Ne_tmp * Vipar + Ni / Ne_tmp * Vipar0 -Zi * Ni / Ne_tmp * Vepar0 ; 
        } else {
          // Vepar = Vipar - B0 * (Jpar) / N_tmp * Vepara;
          Vepar = Vipar - Jpar / Ne_tmp * Vepara; //NOTE(malamast): That's weird. Why was it like that???
        }
      } else {
        // Vepar = Vipar - B0 * (Jpar) / N0 * Vepara; //NOTE(malamast): That's weird. Why was it like that???
        Vepar = Vipar - Jpar / Ne0 * Vepara;  
      }
      Vepar.applyBoundary();
      mesh->communicate(Vepar);
    }

    // Get Delp2(J) from J
    Jpar2 = -Delp2(Jpar); // NOTE(malamast): Do we actually use that?

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

  #if DEBUG_6F>0
    output.write("I see you 1!\n");
  #endif

    /////////////////////////////////////////////////////////////////
    // Parallel electric field
    {
      if (!emass) {
        if (evolve_psi) {
          TRACE("ddt(Psi)");
	  
	        ddt(Psi) = 0.0;
          ddt(Psi) = -Grad_parP(phi) / B0 - eta * Jpar / B0;

          if (diamag && diamag_phi0) {
            if (diamag_er)
              ddt(Psi) -= V_dot_Grad(Ve0, Psi);
            else
              ddt(Psi) -= bracket(phi0, Psi, bm_exb); // Equilibrium flow  //NOTE(malamast): We will need to add this term when taken from 2D and diamag_phi0 is false. 
          }
	  
          if (thermal_force) {
            // grad_par(T_e) correction
            ddt(Psi) += 0.71 * Psipara1 * Grad_parP(Te) / B0;
            ddt(Psi) -= 0.71 * Psipara1 * bracket(Psi, Te0, bm_mag);
          }

          if (eHall) {
            ddt(Psi) += Psipara1 * Grad_parP(Pe) / B0 / Ne0;
            ddt(Psi) -= Psipara1 * bracket(Psi, Pe0, bm_mag) / Ne0;
          }

          if (hyperdiff_par_apar4 > 0.0) {
            tmpA2 = Grad2_par2new(Psi);
            mesh->communicate(tmpA2);
            tmpA2.applyBoundary();
            ddt(Psi) -= hyperdiff_par_apar4 * Grad2_par2new(tmpA2);
          }

	        // Hyper-resistivity
          if (hyperresist > 0.0) {
            ddt(Psi) += hyperresist * Delp2(Jpar / B0);
	        }

	      } else {
          TRACE("ddt(Apar)");

	        ddt(Apar) = 0.0;
          ddt(Apar) = -Grad_parP(phi) - eta * Jpar;

          if (diamag && diamag_phi0) {
            if (diamag_er)
              ddt(Apar) -= V_dot_Grad(Ve0, Apar);
            else
              ddt(Apar) -= bracket(phi0, Apar, bm_exb); // Equilibrium flow //NOTE(malamast): We will need to add this term when taken from 2D and diamag_phi0 is false. 
          }

          if (thermal_force) {
            // grad_par(T_e) correction
            ddt(Apar) += 0.71 * Psipara1 * Grad_parP(Te);
            ddt(Apar) -= 0.71 * Psipara1 * bracket(Apar, Te0, bm_mag);
          }

          if (eHall) {
            ddt(Apar) += Psipara1 * Grad_parP(Pe) / Ne0; // NOTE(malamast): Why is Ne0 and not the total?
            ddt(Apar) -= Psipara1 * bracket(Apar, Pe0, bm_mag) / Ne0;  
          }

          if (hyperdiff_par_apar4 > 0.0) {
            tmpA2 = Grad2_par2new(Apar);
            mesh->communicate(tmpA2);
            tmpA2.applyBoundary();
            ddt(Apar) -= hyperdiff_par_apar4 * Grad2_par2new(tmpA2);
	        }

	        // Hyper-resistivity
          if (hyperresist > 0.0) {
            ddt(Apar) += hyperresist * Delp2(Jpar);
	        }
	      }
      } else { // emass

	      ddt(Ajpar) = 0.0;
      	ddt(Ajpar) = -Grad_parP(phi) / B0 - eta * Jpar / B0;

	      if (diamag && diamag_phi0) {
          if (diamag_er)
            ddt(Ajpar) -= V_dot_Grad(Ve0, Psi);
          else
            ddt(Ajpar) -= bracket(phi0, Psi, bm_exb); // Equilibrium flow
        }

        if (thermal_force) {
          // grad_par(T_e) correction
          ddt(Ajpar) += 0.71 * Psipara1 * Grad_parP(Te) / B0;
          ddt(Ajpar) -= 0.71 * Psipara1 * bracket(Psi, Te0, bm_mag);
	      }

        if (eHall) {
          ddt(Ajpar) += Psipara1 * Grad_parP(Pe) / B0 / Ne0;
          ddt(Ajpar) -= Psipara1 * bracket(Psi, Pe0, bm_mag) / Ne0;
        }

	      // Hyper-resistivity
        if (hyperresist > 0.0) {
          ddt(Ajpar) += hyperresist * Delp2(Jpar / B0);
        }
      }

      if (output_ohm) {
        ohm_phi = -Grad_parP(phi) / B0 - bracket(phi0, Psi, bm_exb);
        mesh->communicate(ohm_phi);
        ohm_phi.applyBoundary();
        ohm_hall = Psipara1 * (Grad_parP(Pe) / (B0 * Ne0) - bracket(Psi, Pe0, bm_mag) / Ne0);
        mesh->communicate(ohm_hall);
        ohm_hall.applyBoundary();
        ohm_thermal = 0.71 * Psipara1 * (Grad_parP(Te) / B0 - bracket(Psi, Te0, bm_mag));
        mesh->communicate(ohm_thermal);
        ohm_thermal.applyBoundary();
      }
    }

  #if DEBUG_6F>0
    output.write("I see you 2!\n");
  #endif

    //////////////////////////////////////////////////////////////////
    // Vorticity equation

    {
      TRACE("ddt(U)");

      ddt(U) = 0.0;

      if (BScurrent) {
        if (evolve_psi)
          ddt(U) = -SQ(B0) * bracket(Psi, Jpar_BS0 / B0, bm_mag) * B0;
	          else  
          ddt(U) = -SQ(B0) * bracket(Apar, Jpar_BS0 / B0, bm_mag);
      } else {
        if (evolve_psi)
	        ddt(U) = -SQ(B0) * bracket(Psi, J0 / B0, bm_mag) * B0; // Grad j term
	      else
          ddt(U) = -SQ(B0) * bracket(Apar, J0 / B0, bm_mag);     // Grad j term
      }

      ddt(U) += 2.0 * Upara1 * b0xcv * Grad(P);  // curvature term

      ddt(U) += SQ(B0) * Grad_parP(Jpar / B0);   // b dot grad j

      if (diamag && diamag_phi0) {
	      if (diamag_er)
          ddt(U) -= V_dot_Grad(Ve0, U);
	      else
          ddt(U) -= bracket(phi0, U, bm_exb);    // Equilibrium flow
      }

      if (experiment_Er && KH_term)
        ddt(U) -= bracket(phi, U0_net, bm_exb); //NOTE(malamast): U0_net is zero in our case.

      if (include_U0) {
        ddt(U) -= bracket(phi, U0, bm_exb); // Advection // Added by malamast
      }

      if (compress0 && include_vipar) {

        ddt(U) -= Vipar * Grad_parP(U0_net);  //NOTE(malamast): U0_net is zero in our case.

        if (include_U0) {
          ddt(U) -= Vipar * Grad_parP(U0); // Added by malamast
        }

        if (include_vpar0) {
          ddt(U) -= Vipar0 * Grad_parP(U);  //NOTE(malamast): U0_net is zero in our case.
        }

        if (include_U0 && include_vpar0) {
          if (evolve_psi)
            ddt(U) += Vipar0 * bracket(Psi, U0, bm_mag) * B0; // Added by malamast
          else
            ddt(U) += Vipar0 * bracket(Apar, U0, bm_mag); // Added by malamast
        }

      }
      
      if (nonlinear) {
        ddt(U) -= bracket(phi, U, bm_exb);       // Advection

	      if (compress0 && include_vipar) {
	        ddt(U) -= Vpar_Grad_par(Vipar, U);
          if (evolve_psi)
            ddt(U) += Vipar * bracket(Psi, U, bm_mag) * B0; // Added by malamast
          else
            ddt(U) += Vipar * bracket(Apar, U, bm_mag); //  Added by malamast
        }
      }

      // parallel hyper-viscous diffusion for vector potential
      if (hyperdiff_par_u4 > 0.0) {
        tmpU2 = Grad2_par2new(U);
        mesh->communicate(tmpU2);
        tmpU2.applyBoundary();
        ddt(U) -= hyperdiff_par_u4 * Grad2_par2new(tmpU2);
      }

      if (hyperdiff_perp_u4 > 0.0) {
        tmpU2 = Delp2(U);
        mesh->communicate(tmpU2);
        tmpU2.applyBoundary();
        ddt(U) -= hyperdiff_perp_u4 * Delp2(tmpU2);
      }

      if (parallel_viscous && compress0) {
        ddt(U) -= 0.666667 * Upara1 * b0xcv * Grad(pi_ci); //NOTE(malamast): Should we exclude the background?
      }

      if (gyroviscous) {
        if (diamag_er && diamag_phi0)
          Dperp2Phi0 = Div(cross(Ve0, B0vec));
        else
          Dperp2Phi0 = Field3D(Delp2(phi0));
        Dperp2Phi0.applyBoundary();
        mesh->communicate(Dperp2Phi0);

	      Dperp2Phi = Delp2(phi);
    	  Dperp2Phi.applyBoundary();
        mesh->communicate(Dperp2Phi);
        
	      if (diamag_er && diamag_phi0) {
          // GradPhi02 = Er0*Er0; // test
          GradPhi02 = (cross(Ve0, B0vec)) * (cross(Ve0, B0vec)) / (B0 * B0);
	      } else {
          GradPhi02 = Field3D(Grad_perp(phi0) * Grad_perp(phi0) / (B0 * B0));
	      }
	      GradPhi02.applyBoundary();
        mesh->communicate(GradPhi02);

        if (diamag_er && diamag_phi0) {
          GradcPhi = (cross(Ve0, B0vec)) * Grad_perp(phi) / (B0 * B0);
        } else {
          GradcPhi = Grad_perp(phi0) * Grad_perp(phi) / (B0 * B0);
        }
        GradcPhi.applyBoundary();
        mesh->communicate(GradcPhi);

	      Dperp2Pi0 = Field3D(Delp2(Pi0));
        Dperp2Pi0.applyBoundary();
        mesh->communicate(Dperp2Pi0);

        Dperp2Pi = Delp2(Pi);
        Dperp2Pi.applyBoundary();
        mesh->communicate(Dperp2Pi);

	      if (diamag_er && diamag_phi0)
	        bracketPhi0P = V_dot_Grad(Ve0, Pi);
        else
          bracketPhi0P = bracket(phi0, Pi, bm_exb);
        bracketPhi0P.applyBoundary();
        mesh->communicate(bracketPhi0P);
    
	      bracketPhiP0 = bracket(phi, Pi0, bm_exb);
    	  bracketPhiP0.applyBoundary();
	      mesh->communicate(bracketPhiP0);

	      ddt(U) -= 0.5 * Upara2 * bracket(Pi, Dperp2Phi0, bm_exb) / B0;
    	  ddt(U) -= 0.5 * Upara2 * bracket(Pi0, Dperp2Phi, bm_exb) / B0;
        ddt(U) += Upara3 * B0 * bracket(N0, GradcPhi, bm_exb);
        ddt(U) += 0.5 * Upara3 * B0 * bracket(Ni, GradPhi02, bm_exb);
        ddt(U) += 0.5 * Upara2 * bracket(phi, Dperp2Pi0, bm_exb) / B0;
        if (diamag_er && diamag_phi0)
          ddt(U) += 0.5 * Upara2 * V_dot_Grad(Ve0, Dperp2Pi) / B0;
        else
          ddt(U) += 0.5 * Upara2 * bracket(phi0, Dperp2Pi, bm_exb) / B0;
        ddt(U) -= 0.5 * Upara2 * Delp2(bracketPhi0P) / B0;
        ddt(U) -= 0.5 * Upara2 * Delp2(bracketPhiP0) / B0;

	      if (impurity_prof && impurity_gyro) {
	        Dperp2Pimp0 = Field3D(Delp2(P_imp0));
    	    Dperp2Pimp0.applyBoundary();
          mesh->communicate(Dperp2Pimp0);

          bracketPhiPimp0 = bracket(phi, P_imp0, bm_exb);
          bracketPhiPimp0.applyBoundary();
          mesh->communicate(bracketPhiPimp0);

          ddt(U) -= 0.5 * Upara_imp * Upara2 * bracket(P_imp0, Dperp2Phi, bm_exb) / B0;
          ddt(U) += Upara_imp * Upara3 * B0 * bracket(N_imp0, GradcPhi, bm_exb);
          ddt(U) += 0.5 * Upara_imp * Upara2 * bracket(phi, Dperp2Pimp0, bm_exb) / B0;
          ddt(U) -= 0.5 * Upara_imp * Upara2 * Delp2(bracketPhiPimp0) / B0;
	      }

	      if (nonlinear) {
          GradPhi2 = Grad_perp(phi) * Grad_perp(phi) / (B0 * B0);
          GradPhi2.applyBoundary();
          mesh->communicate(GradPhi2);

          bracketPhiP = bracket(phi, Pi, bm_exb);
          bracketPhiP.applyBoundary();
          mesh->communicate(bracketPhiP);

          ddt(U) -= 0.5 * Upara2 * bracket(Pi, Dperp2Phi, bm_exb) / B0;
          ddt(U) += 0.5 * Upara3 * B0 * bracket(N0, GradPhi2, bm_exb);
          ddt(U) += Upara3 * B0 * bracket(Ni, GradcPhi, bm_exb);
          ddt(U) += 0.5 * Upara2 * bracket(phi, Dperp2Pi, bm_exb) / B0;
          ddt(U) -= 0.5 * Upara2 * Delp2(bracketPhiP) / B0;

	        if (impurity_prof && impurity_gyro) {
            ddt(U) += 0.5 * Upara_imp * Upara3 * B0 * bracket(N_imp0, GradPhi2, bm_exb);
	        }
	      }
      }

      if (output_transfer) {
        T_R = -bracket(phi, U - Upara0 / B0 * Dperp2Pi, bm_exb);
        mesh->communicate(T_R);
        T_R.applyBoundary();
        T_M = SQ(B0) * Grad_parP(Jpar / B0);
        mesh->communicate(T_M);
        T_M.applyBoundary();
        T_ID = -bracket(phi, Upara0 / B0 * Dperp2Pi, bm_exb);
        mesh->communicate(T_ID);
        T_ID.applyBoundary();
        T_C = 2.0 * Upara1 * b0xcv * Grad(P);
        mesh->communicate(T_C);
        T_C.applyBoundary();
        T_G = -0.5 * Upara2 * bracket(Pi, Dperp2Phi, bm_exb) / B0;
        T_G += 0.5 * Upara3 * B0 * bracket(N0, GradPhi2, bm_exb);
        T_G += Upara3 * B0 * bracket(Ni, GradcPhi, bm_exb);
        T_G += 0.5 * Upara2 * bracket(phi, Dperp2Pi, bm_exb) / B0;
        T_G -= 0.5 * Upara2 * Delp2(bracketPhiP) / B0;
        T_G.applyBoundary();
        mesh->communicate(T_G);
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

  #if DEBUG_6F>0
    output.write("I see you 3!\n");
  #endif
    /////////////////////////////////////////////////////////////////
    // number density equation

    {
      TRACE("ddt(Ni)");

      ddt(Ni) = 0.0;

      if (NiAmp > 0) ddt(Ni) += NiSource;

      ddt(Ni) -= bracket(phi, N0, bm_exb);

      if (continuity) {
        // NOTE(malamast):  Grad(phi) decomposed into parallel and perpendicular parts. 
        // NOTE(malamast):  The contribution of the parallel part is zero here since we multiply with b0xcv. 
        ddt(Ni) -= 2.0 * N0 / B0 * b0xcv * Grad(phi);
        ddt(Ni) -= 2.0 * Ni / B0 * b0xcv * Grad(phi0); //NOTE(malamast): This was commented out for some reason.
        // ddt(Ni) -= 2.0 * Nipara1 * b0xcv * Grad(Pi0) / N0 / B0;    // balanced by above
        // ddt(Ni) -= 2.0 * Ni / B0 * b0xcv * cross(Ve0_net, B0vec); // net flow //NOTE(malamast): Ve0_net is 0 so I have commented it out.
        if (diamag) {
          ddt(Ni) -= 2.0 * Nipara1 * b0xcv * Grad(Pi) / B0;
        }
        if (nonlinear) {
          ddt(Ni) -= 2.0 * Ni / B0 * b0xcv * Grad(phi);
        }
      }

      if (diamag && diamag_phi0) {
	      if (diamag_er)
	        ddt(Ni) -= V_dot_Grad(Ve0, Ni);
	      else 
	        ddt(Ni) -= bracket(phi0, Ni, bm_exb); // Equilibrium flow
      }

      if (nonlinear) {
        ddt(Ni) -= bracket(phi, Ni, bm_exb); // Advection
      }

      if (compress0) {
        // ddt(Ni) -= Vipar * Grad_parP(N0, CELL_YLOW);
        // ddt(Ni) -= Vpar_Grad_par(Vipar, N0);        
        if (include_vipar) {
          ddt(Ni) -= Vipar * Grad_parP(N0);

          if (include_vpar0) {
            ddt(Ni) -= Vipar0 * Grad_parP(Ni);
            if (evolve_psi)
              ddt(Ni) += Vipar0 * bracket(Psi, N0, bm_mag) * B0; 
            else
              ddt(Ni) += Vipar0 * bracket(Apar, N0, bm_mag); 
          }
        }

        if (continuity) {
          ddt(Ni) -= N0 * B0 * Grad_parP(Vipar / B0);

          if (include_vpar0) {
            ddt(Ni) -= Ni * B0 * Grad_parP(Vipar0 / B0);
            if (evolve_psi)
              ddt(Ni) += N0 * B0 * bracket(Psi, Vipar0, bm_mag) * B0; 
            else
              ddt(Ni) += N0 * B0 * bracket(Apar, Vipar0, bm_mag); 
          }
        }

        if (nonlinear) {
          // ddt(Ni) -= Vipar * Grad_par(Ni, CELL_YLOW);

          // ddt(Ni) -= Vpar_Grad_par(Vipar, Ni);
          // ddt(Ni) += Vipar * bracket(Psi, N0, bm_mag)*B0;      

          if (include_vipar) {
            ddt(Ni) -= Vpar_Grad_par(Vipar, Ni); //NOTE(malamast): do we include correction for the perturbed magnetic field like we do in Grad_parP?

            if (evolve_psi)
              ddt(Ni) += Vipar * bracket(Psi, Ni, bm_mag) * B0; //NOTE(malamast): Added by malamast
            else
              ddt(Ni) += Vipar * bracket(Apar, Ni, bm_mag); //NOTE(malamast): Added by malamast
          }

          if (continuity)
            // ddt(Ni) -= Ni * B0 * Grad_par(Vipar / B0); //NOTE(malamast): why not Grad_parP?
            ddt(Ni) -= Ni * B0 * Grad_parP(Vipar / B0);
        }
      }

      if (radial_diffusion)
	      ddt(Ni) += diff_radial * Delp2(Ni);

      /*if (neoclassic_i) {
        tmpddx2 = D2DX2(DC(Ni));
        mesh->communicate(tmpddx2);
        tmpddx2.applyBoundary();
        ddt(Ni) += coord->g11*DC(Dri_neo * tmpddx2);

        partf_neo_i = -Dri_neo * sqrt(coord->g11)*DDX(Ni);
        mesh->communicate(partf_neo_i);
        partf_neo_i.applyBoundary();
      }*/

      // 4th order Parallel diffusion terms
      if (hyperdiff_par_n4 > 0.0) {
        tmpN2 = Grad2_par2new(Ni);
        mesh->communicate(tmpN2);
        tmpN2.applyBoundary();
        ddt(Ni) -= hyperdiff_par_n4 * Grad2_par2new(tmpN2);
      }

      if (hyperdiff_perp_n4 > 0.0) {
        tmpN2 = Delp2(Ni);
        mesh->communicate(tmpN2);
        tmpN2.applyBoundary();
        ddt(Ni) -= hyperdiff_perp_n4 * Delp2(tmpN2);
      }
    }

  #if DEBUG_6F>0
    output.write("I see you 4!\n");
  #endif
    //////////////////////////////////////////////////////////////////
    // ion temperature equation
    {
      TRACE("ddt(Ti)");

      ddt(Ti) = 0.0;

      if (TiAmp > 0) ddt(Ti) += TiSource;

      ddt(Ti) -= bracket(phi, Ti0, bm_exb);

      if (continuity) {
        ddt(Ti) -= 4.0 / 3.0 * Ti0 / B0 * b0xcv * Grad(phi);
        // ddt(Ti) -= 4.0/3.0 * Ti * b0xcv*Grad(phi0*B0) / B0; // ---> ?
        // ddt(Ti) -= 4.0/3.0 * Tipara2 * Ti/N0 * b0xcv*Grad(Pi0) / (N0*B0);   // balanced by above ---> ?        
	      // ddt(Ti) -= 4.0 / 3.0 * Ti / B0 * b0xcv * cross(Ve0_net, B0vec); // net flow // NOTE(malamast): Ve0_net is 0 ---> I have commented it out.

        ddt(Ti) -= 4.0 / 3.0 * Ti / B0 * b0xcv * Grad(phi0); // Added by malamast.
        ddt(Ti) -= 4.0 / 3.0 * Tipara2 * Ti / N0 * b0xcv * Grad(Pi0) / B0;  // Added by malamast.

        if (diamag)
          ddt(Ti) -= 4.0 / 3.0 * Tipara2 * Ti0 / N0 * b0xcv * Grad(Pi) / B0;
        if (nonlinear) {
          ddt(Ti) -= 4.0 / 3.0 * Ti / B0 * b0xcv * Grad(phi);
          if (diamag)
            ddt(Ti) -= 4.0 / 3.0 * Tipara2 * Ti / N0 * b0xcv * Grad(Pi) / B0;
        }
      }

      if (energy_flux) {
        ddt(Ti) -= 10.0 / 3.0 * Tipara2 / B0 * V_dot_Grad(Ti0 * b0xcv, Ti);
        ddt(Ti) -= 10.0 / 3.0 * Tipara2 * Ti / B0 * b0xcv * Grad(Ti0);
        if (nonlinear)
          ddt(Ti) -= 10.0 / 3.0 * Tipara2 / B0 * V_dot_Grad(Ti * b0xcv, Ti);
      }

      if (diamag && diamag_phi0) {
        if (diamag_er)
  	      ddt(Ti) -= V_dot_Grad(Ve0, Ti);
	      else
          ddt(Ti) -= bracket(phi0, Ti, bm_exb); // Equilibrium flow
      }

      if (nonlinear) {
        ddt(Ti) -= bracket(phi, Ti, bm_exb); // Advection
      }

      if (compress0) {
        // ddt(Ti) -= Vipar * Grad_parP(Ti0);
        // ddt(Ti) -= Vpar_Grad_par(Vipar, Ti0);        
        if (include_vipar) {
          ddt(Ti) -= Vipar * Grad_parP(Ti0);

          if (include_vpar0) {
            ddt(Ti) -= Vipar0 * Grad_parP(Ti);
            if (evolve_psi)
              ddt(Ti) += Vipar0 * bracket(Psi, Ti0, bm_mag) * B0; 
            else
              ddt(Ti) += Vipar0 * bracket(Apar, Ti0, bm_mag); 
          }          
        }

        if (continuity) {
          ddt(Ti) -= 2.0 / 3.0 * Ti0 * B0 * Grad_parP(Vipar / B0);

          if (include_vpar0) {
            ddt(Ti) -= 2.0 / 3.0 * Ti * B0 * Grad_parP(Vipar0 / B0);
            if (evolve_psi)
              ddt(Ti) += 2.0 / 3.0 * Ti0 * B0 * bracket(Psi, Vipar0 / B0, bm_mag) * B0; 
            else
              ddt(Ti) += 2.0 / 3.0 * Ti0 * B0 * bracket(Apar, Vipar0 / B0, bm_mag); 
          }
        }

        if (nonlinear) {
          // ddt(Ti) -= Vipar * Grad_par(Ti, CELL_YLOW);

          // ddt(Ti) -= Vpar_Grad_par(Vipar, Ti);
          // ddt(Ti) += Vipar * bracket(Psi, Ti0, bm_mag)*B0;

          if (include_vipar) {
            ddt(Ti) -= Vpar_Grad_par(Vipar, Ti);
            if (evolve_psi)
              ddt(Ti) += Vipar * bracket(Psi, Ti, bm_mag) * B0; // Added by malamast
            else
              ddt(Ti) += Vipar * bracket(Apar, Ti, bm_mag); // Added by malamast            
          }

          if (continuity)
            ddt(Ti) -= 2.0 / 3.0 * Ti * B0 * Grad_parP(Vipar / B0);
        }

      }

      if (energy_exch) {
        ddt(Ti) += 2.0 * Zi * Tbar * nu_e / (ratio_pe * AA) * (Te - Ti);
      }

      if (diffusion_par > 0.0) {

        // if (Landau) {
        //   if (diff_par_flutter)
        //     ddt(Ti) -= Grad_parP(q_par_i) / N0;
        //   else
        //     ddt(Ti) -= Grad_par(q_par_i) / N0;
        // } else {

          ddt(Ti) += kappa_par_i * Grad2_par2(Ti) / N0; // Parallel diffusion
          ddt(Ti) += Grad_par(kappa_par_i) * Grad_par(Ti) / N0;

          if (diff_par_flutter) {
            if (nonlinear) {
              if (evolve_psi)
                bracket1i = -bracket(Psi, Ti + Ti0, bm_mag) * B0;
              else
                bracket1i = -bracket(Apar, Ti + Ti0, bm_mag);
            } else {
              if (evolve_psi)
                bracket1i = -bracket(Psi, Ti0, bm_mag) * B0;
              else
                bracket1i = -bracket(Apar, Ti0, bm_mag);
            }
            mesh->communicate(bracket1i);
            bracket1i.applyBoundary();
            gradpar_ti = Grad_par(Ti); //NOTE(malamast): Do I need to add Ti0 here?
            mesh->communicate(gradpar_ti);
            gradpar_ti.applyBoundary();

            ddt(Ti) += Grad_par(kappa_par_i) * bracket1i / N0;
            ddt(Ti) += kappa_par_i * Grad_par(bracket1i) / N0;
            if (nonlinear) {
              if (evolve_psi) {
                ddt(Ti) -= bracket(Psi, kappa_par_i, bm_mag) * B0 * gradpar_ti / N0;
                ddt(Ti) -= kappa_par_i * bracket(Psi, gradpar_ti, bm_mag) * B0 / N0;
                ddt(Ti) -= bracket(Psi, kappa_par_i, bm_mag) * B0 * bracket1i / N0;
                ddt(Ti) -= kappa_par_i * bracket(Psi, bracket1i, bm_mag) * B0 / N0;
              } else {
                ddt(Ti) -= bracket(Apar, kappa_par_i, bm_mag) * gradpar_ti / N0;
                ddt(Ti) -= kappa_par_i * bracket(Apar, gradpar_ti, bm_mag) / N0;
                ddt(Ti) -= bracket(Apar, kappa_par_i, bm_mag) * bracket1i / N0;
                ddt(Ti) -= kappa_par_i * bracket(Apar, bracket1i, bm_mag) / N0;
              }
            }

            if (output_flux_par) {
              heatf_par_flutter_i = -kappa_par_i * bracket1i;
            }
          }
        // }
      }

      if (diffusion_perp > 0.0) {
        ddt(Ti) += kappa_perp_i * Delp2(Ti) / N0; // Perpendicular diffusion
        ddt(Ti) += Grad_perp(kappa_perp_i) * Grad_perp(Ti) / N0; //NOTE(malamast): We have not included the time variation of kappa_perp_i when we substratced the equilibrium term. 
      }

      if (gyroviscous && compress0 & nonlinear) {
        ddt(Ti) -= 1.333333 * Tipara3 * (Ti0 + Ti) * Vipar * b0xcv * Grad(Vipar) / B0;
      }

      if (parallel_viscous) {
        if (compress0) {
	        ddt(Ti) -= 0.444444 * pi_ci / Tau_ie / N0 / sqrt(B0) * Grad_parP(Vipar * sqrt(B0)); //NOTE(malamast): What happended to kB?
	      }
	      if (nonlinear) {
	        ddt(Ti) += 0.222222 * pi_ci / Tau_ie / N0 / N0 / B0 * bracket(N0+Ni, Ti0+Ti, bm_exb); // NOTE(malamast): 
                                                                                                // 1) What happended to kB?  
                                                                                                // 2) Should we remove the background term here to make it consistent with the rest of the analysis?
                                                                                                // 3) Why do we devide by B0? I think it is already included in the bracket operator?
	      }
      }
      
      if (neoclassic_i) {
        tmpddx2 = D2DX2(DC(Ti));
        mesh->communicate(tmpddx2);
        tmpddx2.applyBoundary();
        ddt(Ti) += coord->g11 * DC((xii_neo * DC(tmpddx2)));

        heatf_neo_i = -xii_neo * sqrt(coord->g11) * DDX(Ti);
        mesh->communicate(heatf_neo_i);
        heatf_neo_i.applyBoundary();
      }

      // 4th order Parallel diffusion terms
      if (hyperdiff_par_ti4 > 0.0) {
        tmpTi2 = Grad2_par2new(Ti);
        mesh->communicate(tmpTi2);
        tmpTi2.applyBoundary();
        ddt(Ti) -= hyperdiff_par_ti4 * Grad2_par2new(tmpTi2);
      }

      if (hyperdiff_perp_ti4 > 0.0) {
        tmpTi2 = Delp2(Ti);
    	  mesh->communicate(tmpTi2);
        tmpTi2.applyBoundary();
        ddt(Ti) -= hyperdiff_perp_ti4 * Delp2(tmpTi2);
      }
    }

  #if DEBUG_6F>0
    output.write("I see you 5!\n");
  #endif
    //////////////////////////////////////////////////////////////////
    // electron temperature equation
    {
      TRACE("ddt(Te)");

      ddt(Te) = 0.0;

      if (TeAmp > 0) ddt(Te) += TeSource;

      ddt(Te) -= bracket(phi, Te0, bm_exb);

      if (continuity) {
        ddt(Te) -= 4.0 / 3.0 * Te0 / B0 * b0xcv * Grad(phi);

        // ddt(Te) -= 4.0 / 3.0 * Te / B0 * b0xcv * cross(Ve0, B0vec); //NOTE(malamast): What is this? I have commented it out.
        ddt(Te) -= 4.0 / 3.0 * Te / B0 * b0xcv * Grad(phi0); // Added by malamast.

	      ddt(Te) += 4.0 / 3.0 * Tepara2 * Te / Ne0 * b0xcv * Grad(Pe0) / B0;

        if (diamag)
          ddt(Te) += 4.0 / 3.0 * Tepara2 * Te0 / Ne0 * b0xcv * Grad(Pe) / B0;
        if (nonlinear) {
          ddt(Te) -= 4.0 / 3.0 * Te / B0 * b0xcv * Grad(phi);
          if (diamag)
            ddt(Te) += 4.0 / 3.0 * Tepara2 * Te / Ne0 * b0xcv * Grad(Pe) / B0;
	      }
      }

      if (energy_flux) {
        // ddt(Te) += 10.0/3.0 * Tepara2 * Te0/B0 * b0xcv*Grad(Te);
        ddt(Te) -= 10.0 / 3.0 * Tepara2 / B0 * V_dot_Grad(-Te0 * b0xcv, Te);
	      ddt(Te) += 10.0 / 3.0 * Tepara2 * Te / B0 * b0xcv * Grad(Te0);

	      if (thermal_force) {
          ddt(Te) += 0.71 * 2.0 / 3.0 * Tepara3 * Te0 * B0 / Ne0 * Grad_parP(Jpar / B0);

          if (BScurrent) {
            if (evolve_psi) {
              ddt(Te) -= 0.71 * 2.0 / 3.0 * Tepara3 * Te0 * B0 / Ne0 * bracket(Psi, Jpar_BS0 / B0, bm_mag) * B0;
            } else {
              ddt(Te) -= 0.71 * 2.0 / 3.0 * Tepara3 * Te0 * B0 / Ne0 * bracket(Apar, Jpar_BS0 / B0, bm_mag);
            }
            ddt(Te) += 0.71 * 2.0 / 3.0 * Tepara3 * Te * B0 / Ne0 * Grad_parP(Jpar_BS0 / B0);
          } else {
            if (evolve_psi) {
              ddt(Te) -= 0.71 * 2.0 / 3.0 * Tepara3 * Te0 * B0 / Ne0 * bracket(Psi, J0 / B0, bm_mag) * B0;
            } else {
              ddt(Te) -= 0.71 * 2.0 / 3.0 * Tepara3 * Te0 * B0 / Ne0 * bracket(Apar, J0 / B0, bm_mag);
            }
            ddt(Te) += 0.71 * 2.0 / 3.0 * Tepara3 * Te * B0 / Ne0 * Grad_parP(J0 / B0);
          }
        }

        if (nonlinear) {
          // ddt(Te) += 10.0 / 3.0 * Tepara2 * Te/B0 * b0xcv*Grad(Te);
          ddt(Te) -= 10.0 / 3.0 * Tepara2 / B0 * V_dot_Grad(-Te * b0xcv, Te);
          if (thermal_force)
            ddt(Te) += 0.71 * 2.0 / 3.0 * Tepara3 * Te * B0 / Ne0 * Grad_parP(Jpar / B0);
        }
      }

      if (diamag && diamag_phi0) {
        if (diamag_er)
          ddt(Te) -= V_dot_Grad(Ve0, Te);
        else
          ddt(Te) -= bracket(phi0, Te, bm_exb); // Equilibrium flow
      }

      if (nonlinear) {
        ddt(Te) -= bracket(phi, Te, bm_exb); // Advection
      }

      if (compress0) {
        // ddt(Te) -= Vepar * Grad_parP(Te0);
        // ddt(Te) -= Vpar_Grad_par(Vepar, Te0);

	      if (include_vipar) {
	        ddt(Te) -= Vepar * Grad_parP(Te0);

          if (include_vpar0) {
            ddt(Te) -= Vepar0 * Grad_parP(Te);
            if (evolve_psi)
              ddt(Te) += Vepar0 * bracket(Psi, Te0, bm_mag) * B0; 
            else
              ddt(Te) += Vepar0 * bracket(Apar, Te0, bm_mag); 
          }                  
        }

        if (continuity) {
          ddt(Te) -= 2.0 / 3.0 * Te0 * B0 * Grad_parP(Vepar / B0);

          if (include_vpar0) {
            ddt(Te) -= 2.0 / 3.0 * Te * B0 * Grad_parP(Vepar0 / B0);
            if (evolve_psi)
              ddt(Te) += 2.0 / 3.0 * Te0 * B0 * bracket(Psi, Vepar0 / B0, bm_mag) * B0; 
            else
              ddt(Te) += 2.0 / 3.0 * Te0 * B0 * bracket(Apar, Vepar0 / B0, bm_mag); 
          }          
        }

        if (nonlinear) {
          // ddt(Te) -= Vepar * Grad_par(Te);

          // ddt(Te) -= Vpar_Grad_par(Vepar, Te);
          // ddt(Te) += Vepar * bracket(Psi, Te0, bm_mag)*B0;   

          if (include_vipar) {
            ddt(Te) -= Vpar_Grad_par(Vepar, Te);
            if (evolve_psi)
              ddt(Te) += Vepar * bracket(Psi, Te, bm_mag) * B0; // Added by malamast
            else
              ddt(Te) += Vepar * bracket(Apar, Te, bm_mag); // Added by malamast
          }

          if (continuity) {
            ddt(Te) -= 2.0 / 3.0 * Te * B0 * Grad_parP(Vepar / B0);
          }
        }
      }

      if (energy_exch) {
        ddt(Te) -= 2.0 * Tbar * nu_e / (ratio_pe * AA) * (Te - Ti);
        if (BScurrent)
          ddt(Te) += 4.0 / 3.0 * Tepara4 * eta * Jpar_BS0 * Jpar / Ne0;
        else
          ddt(Te) += 4.0 / 3.0 * Tepara4 * eta * J0 * Jpar / Ne0;

        if (nonlinear)
          ddt(Te) += 2.0 / 3.0 * Tepara4 * eta * Jpar * Jpar / Ne0;
      }

      if (diffusion_par > 0.0) {
        // if (Landau) {
        //   if (diff_par_flutter)
        //     ddt(Te) -= Grad_parP(q_par_e) / Ne0;
        //   else
        //     ddt(Te) -= Grad_par(q_par_e) / Ne0;

        //   if (output_qparcompare) {
        //     q_par_fl = kappa_par_e * Grad_par(Te);
        //     q_par_landau = q_par_e;
        //   }
        // } else {

          ddt(Te) += kappa_par_e * Grad2_par2(Te) / Ne0; // Parallel diffusion 
          ddt(Te) += Grad_par(kappa_par_e) * Grad_par(Te) / Ne0;

          if (diff_par_flutter) {
            if (nonlinear) {
              if (evolve_psi)
                bracket1e = -bracket(Psi, Te + Te0, bm_mag) * B0;
              else
                bracket1e = -bracket(Apar, Te + Te0, bm_mag);
            } else {
              if (evolve_psi)
                bracket1e = -bracket(Psi, Te0, bm_mag) * B0;
              else
                bracket1e = -bracket(Apar, Te0, bm_mag);
            }
            mesh->communicate(bracket1e);
            bracket1e.applyBoundary();
            gradpar_te = Grad_par(Te); //NOTE(malamast): Should we add Te0 here?
            mesh->communicate(gradpar_te);
            gradpar_te.applyBoundary();

            ddt(Te) += Grad_par(kappa_par_e) * bracket1e / Ne0;
            ddt(Te) += kappa_par_e * Grad_par(bracket1e) / Ne0;

            if (nonlinear) {
              if (evolve_psi) {
                ddt(Te) -= bracket(Psi, kappa_par_e, bm_mag) * B0 * gradpar_te / Ne0;
                ddt(Te) -= kappa_par_e * bracket(Psi, gradpar_te, bm_mag) * B0 / Ne0;
                ddt(Te) -= bracket(Psi, kappa_par_e, bm_mag) * B0 * bracket1e / Ne0;
                ddt(Te) -= kappa_par_e * bracket(Psi, bracket1e, bm_mag) * B0 / Ne0;
              } else {
                ddt(Te) -= bracket(Apar, kappa_par_e, bm_mag) * gradpar_te / Ne0;
                ddt(Te) -= kappa_par_e * bracket(Apar, gradpar_te, bm_mag) / Ne0;
                ddt(Te) -= bracket(Apar, kappa_par_e, bm_mag) * bracket1e / Ne0;
                ddt(Te) -= kappa_par_e * bracket(Apar, bracket1e, bm_mag) / Ne0;
              }
            }

            if (output_flux_par) {
              heatf_par_flutter_e = -kappa_par_e * bracket1e;
            }
          }
        // }
      }

      if (diffusion_perp > 0.0) {
        ddt(Te) += kappa_perp_e * Delp2(Te) / Ne0; // Perpendicular diffusion
        ddt(Te) += Grad_perp(kappa_perp_e) * Grad_perp(Te) / Ne0;
      }

      if (neoclassic_e) {
        tmpddx2 = D2DX2(DC(Te));
        mesh->communicate(tmpddx2);
        tmpddx2.applyBoundary();
        ddt(Te) += coord->g11 * DC(xie_neo * tmpddx2);

        heatf_neo_e = -xie_neo * sqrt(coord->g11) * DDX(Te);
        mesh->communicate(heatf_neo_e);
        heatf_neo_e.applyBoundary();
      }

      if (fix_fraction_imp) {
        ddt(Te) += 2./3. * (Wrad - Wrad0) ;    //impurity radiation 
        // Here, we substracted Wrad0 to make it consistent with the gradiemt driven analysis. 
        // We solve for the perturbed part of Te and we have subtracted the equilibrium terms from the equation. 
      }

      if (hyperdiff_par_te4 > 0.0) {
        tmpTe2 = Grad2_par2new(Te);
        mesh->communicate(tmpTe2);
        tmpTe2.applyBoundary();
        ddt(Te) -= hyperdiff_par_te4 * Grad2_par2new(tmpTe2);
      }

      if (hyperdiff_perp_te4 > 0.0) {
        tmpTe2 = Delp2(Te);
        mesh->communicate(tmpTe2);
        tmpTe2.applyBoundary();
        ddt(Te) -= hyperdiff_perp_te4 * Delp2(tmpTe2);
      }

      // right edge sink terms
      if (sink_Ter > 0.0) {
        ddt(Te) -= sink_Ter * sink_tanhxr(Te0, Te, ste_widthr, ste_lengthr); // sol sink
      }

      if (output_flux_par) {
        // gamma_par_i = (N0 + Ni) * Vipar;
        heatf_par_i = -kappa_par_i * Grad_par(Ti);
        mesh->communicate(heatf_par_i);
        heatf_par_i.applyBoundary();
        heatf_par_e = -kappa_par_e * Grad_par(Te);
        mesh->communicate(heatf_par_e);
        heatf_par_e.applyBoundary();
      }

      if (output_vradial) {
        if (diamag_er)
          Vexb = Ve0 + (cross(B0vec, Grad(phi))) / (B0 * B0);
        else
          Vexb = (cross(B0vec, Grad(phi + phi0))) / (B0 * B0);
        mesh->communicate(Vexb);
        Vexb.applyBoundary();
        if (evolve_psi)
          Vbtilde = -cross(B0vec, Grad(Psi));
        else
          Vbtilde = -cross(B0vec, Grad(Apar)) / B0;
        mesh->communicate(Vbtilde);
        Vbtilde.applyBoundary();
        // Vbti_par = Vipar*Vbtilde.x;
        // Vbte_par = Vepar*Vbtilde.x;
        // mesh->communicate(Vbt_par);
      }
    }

  #if DEBUG_6F>0
    output.write("I see you 6!\n");
  #endif

    //////////////////////////////////////////////////////////////////
    if (compress0) { // parallel velocity equation
      TRACE("ddt(Vipar)");

      ddt(Vipar) = 0.0;

      ddt(Vipar) -= Vipara * Grad_parP(P) / N0;
      if (evolve_psi)
        ddt(Vipar) += Vipara * bracket(Psi, P0, bm_mag) * B0 / N0;
      else
	      ddt(Vipar) += Vipara * bracket(Apar, P0, bm_mag) / N0;

      if (diamag && diamag_phi0) {
	      if (diamag_er)
          ddt(Vipar) -= V_dot_Grad(Ve0, Vipar);
	      else
	        ddt(Vipar) -= bracket(phi0, Vipar, bm_exb);
      }

      if (include_vipar && include_vpar0) {

        ddt(Vipar) -= Vpar_Grad_par(Vipar0, Vipar); 
        if (evolve_psi)
          ddt(Vipar) += Vipar0 * bracket(Psi, Vipar, bm_mag) * B0; 
        else
          ddt(Vipar) += Vipar0 * bracket(Apar, Vipar, bm_mag);

        ddt(Vipar) -= Vpar_Grad_par(Vipar, Vipar0); 
        if (evolve_psi)
          ddt(Vipar) += Vipar * bracket(Psi, Vipar0, bm_mag) * B0; 
        else
          ddt(Vipar) += Vipar * bracket(Apar, Vipar0, bm_mag);        

        if (evolve_psi)
          ddt(Vipar) += Vipar0 * bracket(Psi, Vipar0, bm_mag) * B0; 
        else
          ddt(Vipar) += Vipar0 * bracket(Apar, Vipar0, bm_mag);        

      }

      if (nonlinear) {
        ddt(Vipar) -= bracket(phi, Vipar, bm_exb);

        // ddt(Vipar) -= Vipar * Grad_par(Vipar);
        // ddt(Vipar) -= Vpar_Grad_par(Vipar, Vipar);

        if (include_vipar) {
          ddt(Vipar) -= Vpar_Grad_par(Vipar, Vipar); 
          if (evolve_psi)
            ddt(Vipar) += Vipar * bracket(Psi, Vipar, bm_mag) * B0; //NOTE(malamast): Added by malamast
          else
            ddt(Vipar) += Vipar * bracket(Apar, Vipar, bm_mag); //NOTE(malamast): Added by malamast
        }
      }

      if (parallel_viscous && compress0) {
        // Field3D temp_pi;
        // temp_pi = pi_ci / (B0*sqrt(B0));
        ddt(Vipar) -= 0.666667 * Vipara * (B0 * sqrt(B0)) * Grad_parP(pi_ci / (B0 * sqrt(B0))) / N0;
      }
      
      if (gyroviscous) {
        ddt(Vipar) -= Upara0 * bracket(Pi0, Vipar, bm_exb) / N0; // NOTE(malamast): Should we use bm_exb or bm_mag here?
        if (include_vpar0)
          ddt(Vipar) -= Upara0 * bracket(Pi, Vipar0, bm_exb) / N0;
        if (nonlinear)
          ddt(Vipar) -= Upara0 * bracket(Pi, Vipar, bm_exb) / N0;
      }

      // parallel hyper-viscous diffusion for vector potential
      if (hyperdiff_par_v4 > 0.0) {
        tmpVp2 = Grad2_par2new(Vipar);
        mesh->communicate(tmpVp2);
        tmpVp2.applyBoundary();
        ddt(Vipar) -= hyperdiff_par_v4 * Grad2_par2new(tmpVp2);
      }

      if (hyperdiff_perp_v4 > 0.0) {
        tmpVp2 = Delp2(Vipar);
        mesh->communicate(tmpVp2);
        tmpVp2.applyBoundary();
        ddt(Vipar) -= hyperdiff_perp_v4 * Delp2(tmpVp2);
      }

      if (sink_vp > 0.0) {
        if (include_vpar0){
          ddt(Vipar) -= sink_vp * sink_tanhxl(Vipar0, Vipar, sp_width, sp_length); // sink
        } else {
          //Field2D V0tmp = 0.;
          ddt(Vipar) -= sink_vp * sink_tanhxl(F2D_tmp, Vipar, sp_width, sp_length); // sink
        }
      }
    }
  #if DEBUG_6F>0
    output.write("I see you 7!\n");
  #endif

    //////////////////////////////////////////////////////////////////
    if (neutral) { // neutral model
    
      N_tmp = field_larger(N0 + Ni, Low_limit);
      Te_tmp = field_larger(Te0 + Te, Low_limit);
      Ti_tmp = field_larger(Ti0 + Ti, Low_limit);
      if (impurity_prof)
        Ne_tmp = field_larger(Ne0 + Zi * Ni, Low_limit);
      else
        Ne_tmp = Zi * N_tmp;

      //**************************************************************
      // Molecule density Nm and Perpendicular Velocity in X ---Vmx---
      //**************************************************************
      if (with_fueling) {
        ddt(Nm) = 0.0;
	      ddt(Vm) = 0.0;
	
        //******---Nm---******
        Nm_tmp = Nm;
        nu_diss = 3.e4 * Nm_tmp * Nbar * Te_tmp * Te_tmp * Tebar * Tebar/(3.+0.01 * Te_tmp * Te_tmp * Tebar * Tebar);
        S_diss = N_tmp * nu_diss * Tbar;
              mesh->communicate(S_diss);
              S_diss.applyBoundary("neumann");
              
        ddt(Nm) -= V_dot_Grad(Vm, Nm)-Nm * Div(Vm) + S_diss;
              
        if (gas_puffing) {
                ddt(Nm) += Sgas;
        }

        //******---Vmx---******
        pm = Nm * Tm_x / Tebar;
        //Vm.x=Vmx;
        //Vm.y=0.;
        //Vm.z=0.;
        ddt(Vm) -= V_dot_Grad(Vm,Vm) + Grad(pm) / Nm_tmp / Mm;
      }
      
      //**************************************************************
      // Atom equations
      //**************************************************************
      Pn = Nn * Ti_tmp;
      BoutReal minimum_val = 1.e-10;
      Pn=field_larger(Pn,minimum_val);

      ddt(Nn) = 0.;
      ddt(Vn) = 0.;      
	
      // neutral density equations
      if (Solving_Eq_Nn) {
        ddt(Nn) += Dn * Delp2(Nn);
        term2 = Dn * Delp2(Nn);

        if (Gradperp_Dn) {
          ddt(Nn) += Grad_perp(Dn) * Grad_perp(Nn);
          term2 += Grad_perp(Dn) * Grad_perp(Nn);
        }
        mesh->communicate(term2);
        term2.applyBoundary("neumann");

        if (with_vipar) {
          Field3D N_eff = Nn * (nu_cx + nu_rc)/(nu_cx + nu_iz);
          ddt(Nn) -= Grad_par(Vipar) * N_eff + Vpar_Grad_par(Vipar, N_eff); //NOTE(malamast): Inconsistency: Vipar is the perturbed part. Not the total. 
          term1 = -Grad_par(Vipar) * N_eff - Vpar_Grad_par(Vipar, N_eff);
          mesh->communicate(term1);
          term1.applyBoundary("neumann");
        }
	  
        if (!Solving_Eq_Vn){
          ddt(Nn) += Dn * Grad2_par2(Nn);
          term3 = Dn * Grad2_par2(Nn);
          if (Gradpar_Dn) {
            ddt(Nn) += Grad_par(Dn) * Grad_par(Nn);
            term3 += Grad_par(Dn) * Grad_par(Nn);
          }
          mesh->communicate(term3);
          term3.applyBoundary("neumann");
        } else {
          ddt(Nn) -= Grad_par(Vn) * Nn + Vpar_Grad_par(Vn, Nn);
        }

        // source/sink terms
        if (external_source) {
                ddt(Nn) += Sn_ext;
        }
	  
        if (with_fueling) {
                ddt(Nn) += 2.0 * S_diss;
        }
        
        Sn = nu_rc * Ne_tmp * N_tmp - nu_iz * Nn * N_tmp;
        mesh->communicate(Sn);
        Sn.applyBoundary("neumann");
        ddt(Nn) += Sn;
        ddt(Ni) -= Sn;
      }

      // neutral parallel velcocity equations
      if (Solving_Eq_Vn) {
        ddt(Vn) = -Vn * Grad_par(Vn);
        ddt(Vn) -= Upara1 * Tau_ie * Grad_par(Pn);	  
        //ddt(Vn) += (Grad_par(etan) * Grad_par(Vn) + etan * Grad2_par2new(Vn)) / Nn;
        ddt(Vn) += etan * Grad2_par2(Vn);

        if (Gradpar_etan) {
          ddt(Vn) += (Grad_par(etan) * Grad_par(Vn));
        }

        if (Gradperp_etan) {
          ddt(Vn) += etan_perp * Delp2(Vn);
          ddt(Vn) += Grad_perp(etan_perp) * Grad_perp(Vn);
	      }
	
      	if (compress0) {
          Field3D Nn_tmp = Nn;
          Nn_tmp=field_larger(Nn_tmp, minimum_val);
          Sv = N_tmp * (Ne_tmp / Nn_tmp * nu_rc + nu_cx) * (Vipar - Vn);
          mesh->communicate(Sv);
          Sv.applyBoundary("neumann");
          ddt(Vn) += Sv;
	      }
      }

      // NOTE(malamast): The below are all wrong/inconsistent. Te, Ti, U, Vipar are all the perturbed parts. You cannot mix up terms.  
      ddt(Te) -= nu_iz * Nn * (Te0 + Te + 2./3.* Wiz /Tebar); // ionization effect
      ddt(Te) += nu_rc * Nn * Wrc / Tebar; // recombanation effect
      ddt(Ti) -= (nu_iz-nu_rc) * Nn * (Ti0 + Ti);
      ddt(U) -= (nu_iz+nu_rc) * Nn * U;
      ddt(Vipar) -= Nn * (nu_iz + nu_cx) * (Vipar - Vn);
    }
    
  #if DEBUG_6F>0
    output.write("I see you 8!\n");
  #endif
    //////////////////////////////////////////////////////////////////
    
    if (PF_limit) {
      phi = PF_filter(phi, PF_limit_range);
      Jpar = PF_filter(Jpar, PF_limit_range);
      if (evolve_psi)
        Psi = PF_filter(Psi, PF_limit_range);
      else
        Apar = PF_filter(Apar, PF_limit_range);
      Ni = PF_filter(Ni, PF_limit_range);
      Ti = PF_filter(Ti, PF_limit_range);
      Te = PF_filter(Te, PF_limit_range);
      U = PF_filter(U, PF_limit_range);
      P = PF_filter(P, PF_limit_range);
      if (compress0) {
        Vipar = PF_filter(Vipar, PF_limit_range);
      }

      // phi = PF_filter_2(phi, PF_limit_range);
      // Jpar = PF_filter_2(Jpar, PF_limit_range);
      // if (evolve_psi)
      //   Psi = PF_filter_2(Psi, PF_limit_range);
      // else
      //   Apar = PF_filter_2(Apar, PF_limit_range);
      // Ni = PF_filter_2(Ni, PF_limit_range);
      // Ti = PF_filter_2(Ti, PF_limit_range);
      // Te = PF_filter_2(Te, PF_limit_range);
      // U = PF_filter_2(U, PF_limit_range);
      // P = PF_filter_2(P, PF_limit_range);
      // if (compress0) {
      //   Vipar = PF_filter_2(Vipar, PF_limit_range);
      // }      

    }
    
    if (filter_z) {
      // Filter out all except filter_z_mode
      TRACE("filter_z");

      if (!emass) {
	      if (evolve_psi)
          ddt(Psi) = filter(ddt(Psi), filter_z_mode);
        else
	        ddt(Apar) = filter(ddt(Apar), filter_z_mode);
      } else {
        ddt(Ajpar) = filter(ddt(Ajpar), filter_z_mode);
      }

      ddt(U) = filter(ddt(U), filter_z_mode);

      ddt(Ni) = filter(ddt(Ni), filter_z_mode);

      ddt(Ti) = filter(ddt(Ti), filter_z_mode);

      ddt(Te) = filter(ddt(Te), filter_z_mode);

      if (compress0) {
        ddt(Vipar) = filter(ddt(Vipar), filter_z_mode);
      }
    } else if (filter_z_nonlinear) {
      if (!emass) {
        if (evolve_psi)
	        ddt(Psi) = filter_z_non(ddt(Psi), 0, filter_z_mode);
	      else 
	        ddt(Apar) = filter_z_non(ddt(Apar), 0, filter_z_mode);
      } else {
        ddt(Ajpar) = filter_z_non(ddt(Ajpar), 0, filter_z_mode);
      }

      ddt(U) = filter_z_non(ddt(U), 0, filter_z_mode);
  
      ddt(Ni) = filter_z_non(ddt(Ni), 0, filter_z_mode);
      
      ddt(Ti) = filter_z_non(ddt(Ti), 0, filter_z_mode);

      ddt(Te) = filter_z_non(ddt(Te), 0, filter_z_mode);

      if (compress0) {
        ddt(Vipar) = filter_z_non(ddt(Vipar), 0, filter_z_mode);
      }
    }

    if (PF_sink > 0.) {
      if (evolve_psi)
        ddt(Psi) -= PF_sink * sink_PF(F2D_tmp, Psi, PFs_width, PFs_length);
      else
        ddt(Apar) -= PF_sink * sink_PF(F2D_tmp, Apar, PFs_width, PFs_length);
     
      ddt(Ni) -= PF_sink * sink_PF(N0, Ni, PFs_width, PFs_length);
     
      ddt(Ti) -= PF_sink * sink_PF(Ti0, Ti, PFs_width, PFs_length);
   
      ddt(Te) -= PF_sink * sink_PF(Te0, Te, PFs_width, PFs_length);
    
      ddt(U) -= PF_sink * sink_PF(F2D_tmp, U, PFs_width, PFs_length);
    
      //  sink_PFtmp = sink_PF(Te0,Te,PFs_width,PFs_length);
    }

    //////////////////////////////////////////////////////////////////

    if (low_pass_z > 0) {
      // Low-pass filter, keeping n up to low_pass_z
      TRACE("low_pass_z");

      if (!emass) {
        if (evolve_psi)
     	    ddt(Psi) = lowPass(ddt(Psi), low_pass_z, zonal_field);
	      else 
	        ddt(Apar) = lowPass(ddt(Apar), low_pass_z, zonal_field);
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

      if (neutral) {
        if (Solving_Eq_Nn){
          ddt(Nn) = lowPass(ddt(Nn), low_pass_z, zonal_bkgd);
        }
        if (Solving_Eq_Vn){
          ddt(Vn) = lowPass(ddt(Vn), low_pass_z, zonal_bkgd);
        }
      }

      if (pos_filter) {
        Ti = lowPass_pos2(Ti, Ti);
        Te = lowPass_pos2(Te, Te);
        Ni = lowPass_pos2(Ni, Ni);
      }

      if (pos_filter2) {
        ddt(Ti) = lowPass_pos(ddt(Ti), filter_position_ti);
        ddt(Te) = lowPass_pos(ddt(Te), filter_position_te);
        ddt(Ni) = lowPass_pos(ddt(Ni), filter_position_ni);
      }

      if (pos_filter_zf) {
        ddt(Ti) = sink_zonal_core(ddt(Ti), filter_position_ti);
        ddt(Te) = sink_zonal_core(ddt(Te), filter_position_te);
        ddt(Ni) = sink_zonal_core(ddt(Ni), filter_position_ni);
      }

      if (pos_sink_zf > 0.) {
        ddt(Ti) -= pos_sink_zf * DC(sink_tanhxl(Ti0, Ti, filter_position_ti * pos_filter_width / int(Grid_NX), filter_position_ti *pos_filter_length / int(Grid_NX)));
        ddt(Te) -= pos_sink_zf * DC(sink_tanhxl(Te0, Te, filter_position_te * pos_filter_width / int(Grid_NX), filter_position_te *pos_filter_length / int(Grid_NX)));
        ddt(Ni) -= pos_sink_zf * DC(sink_tanhxl(N0, Ni, filter_position_ni * pos_filter_width / int(Grid_NX), filter_position_ni *pos_filter_length / int(Grid_NX)));
      }

      if (nonlinear && pos_filter) {
        Pi = lowPass_pos2(Pi, Pi);
        Pe = lowPass_pos2(Pe, Pe);
        P = lowPass_pos2(P, P);
      }

      if (nonlinear && pos_filter2) {
        Pi = lowPass_pos(Pi, position_tmpi);
        Pe = lowPass_pos(Pe, position_tmpe);
        P = lowPass_pos(P, position_tmp);
      }

      if (nonlinear && (pos_filter2 || pos_filter_zf)) {
        Pi = sink_zonal_core(Pi, position_tmpi);
        Pe = sink_zonal_core(Pe, position_tmpe);
        P = sink_zonal_core(P, position_tmp);
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

  //****************BOUNDARY FUNCTIONS******************************************
  // Sheath Boundary Conditions on Phi
  // Linearized
  void SBC_Dirichlet(Field3D &var, const Field3D &value, bool PF_limit, BoutReal PF_limit_range) {
    // let the boundary equall to the value next to the boundary
    SBC_yup_eq(var, value, PF_limit, PF_limit_range);
    SBC_ydown_eq(var, -value, PF_limit, PF_limit_range);
  } 

  void SBC_Gradpar(Field3D &var, const Field3D &value, bool PF_limit, BoutReal PF_limit_range) {
    SBC_yup_Grad_par(var, value, PF_limit, PF_limit_range);
    SBC_ydown_Grad_par(var, -value, PF_limit, PF_limit_range);
  }

  // Boundary to specified Field3D object
  void SBC_yup_eq(Field3D &var, const Field3D &value, bool PF_limit, BoutReal PF_limit_range) {
    RangeIterator xrup = mesh->iterateBndryUpperY();
    
    // for(xrup->first(); !xrup->isDone(); xrup->next())
    for (; !xrup.isDone(); xrup++) {
      xind = xrup.ind;
      indx = mesh->getGlobalXIndex(xind);
      if (PF_limit && (BoutReal(indx) > ixsep * PF_limit_range)) {
        for (int jy = mesh->yend + 1 - Sheath_width; jy < mesh->LocalNy; jy++)
          for (int jz = 0; jz < mesh->LocalNz; jz++) {
            var(xind,jy,jz) = value(xind,jy,jz);
	  }
      } else if (PF_limit && (BoutReal(indx) <= ixsep * PF_limit_range)) {
        for (int jy = mesh->yend + 1 - Sheath_width; jy < mesh->LocalNy; jy++)
          for (int jz = 0; jz < mesh->LocalNz; jz++) {
            var(xind,jy,jz) = 0.;
          }
      } else if (!PF_limit) {
        for (int jy = mesh->yend + 1 - Sheath_width; jy < mesh->LocalNy; jy++) {
          for (int jz = 0; jz < mesh->LocalNz; jz++) {
            var(xind,jy,jz) = value(xind,jy,jz);
	  }
	}
      }
    }
  }

  void SBC_ydown_eq(Field3D &var, const Field3D &value, bool PF_limit, BoutReal PF_limit_range) {
    RangeIterator xrdn = mesh->iterateBndryLowerY();

    // for(xrdn->first(); !xrdn->isDone(); xrdn->next())
    for (; !xrdn.isDone(); xrdn++) {
      xind = xrdn.ind;
      indx = mesh->getGlobalXIndex(xind);
      if (PF_limit && (BoutReal(indx) > ixsep * PF_limit_range)) {
        for (int jy = mesh->ystart - 1 + Sheath_width; jy >= 0; jy--)
          for (int jz = 0; jz < mesh->LocalNz; jz++) {
            var(xind,jy,jz) = value(xind,jy,jz);
          }
      } else if (PF_limit && (BoutReal(indx) <= ixsep * PF_limit_range)) {
        for (int jy = mesh->ystart - 1 + Sheath_width; jy >= 0; jy--)
          for (int jz = 0; jz < mesh->LocalNz; jz++) {
            var(xind,jy,jz) = 0.;
          }
      } else if (!PF_limit) {
        for (int jy = mesh->ystart - 1 + Sheath_width; jy >= 0; jy--) {
          for (int jz = 0; jz < mesh->LocalNz; jz++) {
            var(xind,jy,jz) = value(xind,jy,jz);
          }
	}
      }
    }
  }

  // Boundary gradient to specified Field3D object
  void SBC_yup_Grad_par(Field3D &var, const Field3D &value, bool PF_limit, BoutReal PF_limit_range) {
    RangeIterator xrup = mesh->iterateBndryUpperY();

    Coordinates* coord = mesh->getCoordinates();
    
    // for(xrup->first(); !xrup->isDone(); xrup->next())
    for (; !xrup.isDone(); xrup++) {
      xind = xrup.ind;
      indx = mesh->getGlobalXIndex(xind);
      if (PF_limit && (BoutReal(indx) > ixsep * PF_limit_range)) {
        for (int jy = mesh->yend + 1 - Sheath_width; jy < mesh->LocalNy; jy++)
          for (int jz = 0; jz < mesh->LocalNz; jz++) {
            var(xind,jy,jz) = var(xind,jy - 1,jz) + coord->dy(xind,jy) * sqrt(coord->g_22(xind,jy)) * value(xind,jy,jz);
          }
      } else if (PF_limit && (BoutReal(indx) <= ixsep * PF_limit_range)) {
        for (int jy = mesh->yend + 1 - Sheath_width; jy < mesh->LocalNy; jy++)
          for (int jz = 0; jz < mesh->LocalNz; jz++) {
            var(xind,jy,jz) = var(xind,jy - 1,jz);
	  }
      } else if (!PF_limit) {
        for (int jy = mesh->yend + 1 - Sheath_width; jy < mesh->LocalNy; jy++)
          for (int jz = 0; jz < mesh->LocalNz; jz++) {
            var(xind,jy,jz) = var(xind,jy - 1,jz) + coord->dy(xind,jy) * sqrt(coord->g_22(xind,jy)) * value(xind,jy,jz);
          }
      }
    }
  }

  void SBC_ydown_Grad_par(Field3D &var, const Field3D &value, bool PF_limit, BoutReal PF_limit_range) {
    RangeIterator xrdn = mesh->iterateBndryLowerY();
    Coordinates* coord = mesh->getCoordinates();

    for (; !xrdn.isDone(); xrdn++) {
      xind = xrdn.ind;
      indx = mesh->getGlobalXIndex(xind);
        if (PF_limit && (BoutReal(indx) > ixsep * PF_limit_range)) {
          for (int jy = mesh->ystart - 1 + Sheath_width; jy >= 0; jy--)
            for (int jz = 0; jz < mesh->LocalNz; jz++) {
              var(xind,jy,jz) = var(xind,jy + 1,jz) + coord->dy(xind,jy) * sqrt(coord->g_22(xind,jy)) * value(xind,jy,jz);
	    }
      }
      if (PF_limit && (BoutReal(indx) <= ixsep * PF_limit_range)) {
        for (int jy = mesh->ystart - 1 + Sheath_width; jy >= 0; jy--)
          for (int jz = 0; jz < mesh->LocalNz; jz++) {
            var(xind,jy,jz) = var(xind,jy + 1,jz);
	  }
      } else if (!PF_limit) {
        for (int jy = mesh->ystart - 1 + Sheath_width; jy >= 0; jy--)
          for (int jz = 0; jz < mesh->LocalNz; jz++) {
            var(xind,jy,jz) = var(xind,jy + 1,jz) + coord->dy(xind,jy) * sqrt(coord->g_22(xind,jy)) * value(xind,jy,jz);
	  }
      }
    }
  }

  /*****************************************************************************
   * Preconditioner
   *
   * o System state in variables (as in rhs function)
   * o Values to be inverted in time derivatives
   *
   * o Return values should be in time derivatives
   *
   * enable by setting solver / use_precon = true in BOUT.inp
   *****************************************************************************/
  
  int precon(BoutReal UNUSED(t), BoutReal gamma, BoutReal UNUSED(delta)) {
    ni_tmp = ddt(Ni);
    ti_tmp = ddt(Ti);
    te_tmp = ddt(Te);
    vi_tmp = ddt(Vipar);
    if (evolve_psi)
      psi_tmp = ddt(Psi);
    else 
      psi_tmp = ddt(Apar);
    u_tmp = ddt(U);
    ni_tmp.applyBoundary("neumann");
    ti_tmp.applyBoundary("neumann");
    te_tmp.applyBoundary("neumann");
    vi_tmp.applyBoundary("neumann");
    psi_tmp.applyBoundary("neumann");
    u_tmp.applyBoundary("neumann");
    
    // First matrix, applying L
    mesh->communicate(ni_tmp, ti_tmp, te_tmp, psi_tmp);
    if (evolve_psi)
      jpar1 = -B0 * Delp2(psi_tmp);
    else
      jpar1 = - Delp2(psi_tmp);
    mesh->communicate(vi_tmp, u_tmp, jpar1);
    jpar1.applyBoundary();

    if (smooth_j_x) {
      jpar1 = smooth_x(jpar1);
      jpar1 = smooth_x(jpar1);
      jpar1 = smooth_x(jpar1);
      jpar1 = smooth_x(jpar1);
      jpar1 = smooth_x(jpar1);
      jpar1 = smooth_y(jpar1);
      jpar1 = smooth_y(jpar1);
      jpar1 = smooth_y(jpar1);
      jpar1 = smooth_y(jpar1);
      jpar1 = smooth_y(jpar1);
      mesh->communicate(jpar1);
      jpar1.applyBoundary();
    }

    if (jpar_bndry_width > 0) {
      // Boundary in jpar
      if (mesh->firstX()) {
        for (int i = jpar_bndry_width; i >= 0; i--) {
          for (int j = 0; j < mesh->LocalNy; j++) {
            for (int k = 0; k < mesh->LocalNz; k++) {
              jpar1(i, j, k) = 0.5 * jpar1(i + 1, j, k);
            }
          }
        }
      }
      if (mesh->lastX()) {
        for (int i = mesh->LocalNx - jpar_bndry_width - 1; i < mesh->LocalNx; i++) {
          for (int j = 0; j < mesh->LocalNy; j++) {
            for (int k = 0; k < mesh->LocalNz; k++) {
              jpar1(i, j, k) = 0.5 * jpar1(i - 1, j, k);
            }
          }
        }
      }
    }

    p_tmp = N0 * (Tau_ie * ti_tmp + te_tmp) + ni_tmp * (Tau_ie * Ti0 + Te0);
    u_tmp1 = u_tmp + gamma * ((B0 * B0) * Grad_par(jpar1 / B0, CELL_CENTRE) + 2.0 * Upara1 * b0xcv * Grad(p_tmp)); // curvature term
    mesh->communicate(u_tmp1);
    u_tmp1.applyBoundary();
    Vipar = vi_tmp;
    mesh->communicate(Vipar);
    Vipar.applyBoundary();

    // second metrix
    if (!diffusion_par) {
      kappa_par_i_lin = 0.;
      kappa_par_e_lin = 0.;
    }

    //static InvertPar *invi = 0;
    static std::unique_ptr<InvertPar> invi{nullptr};
    if (!invi) {
      invi = InvertPar::create(); // Create parallel solver
      invi->setCoefA(1.0);
    }
    invi->setCoefB(-gamma * kappa_par_i_lin);
    ti_tmp2 = invi->solve(ti_tmp); // Solve Pshur
    mesh->communicate(ti_tmp2);
    ti_tmp2.applyBoundary();

    //static InvertPar *inve = 0;
    static std::unique_ptr<InvertPar> inve{nullptr};
    if (!inve) {
      inve = InvertPar::create(); // Create parallel solver
      inve->setCoefA(1.0);
    }
    inve->setCoefB(-gamma * kappa_par_e_lin);
    te_tmp2 = inve->solve(te_tmp); // Solve Pshur
    mesh->communicate(te_tmp2);
    te_tmp2.applyBoundary();

    //static InvertPar *invu = 0;
    static std::unique_ptr<InvertPar> invu{nullptr};
    if (!invu) {
      invu = InvertPar::create(); // Create parallel solver
      Field2D rb_tmp = Grad_par(Rxy * Bpxy);
      mesh->communicate(rb_tmp);
      rb_tmp.applyBoundary("dirichlet");
      invu->setCoefA(1.0 + 2. * gamma * gamma * Grad_par(rb_tmp) / (Rxy * Bpxy * SQ(B0)));
    }
    invu->setCoefB(-SQ(gamma * B0));
    U = invu->solve(u_tmp1); // Solve Pshur
    mesh->communicate(U);
    U.applyBoundary();

    // third metrix
    BoutReal Ntemp = max(N0, true);
    // TODO: diamag in U?
    phi_tmp = phiSolver->solve(U * B0 / Ntemp);
    mesh->communicate(phi_tmp);
    phi_tmp.applyBoundary();
    Ni = ni_tmp - gamma * bracket(phi_tmp, N0, bm_exb);
    Ti = ti_tmp2 - gamma * bracket(phi_tmp, Ti0, bm_exb);
    Te = te_tmp2 - gamma * bracket(phi_tmp, Te0, bm_exb);
    if (evolve_psi)
      Psi = psi_tmp - gamma * Grad_par(phi_tmp) / B0;
    else
      Apar = psi_tmp - gamma * Grad_par(phi_tmp);
    mesh->communicate(Ni, Ti, Te);
    Ni.applyBoundary();
    Ti.applyBoundary();
    Te.applyBoundary();
    
    if (evolve_psi) {
      mesh->communicate(Psi);
      Psi.applyBoundary();
    } else {
      mesh->communicate(Apar);
      Apar.applyBoundary();
    }
    
    return 0;
  }

  /*****************************************************************************
   * Preconditioner for when phi solved as a constraint
   * Currently only possible with the IDA solver
   *
   * o System state in variables (as in rhs function)
   * o Values to be inverted in F_vars
   *
   * o Return values should be in vars (overwriting system state)
   *****************************************************************************/

  int precon_phi(BoutReal UNUSED(t), BoutReal UNUSED(cj), BoutReal UNUSED(delta)) {
    ddt(phi) = phiSolver->solve(C_phi - ddt(U));
    return 0;
  }
};

BOUTMAIN(Elm_6f)