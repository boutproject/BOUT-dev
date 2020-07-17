/****************************************************************
 * GEM Gyro-fluid model
 *
 * 6 moments for each species
 *
 * "GEM - An Energy Conserving Electromagnetic Gyrofluid Model"
 *  by Bruce D Scott. arXiv:physics/0501124v1 23 Jan 2005 
 *
 * This version uses global parameters for collisionality etc.
 ****************************************************************/

#include <bout/physicsmodel.hxx>
#include <bout/constants.hxx>

#include <gyro_average.hxx>
#include <interpolation.hxx>
#include <invert_laplace.hxx>

/// Fundamental constants

const BoutReal e0 = 8.854e-12; // Permittivity of free space
const BoutReal qe = 1.602e-19; // Electron charge
const BoutReal Me = 9.109e-31; // Electron mass
const BoutReal Mp = 1.67262158e-27; // Proton mass

class GEM : public PhysicsModel {
  //////////////////////////////////////////
  // Evolving quantities
  
  // Ion species
  Field3D Ni;   // Gyro-center density
  Field3D ApUi; // beta_e*Apar + mu_i * Ui
  Field3D Tipar, Tiperp; // Parallel or perpendicular temp
  Field3D qipar, qiperp; // Parallel and perpendicular heat flux
  
  // Electron species
  Field3D Ne;
  Field3D ApUe; // beta_e*Apar + mu_e * Ue
  Field3D Tepar, Teperp;
  Field3D qepar, qeperp;

  //////////////////////////////////////////
  // Derived quantities
  
  Field3D phi;    // Electrostatic potential
  Field3D Apar;   // Parallel vector potential
  
  Field3D Ui, Ue; // Ion and electron parallel velocity
  
  Field3D Jpar;   // Parallel current
  
  Field3D phi_G, Phi_G; // Gyro-reduced potential
  
  //////////////////////////////////////////
  // Equilibrium

  Coordinates *coord;
  
  Vector3D B0vec; // Equilibrium B field vector
  Field2D logB;   // For curvature
  Coordinates::FieldMetric Grad_par_logB; // Grad_par(log(B))

  Field2D Ni0, Ne0; // Gyro-center densities
  Field2D Ti0, Te0; // Starting isotropic temperatures

  BoutReal tau_e, tau_i; // T_z / (Te * Z)
  BoutReal mu_e, mu_i;   // M_z / (M_i * Z)

  BoutReal beta_e; // Electron dynamical beta

  BoutReal rho_e, rho_i; // Electron, ion gyroradius
  BoutReal rho_s, delta;

  // Collisional transport coefficients
  const BoutReal eta     = 0.51;

  const BoutReal alpha_e = 0.71;
  const BoutReal kappa_e = 3.2;
  const BoutReal pi_e    = 0.73;

  // alpha_i = 0
  const BoutReal kappa_i = 3.9;
  const BoutReal pi_i    = 0.73;

  Field3D Rei;

  //////////////////////////////////////////
  // Options

  bool adiabatic_electrons;  // Solve adiabatic electrons
  bool small_rho_e;          // Neglect electron gyro-radius
  bool include_grad_par_B;   // Include terms like Grad_par(log(B))

  bool curv_logB;

  BoutReal Landau; // Multiplier for Landau damping terms

  BoutReal nu_e, nu_i; // Collisional dissipation

  BoutReal nu_perp, nu_par;  // Artificial dissipation

  bool fix_profiles; // Subtract toroidal averages

  int  jpar_bndry_width; // Set jpar = 0 in a boundary region

  bool nonlinear; // Include nonlinear terms

  BoutReal flat_temp; // Flat temperature profile
  BoutReal flat_dens; // Flat density profile

  // Method to use for brackets: BRACKET_ARAKAWA, BRACKET_STD or BRACKET_SIMPLE
  const BRACKET_METHOD bm = BRACKET_SIMPLE;

  int low_pass_z; // Toroidal (Z) filtering of all variables

  //////////////////////////////////////////
  // Terms in the equations

  bool ne_ddt, ne_ne1, ne_te0, ne_te1, ne_ue, ne_curv;
  bool apue_ddt, apue_phi1, apue_pet, apue_curv, apue_gradB, apue_Rei; // Linear terms
  bool apue_ue1_phi1, apue_qe1_phi1, apue_apar1_phi1, apue_apar1_pe1;  // Nonlinear terms
  bool tepar_ddt;
  bool teperp_ddt;
  bool qepar_ddt;
  bool qeperp_ddt;

  bool ni_ddt, ni_ni1, ni_ti0, ni_ti1, ni_ui, ni_curv;
  bool apui_ddt, apui_phi1, apui_pit, apui_curv, apui_gradB, apui_Rei; // Linear terms
  bool apui_ui1_phi1, apui_qi1_phi1, apui_apar1_phi1, apui_apar1_pi1;  // Nonlinear terms
  bool tipar_ddt;
  bool tiperp_ddt;
  bool qipar_ddt;
  bool qiperp_ddt;

  //////////////////////////////////////////
  // Normalisation factors

  BoutReal Lbar;   // Perpendicular scale length
  BoutReal Tenorm; // Typical value of Te for normalisation
  BoutReal Ninorm; // Typical density value for normalisation
  BoutReal Bbar;   // Magnetic field  

  BoutReal Cs;    // Sound speed sqrt(Tenorm / Mi)
  BoutReal Tbar;  // Timescale Lbar / Cs

  FieldGroup comms; // Communications

  /// Solver for inverting Laplacian
  std::unique_ptr<Laplacian> phiSolver{nullptr};
  std::unique_ptr<Laplacian> aparSolver{nullptr};
  
  ////////////////////////////////////////////////////////////////////////
  // Initialisation
  
  //Field3D xzdamp;
  
  int init(bool restarting) override {
    //////////////////////////////////
    // Read options

    auto globalOptions = Options::root();
    auto options = globalOptions["gem"];

    adiabatic_electrons = options["adiabatic_electrons"].withDefault(false);
    small_rho_e = options["small_rho_e"].withDefault(true);
    include_grad_par_B = options["include_grad_par_B"].withDefault(true);

    Landau = options["Landau"].withDefault(1.0);

    nu_perp =
        options["nu_perp"].withDefault(0.01);     // Artificial perpendicular dissipation
    nu_par = options["nu_par"].withDefault(3e-3); // Artificial parallel dissipation

    low_pass_z = options["low_pass_z"].withDefault(-1); // Default is no filtering

    curv_logB =
        options["curv_logB"].withDefault(false); // Read in a separate logB variable

    fix_profiles = options["fix_profiles"].withDefault(false); // Subtract DC components

    nonlinear = options["nonlinear"].withDefault(true);
    flat_temp = options["flat_temp"].withDefault(-1.0);
    flat_dens = options["flat_dens"].withDefault(-1.0);

    //////////////////////////////////
    // Read profiles
    
    // Mesh
    Field2D Rxy, Bpxy, Btxy, Bxy, hthe;
    GRID_LOAD(Rxy);    // Major radius [m]
    GRID_LOAD(Bpxy);   // Poloidal B field [T]
    GRID_LOAD(Btxy);   // Toroidal B field [T]
    GRID_LOAD(Bxy);    // Total B field [T]
    GRID_LOAD(hthe);   // Poloidal arc length [m / radian]
    
    GRID_LOAD(Te0); // Electron temperature in eV
    GRID_LOAD(Ni0); // Ion number density in 10^20 m^-3
    
    Ni0 *= 1.e20; // Convert to m^-3
    
    Ti0 = Te0;
    Ne0 = Ni0;
    
    Field2D p_e = qe * Te0 * Ne0; // Electron pressure in Pascals
    
    if (flat_temp > 0.) {
      Te0 = flat_temp;
      Ti0 = flat_temp;
      
      Ne0 = Ni0 = p_e / (qe*Te0);
    } else if (flat_dens > 0.) {
      Ne0 = flat_dens;
      Ni0 = flat_dens;
      
      Te0 = Ti0 = p_e / (qe*Ne0);
    }
    
    if (curv_logB) {
      GRID_LOAD(logB);
    }

    // Set location of staggered fields
    ApUe.setLocation(CELL_YLOW);
    ApUi.setLocation(CELL_YLOW);
    qepar.setLocation(CELL_YLOW);
    qipar.setLocation(CELL_YLOW);
    
    //////////////////////////////////
    // Pick normalisation factors
    
    if (mesh->get(Lbar, "Lbar")) {// Try to read from grid file
      if (mesh->get(Lbar, "rmag")) {
        Lbar = 1.0;
      }
    }
    Lbar = options["Lbar"].withDefault(Lbar); // Override in options file
    SAVE_ONCE(Lbar);   // Save in output file

    BoutReal AA; // Ion atomic mass
    BoutReal ZZ; // Ion charge
    AA = options["AA"].withDefault(2.0); // Deuterium by default
    ZZ = options["ZZ"].withDefault(1.0);

    Tenorm = max(Te0,true); SAVE_ONCE(Tenorm); // Maximum value over the grid
    Ninorm = max(Ni0, true); SAVE_ONCE(Ninorm);
    
    Cs = sqrt(qe*Tenorm / (AA*Mp)); SAVE_ONCE(Cs); // Sound speed in m/s
    
    Tbar = Lbar / Cs;
    Tbar = options["Tbar"].withDefault(Tbar); // Override in options file
    SAVE_ONCE(Tbar); // Timescale in seconds
    
    if (mesh->get(Bbar, "Bbar")) {
      if (mesh->get(Bbar, "bmag")) {
        Bbar = max(Bxy, true);
      }
    }
    Bbar = options["Bbar"].withDefault(Bbar); // Override in options file
    SAVE_ONCE(Bbar);
    
    beta_e =  4.e-7*PI * max(p_e,true) / (Bbar*Bbar); SAVE_ONCE(beta_e); 
    
    // Mass to charge ratios
    mu_i = 1. / ZZ;
    mu_e = -1. / (AA * 1860.);
    
    tau_e = -1;
    tau_i = 1. /ZZ;
    
    // Gyro-radii (SI units)
    rho_s = Cs * AA * Mp / (qe * Bbar); 
    rho_e = rho_s * sqrt(fabs(mu_e * tau_e));
    rho_i = rho_s * sqrt(fabs(mu_i * tau_i));
    SAVE_ONCE3(rho_s, rho_e, rho_i);
    
    delta = rho_s / Lbar; SAVE_ONCE(delta); // This should be small
    
    ////////////////////////////////////////////////////
    // Terms in equations

    jpar_bndry_width = options["jpar_bndry_width"].withDefault(-1);

    OPTION4(options, ne_ddt, ne_te0, ne_ue, ne_curv, true); // Linear
    OPTION2(options, ne_ne1, ne_te1, nonlinear); // Nonlinear
    OPTION6(options, apue_ddt, apue_phi1, apue_pet, apue_curv, apue_gradB, apue_Rei, true);
    OPTION4(options, apue_ue1_phi1, apue_qe1_phi1, apue_apar1_phi1, apue_apar1_pe1, nonlinear);

    tepar_ddt = options["tepar_ddt"].withDefault(true);
    teperp_ddt = options["teperp_ddt"].withDefault(true);
    qepar_ddt = options["qepar_ddt"].withDefault(true);
    qeperp_ddt = options["qeperp_ddt"].withDefault(true);

    OPTION4(options, ni_ddt, ni_ti0, ni_ui, ni_curv, true); // Linear
    OPTION2(options, ni_ni1, ni_ti1, nonlinear); // Nonlinear
    OPTION6(options, apui_ddt, apui_phi1, apui_pit, apui_curv, apui_gradB, apui_Rei, true);
    OPTION4(options, apui_ui1_phi1, apui_qi1_phi1, apui_apar1_phi1, apui_apar1_pi1, nonlinear);
    tipar_ddt = options["tipar_ddt"].withDefault(true);
    tiperp_ddt = options["tiperp_ddt"].withDefault(true);
    qipar_ddt = options["qipar_ddt"].withDefault(true);
    qiperp_ddt = options["qiperp_ddt"].withDefault(true);

    ////////////////////////////////////////////////////
    // Collisional parameters
    
    /// Coulomb logarithm
    BoutReal Coulomb = 6.6 - 0.5*log(Ninorm * 1e-20) + 1.5*log(Tenorm);
    
    BoutReal t_e, t_i; // Braginskii collision times
    
    t_e = 1. / (2.91e-6 * (Ninorm / 1e6) * Coulomb * pow(Tenorm, -3./2));
    
    t_i = pow(ZZ, -4.) * sqrt(AA) / (4.80e-8 * (Ninorm / 1e6) * Coulomb * pow(Tenorm, -3./2));
    
    output << "\n\tParameters\n";
    output.write("\tt_e = {:e} [s], t_i = {:e} [s]\n", t_e, t_i);
    output.write("\tLbar = {:e} [m], Cs = {:e} [m/s]\n",
                 Lbar, Cs);
    output.write("\tTbar = {:e} [s]\n", Tbar);
    nu_e = Lbar / (Cs*t_e); SAVE_ONCE(nu_e);
    nu_i = Lbar / (Cs*t_i); SAVE_ONCE(nu_i);
    output.write("\tNormalised nu_e = {:e}, nu_i = {:e}\n", nu_e, nu_i);
    output << "\tbeta_e = " << beta_e << endl;
    output << "\tdelta = " << delta << endl;
    
    ////////////////////////////////////////////////////
    // Normalise
    
    Te0 /= Tenorm * delta; SAVE_ONCE(Te0);
    Ti0 /= Tenorm * delta; SAVE_ONCE(Ti0);
    
    Ni0 /= Ninorm * delta; SAVE_ONCE(Ni0);
    Ne0 /= Ninorm * delta; SAVE_ONCE(Ne0);
    
    rho_e /= rho_s;
    rho_i /= rho_s;
    
    output << "\tNormalised rho_e = " << rho_e << endl;
    output << "\tNormalised rho_i = " << rho_i << endl;
    
    //////////////////////////////////
    // Metric tensor components

    coord = mesh->getCoordinates();
    
    // Normalise
    hthe /= Lbar; // parallel derivatives normalised to Lperp
    
    Bpxy /= Bbar;
    Btxy /= Bbar;
    Bxy  /= Bbar;
  
    Rxy  /= rho_s; // Perpendicular derivatives normalised to rho_s
    coord->dx /= rho_s*rho_s*Bbar;
  
    // Metric components
    
    coord->g11 = SQ(Rxy*Bpxy);
    coord->g22 = 1.0 / SQ(hthe);
    coord->g33 = SQ(Bxy)/coord->g11;
    coord->g12 = 0.0;
    coord->g13 = 0.;
    coord->g23 = -Btxy/(hthe*Bpxy*Rxy);
  
    coord->J = hthe / Bpxy;
    coord->Bxy = Bxy;
  
    coord->g_11 = 1.0/coord->g11;
    coord->g_22 = SQ(Bxy*hthe/Bpxy);
    coord->g_33 = Rxy*Rxy;
    coord->g_12 = 0.;
    coord->g_13 = 0.;
    coord->g_23 = Btxy*hthe*Rxy/Bpxy;
  
    coord->geometry();
  
    // Set B field vector
  
    B0vec.covariant = false;
    B0vec.x = 0.;
    B0vec.y = Bpxy / hthe;
    B0vec.z = 0.;

    // Precompute this for use in RHS
    if (include_grad_par_B) {
      if (curv_logB) {
        Grad_par_logB = Grad_par(logB);
      } else {
        Grad_par_logB = Grad_par(log(coord->Bxy));
      }
    } else {
      Grad_par_logB = 0.;
    }
    
    //////////////////////////////////
    
    // Add ion equations
    if (ni_ddt) {
      SOLVE_FOR(Ni);
      comms.add(Ni);
    } else Ni = 0.;
  
    if (apui_ddt) {
      SOLVE_FOR(ApUi);
      comms.add(ApUi);
    } else ApUi = 0.;
    
    if (tipar_ddt) {
      SOLVE_FOR(Tipar);
      comms.add(Tipar);
    } else Tipar = 0.;
  
    if (tiperp_ddt) {
      SOLVE_FOR(Tiperp);
      comms.add(Tiperp);
    } else Tiperp = 0.;
  
    if (qipar_ddt) {
      SOLVE_FOR(qipar);
      comms.add(qipar);
    } else qipar = 0.;
  
    if (qiperp_ddt) {
      SOLVE_FOR(qiperp);
      comms.add(qiperp);
    } else qiperp = 0.;
  
    /// Split operator, with artificial dissipation in second function
    setSplitOperator(); // Split into convective and diffusive (stiff)
    
    if (adiabatic_electrons) {
      // Solving with adiabatic electrons
      output << "Error: no adiabatic electrons yet\n";
      return 1;
    } else {
      // Add electron equations
      
      if (ne_ddt) {
        SOLVE_FOR(Ne);
        comms.add(Ne);
      } else Ne = 0.;
    
      if (apue_ddt) {
        SOLVE_FOR(ApUe);
        comms.add(ApUe);
      } else ApUe = 0.;
    
      if (tepar_ddt) {
        SOLVE_FOR(Tepar);
        comms.add(Tepar);
      } else Tepar = 0.;
    
      if (teperp_ddt) {
        SOLVE_FOR(Teperp);
        comms.add(Teperp);
      } else Teperp = 0.;
    
      if (qepar_ddt) {
        SOLVE_FOR(qepar);
        comms.add(qepar);
      } else qepar = 0.;
    
      if (qeperp_ddt) {
        SOLVE_FOR(qeperp);
        comms.add(qeperp);
      } else qeperp = 0.;
    }
    
    bool output_ddt;
    output_ddt = options["output_ddt"].withDefault(false);
    if (output_ddt) {
      // Output the time derivatives
      
      if (ni_ddt)
        dump.add(ddt(Ni),     "F_Ni", 1);
      if (apui_ddt)
        dump.add(ddt(ApUi),   "F_ApUi", 1);
      if (tipar_ddt)
        dump.add(ddt(Tipar),  "F_Tipar", 1);
      if (tiperp_ddt)
        dump.add(ddt(Tiperp), "F_Tiperp", 1);
      if (qipar_ddt)
        dump.add(ddt(qipar),  "F_qipar", 1);
      if (qiperp_ddt)
        dump.add(ddt(qiperp), "F_qiperp", 1);
    
      if (!adiabatic_electrons) {
        if (ne_ddt)
          dump.add(ddt(Ne),     "F_Ne", 1);
        if (apue_ddt)
          dump.add(ddt(ApUe),   "F_ApUe", 1);
        if (tepar_ddt)
          dump.add(ddt(Tepar),  "F_Tepar", 1);
        if (teperp_ddt)
          dump.add(ddt(Teperp), "F_Teperp", 1);
        if (qepar_ddt)
          dump.add(ddt(qepar),  "F_qepar", 1);
        if (qeperp_ddt)
          dump.add(ddt(qeperp), "F_qeperp", 1);
      }
    }
    
    dump.add(phi, "phi", 1);
    dump.add(Apar, "Apar", 1);
    dump.add(Ui, "Ui", 1);
    dump.add(Ue, "Ue", 1);
    dump.add(Jpar, "Jpar", 1);
    
    // To ensure X periodicity
    mesh->communicate(comms);
    
    comms.add(phi, Apar, Ui, Ue, Jpar);
    
    dump.add(phi_G, "phi_G", 1);
    
    //////////////////////////////////
    
    if (!restarting) {
      // Initial current
      
      Field2D Jpar0;
      if (mesh->get(Jpar0, "Jpar0") == 0) {
        // Initial current specified. Set parallel electron velocity
        
      }
      
      // Initial potential
      
      Field2D phi0;
      if (mesh->get(phi0, "phi0") == 0) {
        
      }
    }
    
    phi.setBoundary("phi");
    Apar.setBoundary("Apar");
    
    // Create a solver for the Laplacian
    phiSolver = Laplacian::create(&options["phiSolver"]);

    aparSolver = Laplacian::create(&options["aparSolver"]);
    aparSolver->setCoefA(beta_e * (1./mu_e - 1./mu_i));
    
    return 0;
  }

  ////////////////////////////////////////////////////////////////////////
  // Calculate auxilliary quantities
  
  void calc_aux() {
    ////////////////////////////////////////////
    // Adiabatic electrons
    
    if (adiabatic_electrons) {
      // Solve adiabatic electrons using surface-averaged phi
      
      Field2D phi_zonal = averageY(DC(phi)); // Average over Y and Z
      Ne = phi - phi_zonal;
      
      // Need to solve with polarisation!
    }
    
    ////////////////////////////////////////////
    // Polarisation equation (quasi-neutrality)
    
    if (small_rho_e) {
      // Neglect electron Larmor radius
      
      Field3D dn = Ne - gyroPade1(Ni, rho_i) - gyroPade2(Tiperp, rho_i);
      
      phi = phiSolver->solve(tau_i * dn / SQ(rho_i));
      phi -= tau_i * dn;
    } else {
      Field3D dn = gyroPade1(Ne, rho_e) + gyroPade2(Teperp, rho_e)
        - gyroPade1(Ni, rho_i) - gyroPade2(Tiperp, rho_i);
      
      // Neglect electron gyroscreening
      phi = phiSolver->solve(tau_i * dn / (rho_i * rho_i));
      phi -= tau_i * dn;
    }
    
    phi.applyBoundary();
    
    ////////////////////////////////////////////
    // Helmholtz equation for Apar
    
    Field2D a = beta_e * (1./mu_e - 1./mu_i);
    Apar = aparSolver->solve(ApUe/mu_e - ApUi/mu_i);
    
    Apar.applyBoundary();
    
    Ui = (ApUi - beta_e*Apar) / mu_i;
    Ue = (ApUe - beta_e*Apar) / mu_e;
    
    if (jpar_bndry_width > 0) {
      // Zero j in boundary regions. Prevents vorticity drive
      // at the boundary
      
      for (int i=0;i<jpar_bndry_width;i++)
        for (int j=0;j<mesh->LocalNy;j++)
          for (int k=0;k<mesh->LocalNz;k++) {
            if (mesh->firstX()) {
              Ui(i,j,k) = 0.0;
              Ue(i,j,k) = 0.0;
            }
            if (mesh->lastX()) {
              Ui(mesh->LocalNx-1-i,j,k) = 0.0;
              Ue(mesh->LocalNx-1-i,j,k) = 0.0;
            }
          }
    }
    
    Jpar = Ui - Ue;
  }
  
  ////////////////////////////////////////////////////////////////////////
  // Non-stiff part of the RHS function
  
  int convective(BoutReal UNUSED(time)) override {
    calc_aux();
    
    ////////////////////////////////////////////
    // Communicate
    
    mesh->communicate(comms);
    
    ////////////////////////////////////////////
    // Electron equations
    
    if (!adiabatic_electrons) {
      // Electron equations
      
      if (small_rho_e) {
        // No gyro-averaging for small rho_e
        phi_G = phi;
        Phi_G = 0.0;
      } else {
        // Gyro-reduced potentials
        phi_G = gyroPade1(phi, rho_e, INVERT_RHS, INVERT_RHS);
        Phi_G = gyroPade2(phi, rho_e, INVERT_RHS, INVERT_RHS);
        
        mesh->communicate(phi_G, Phi_G);
      }
      
      if (ne_ddt) {
        ddt(Ne) = -UE_Grad(Ne0, phi_G);
        if (ne_ne1)
        ddt(Ne) -= UE_Grad(Ne, phi_G);
        if (ne_te0)
          ddt(Ne) -= WE_Grad(Te0, Phi_G);
        if (ne_te1)
          ddt(Ne) -= WE_Grad(Teperp, Phi_G);
        
        if (ne_ue)
          ddt(Ne) -= Div_parP(Ue, CELL_CENTRE);
        
        if (ne_curv)
          ddt(Ne) += curvature(phi_G + tau_e*Ne + 0.5*(tau_e*Tepar + tau_e*Teperp + Phi_G));
        
        if (low_pass_z > 0)
          ddt(Ne) = lowPass(ddt(Ne), low_pass_z);
        
        if (fix_profiles)
          ddt(Ne) -= DC(ddt(Ne));
      }
      
      if (apue_ddt) {
        if (apue_ue1_phi1) {
          ddt(ApUe) = -mu_e*UE_Grad(Ue, phi_G);
        } else
          ddt(ApUe) = 0.0;
        
        if (apue_qe1_phi1) 
          ddt(ApUe) -= mu_e*WE_Grad(qeperp, Phi_G);
      
        if (apue_phi1) // Linear term
          ddt(ApUe) -= Grad_par(phi_G, CELL_YLOW);
        if (apue_apar1_phi1) // Nonlinear term
          ddt(ApUe) += beta_e*bracket(Apar, phi_G, BRACKET_ARAKAWA);
        
        if (apue_pet) // Linear terms
          ddt(ApUe) -= tau_i*Grad_parP(Ne0 + Te0, CELL_YLOW)
                       + tau_i*Grad_par(Ne+Tepar, CELL_YLOW);
        if (apue_apar1_pe1) // Nonlinear terms
          ddt(ApUe) +=  tau_i*beta_e*bracket(Apar, Ne+Tepar, BRACKET_ARAKAWA);
        
        if (apue_curv)
          ddt(ApUe) += mu_e * tau_e * curvature(2.*Ue + qepar + 0.5*qeperp);
        
        if (apue_gradB)
          ddt(ApUe) -= tau_e * (Phi_G + tau_e*Teperp - tau_e*Tepar)*Grad_par_logB;
        
        if (low_pass_z > 0)
          ddt(ApUe) = lowPass(ddt(ApUe), low_pass_z);
        
        if (fix_profiles)
          ddt(ApUe) -= DC(ddt(ApUe));
      }
      
      if (tepar_ddt) {
        ddt(Tepar) = curvature(phi_G + tau_e*(Ne+Tepar) + 2.*tau_e*Tepar)
          - (Ue + qeperp)*Grad_par_logB
          ;
        
        if (nonlinear) {
          ddt(Tepar) += -UE_Grad(Te0 + Tepar, phi_G) 
            - 2.*Div_parP(Ue + qepar, CELL_CENTRE);
        } else {
          ddt(Tepar) += -UE_Grad(Te0, phi_G)
            - 2.*Div_par(Ue + qepar, CELL_CENTRE);
        }
        
        if (low_pass_z > 0)
          ddt(Tepar) = lowPass(ddt(Tepar), low_pass_z);
        
        if (fix_profiles)
          ddt(Tepar) -= DC(ddt(Tepar));
      }
      
      if (teperp_ddt) {
        ddt(Teperp) = 
          + 0.5*curvature(phi_G + Phi_G + tau_e*(Ne + Teperp) 
                          + 3.*(Phi_G + tau_e*Teperp))
          + (Ue + qeperp)*Grad_par_logB
          ;
        
        if (nonlinear) {
          ddt(Teperp) +=
            - UE_Grad(Te0 + Teperp, phi_G)
            - WE_Grad(Ne0 + Ne + 2.*(Te0 + Teperp), Phi_G)
            - Div_parP(qeperp, CELL_CENTRE);
        } else {
          ddt(Teperp) +=
            - UE_Grad(Te0, phi_G)
            - WE_Grad(Ne0 + 2.*Te0, Phi_G)
            - Div_par(qeperp, CELL_CENTRE);
        }
        
        if (low_pass_z > 0)
          ddt(Teperp) = lowPass(ddt(Teperp), low_pass_z);
        
        if (fix_profiles)
          ddt(Teperp) -= DC(ddt(Teperp));
      }
      
      if (qepar_ddt) {
        ddt(qepar) = 
          - 1.5*(1./mu_e)*Grad_parP(tau_e*Te0, CELL_YLOW)
          + 0.5*mu_e*tau_e*curvature(3.*Ue + 8.*qepar)
          ;
        
        if (nonlinear) {
          ddt(qepar) += 
            - UE_Grad(qepar, phi_G)
            - 1.5*(1./mu_e)*Grad_parP(tau_e*Tepar, CELL_YLOW);
        } else {
          ddt(qepar) += 
            - 1.5*(1./mu_e)*Grad_par(tau_e*Tepar, CELL_YLOW);
        }
        
        if (low_pass_z > 0)
          ddt(qepar) = lowPass(ddt(qepar), low_pass_z);
        
        if (fix_profiles)
          ddt(qepar) -= DC(ddt(qepar));
      }
      
      if (qeperp_ddt) {
        ddt(qeperp) = 
          + (1./mu_e)*beta_e*bracket(Apar, tau_e*Te0, BRACKET_ARAKAWA)
          + 0.5*tau_e*curvature(Ue + 6.*qeperp)
          - (tau_e/mu_e)*(Phi_G + tau_e*Teperp - tau_e*Tepar)*Grad_par_logB
          ;
        
        if (nonlinear) {
          ddt(qeperp) +=
            - UE_Grad(qeperp, phi_G)
            - WE_Grad(Ue + 2.*qeperp, Phi_G)
            - (1./mu_e)*Grad_parP(Phi_G + tau_e*Teperp, CELL_YLOW);
        } else {
          ddt(qeperp) +=
            -(1./mu_e)*Grad_par(Phi_G + tau_e*Teperp, CELL_YLOW);
        }
        
        if (low_pass_z > 0)
          ddt(qeperp) = lowPass(ddt(qeperp), low_pass_z);
        
        if (fix_profiles)
          ddt(qeperp) -= DC(ddt(qeperp));
      }
    }
    
    ////////////////////////////////////////////
    // Ion equations
  
    // Calculate gyroreduced potentials
    phi_G = gyroPade1(phi, rho_i, INVERT_RHS, INVERT_RHS);
    Phi_G = gyroPade2(phi, rho_i, INVERT_RHS, INVERT_RHS);
    
    mesh->communicate(phi_G, Phi_G);
    
    if (ni_ddt) {
      ddt(Ni) = -UE_Grad(Ni0, phi_G);
      if (ni_ni1)
        ddt(Ni) -= UE_Grad(Ni, phi_G);
      if (ni_ti0)
        ddt(Ni) -= WE_Grad(Ti0, Phi_G);
      if (ni_ti1)
        ddt(Ni) -= WE_Grad(Tiperp, Phi_G);
      
      if (ni_ui)
        ddt(Ni) -= Div_parP(Ui, CELL_CENTRE);
      
      if (ni_curv)
        ddt(Ni) += curvature(phi_G + tau_i*Ni + 0.5*(tau_i*Tipar + tau_i*Tiperp + Phi_G));
      
      if (low_pass_z > 0)
        ddt(Ni) = lowPass(ddt(Ni), low_pass_z);
      
      if (fix_profiles)
        ddt(Ni) -= DC(ddt(Ni));
    }
    
    if (apui_ddt) {
      if (apui_ui1_phi1) {
        ddt(ApUi) = -mu_i*UE_Grad(Ui, phi_G);
      } else
        ddt(ApUi) = 0.0;
      
      if (apui_qi1_phi1) 
        ddt(ApUi) -= mu_i*WE_Grad(qiperp, Phi_G);
    
      if (apui_phi1)
        ddt(ApUi) -= Grad_par(phi_G, CELL_YLOW);
      if (apui_apar1_phi1) // Nonlinear term
        ddt(ApUi) += beta_e*bracket(Apar, phi_G, BRACKET_ARAKAWA);
      
      if (apui_pit) // Linear terms
        ddt(ApUi) -= tau_i*Grad_parP(Ni0 + Ti0, CELL_YLOW)
                     + tau_i*Grad_par(Ni+Tipar, CELL_YLOW);
      if (apui_apar1_pi1) // Nonlinear terms
        ddt(ApUi) +=  tau_i*beta_e*bracket(Apar, Ni+Tipar, BRACKET_ARAKAWA);
      
      if (apui_curv)
        ddt(ApUi) += mu_i * tau_i * curvature(2.*Ui + qipar + 0.5*qiperp);
      
      if (apui_gradB)
        ddt(ApUi) -= tau_i * (Phi_G + tau_i*Tiperp - tau_i*Tipar)*Grad_par_logB;
      
      if (low_pass_z > 0)
        ddt(ApUi) = lowPass(ddt(ApUi), low_pass_z);
      
      if (fix_profiles)
        ddt(ApUi) -= DC(ddt(ApUi));
    }
    
    if (tipar_ddt) {
      ddt(Tipar) = 
        + curvature(phi_G + tau_i*(Ni+Tipar) + 2.*tau_i*Tipar)
        - (Ui + qiperp)*Grad_par_logB
        ;
      
      if (nonlinear) {
        ddt(Tipar) += 
          -UE_Grad(Ti0 + Tipar, phi_G)
          - 2.*Div_parP(Ui + qipar, CELL_CENTRE);
      } else {
        ddt(Tipar) += 
          -UE_Grad(Ti0, phi_G)
          - 2.*Div_par(Ui + qipar, CELL_CENTRE);
      }
      
      if (low_pass_z > 0)
        ddt(Tipar) = lowPass(ddt(Tipar), low_pass_z);
      
      if (fix_profiles)
        ddt(Tipar) -= DC(ddt(Tipar));
    }
    
    if (tiperp_ddt) {
      ddt(Tiperp) = 
        + 0.5*curvature(phi_G + Phi_G + tau_i*(Ni + Tiperp) 
                        + 3.*(Phi_G + tau_i*Tiperp))
        + (Ui + qiperp)*Grad_par_logB
        ;
      
      if (nonlinear) {
        ddt(Tiperp) +=
          - UE_Grad(Ti0 + Tiperp, phi_G)
          - WE_Grad(Ni0 + Ni + 2.*(Ti0 + Tiperp), Phi_G)
          - Div_parP(qiperp, CELL_CENTRE);
      } else {
        ddt(Tiperp) +=
          - UE_Grad(Ti0, phi_G)
          - WE_Grad(Ni0 + 2.*Ti0, Phi_G)
          - Div_par(qiperp, CELL_CENTRE);
      }
      
      if (low_pass_z > 0)
        ddt(Tiperp) = lowPass(ddt(Tiperp), low_pass_z);
      
      if (fix_profiles)
        ddt(Tiperp) -= DC(ddt(Tiperp));
    }
    
    if (qipar_ddt) {
      ddt(qipar) = 
        + 1.5*(1./mu_i)*beta_e*bracket(Apar, tau_i*Ti0, BRACKET_ARAKAWA)
        + 0.5*tau_i*curvature(3.*Ui + 8.*qipar);
      
      if (nonlinear) {
        ddt(qipar) +=
          - UE_Grad(qipar, phi_G)
          - 1.5*(1./mu_i)*Grad_parP(tau_i*Tipar, CELL_YLOW);
      } else {
        ddt(qipar) +=
          - 1.5*(1./mu_i)*Grad_par(tau_i*Tipar, CELL_YLOW);
      }
      
      if (low_pass_z > 0)
        ddt(qipar) = lowPass(ddt(qipar), low_pass_z);
      
      if (fix_profiles)
        ddt(qipar) -= DC(ddt(qipar));
    }
    
    if (qiperp_ddt) {
      ddt(qiperp) = 
        + (1./mu_i)*beta_e*bracket(Apar, tau_i*Ti0, BRACKET_ARAKAWA)
        + 0.5*tau_i*curvature(Ui + 6.*qiperp)
        - (tau_i/mu_i)*(Phi_G + tau_i*Tiperp - tau_i*Tipar)*Grad_par_logB;
      
      if (nonlinear) {
        ddt(qiperp) +=
          - UE_Grad(qiperp, phi_G)
          - WE_Grad(Ui + 2.*qiperp, Phi_G)
          - (1./mu_i)*Grad_parP(Phi_G + tau_i*Tiperp, CELL_YLOW);
      } else {
        ddt(qiperp) +=
          - (1./mu_i)*Grad_par(Phi_G + tau_i*Tiperp, CELL_YLOW);
      }
      
      if (low_pass_z > 0)
        ddt(qiperp) = lowPass(ddt(qiperp), low_pass_z);
      
      if (fix_profiles)
        ddt(qiperp) -= DC(ddt(qiperp));
    }
    
    return 0;
  }

  // Stiff part of the RHS function
  // Artificial dissipation terms
  int diffusive(BoutReal UNUSED(time)) override {
    Field3D S_D, K_par, K_perp, K_D; // Collisional dissipation terms
    
    calc_aux();
    
    ////////////////////////////////////////////
    // Communicate
    
    mesh->communicate(comms);
    
    ////////////////////////////////////////////
    // Resistivity
    
    Rei = mu_e*nu_e*(eta*Jpar + 
                     (alpha_e/kappa_e)*(qepar + qeperp + alpha_e*Jpar));
    
    ////////////////////////////////////////////
    // Electron equations
    
    if (!adiabatic_electrons) {
      if (small_rho_e) {
        // No gyro-averaging for small rho_e
        phi_G = phi;
        Phi_G = 0.0;
      } else {
        // Gyro-reduced potentials
        phi_G = gyroPade1(phi, rho_e, INVERT_RHS, INVERT_RHS);
        Phi_G = gyroPade2(phi, rho_e, INVERT_RHS, INVERT_RHS);
        
        mesh->communicate(phi_G, Phi_G);
      }
      
      // Collisional dissipation
      S_D = (nu_e / (3.*pi_e)) * (Tepar - Teperp);
      K_par = mu_e*tau_e*nu_e*((5./2.)/kappa_e)*(qepar + 0.6*alpha_e*Jpar);
      K_perp = mu_e*tau_e*nu_e*((5./2.)/kappa_e)*(qeperp + 0.4*alpha_e*Jpar);
      K_D = 1.28*mu_e*tau_e*nu_e*((5./2.)/kappa_e)*(qepar - 1.5*qeperp);
      
      if (ne_ddt) {
        ddt(Ne) = 0.;
        if (ne_ne1) {
          ddt(Ne) -= UE_Grad_D(Ne, phi_G);
        }
      }
      if (apue_ddt) {
        ddt(ApUe) = 0.0;
        if (apue_ue1_phi1)
          ddt(ApUe) -= mu_e*UE_Grad_D(Ue, phi_G);
        if (apue_Rei)
          ddt(ApUe) -= Rei;
      }
      if (tepar_ddt) {
        ddt(Tepar) = -UE_Grad_D(Tepar, phi_G) - 2.*S_D;
        
      }
      if (teperp_ddt) {
        ddt(Teperp) = -UE_Grad_D(Teperp, phi_G) + S_D;
      }
      if (qepar_ddt) {
        ddt(qepar) = -UE_Grad_D(qepar, phi_G)
          //- Landau*(tau_e/mu_e)*(1. - 0.125*Grad2_par2(qepar))
          - (1./mu_e)*K_par
          - (1./mu_e)*K_D;
      }
      if (qeperp_ddt) {
        ddt(qeperp) = -UE_Grad_D(qeperp, phi_G)
          - (1./mu_e)*K_perp
          + (1./mu_e)*K_D;
      }
    }
    
    ////////////////////////////////////////////
    // Ion equations
    
    // Calculate gyroreduced potentials
    phi_G = gyroPade1(phi, rho_i, INVERT_RHS, INVERT_RHS);
    Phi_G = gyroPade2(phi, rho_i, INVERT_RHS, INVERT_RHS);
    
    mesh->communicate(phi_G, Phi_G);
    
    // Collisional dissipation
    S_D = (nu_i / (3.*pi_i)) * (Tipar - Tiperp);
    K_par = mu_i*tau_i*nu_i*((5./2.)/kappa_i)*qipar;
    K_perp = mu_i*tau_i*nu_i*((5./2.)/kappa_i)*qiperp;
    K_D = 1.28*mu_i*tau_i*nu_i*((5./2.)/kappa_i)*(qipar - 1.5*qiperp);
    
    if (ni_ddt) {
      ddt(Ni) = 0.;
      if (ni_ni1)
        ddt(Ni) -= UE_Grad_D(Ni, phi_G);
    }
    if (apui_ddt) {
      ddt(ApUi) = 0.0;
      if (apui_ui1_phi1)
        ddt(ApUi) = -mu_i*UE_Grad_D(Ui, phi_G);
      if (apui_Rei)
        ddt(ApUi) -= Rei;
    }
    if (tipar_ddt) {
      ddt(Tipar) = - UE_Grad_D(Tipar, phi_G) - 2.*S_D;
    }
    if (tiperp_ddt) {
      ddt(Tiperp) = - UE_Grad_D(Tiperp, phi_G) + S_D;
    }
    if (qipar_ddt) {
      ddt(qipar) = - UE_Grad_D(qipar, phi_G)
        - (1./mu_e)*K_par
        - (1./mu_e)*K_D;
    }
    if (qiperp_ddt) {
      ddt(qiperp) = - UE_Grad_D(qiperp, phi_G)
        - (1./mu_e)*K_perp
        + (1./mu_e)*K_D;
    }
    
    return 0;
  }
  
  ////////////////////////////////////////////////////////////////////////
  // Curvature operator
  
  // K(f) = Div((c/B^2) B x Grad(f))
  // Simple implementation. Could be improved to eliminate the communication
  const Field3D curvature(const Field3D &f) {
    if (curv_logB) {
      return -bracket(2.*logB, f, BRACKET_ARAKAWA);
    }
    return -bracket(2.*log(coord->Bxy), f, BRACKET_ARAKAWA);
  }
  
  ////////////////////////////////////////////////////////////////////////
  // Advection terms
  
  /// ExB advection
  const Field3D UE_Grad(const Field3D &f, const Field3D &p) {
    return bracket(p, f, bm);
  }
  
  /// Artificial dissipation terms in advection
  const Field3D UE_Grad_D(const Field3D &f, const Field3D &UNUSED(p)) {
    Field3D delp2 = Delp2(f);
    delp2.applyBoundary("neumann");
    mesh->communicate(delp2);
    
    return nu_perp*Delp2( delp2 * SQ(SQ(1./coord->Bxy)) )
      - nu_par*Grad2_par2(f) // NB: This should be changed for variable B
      ;
  }
  
  const Field3D WE_Grad(const Field3D &f, const Field3D &p) {
    return bracket(p, f, bm);
  }
  
  ////////////////////////////////////////////////////////////////////////
  // Parallel derivative
  
  const Field3D Grad_parP(const Field3D &f, CELL_LOC loc = CELL_DEFAULT) {
    return Grad_par(f, loc)
           - beta_e*bracket(interp_to(Apar, loc), interp_to(f, loc), BRACKET_ARAKAWA);
  }
  
  const Field3D Div_parP(const Field3D &f, CELL_LOC loc = CELL_DEFAULT) {
    return interp_to(coord->Bxy, loc)
           *Grad_parP(f/interp_to(coord->Bxy, f.getLocation()), loc);
  }

};

BOUTMAIN(GEM);

