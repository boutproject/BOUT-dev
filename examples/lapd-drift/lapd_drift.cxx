/*******************************************************************************
 * 2-fluid equations
 * Same as Maxim's version of BOUT - simplified 2-fluid for benchmarking
 * LAPD standard case with simulation results published in PoP in Popovich et. al. 2010
 *******************************************************************************/
#include <bout/physicsmodel.hxx>

#include <initialprofiles.hxx>
#include <derivs.hxx>
#include <interpolation.hxx>
#include <invert_laplace.hxx>
#include <boutexception.hxx>

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/// Solves 2-fluid equations for turbulence in a linear device
/// 
class LAPDdrift : public PhysicsModel {
private:
  // 2D initial profiles
  Field2D Ni0, Ti0, Te0, Vi0, phi0, Ve0, rho0, Ajpar0, src_ni0;
  Vector2D b0xcv; // for curvature terms
  
  // 3D evolving fields
  Field3D rho, ni, ajpar, te;
  
  // Derived 3D variables
  Field3D phi, Apar, Ve, jpar;
  
  // Non-linear coefficients
  Field3D nu, mu_i, kapa_Te, kapa_Ti;
  
  // 3D total values
  Field3D Nit, Tit, Tet, Vit, phit, VEt, dphi_bc_ydown, dphi_bc_yup;

  // pressures
  Field3D pei, pe;
  Field2D pei0, pe0;
  
  // Metric coefficients
  Field2D Rxy, Bpxy, Btxy, hthe,Zxy;
  
  // parameters
  BoutReal Te_x, Ti_x, Ni_x, Vi_x, bmag, rho_s, fmei, AA, ZZ;
  BoutReal lambda_ei, lambda_ii;
  BoutReal nu_hat, mui_hat, wci, nueix, nuiix;
  BoutReal ni_perpdiff, rho_perpdiff, te_perpdiff;
  BoutReal beta_p;
  BoutReal nuIonNeutral; // Ion-neutral collision rate (normalised by wci)
  
  // settings
  bool estatic, ZeroElMass; // Switch for electrostatic operation (true = no Apar)
  
  bool arakawa;   // Use Arakawa scheme for ExB advection
  bool bout_exb;  // Use BOUT-06 expression for ExB velocity
  
  BoutReal zeff, nu_perp;
  bool evolve_rho, evolve_ni, evolve_ajpar, evolve_te;
  BoutReal ShearFactor;
  BoutReal time_step;
  
  bool nonlinear;
  bool neg_fix;
  
  BoutReal ni_floor, minNit;
  
  bool filter_z;
  int filter_z_mode;
  
  bool log_density;  // Evolve logarithm of the density
  
  int phi_flags, apar_flags; // Inversion flags
  
  bool niprofile;
  
  bool evolve_source_ni, evolve_source_te; // If true, evolve a source/sink profile
  BoutReal source_response;  // Initial source response (inverse timescale) 
  BoutReal source_converge;  // Timescale for convergence
  Field2D Sn,St; // Density source (inverse timescale)
  bool input_source; // Read Sn from the input file
  bool remove_tor_av_ni, remove_tor_av_te; // Subtract the toroidal averages
  
  // Switches for terms in the ni equation
  bool ni_jpar1, ni_ni0_phi1, ni_ni1_phi0, ni_ni1_phi1, ni_src_ni0, ni_diff;
  
  // Switches for terms in the rho equation
  bool rho_jpar1, rho_nuin_rho1, rho_rho1, rho_rho1_phi1, rho_ve2t, rho_diff;
  bool rho_rho0_phi1, rho_rho1_phi0, rho_ve2lin;
  
  // Switches for terms in the ajpar equation
  bool ajpar_phi1, ajpar_jpar1, ajpar_te_ni, ajpar_te; 
  bool ajpar_ajpar1_phi1, ajpar_ajpar1_phi0, ajpar_ve1_ve1;
  
  // Switches for terms in the te equation
  bool te_te1_phi0, te_te0_phi1, te_te1_phi1, te_ajpar_te;
  bool te_te_ajpar, te_nu_te1, te_nu_tet, te_jpar, te_diff;
  
  // Coefficients for linear sheath problem
  Field2D LAMBDA1, LAMBDA2;
  
  // My ixseps variables
  int my_ixseps;
  
  // Communication object
  FieldGroup comms;

  // Laplacian inversion object
  Laplacian *phiSolver;
protected:
  
  /// Function called once at the start of the simulation
  ///
  /// @param[in] restarting  True if simulation is restarting
  int init(bool UNUSED(restarting)) {
    Field2D I; // Shear factor 
    
    output.write("Solving LAPD drift test case\n");
    
    /************* LOAD DATA FROM GRID FILE ****************/
    
    // Load 2D profiles (set to zero if not found)
    mesh->get(Ni0,    "Ni0");
    mesh->get(Ti0,    "Ti0");
    mesh->get(Te0,    "Te0");
    mesh->get(Vi0,    "Vi0");
    mesh->get(Ve0,    "Ve0");
    mesh->get(phi0,   "phi0");
    mesh->get(rho0,   "rho0");
    mesh->get(Ajpar0, "Ajpar0");
    mesh->get(src_ni0, "src_ni0");
    
    // Load magnetic curvature term
    b0xcv.covariant = false; // Read contravariant components
    mesh->get(b0xcv, "bxcv"); // b0xkappa terms

    Coordinates *coord = mesh->getCoordinates();
    
    // Load metrics
    mesh->get(Rxy,  "Rxy");
    mesh->get(Zxy,  "Zxy");
    mesh->get(Bpxy, "Bpxy");
    mesh->get(Btxy, "Btxy");
    mesh->get(hthe, "hthe");
    mesh->get(coord->dx,   "dpsi");
    mesh->get(I,    "sinty");
    
    // Load normalisation values
    mesh->get(Te_x, "Te_x");
    mesh->get(Ti_x, "Ti_x");
    mesh->get(Ni_x, "Ni_x");
    mesh->get(bmag, "bmag");
    
    // Get separatrix location
    mesh->get(my_ixseps, "ixseps1");
    
    Ni_x *= 1.0e14;
    bmag *= 1.0e4;

    /*************** READ OPTIONS *************************/
    // Read some parameters
    
    auto globalOptions = Options::root();
    
    time_step = globalOptions["TIMESTEP"].withDefault(1.0);

    auto options = globalOptions["2fluid"];
    AA = options["AA"].withDefault(4.0); // <=> AA = options["AA"].withDefault(1.0);
    ZZ = options["ZZ"].withDefault(1.0);

    estatic = options["estatic"].withDefault(false);
    ZeroElMass = options["ZeroElMass"].withDefault(false);
    zeff = options["zeff"].withDefault(1.0);
    nu_perp = options["nu_perp"].withDefault(0.0);
    ShearFactor = options["ShearFactor"].withDefault(1.0);
    nuIonNeutral = options["nuIonNeutral"].withDefault(-1.);
    arakawa = options["arakawa"].withDefault(false);
    bout_exb = options["bout_exb"].withDefault(false);

    niprofile = options["niprofile"].withDefault(false);
    evolve_source_ni = options["evolve_source_ni"].withDefault(false);
    evolve_source_te = options["evolve_source_te"].withDefault(false);
    source_response = options["source_response"].withDefault(1.0);
    source_converge = options["source_converge"].withDefault(-1);

    ni_perpdiff = options["ni_perpdiff"].withDefault(0.0);
    rho_perpdiff = options["rho_perpdiff"].withDefault(0.0);
    te_perpdiff = options["te_perpdiff"].withDefault(0.0);

    input_source = options["input_source"].withDefault(false);
    remove_tor_av_ni = options["remove_tor_av_ni"].withDefault(false);
    remove_tor_av_te = options["remove_tor_av_te"].withDefault(false);

    phi_flags = options["phi_flags"].withDefault(0);
    apar_flags = options["apar_flags"].withDefault(0);

    nonlinear = options["nonlinear"].withDefault(true);

    // Toroidal filtering
    filter_z = options["filter_z"].withDefault(false); // Filter a single n
    filter_z_mode = options["filter_z_mode"].withDefault(1);

    // Set default values for terms in each equation
    // Allows default to be overridden in BOUT.inp file
    auto option_rho = globalOptions["rho"];
    evolve_rho = option_rho["evolve_rho"].withDefault(true);
    rho_jpar1 = option_rho["rho_jpar1"].withDefault(false);
    rho_nuin_rho1 = option_rho["rho_nuin_rho1"].withDefault(false);
    rho_rho1 = option_rho["rho_rho1"].withDefault(false);
    rho_rho0_phi1 = option_rho["rho_rho0_phi1"].withDefault(false);
    rho_rho1_phi0 = option_rho["rho_rho1_phi0"].withDefault(false);
    rho_ve2lin = option_rho["rho_ve2lin"].withDefault(false);
    rho_rho1_phi1 = option_rho["rho_rho1_phi1"].withDefault(false);
    rho_ve2t = option_rho["rho_ve2t"].withDefault(false);
    rho_diff = option_rho["rho_diff"].withDefault(false);

    auto option_ni = globalOptions["ni"];
    evolve_ni = option_ni["evolve_ni"].withDefault(true);
    ni_jpar1 = option_ni["ni_jpar1"].withDefault(false);
    ni_ni0_phi1 = option_ni["ni_ni0_phi1"].withDefault(false);
    ni_ni1_phi0 = option_ni["ni_ni1_phi0"].withDefault(false);
    ni_ni1_phi1 = option_ni["ni_ni1_phi1"].withDefault(false);
    ni_src_ni0 = option_ni["ni_src_ni0"].withDefault(false);
    ni_diff = option_ni["ni_diff"].withDefault(false);

    auto option_ajpar = globalOptions["ajpar"];
    evolve_ajpar = option_ajpar["evolve_ajpar"].withDefault(true);
    ajpar_phi1 = option_ajpar["ajpar_phi1"].withDefault(false);
    ajpar_jpar1 = option_ajpar["ajpar_jpar1"].withDefault(false);
    ajpar_te_ni = option_ajpar["ajpar_te_ni"].withDefault(false);
    ajpar_te = option_ajpar["ajpar_te"].withDefault(false);
    ajpar_ajpar1_phi0 = option_ajpar["ajpar_ajpar1_phi0"].withDefault(false);
    ajpar_ajpar1_phi1 = option_ajpar["ajpar_ajpar1_phi1"].withDefault(false);
    ajpar_ve1_ve1 = option_ajpar["ajpar_ve1_ve1"].withDefault(false);

    auto option_te = globalOptions["te"];
    evolve_te = option_te["evolve_te"].withDefault(true);
    te_te1_phi0 = option_te["te_te1_phi0"].withDefault(false);
    te_te0_phi1 = option_te["te_te0_phi1"].withDefault(false);
    te_te1_phi1 = option_te["te_te1_phi1"].withDefault(false);
    te_ajpar_te = option_te["te_ajpar_te"].withDefault(false);
    te_te_ajpar = option_te["te_te_ajpar"].withDefault(false);
    te_nu_te1 = option_te["te_nu_te1"].withDefault(false);
    te_nu_tet = option_te["te_nu_tet"].withDefault(false);
    te_jpar = option_te["te_jpar"].withDefault(false);
    te_diff = option_te["te_diff"].withDefault(false);

    if (ZeroElMass) {
      evolve_ajpar = false; // Don't need ajpar - calculated from ohm's law
    }
      
    /************* SHIFTED RADIAL COORDINATES ************/
    
    // Check type of parallel transform
    std::string ptstr = Options::root()["mesh"]["paralleltransform"].withDefault<std::string>("identity");

    if (lowercase(ptstr) == "shifted") {
      ShearFactor = 0.0;  // I disappears from metric
      b0xcv.z += I*b0xcv.x;
    }
    
    /************** CALCULATE PARAMETERS *****************/
    
    rho_s = 1.02*sqrt(AA*Te_x)/ZZ/bmag;
    fmei  = 1./1836.2/AA;
    
    lambda_ei = 24.-log(sqrt(Ni_x)/Te_x);
    lambda_ii = 23.-log(ZZ*ZZ*ZZ*sqrt(2.*Ni_x)/pow(Ti_x, 1.5));
    wci       = 9.58e3*ZZ*bmag/AA;
    nueix     = 2.91e-6*Ni_x*lambda_ei/pow(Te_x, 1.5);
    nuiix     = 4.78e-8*pow(ZZ,4.)*Ni_x*lambda_ii/pow(Ti_x, 1.5)/sqrt(AA);
    nu_hat    = zeff*nueix/wci;
    mui_hat   = 0.96*wci/nuiix*pow(Ti_x/Te_x, -1.5);
    
    
    if (estatic) {
      beta_p    = 1.e-29;
    } else {
      beta_p    = 4.03e-11*Ni_x*Te_x/bmag/bmag;
    }

    Vi_x = wci * rho_s;
    
    output.write("Collisions: nueix = %e, nu_hat = %e\n", nueix, nu_hat);
    
    /************** PRINT Z INFORMATION ******************/
    
    BoutReal hthe0;
    if(mesh->get(hthe0, "hthe0") == 0) {
      output.write("    ****NOTE: input from BOUT, Z length needs to be divided by %e\n", hthe0/rho_s);
    }
    
    /************** SHIFTED GRIDS LOCATION ***************/
    
    // Velocities defined on cell boundaries
    ajpar.setLocation(CELL_YLOW);

    // Apar and jpar too
    Apar.setLocation(CELL_YLOW); 
    jpar.setLocation(CELL_YLOW);
    
    /************** NORMALISE QUANTITIES *****************/
    
    output.write("\tNormalising to rho_s = %e\n", rho_s);

    // Normalise profiles
    Ni0  /= Ni_x/1.0e14;
    Ti0  /= Te_x;
    Te0  /= Te_x;
    phi0 /= Te_x;
    Vi0  /= Vi_x;
    
    // Normalise curvature term
    b0xcv.x /= (bmag/1e4);
    b0xcv.y *= rho_s*rho_s;
    b0xcv.z *= rho_s*rho_s;
    
    // Normalise geometry 
    Rxy /= rho_s;
    hthe /= rho_s;
    I *= rho_s*rho_s*(bmag/1e4)*ShearFactor;
    coord->dx /= rho_s*rho_s*(bmag/1e4);
    
    // Normalise magnetic field
    Bpxy /= (bmag/1.e4);
    Btxy /= (bmag/1.e4);
    coord->Bxy  /= (bmag/1.e4);
    
    // calculate pressures
    pei0 = (Ti0 + Te0)*Ni0;
    pe0 = Te0*Ni0;
    
    /**************** CALCULATE METRICS ******************/
    
    coord->g11 = SQ(Rxy*Bpxy);
    coord->g22 = 1.0 / SQ(hthe);
    coord->g33 = SQ(I)*coord->g11 + SQ(coord->Bxy)/coord->g11;
    coord->g12 = 0.0;
    coord->g13 = -I*coord->g11;
    coord->g23 = -Btxy/(hthe*Bpxy*Rxy);
    
    coord->J = hthe / Bpxy;
    
    coord->g_11 = 1.0/coord->g11 + SQ(I*Rxy);
    coord->g_22 = SQ(coord->Bxy*hthe/Bpxy);
    coord->g_33 = Rxy*Rxy;
    coord->g_12 = Btxy*hthe*I*Rxy/Bpxy;
    coord->g_13 = I*Rxy*Rxy;
    coord->g_23 = Btxy*hthe*Rxy/Bpxy;
    
    coord->geometry();
    
    rho0 = Ni0*Delp2(phi0) + Perp_Grad_dot_Grad(phi0,Ni0);
    
    /**************** SET EVOLVING VARIABLES *************/
    
    // Tell BOUT++ which variables to evolve
    // add evolving variables to the communication object
    if (evolve_rho) {
      SOLVE_FOR(rho);
      comms.add(rho);
    } else
      initial_profile("rho", rho);
    
    if (evolve_ni) {
      SOLVE_FOR(ni);
      comms.add(ni);
    } else
      initial_profile("ni", ni);
    
    if (evolve_ajpar) {
      SOLVE_FOR(ajpar);
      comms.add(ajpar);
    } else {
      initial_profile("ajpar", ajpar);
      if (ZeroElMass)
        dump.add(ajpar, "ajpar", 1); // output calculated Ajpar
    }
    
    if (evolve_te) {
      SOLVE_FOR(te);
      comms.add(te);
    } else
      initial_profile("te", te);
    
    // Set boundary conditions on jpar and VEt
    jpar.setBoundary("jpar");
    VEt.setBoundary("VEt");
    
    if (evolve_source_ni) {
      SOLVE_FOR(Sn);
    }
    if (input_source) {
      mesh->get(Sn, "Sn");
      dump.add(Sn, "Sn");
    }
    
    if (evolve_source_te) {
      SOLVE_FOR(St);
    }

    /************** SETUP COMMUNICATIONS **************/
    
    // add extra variables to communication
    comms.add(phi);
    
    /*************** DUMP VARIABLES TO OUTPUT**********/
    dump.add(phi,  "phi",  1);  dump.add(jpar, "jpar", 1);
    
    SAVE_ONCE(Ni0,Te0,phi0,rho0);
    SAVE_ONCE(Rxy,Bpxy,Btxy,Zxy,hthe);
    dump.addOnce(coord->Bxy, "Bxy");
    dump.addOnce(my_ixseps, "ixseps");
    
    SAVE_ONCE(Te_x,Ti_x,Ni_x);
    SAVE_ONCE(AA,ZZ,zeff,rho_s,wci,bmag);
    dump.addOnce(mesh->LocalNx, "ngx");
    dump.addOnce(mesh->LocalNy, "ngy"); 
    dump.addOnce(mesh->LocalNz, "ngz");
    SAVE_ONCE(mui_hat,nu_hat,nuIonNeutral,beta_p,time_step,hthe0);
    SAVE_ONCE(ni_perpdiff,rho_perpdiff,te_perpdiff);

    // Laplacian inversion solver
    phiSolver = Laplacian::create();
    phiSolver->setFlags(phi_flags);
    phiSolver->setCoefC(Ni0);
    
    return 0;
  }
  // End of physics_init()
  //////////////////////////////////////
  

  //////////////////////////////////////
  /// Function called at each time step
  /// Time derivatives calculated here
  int rhs(BoutReal t) {

    Coordinates *coord = mesh->getCoordinates();
    
    // Invert vorticity to get phi
    
    // Solves \nabla^2_\perp x + (1./c)*\nabla_perp c\cdot\nabla_\perp x + a x = b
    // Arguments are:   (b,   bit-field, a,    c)
    // Passing NULL -> missing term
    if (nonlinear) {
      phi = phiSolver->solve(rho/(Ni0+ni));
    } else {
      phi = phiSolver->solve(rho/Ni0);
    }
    
    // Communicate variables
    mesh->communicate(comms);
    
    // Update profiles
    if (nonlinear) {
      Nit = Ni0 + ni;
      phit = phi0 + phi;
      Tit = Ti0;
      Tet = Te0 + te;
    } else {
      Nit = Ni0;
      phit = phi0;
      Tit = Ti0;
      Tet = Te0;
    }
    
    BoutReal source_alpha;
    
    // Calculate source response
    if (source_converge > 0.) {
      source_alpha = source_response * exp(-1.*t/source_converge);
    } else {
      source_alpha = source_response;
    }

    // Exit if the density goes negative
    if (nonlinear && min(Nit) < 0.0) {
      output.enable();  // Use stdout for the next line
      throw BoutException("Unphysical negative density encountered. Exiting...\n");
    }


    // Exit if the temperature goes negative or make negatives zero
    if (nonlinear && evolve_te && min(Tet) < 0.0) {
      output.enable();  // Use stdout for the next line
      throw BoutException("Unphysical negative temperature encountered. Exiting...\n");
    }
  
    // Update non-linear coefficients
    nu      = nu_hat * Nit / pow(Tet,1.5);
    mu_i    = mui_hat * pow(Tit,2.5)/Nit;
    //kapa_Te = 3.2*(1./fmei)*(wci/nueix)*pow(Tet,2.5);
    //kapa_Ti = 3.9*(wci/nuiix)*(Tit^2.5);
    
    // Calculate pressures
    //pei = (Tet+Tit)*Nit;
    //pe  = Tet*Nit;
    
    // Calculate E cross B velocity
    if (nonlinear) {
      VEt = sqrt(coord->g11*DDX(phit)*DDX(phit) + coord->g33*DDZ(phit)*DDZ(phit));
      
      // Set boundary condition on VEt
      VEt.applyBoundary();
      
      // Communicate VEt
      mesh->communicate(VEt);
    }
    
    if (ZeroElMass) {
      // Set jpar,Ve,Ajpar neglecting the electron inertia term
      jpar = ((Tet*Grad_par_LtoC(ni)) - (Nit*Grad_par_LtoC(phi)))/(fmei*0.51*nu);
      
      // Set boundary condition on jpar
      jpar.applyBoundary();
      
      // Need to communicate jpar
      mesh->communicate(jpar);
      
      Ve = -jpar/Nit;
      ajpar = Ve;
    } else {
    
      Ve = ajpar;
      jpar = -Nit*Ve;
      //jpar = -Ni0*Ve; //Linearize as in BOUT06
    }
    //////////////////////////////////////////////////////////////////////////////////////////
    // DENSITY EQUATION
    
    ddt(ni) = 0.0;
    if (evolve_ni) {
      
      if (ni_ni0_phi1) {
        ddt(ni) -= DDX(Ni0)*DDZ(phi);
      }
      
      if (ni_ni1_phi0) {
        ddt(ni) -= vE_Grad(ni, phi0);
      }
      
      if (ni_ni1_phi1) {
        ddt(ni) -= vE_Grad(ni, phi);
      }
      
      if (ni_jpar1) {
        ddt(ni) += Grad_par_CtoL(jpar); // Left hand differencing
      }
      
      if (ni_src_ni0) {
        ddt(ni) += src_ni0;
      }
      
      if (ni_diff) {
        ddt(ni) += ni_perpdiff * Delp2(ni);
      }

      if (evolve_source_ni) {
        ddt(Sn) = averageY(-1. * source_alpha * DC(ni) / Ni0);
        
        // Add density source/sink
        ddt(ni) += Sn*where(Sn, Ni0, Nit); // Sn*Ni0 if Sn > 0, Sn*Nit if Sn < 0
      }
      
      if(remove_tor_av_ni) {
        ddt(ni) -= DC(ddt(ni)); // REMOVE TOROIDAL AVERAGE DENSITY
      }
    }
    ///////////////////////////////////////////////////////////////////////////////////////////
    // VORTICITY
    
    ddt(rho) = 0.0;
    if (evolve_rho) {
      
      if (rho_jpar1) { 
        ddt(rho) += Grad_par_CtoL(jpar); // Left hand differencing
      }
      
      if (rho_nuin_rho1) {
        ddt(rho) -= nuIonNeutral * rho;
      }
      
      if (rho_rho1) {
        ddt(rho) += mu_i * Delp2(rho);
      }
      
      if (rho_diff) {
        ddt(rho) += rho_perpdiff * Delp2(rho);
      }
      
      if (rho_rho1_phi1) {
        ddt(rho) -= vE_Grad(rho, phi);
      }

      if (rho_rho0_phi1) {
        ddt(rho) -= vE_Grad(rho0, phi);
      }
      
      if (rho_rho1_phi0) {
        ddt(rho) -= vE_Grad(rho, phi0);
      }
      
      if (rho_ve2lin) {
        ddt(rho) -= coord->g11*coord->g33 * DDX(phi0)*(DDX(Ni0)*D2DXDZ(phi) - D2DX2(phi0)*DDZ(ni));
      }
      
      if (rho_ve2t) {
        ddt(rho) += VEt * vE_Grad(VEt,Nit);
      }
      
    }
  
    /////////////////////////////////////////////////////////////////////////////////////
    // AJPAR
    
    ddt(ajpar) = 0.0;
    if (evolve_ajpar) {

      if (ajpar_phi1) {
        ddt(ajpar) += (1./fmei)*Grad_par_LtoC(phi); // Right-hand deriv with b.c. Necessary for sheath mode
      }
      
      if (ajpar_jpar1) {
        ddt(ajpar) -= 0.51*nu*ajpar;
      }
      

      if (ajpar_te_ni) {
        ddt(ajpar) -= (1./fmei)*(Tet/Nit)*Grad_par_LtoC(ni);
      }
      
      if (ajpar_te) {
        ddt(ajpar) -= (1.71/fmei)*Grad_par_LtoC(te);
      }
      
      if (ajpar_ajpar1_phi0) {
        ddt(ajpar) -= vE_Grad(ajpar,phi0);
      }
      
      if (ajpar_ajpar1_phi1) {
        ddt(ajpar) -= vE_Grad(ajpar,phi);
      }
      
      if (ajpar_ve1_ve1) {
        ddt(ajpar) -= Vpar_Grad_par(ajpar,ajpar);
      }
      
    }
    
    
    /////////////////////////////////////////////////////////////////////////////////////////
    // TEMPERATURE EQUATION
    
    ddt(te) = 0.0;
    if (evolve_te) {

      if (te_te0_phi1) {
        ddt(te) -= vE_Grad(Te0, phi);
      }
      
      
      if (te_te1_phi0) {
        ddt(te) -= vE_Grad(te, phi0);
      }
      
      if (te_te1_phi1) {
        ddt(te) -= vE_Grad(te, phi);
      }
      
      if (te_ajpar_te) {
        ddt(te) -= ajpar * Grad_par_LtoC(te);
      }
      
      if (te_te_ajpar) {
        ddt(te) -= 2./3. * Tet * Grad_par_CtoL(ajpar);
      }
      
      if (te_nu_te1) {
        ddt(te) -= 2.*fmei*nu_hat/sqrt(Te0)*(ni - 1./2.*Ni0/Te0*te);  // Explicitly linear
      }
      
      if (te_nu_tet) {
        ddt(te) -= 2.*fmei*(nu*Tet -  nu_hat*Ni0/sqrt(Te0)); // All Nonlinear
      }

      if (te_jpar) {
        ddt(te) += 0.71*2./3. * Tet/Nit * Grad_par_CtoL(jpar);
      }
      
      if (te_diff) {
        ddt(te) += te_perpdiff * Delp2(te);
      }
      
      if (remove_tor_av_te) {
        ddt(te) -= DC(ddt(te)); // REMOVE TOROIDAL AVERAGE TEMPERATURE
      }
      
      if (evolve_source_te) {
        // Evolve source
        ddt(St) = averageY(-1. * source_alpha * DC(te) / Te0);
        
        // Add heat source/sink
        ddt(te) += St*where(St, Te0, Tet);
        
      }
      
      
      // There is an ion collision term that can be added with finite T_i
      
    }
    
    /////////////////////////////////////////////////////////////////////////////////////////
    // Z filtering
    if (filter_z) {
      // Filter out all except filter_z_mode
      
      ddt(rho) = filter(ddt(rho), filter_z_mode);
      ddt(ni) = filter(ddt(ni), filter_z_mode);
      ddt(ajpar) = filter(ddt(ajpar), filter_z_mode);
      ddt(te) = filter(ddt(te), filter_z_mode);
    }
    
    
    return 0;
  }
  //End of physics_run
  /////////////////////////////////////////////////////////////////
  
  
  
  /****************SPECIAL DIFFERENTIAL OPERATORS******************/
  const Field2D Perp_Grad_dot_Grad(const Field2D &p, const Field2D &f) {
    
    return DDX(p)*DDX(f)*mesh->getCoordinates()->g11;
  }
  
  
  /////////////////////////////////////////////////////////////////
  // ExB terms. These routines allow comparisons with BOUT-06
  // if bout_exb=true is set in BOUT.inp
  /////////////////////////////////////////////////////////////////
  const Field2D vE_Grad(const Field2D &f, const Field2D &p) {
    Field2D result;
    if (bout_exb) {
      // Use a subset of terms for comparison to BOUT-06
      result = 0.0;
    } else {
      // Use full expression with all terms
      
      result = b0xGrad_dot_Grad(p, f) / mesh->getCoordinates()->Bxy;
    }
    return result;
  }

  const Field3D vE_Grad(const Field2D &f, const Field3D &p) {
    Coordinates *coord = mesh->getCoordinates();
    Field3D result;
    if (arakawa) {
      // Arakawa scheme for perpendicular flow. Here as a test
      result = bracket(p, f, BRACKET_ARAKAWA);
    }else if(bout_exb) {
      // Use a subset of terms for comparison to BOUT-06
      result = bracket(p, f, BRACKET_SIMPLE);
    }else {
      // Use full expression with all terms
      result = bracket(p, f, BRACKET_STD);
    }
    return result;
  }
  
  const Field3D vE_Grad(const Field3D &f, const Field2D &p) {
    Field3D result;
    if (bout_exb) {
      // Use a subset of terms for comparison to BOUT-06
      result = VDDZ(-DDX(p), f);
    } else {
      // Use full expression with all terms
      result = b0xGrad_dot_Grad(p, f) / mesh->getCoordinates()->Bxy;
    }
    return result;
  }
  
  const Field3D vE_Grad(const Field3D &f, const Field3D &p) {
    Field3D result;
    
    Coordinates *coord = mesh->getCoordinates();
    if (arakawa) {
      // Arakawa scheme for perpendicular flow. Here as a test
      result = bracket(p, f, BRACKET_ARAKAWA);
    }else if(bout_exb) {
      // Use a subset of terms for comparison to BOUT-06
      result = VDDX(DDZ(p), f) + VDDZ(-DDX(p), f);
    }else {
      // Use full expression with all terms
      result = b0xGrad_dot_Grad(p, f) / coord->Bxy;
    }
    return result;
  }

};

BOUTMAIN(LAPDdrift);
