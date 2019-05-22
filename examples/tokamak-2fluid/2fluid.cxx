/*******************************************************************************
 * 2-fluid turbulence model
 * This version intended to have inputs as similar to BOUT-06 as possible
 * for cross-benchmarking etc.
 *******************************************************************************/

#include <bout/physicsmodel.hxx>

#include <initialprofiles.hxx>
#include <derivs.hxx>
#include <interpolation.hxx>
#include <invert_laplace.hxx>
#include <msg_stack.hxx>

class TwoFluid : public PhysicsModel {
private:
  // 2D initial profiles
  Field2D Ni0, Ti0, Te0, Vi0, phi0, Ve0, rho0, Ajpar0;
  Vector2D b0xcv; // for curvature terms
  
  // 3D evolving fields
  Field3D rho, Te, Ni, Ajpar, Vi, Ti;
  
  // Derived 3D variables
  Field3D phi, Apar, Ve, jpar;
  
  // Non-linear coefficients
  Field3D nu, mu_i, kapa_Te, kapa_Ti;
  
  // 3D total values
  Field3D Nit, Tit, Tet, Vit;
  
  // pressures
  Field3D pei, pe;
  Field2D pei0, pe0;
  
  // Metric coefficients
  Field2D Rxy, Bpxy, Btxy, hthe;
  
  // parameters
  BoutReal Te_x, Ti_x, Ni_x, Vi_x, bmag, rho_s, fmei, AA, ZZ;
  BoutReal lambda_ei, lambda_ii;
  BoutReal nu_hat, mui_hat, wci, nueix, nuiix;
  BoutReal beta_p;
  BoutReal nuIonNeutral; // Ion-neutral collision rate (normalised by wci)
  
  // settings
  bool estatic, ZeroElMass; // Switch for electrostatic operation (true = no Apar)
  BoutReal zeff, nu_perp;
  bool evolve_rho, evolve_te, evolve_ni, evolve_ajpar, evolve_vi, evolve_ti;
  BoutReal ShearFactor;
  
  bool curv_upwind; // Use upwinding methods for curvature terms
  
  bool laplace_extra_rho_term; // An extra first-order term in the vorticity inversion
  bool vort_include_pi;    // Include Pi in vorticity
  
  bool bout_jpar; // Use BOUT-06 method for Jpar
  bool OhmPe;     // Include the Pe term in Ohm's law
  
  int bkgd;   // Profile options for coefficients (same options as BOUT-06)
  int iTe_dc; // Profile evolution options
  
  bool stagger; // Use CtoL and LtoC for parallel derivs
  
  // Poisson brackets: b0 x Grad(f) dot Grad(g) / B = [f, g]
  // Method to use: BRACKET_ARAKAWA, BRACKET_STD or BRACKET_SIMPLE
  BRACKET_METHOD bm; // Bracket method for advection terms
  
  // Switches for the equation terms
  bool ni_ni1_phi0, ni_ni0_phi1, ni_ni1_phi1, ni_nit_phit;
  bool ni_vi1_ni0,  ni_vi0_ni1,  ni_vi1_ni1,  ni_vit_nit;
  bool ni_jpar1, ni_pe1, ni_ni1;
  bool ni_ni0_curv_phi1, ni_ni1_curv_phi0, ni_ni1_curv_phi1, ni_nit_curv_phit;
  
  bool rho_rho0_phi1, rho_rho1_phi0, rho_rho1_phi1;
  bool rho_vi1_rho0, rho_vi0_rho1, rho_vi1_rho1;
  bool rho_pei1, rho_jpar1, rho_rho1;
  
  bool vi_vi0_phi1, vi_vi1_phi0, vi_vi1_phi1, vi_vit_phit;
  bool vi_vi1_vi0, vi_vi0_vi1, vi_vi1_vi1, vi_vit_vit;
  bool vi_pei1, vi_peit, vi_vi1;
  
  bool te_te1_phi0, te_te0_phi1, te_te1_phi1;
  
  bool ti_ti1_phi0, ti_ti0_phi1, ti_ti1_phi1;
  
  int lowPass_z; // Low-pass filter result
  
  int phi_flags, apar_flags; // Inversion flags
  
  // Group of objects for communications
  FieldGroup comms;

  // Coordinate system metrics
  Coordinates *coord;

  // Inverts a Laplacian to get potential
  Laplacian *phiSolver;

  // Solves the electromagnetic potential
  Laplacian *aparSolver;
  Field2D acoef; // Coefficient in the Helmholtz equation
  
  int init(bool UNUSED(restarting)) override {
    TRACE("int init(bool) ");
    
    Field2D I; // Shear factor 

    // Get the coordinate system
    coord = mesh->getCoordinates();
    
    output.write("Solving 6-variable 2-fluid equations\n");
    
    ////////////////////////////////////////////////////////
    // LOAD DATA FROM GRID FILE
    
    // Load 2D profiles (set to zero if not found)
    GRID_LOAD(Ni0);
    GRID_LOAD(Ti0);
    GRID_LOAD(Te0);
    GRID_LOAD(Vi0);
    GRID_LOAD(Ve0);
    GRID_LOAD(phi0);
    GRID_LOAD(rho0);
    GRID_LOAD(Ajpar0);
    
    // Load magnetic curvature term
    b0xcv.covariant = false; // Read contravariant components
    mesh->get(b0xcv, "bxcv"); // b0xkappa terms

    // Load metrics
    GRID_LOAD(Rxy);
    GRID_LOAD(Bpxy);
    GRID_LOAD(Btxy);
    GRID_LOAD(hthe);
    
    Field2D dx;
    if(!mesh->get(dx,   "dpsi")) {
      output << "\tUsing dpsi as the x grid spacing\n";
      coord->dx = dx; // Only use dpsi if found
    }else {
      // dx will have been read already from the grid
      output << "\tUsing dx as the x grid spacing\n";
    }
    mesh->get(I,    "sinty");
    
    // Load normalisation values
    GRID_LOAD(Te_x);
    GRID_LOAD(Ti_x);
    GRID_LOAD(Ni_x);
    GRID_LOAD(bmag);
    
    Ni_x *= 1.0e14;
    bmag *= 1.0e4;
    
    ////////////////////////////////////////////////////////
    // READ OPTIONS
    
    // Read some parameters
    auto globalOptions = Options::root();
    auto options = globalOptions["2fluid"];

    AA = options["AA"].withDefault(2.0);
    ZZ = options["ZZ"].withDefault(1.0);

    estatic = options["estatic"].withDefault(false);
    ZeroElMass = options["ZeroElMass"].withDefault(false);
    zeff = options["zeff"].withDefault(1.0);
    nu_perp = options["nu_perp"].withDefault(0.0);
    ShearFactor = options["ShearFactor"].withDefault(1.0);
    OhmPe = options["OhmPe"].withDefault(true);
    bout_jpar = options["bout_jpar"].withDefault(false);
    curv_upwind = options["curv_upwind"].withDefault(false);

    // Choose method to use for Poisson bracket advection terms
    int bracket_method;
    bracket_method = options["bracket_method"].withDefault(0);
    switch(bracket_method) {
    case 0: {
      bm = BRACKET_STD; 
      output << "\tBrackets: default differencing\n";
      break;
    }
    case 1: {
      bm = BRACKET_SIMPLE; 
      output << "\tBrackets: simplified operator\n";
      break;
    }
    case 2: {
      bm = BRACKET_ARAKAWA; 
      output << "\tBrackets: Arakawa scheme\n";
      break;
    }
    case 3: {
      bm = BRACKET_CTU; 
      output << "\tBrackets: Corner Transport Upwind method\n";
      break;
    }
    default:
      output << "ERROR: Invalid choice of bracket method. Must be 0 - 3\n";
      return 1;
    }

    nuIonNeutral = options["nuIonNeutral"].withDefault(-1.);

    bkgd = options["bkgd"].withDefault(2);
    iTe_dc = options["iTe_dc"].withDefault(2);

    stagger = options["stagger"].withDefault(false);

    laplace_extra_rho_term = options["laplace_extra_rho_term"].withDefault(false);
    vort_include_pi = options["vort_include_pi"].withDefault(false);

    lowPass_z = options["lowPass_z"].withDefault(-1);

    phi_flags = options["phi_flags"].withDefault(0);
    apar_flags = options["apar_flags"].withDefault(0);

    evolve_ni = globalOptions["Ni"]["evolve"].withDefault(true);
    evolve_rho = globalOptions["rho"]["evolve"].withDefault(true);
    evolve_vi = globalOptions["vi"]["evolve"].withDefault(true);
    evolve_ti = globalOptions["ti"]["evolve"].withDefault(true);
    evolve_ajpar = globalOptions["Ajpar"]["evolve"].withDefault(true);
    
    if (ZeroElMass) {
      evolve_ajpar = 0; // Don't need ajpar - calculated from ohm's law
    }
      
    ////////////////////////////////////////////////////////
    // Equation terms

    if (evolve_ni) {
      Options* options = &globalOptions["Ni"];
      options->get("ni1_phi0", ni_ni1_phi0, false);
      options->get("ni0_phi1", ni_ni0_phi1, false);
      options->get("ni1_phi1", ni_ni1_phi1, false);
      options->get("nit_phit", ni_nit_phit, false);
      options->get("vi1_ni0",  ni_vi1_ni0, false);
      options->get("vi0_ni1",  ni_vi0_ni1, false);
      options->get("vi1_ni1",  ni_vi1_ni1, false);
      options->get("vit_nit",  ni_vit_nit, false);
      options->get("jpar1",    ni_jpar1,  false);
      options->get("pe1",      ni_pe1,    false);
      options->get("ni0_curv_phi1", ni_ni0_curv_phi1, false);
      options->get("ni1_curv_phi0", ni_ni1_curv_phi0, false);
      options->get("ni1_curv_phi1", ni_ni1_curv_phi1, false);
      options->get("nit_curv_phit", ni_nit_curv_phit, false);
    }    
    
    if (evolve_rho) {
      Options* options = &globalOptions["rho"];
      options->get("rho0_phi1", rho_rho0_phi1, false);
      options->get("rho1_phi0", rho_rho1_phi0, false);
      options->get("rho1_phi1", rho_rho1_phi1, false);
      options->get("vi1_rho0",  rho_vi1_rho0, false);
      options->get("vi0_rho1",  rho_vi0_rho1, false);
      options->get("vi1_rho1",  rho_vi1_rho1, false);
      options->get("pei1",   rho_pei1, false);
      options->get("jpar1",  rho_jpar1, false);
      options->get("rho1",   rho_rho1, false);
    }
    
    if (evolve_vi) {
      Options* options = &globalOptions["vi"];
      options->get("vi0_phi1", vi_vi0_phi1, false);
      options->get("vi1_phi0", vi_vi1_phi0, false);
      options->get("vi1_phi1", vi_vi1_phi1, false);
      options->get("vit_phit", vi_vit_phit, false);
      options->get("vi1_vi0", vi_vi1_vi0, false);
      options->get("vi0_vi1", vi_vi0_vi1, false);
      options->get("vi1_vi1", vi_vi1_vi1, false);
      options->get("vit_vit", vi_vit_vit, false);
      options->get("pei1", vi_pei1, false);
      options->get("peit", vi_peit, false);
      options->get("vi1", vi_vi1, false);
    }
    
    if (evolve_te) {
      Options* options = &globalOptions["te"];
      options->get("te1_phi0", te_te1_phi0, false);
      options->get("te0_phi1", te_te0_phi1, false);
      options->get("te1_phi1", te_te1_phi1, false);
    }

    if (evolve_ti) {
      Options* options = &globalOptions["ti"];
      options->get("ti1_phi0", ti_ti1_phi0, false);
      options->get("ti0_phi1", ti_ti0_phi1, false);
      options->get("ti1_phi1", ti_ti1_phi1, false);
    }
    
    ////////////////////////////////////////////////////////
    // SHIFTED RADIAL COORDINATES

    // Check type of parallel transform
    std::string ptstr = Options::root()["mesh"]["paralleltransform"].withDefault<std::string>("identity");

    if (lowercase(ptstr) == "shifted") {
      ShearFactor = 0.0;  // I disappears from metric
      b0xcv.z += I*b0xcv.x;
    }

    ////////////////////////////////////////////////////////
    // CALCULATE PARAMETERS
    
    rho_s = 1.02*sqrt(AA*Te_x)/ZZ/bmag;
    fmei  = 1./1836.2/AA;
    
    lambda_ei = 24.-log(sqrt(Ni_x)/Te_x);
    lambda_ii = 23.-log(ZZ*ZZ*ZZ*sqrt(2.*Ni_x)/pow(Ti_x, 1.5));
    wci       = 9.58e3*ZZ*bmag/AA;
    nueix     = 2.91e-6*Ni_x*lambda_ei/pow(Te_x, 1.5);
    nuiix     = 4.78e-8*pow(ZZ,4.)*Ni_x*lambda_ii/pow(Ti_x, 1.5)/sqrt(AA);
    nu_hat    = zeff*nueix/wci;
    
    if (nu_perp < 1.e-10) {
      mui_hat      = (3./10.)*nuiix/wci;
    } else {
      mui_hat      = nu_perp;
    }
    
    if (estatic) {
      beta_p    = 1.e-29;
    } else {
      beta_p    = 4.03e-11*Ni_x*Te_x/bmag/bmag;
    }
    
    Vi_x = wci * rho_s;
    
    output.write("Collisions: nueix = %e, nu_hat = %e\n", nueix, nu_hat);
    
    ////////////////////////////////////////////////////////
    // PRINT Z INFORMATION
    
    BoutReal hthe0;
    if(mesh->get(hthe0, "hthe0") == 0) {
      output.write("    ****NOTE: input from BOUT, Z length needs to be divided by %e\n", hthe0/rho_s);
    }
    
    ////////////////////////////////////////////////////////
    // SHIFTED GRIDS LOCATION
    
    // Velocities defined on cell boundaries
    Vi.setLocation(CELL_YLOW);
    Ajpar.setLocation(CELL_YLOW);
    
    // Apar and jpar too
    Apar.setLocation(CELL_YLOW); 
    jpar.setLocation(CELL_YLOW);
    
    ////////////////////////////////////////////////////////
    // NORMALISE QUANTITIES
    
    output.write("\tNormalising to rho_s = %e\n", rho_s);
    
    // Normalise profiles
    Ni0 /= Ni_x/1.0e14;
    Ti0 /= Te_x;
    Te0 /= Te_x;
    phi0 /= Te_x;
    Vi0 /= Vi_x;
    
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
    
    ////////////////////////////////////////////////////////
    // CALCULATE METRICS
    
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
    
    ////////////////////////////////////////////////////////
    // SET EVOLVING VARIABLES
    
    // Tell BOUT++ which variables to evolve
    // add evolving variables to the communication object
    if (evolve_rho) {
      SOLVE_FOR(rho);
      comms.add(rho);
      output.write("rho\n");
    } else {
      initial_profile("rho", rho);
      rho.setBoundary("rho");
      rho.applyBoundary();
    }
    
    if (evolve_ni) {
      SOLVE_FOR(Ni);
      comms.add(Ni);
      output.write("ni\n");
    } else {
      initial_profile("Ni", Ni);
      Ni.setBoundary("Ni");
      Ni.applyBoundary();
    }
    
    if (evolve_te) {
      SOLVE_FOR(Te);
      comms.add(Te);
      output.write("te\n");
    } else {
      initial_profile("Te", Te);
      Te.setBoundary("Te");
      Te.applyBoundary();
    }
    
    if (evolve_ajpar) {
      SOLVE_FOR(Ajpar);
      comms.add(Ajpar);
      output.write("ajpar\n");
    }else {
      initial_profile("Ajpar", Ajpar);
      if (ZeroElMass) {
        dump.add(Ajpar, "Ajpar", 1); // output calculated Ajpar
      }
      Ajpar.setBoundary("Ajpar");
      Ajpar.applyBoundary();
    }
    
    if (evolve_vi) {
      SOLVE_FOR(Vi);
      comms.add(Vi);
      output.write("vi\n");
    } else {
      initial_profile("Vi", Vi);
      Vi.setBoundary("Vi");
      Vi.applyBoundary();
    }
    
    if (evolve_ti) {
      SOLVE_FOR(Ti);
      comms.add(Ti);
      output.write("ti\n");
    } else {
      initial_profile("Ti", Ti);
      Ti.setBoundary("Ti");
      Ti.applyBoundary();
    }
    
    jpar.setBoundary("jpar");
    
    ////////////////////////////////////////////////////////
    // SETUP COMMUNICATIONS
    
    // add extra variables to communication
    comms.add(phi);
    comms.add(Apar);
    
    // Add any other variables to be dumped to file
    dump.addRepeat(phi,  "phi");
    dump.addRepeat(Apar, "Apar");
    dump.addRepeat(jpar, "jpar");
    
    dump.addOnce(Ni0, "Ni0");
    dump.addOnce(Te0, "Te0");
    dump.addOnce(Ti0, "Ti0");
    
    dump.addOnce(Te_x,  "Te_x");
    dump.addOnce(Ti_x,  "Ti_x");
    dump.addOnce(Ni_x,  "Ni_x");
    dump.addOnce(rho_s, "rho_s");
    dump.addOnce(wci,   "wci");
    
    // Create a solver for the Laplacian
    phiSolver = Laplacian::create();
    phiSolver->setFlags(phi_flags);
    if (laplace_extra_rho_term) {
      // Include the first order term Grad_perp Ni dot Grad_perp phi
      phiSolver->setCoefC(Ni0);
    }

    if (! (estatic || ZeroElMass)) {
      // Create a solver for the electromagnetic potential
      aparSolver = Laplacian::create();
      aparSolver->setFlags(apar_flags);
      acoef = (-0.5 * beta_p / fmei) * Ni0;
      aparSolver->setCoefA(acoef);
    }
    
    return 0;
  }
  
  // ExB terms using Poisson bracket
#define vE_Grad(f, p) ( bracket(p, f, bm) )

  int rhs(BoutReal UNUSED(t)) override {
  
    ////////////////////////////////////////////////////////
    // Invert vorticity to get phi
    //
    // Solves \nabla^2_\perp x + \nabla_perp c\cdot\nabla_\perp x + a x = b
    // Arguments are:   (b,   bit-field, a,    c)
    // Passing NULL -> missing term
    
    {
      TRACE("Solving for phi");
      
      phi = phiSolver->solve(rho/Ni0);
      
      if (vort_include_pi) {
        // Include Pi term in vorticity
        phi -= (Ti*Ni0 + Ni*Te0) / Ni0;
      }
    }

    ////////////////////////////////////////////////////////
    // Invert Ajpar to get Apar

    {
      TRACE("Solving for Apar");
      if (estatic || ZeroElMass) {
        // Electrostatic operation
        Apar = 0.0;
      } else {
        Apar = aparSolver->solve(acoef*(Vi - Ajpar));
      }
    }
    
    ////////////////////////////////////////////////////////
    // Communicate variables
    mesh->communicate(comms);

    ////////////////////////////////////////////////////////
    // Update profiles for calculating nu, mu_i, kapa_Te,i
    switch(bkgd) {
    case 0: { // Toroidal averages
      Nit = Ni0 + DC(Ni);
      Tit = Ti0 + DC(Ti);
      Tet = Te0 + DC(Te);
      Vit = Vi0 + DC(Vi);
      break;
    }
    case 1: { // Full perturbed values
      Nit = Ni0 + Ni;
      Tit = Ti0 + Ti;
      Tet = Te0 + Te;
      Vit = Vi0 + Vi;
      break; 
    }
    case 2: { // Unperturbed values
      Nit = Ni0;
      Tit = Ti0;
      Tet = Te0;
      Vit = Vi0;
      break;
    }
    default: {
      throw BoutException("ERROR: Invalid bkgd option\n");
    }
    }
    
    ////////////////////////////////////////////////////////
    // Update non-linear coefficients on the mesh
    nu      = nu_hat * Nit / pow(Tet,1.5);
    mu_i    = mui_hat * Nit / sqrt(Tit);
    kapa_Te = 3.2*(1./fmei)*(wci/nueix)*pow(Tet,2.5);
    kapa_Ti = 3.9*(wci/nuiix)*pow(Tit,2.5);
  
    // note: nonlinear terms are not here
    pei = (Te0+Ti0)*Ni + (Te + Ti)*Ni0;
    pe  = Te0*Ni + Te*Ni0;
    
    ////////////////////////////////////////////////////////
    if (ZeroElMass) {
      // Set jpar,Ve,Ajpar neglecting the electron inertia term
      // Calculate Jpar, communicating across processors
      if (!stagger) {
        jpar = -(Ni0*Grad_par(phi, CELL_YLOW)) / (fmei*0.51*nu);
        
        if (OhmPe) {
          jpar += (Te0*Grad_par(Ni, CELL_YLOW)) / (fmei*0.51*nu);
        }
      } else {
        jpar = -(Ni0*Grad_par_LtoC(phi))/(fmei*0.51*nu);
      
        if (OhmPe) {
          jpar += (Te0*Grad_par_LtoC(Ni)) / (fmei*0.51*nu);
        }
      }
      
      // Need to communicate jpar
      mesh->communicate(jpar);
      jpar.applyBoundary();
      
      Ve = Vi - jpar/Ni0;
      Ajpar = Ve;
    } else {
    
      Ve = Ajpar + Apar;
      jpar = Ni0*(Vi - Ve);
    }
    
    ////////////////////////////////////////////////////////
    // DENSITY EQUATION
    
    ddt(Ni) = 0.0;
    if (evolve_ni) {
      TRACE("Density equation");
      
      if (ni_ni1_phi0)
        ddt(Ni) -= vE_Grad(Ni, phi0);
      
      if (ni_ni0_phi1) 
        ddt(Ni) -= vE_Grad(Ni0, phi);
      
      if (ni_ni1_phi1)
        ddt(Ni) -= vE_Grad(Ni, phi);
      
      if (ni_nit_phit)
        ddt(Ni) -= vE_Grad(Nit, phi0 + phi) - vE_Grad(Ni0, phi0);
      
      if (ni_vi1_ni0)
        ddt(Ni) -= Vpar_Grad_par(Vi, Ni0);
      
      if (ni_vi0_ni1)
        ddt(Ni) -= Vpar_Grad_par(Vi0, Ni);
      
      if (ni_vi1_ni1)
        ddt(Ni) -= Vpar_Grad_par(Vi, Ni);
      
      if (ni_vit_nit)
        ddt(Ni) -= Vpar_Grad_par(Vit, Nit) - Vpar_Grad_par(Vi0, Ni0);
      
      if (ni_jpar1) {
        if (stagger) {
          ddt(Ni) += Div_par_CtoL(jpar);
        } else {
          ddt(Ni) += Div_par(jpar);
        }
      }
      
      if (ni_pe1)
        ddt(Ni) += 2.0*V_dot_Grad(b0xcv, pe);
      
      if (ni_ni0_curv_phi1)
        ddt(Ni) -= 2.0*Ni0*V_dot_Grad(b0xcv, phi);
    
      if (ni_ni1_curv_phi0)
        ddt(Ni) -= 2.0*Ni*V_dot_Grad(b0xcv, phi0);
    
      if (ni_ni1_curv_phi1)
        ddt(Ni) -= 2.0*Ni*V_dot_Grad(b0xcv, phi);

      if (ni_nit_curv_phit)
        ddt(Ni) -= 2.0*Nit*V_dot_Grad(b0xcv, phi+phi0) - 2.0*Ni0*V_dot_Grad(b0xcv, phi0);
    
      if (ni_ni1)
        ddt(Ni) += mu_i * Delp2(Ni);
      
      //ddt(Ni) -= Ni0*Div_par(Vi) + Ni*Div_par(Vi0) + Ni*Div_par(Vi);
      
      if (lowPass_z > 0)
        ddt(Ni) = lowPass(ddt(Ni), lowPass_z);
    }

    ////////////////////////////////////////////////////////
    // ION VELOCITY
    
    ddt(Vi) = 0.0;
    if (evolve_vi) {
      TRACE("Ion velocity equation");
      
      if (vi_vi0_phi1)
        ddt(Vi) -= vE_Grad(Vi0, phi);
      
      if (vi_vi1_phi0)
        ddt(Vi) -= vE_Grad(Vi, phi0);

      if (vi_vi1_phi1)
        ddt(Vi) -= vE_Grad(Vi, phi);
    
      if (vi_vit_phit)
        ddt(Vi) -= vE_Grad(Vit, phi+phi0) - vE_Grad(Vi0, phi+phi0);
    
      if (vi_vi1_vi0)
        ddt(Vi) -= Vpar_Grad_par(Vi0, Vi);

      if (vi_vi0_vi1)
        ddt(Vi) -= Vpar_Grad_par(Vi, Vi0);
    
      if (vi_vi1_vi1)
        ddt(Vi) -= Vpar_Grad_par(Vi, Vi);
      
      if (vi_vit_vit)
        ddt(Vi) -= Vpar_Grad_par(Vit, Vit) - Vpar_Grad_par(Vi0, Vi0);
      
      if (vi_pei1)
        ddt(Vi) -= Grad_par(pei)/Ni0;

      if (vi_peit)
        ddt(Vi) -= Grad_par(pei)/Nit;
      
      if (vi_vi1)
        ddt(Vi) -= mu_i*Delp2(Vi);

      if (lowPass_z > 0)
        ddt(Vi) = lowPass(ddt(Vi), lowPass_z);
    }
  
    ////////////////////////////////////////////////////////
    // ELECTRON TEMPERATURE
    
    ddt(Te) = 0.0;
    if (evolve_te) {
      TRACE("Electron temperature equation");
      
      if (te_te1_phi0)
        ddt(Te) -= vE_Grad(Te, phi0);
      if (te_te0_phi1)
        ddt(Te) -= vE_Grad(Te0, phi);
      if (te_te1_phi1)
        ddt(Te) -= vE_Grad(Te, phi);
    
      /*
        ddt(Te) -= vE_Grad(Te0, phi) + vE_Grad(Te, phi0) + vE_Grad(Te, phi);
        ddt(Te) -= Vpar_Grad_par(Ve, Te0) + Vpar_Grad_par(Ve0, Te) + Vpar_Grad_par(Ve, Te);
        ddt(Te) += 1.333*Te0*( V_dot_Grad(b0xcv, pe)/Ni0 - V_dot_Grad(b0xcv, phi) );
        ddt(Te) += 3.333*Te0*V_dot_Grad(b0xcv, Te);
        ddt(Te) += (0.6666667/Ni0)*Div_par_K_Grad_par(kapa_Te, Te);
        
      */
      if (lowPass_z > 0)
        ddt(Te) = lowPass(ddt(Te), lowPass_z);
    }
    
    ////////////////////////////////////////////////////////
    // ION TEMPERATURE
    
    ddt(Ti) = 0.0;
    if (evolve_ti) {
      TRACE("Ion temperature equation");
      
      if (ti_ti1_phi0)
        ddt(Ti) -= vE_Grad(Ti, phi0);
      if (ti_ti0_phi1)
        ddt(Ti) -= vE_Grad(Ti0, phi);
      if (ti_ti1_phi1)
        ddt(Ti) -= vE_Grad(Ti, phi);
      
      /*
        ddt(Ti) -= vE_Grad(Ti0, phi) + vE_Grad(Ti, phi0) + vE_Grad(Ti, phi);
        ddt(Ti) -= Vpar_Grad_par(Vi, Ti0) + Vpar_Grad_par(Vi0, Ti) + Vpar_Grad_par(Vi, Ti);
        ddt(Ti) += 1.333*( Ti0*V_dot_Grad(b0xcv, pe)/Ni0 - Ti*V_dot_Grad(b0xcv, phi) );
        ddt(Ti) -= 3.333*Ti0*V_dot_Grad(b0xcv, Ti);
        ddt(Ti) += (0.6666667/Ni0)*Div_par_K_Grad_par(kapa_Ti, Ti);
      */
      
      if (lowPass_z > 0)
        ddt(Ti) = lowPass(ddt(Ti), lowPass_z);
    }
  
    ////////////////////////////////////////////////////////
    // VORTICITY
    
    ddt(rho) = 0.0;
    if (evolve_rho) {
      TRACE("Vorticity equation");
      
      if (rho_rho0_phi1)
        ddt(rho) -= vE_Grad(rho0, phi);
      
      if (rho_rho1_phi0)
        ddt(rho) -= vE_Grad(rho, phi0);
      
      if (rho_rho1_phi1)
        ddt(rho) -= vE_Grad(rho, phi);
      
      if (rho_vi1_rho0)
        ddt(rho) -= Vpar_Grad_par(Vi, rho0);
      
      if (rho_vi0_rho1)
        ddt(rho) -= Vpar_Grad_par(Vi0, rho);
      
      if (rho_vi1_rho1)
        ddt(rho) -= Vpar_Grad_par(Vi, rho);
      
      if (rho_pei1) {
        if (curv_upwind) {
          ddt(rho) += 2.0*coord->Bxy*V_dot_Grad(b0xcv, pei);  // Use upwinding
        } else {
          ddt(rho) += 2.0*coord->Bxy*b0xcv*Grad(pei);     // Use central differencing
        }
      }    
      
      if (rho_jpar1) {
        if (stagger) {
          ddt(rho) += SQ(coord->Bxy)*Div_par_CtoL(jpar);
        } else  {
          ddt(rho) += SQ(coord->Bxy)*Div_par(jpar, CELL_CENTRE);
        }
      }
      
      if (rho_rho1)
        ddt(rho) += mu_i * Delp2(rho);
      
      if (lowPass_z > 0)
        ddt(rho) = lowPass(ddt(rho), lowPass_z);
    }
  
    ////////////////////////////////////////////////////////
    // AJPAR
    
    ddt(Ajpar) = 0.0;
    if (evolve_ajpar) {
      TRACE("Ajpar equation");
      
      //ddt(Ajpar) -= vE_Grad(Ajpar0, phi) + vE_Grad(Ajpar, phi0) + vE_Grad(Ajpar, phi);
      //ddt(Ajpar) -= (1./fmei)*1.71*Grad_par(Te);
      
      if (stagger) {
        ddt(Ajpar) += (1./fmei)*Grad_par_LtoC(phi); // Right-hand differencing
      } else {
        ddt(Ajpar) += (1./fmei)*Grad_par(phi, CELL_YLOW);
      }
      
      if (OhmPe) {
        if (stagger) {
          ddt(Ajpar) -= (1./fmei)*(Tet/Nit)*Grad_par_LtoC(Ni);
        } else {
          ddt(Ajpar) -= (1./fmei)*(Te0/Ni0)*Grad_par(Ni, CELL_YLOW);
        }
      }
      
      ddt(Ajpar) += 0.51*interp_to(nu, CELL_YLOW)*jpar/Ni0;
      
      if(lowPass_z > 0)
        ddt(Ajpar) = lowPass(ddt(Ajpar), lowPass_z);
    }
    
    ////////////////////////////////////////////////////////
    // Profile evolution options
    
    switch(iTe_dc) {
    case 1: { // subtacting out toroidal averages for all fields
      if (evolve_ni)
        ddt(Ni) -= DC(ddt(Ni));
      if (evolve_rho)
        ddt(rho) -= DC(ddt(rho));
      if (evolve_te)
        ddt(Te) -= DC(ddt(Te));
      if (evolve_ti)
        ddt(Ti) -= DC(ddt(Ti));
      if (evolve_ajpar)
        ddt(Ajpar) -= DC(ddt(Ajpar));
      break;
    }
    case 2: { // not subtacting out toroidal averages for any field
      break;
    }
    case 4: { // using toroidal averages in right-hand sides, e.g., axisymmetric mode
      if (evolve_ni)
        ddt(Ni) = DC(ddt(Ni));
      if (evolve_rho)
        ddt(rho) = DC(ddt(rho));
      if (evolve_te)
        ddt(Te) = DC(ddt(Te));
      if (evolve_ti)
        ddt(Ti) = DC(ddt(Ti));
      if (evolve_ajpar)
        ddt(Ajpar) = DC(ddt(Ajpar));
      break;
    }
    default: {
      throw BoutException("ERROR: invalid option for iTe_dc\n");
    }
    }
    
    return(0);
  }
};

BOUTMAIN(TwoFluid);

