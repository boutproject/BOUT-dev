/*
  Global Braginskii Solver (GBS) equations
  
  http://infoscience.epfl.ch/record/182434/files/PPCF2012.pdf

  NOTES:
  
  1. A slightly different definition of Poisson brackets is used in BOUT++:
     [f,g] = b.(Grad f cross Grad g) / B
     
  2. This version is collocated, and uses a
     4th-order smoother to stabilise the method
 */

#include "gbs.hxx"

#include <initialprofiles.hxx>
#include <derivs.hxx>
#include <field_factory.hxx>

int GBS::init(bool restarting) {
  Options *opt = Options::getRoot();
  coords = mesh->getCoordinates();

  // Switches in model section
  Options *optgbs = opt->getSection("GBS");
  OPTION(optgbs, ionvis, false);
  if(ionvis)
    OPTION(optgbs, Ti, 10); // Ion temperature [eV]
  OPTION(optgbs, elecvis, false);    // Include electron viscosity?
  OPTION(optgbs, resistivity, true); // Include resistivity?
  OPTION(optgbs, estatic, true);     // Electrostatic?
  
  // option for ExB Poisson Bracket 
  int bm_exb_flag;
  OPTION(optgbs, bm_exb_flag,         2); // Arakawa default
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

  OPTION(opt->getSection("solver"), mms, false);

  if(ionvis)
    SAVE_REPEAT(Gi);
  if(elecvis)
    SAVE_REPEAT(Ge);
  if(resistivity)
    SAVE_REPEAT(nu);

  // Normalisation
  OPTION(optgbs, Tnorm, 100);  // Reference temperature [eV]
  OPTION(optgbs, Nnorm, 1e19); // Reference density [m^-3]
  OPTION(optgbs, Bnorm, 1.0);  // Reference magnetic field [T]
  OPTION(optgbs, AA, 2.0);     // Ion mass

  output.write("Normalisation Te=%e, Ne=%e, B=%e\n", Tnorm, Nnorm, Bnorm);
  SAVE_ONCE4(Tnorm, Nnorm, Bnorm, AA); // Save

  Cs0      = sqrt(SI::qe*Tnorm / (AA*SI::Mp)); // Reference sound speed [m/s]
  Omega_ci = SI::qe*Bnorm / (AA*SI::Mp);       // Ion cyclotron frequency [1/s]
  rho_s0   = Cs0 / Omega_ci;
  
  mi_me  = AA*SI::Mp/SI::Me;
  beta_e = SI::qe*Tnorm*Nnorm / (SQ(Bnorm)/SI::mu0);
  
  output.write("\tmi_me=%e, beta_e=%e\n", mi_me, beta_e);
  SAVE_ONCE2(mi_me, beta_e);

  output.write("\t Cs=%e, rho_s=%e, Omega_ci=%e\n", Cs0, rho_s0, Omega_ci);
  SAVE_ONCE3(Cs0, rho_s0, Omega_ci);
  
  // Collision times
  BoutReal Coulomb = 6.6 - 0.5*log(Nnorm * 1e-20) + 1.5*log(Tnorm);
  tau_e0 = 1. / (2.91e-6 * (Nnorm / 1e6) * Coulomb * pow(Tnorm, -3./2));
  tau_i0 = sqrt(AA) / (4.80e-8 * (Nnorm / 1e6) * Coulomb * pow(Tnorm, -3./2));
  
  output.write("\ttau_e0=%e, tau_i0=%e\n", tau_e0, tau_i0);
  
  // Get switches from each variable section
  Options *optne = opt->getSection("Ne");
  optne->get("evolve", evolve_Ne, true);
  optne->get("D", Dn, 1e-3);
  optne->get("H", Hn, -1.0);
  if(mms) {
    Sn = 0.0;
  }else {
    string source;
    optne->get("source", source, "0.0");
    Sn = FieldFactory::get()->create3D(source, NULL, mesh);
    Sn /= Omega_ci;
    SAVE_ONCE(Sn);
  }
  
  Options *optte = opt->getSection("Te");
  optte->get("evolve", evolve_Te, true);
  optte->get("D", Dte, 1e-3);
  optte->get("H", Hte, -1.0);
  if(mms) {
    Sp = 0.0;
  }else {
    string source;
    optte->get("source", source, "0.0");
    Sp = FieldFactory::get()->create3D(source, NULL, mesh);
    Sp /= Omega_ci;
    SAVE_ONCE(Sp);
  }
  
  Options *optvort = opt->getSection("Vort");
  optvort->get("evolve", evolve_Vort, true);
  optvort->get("D", Dvort, 1e-3);
  optvort->get("H", Hvort, -1.0);

  Options *optve = opt->getSection("VePsi");
  optve->get("evolve", evolve_Ve, true);
  optve->get("D", Dve, 1e-3);
  optve->get("H", Hve, -1.0);

  Options *optvi = opt->getSection("Vi");
  optvi->get("evolve", evolve_Vi, true);
  optvi->get("D", Dvi, 1e-3);
  optvi->get("H", Hvi, -1.0);
  
  OPTION(optgbs, parallel, true);
  if(!parallel) {
    // No parallel dynamics
    evolve_Ve = false;
    evolve_Vi = false;
  }
  
  // If evolving, add to solver. Otherwise, set to an initial (fixed) value.
  if(evolve_Ne)    { solver->add(Ne, "Ne"); evars.add(Ne); }
       else { initial_profile("Ne", Ne); SAVE_ONCE(Ne); }
  if(evolve_Vort)  { solver->add(Vort, "Vort"); evars.add(Vort); }
       else { initial_profile("Vort", Vort); SAVE_ONCE(Vort); }
  if(evolve_Ve) { solver->add(VePsi, "VePsi"); evars.add(VePsi); }
       else { initial_profile("VePsi", VePsi); SAVE_ONCE(VePsi); }
  if(evolve_Vi)    { solver->add(Vi, "Vi"); evars.add(Vi); }
       else { initial_profile("Vi", Vi); SAVE_ONCE(Vi);}
  if(evolve_Te)    { solver->add(Te, "Te"); evars.add(Te); }
       else { initial_profile("Te", Te); SAVE_ONCE(Te);}
  
  if(evolve_Ve && (!estatic)) {
    SAVE_REPEAT2(psi, Ve);
  }
    
  // Load metric tensor from the mesh, passing length and B field normalisations
  LoadMetric(rho_s0, Bnorm);
  
  if(!restarting) {
    bool startprofiles;
    OPTION(optgbs, startprofiles, true);
    if(startprofiles) {
      // Read profiles from the mesh
      Field2D NeMesh; 
      if(mesh->get(NeMesh, "Ne0")) {
        output << "\nHERE\n";
        // No Ne0. Try Ni0
        if(mesh->get(NeMesh, "Ni0")) {
          output << "WARNING: Neither Ne0 nor Ni0 found in mesh input\n";
        }
      }
      NeMesh *= 1e20; // Convert to m^-3
      
      Field2D TeMesh;
      if(mesh->get(TeMesh, "Te0")) {
        // No Te0
        output << "WARNING: Te0 not found in mesh\n";
      }
      Field3D ViMesh; mesh->get(ViMesh, "Vi0");
      
      // Normalise
      NeMesh /= Nnorm;
      TeMesh /= Tnorm;
      ViMesh /= Cs0;
      
      // Add profiles in the mesh file
      Ne += NeMesh;
      Te += TeMesh;
      Vi += ViMesh;
    }
    
    // Check for negatives
    if(min(Ne, true) < 0.0) {
      throw BoutException("Starting density is negative");
    }
    if(max(Ne, true) < 1e-5) {
      throw BoutException("Starting density is too small");
    }
    if(min(Te, true) < 0.0) {
      throw BoutException("Starting temperature is negative");
    }
    if(max(Te, true) < 1e-5) {
      throw BoutException("Starting temperature is too small");
    }
  }
  
  SAVE_REPEAT(phi);
  
  phi.setBoundary("phi"); // For Y boundaries (if any)

  // Curvature
  OPTION(optgbs, curv_method, 1); // Get the method to use
  
  switch(curv_method) {
  case 0:   // bxcv vector, upwinding
  case 1: { // bxcv vector, central differencing
    bxcv.covariant = false; // Read contravariant components
    GRID_LOAD(bxcv);        // Specified components of b0 x kappa
    
    bool ShiftXderivs;
    Options::getRoot()->get("shiftXderivs", ShiftXderivs, false); // Read global flag
    if(ShiftXderivs) {
      Field2D sinty; GRID_LOAD(sinty);
      bxcv.z += sinty*bxcv.x;
    }
    // Normalise curvature
    bxcv.x /= Bnorm;
    bxcv.y /= rho_s0*rho_s0;
    bxcv.z *= rho_s0*rho_s0;
    break;
  }
  case 2: { // logB, read from input
    GRID_LOAD(logB);
    mesh->communicate(logB);
    SAVE_ONCE(logB);
    break;
  }
  case 3: { // logB, taken from mesh
    logB = log(coords->Bxy);
  }
  default:
    throw BoutException("Invalid value for curv_method");
  };
  
  
  // Phi solver
  phiSolver  = Laplacian::create(opt->getSection("phiSolver"));
  aparSolver = Laplacian::create(opt->getSection("aparSolver"));

  dx4 = SQ(SQ(coords->dx));
  dy4 = SQ(SQ(coords->dy));
  dz4 = SQ(SQ(coords->dz));

  SAVE_REPEAT(Ve);
  
  output.write("dx = %e, dy = %e, dz = %e\n", coords->dx(2,2), coords->dy(2,2), coords->dz);
  output.write("g11 = %e, g22 = %e, g33 = %e\n", coords->g11(2,2), coords->g22(2,2), coords->g33(2,2));
  output.write("g12 = %e, g23 = %e\n", coords->g12(2,2), coords->g23(2,2));
  output.write("g_11 = %e, g_22 = %e, g_33 = %e\n", coords->g_11(2,2), coords->g_22(2,2), coords->g_33(2,2));
  output.write("g_12 = %e, g_23 = %e\n", coords->g_12(2,2), coords->g_23(2,2));


  std::shared_ptr<FieldGenerator> gen = FieldFactory::get()->parse("source", Options::getRoot()->getSection("ne"));
  output << "Ne::source = " << gen->str() << endl;
  
  return 0;
}

void GBS::LoadMetric(BoutReal Lnorm, BoutReal Bnorm) {
  // Load metric coefficients from the mesh
  Field2D Rxy, Bpxy, Btxy, hthe, sinty;
  GRID_LOAD5(Rxy, Bpxy, Btxy, hthe, sinty); // Load metrics
  
  // Checking for dpsi and qinty used in BOUT grids
  Field2D dx;
  if(!mesh->get(dx,   "dpsi")) {
    output << "\tUsing dpsi as the x grid spacing\n";
    coords->dx = dx; // Only use dpsi if found
  }else {
    // dx will have been read already from the grid
    output << "\tUsing dx as the x grid spacing\n";
  }

  Rxy      /= Lnorm;
  hthe     /= Lnorm;
  sinty    *= SQ(Lnorm)*Bnorm;
  coords->dx /= SQ(Lnorm)*Bnorm;
  
  Bpxy /= Bnorm;
  Btxy /= Bnorm;
  coords->Bxy  /= Bnorm;
  
  // Calculate metric components
  bool ShiftXderivs;
  Options::getRoot()->get("shiftXderivs", ShiftXderivs, false); // Read global flag
  if(ShiftXderivs) {
    sinty = 0.0;  // I disappears from metric
  }
  
  BoutReal sbp = 1.0; // Sign of Bp
  if(min(Bpxy, true) < 0.0)
    sbp = -1.0;
  
  coords->g11 = SQ(Rxy*Bpxy);
  coords->g22 = 1.0 / SQ(hthe);
  coords->g33 = SQ(sinty)*coords->g11 + SQ(coords->Bxy)/coords->g11;
  coords->g12 = 0.0;
  coords->g13 = -sinty*coords->g11;
  coords->g23 = -sbp*Btxy/(hthe*Bpxy*Rxy);
  
  coords->J = hthe / Bpxy;
  
  coords->g_11 = 1.0/coords->g11 + SQ(sinty*Rxy);
  coords->g_22 = SQ(coords->Bxy*hthe/Bpxy);
  coords->g_33 = Rxy*Rxy;
  coords->g_12 = sbp*Btxy*hthe*sinty*Rxy/Bpxy;
  coords->g_13 = sinty*Rxy*Rxy;
  coords->g_23 = sbp*Btxy*hthe*Rxy/Bpxy;
  
  coords->geometry();
}

// just define a macro for V_E dot Grad
#define vE_Grad(f, p) ( bracket(p, f, bm_exb) )

int GBS::rhs(BoutReal t) {
  
  printf("TIME = %e\r", t); // Bypass logging, only to stdout
  
  // Communicate evolving variables
  mesh->communicate(evars);
  
  // Floor small values
  Te = floor(Te, 1e-3);
  Ne = floor(Ne, 1e-3);

  // Solve phi from Vorticity
  if(mms) {
    // Solve for potential, adding a source term
    Field3D phiS = FieldFactory::get()->create3D("phi:source", Options::getRoot(), mesh, CELL_CENTRE, t);
    phi = phiSolver->solve(Vort + phiS);
  }else {
    phi = phiSolver->solve(Vort);
  }

  if(estatic) {
    // Electrostatic
    Ve = VePsi;
    mesh->communicate(Ve);
  }else {
    aparSolver->setCoefA(-Ne*0.5*mi_me*beta_e);
    psi = aparSolver->solve(Ne*(Vi - VePsi));

    Ve = VePsi - 0.5*mi_me*beta_e*psi;
    mesh->communicate(psi, Ve);
  }

  // Communicate auxilliary variables
  mesh->communicate(phi);

  // Y boundary condition on phi
  phi.applyBoundary();

  // Stress tensor
  
  Field3D tau_e = Omega_ci*tau_e0 * pow(Te, 1.5)/Ne; // Normalised collision time

  Gi = 0.0;
  if(ionvis) {
    Field3D tau_i = Omega_ci*tau_i0 * pow(Ti, 1.5) / Ne;
    Gi = -(0.96*Ti*Ne*tau_i) * ( 2.*Grad_par(Vi) + C(phi)/coords->Bxy );
    mesh->communicate(Gi);
    Gi.applyBoundary("neumann");
  }else{
    mesh->communicate(Gi);
  }

  Field3D logNe = log(Ne);
  mesh->communicate(logNe);

  Ge = 0.0;
  if(elecvis) {
    Ge = -(0.73*Te*Ne*tau_e) * (2.*Grad_par(Ve) + (5.*C(Te) + 5.*Te*C(logNe) + C(phi))/coords->Bxy);
    mesh->communicate(Ge);
    Ge.applyBoundary("neumann");
  }else{
    mesh->communicate(Ge);
  }
  
  // Collisional damping (normalised)
  nu = 1. / (1.96 * Ne * tau_e * mi_me);
  
  Field3D Pe = Ne * Te;

  if(evolve_Ne) { 
    // Density
    ddt(Ne) = 
      - vE_Grad(Ne, phi)    // ExB term
      + (2./coords->Bxy) * (C(Pe) - Ne*C(phi)) // Perpendicular compression
      + D(Ne, Dn)
      + H(Ne, Hn)
      ;
    
    if(parallel) {
      ddt(Ne) -= Ne*Grad_par(Ve) + Vpar_Grad_par(Ve, Ne); // Parallel compression, advection
    }
    
    if(!mms) {
      // Source term
      ddt(Ne) += Sn*where(Sn, 1.0, Ne);
    }
  }
  
  if(evolve_Te) {
    // Electron temperature
    ddt(Te) = 
      - vE_Grad(Te, phi)
      + (4./3.)*(Te/coords->Bxy)*( (7./2.)*C(Te) + (Te/Ne)*C(Ne) - C(phi) )
      + D(Te, Dte)
      + H(Te, Hte)
      ;
    
    if(parallel) {
      ddt(Te) -= Vpar_Grad_par(Ve, Te);
      ddt(Te) += (2./3.)*Te*( 0.71*Grad_par(Vi) - 1.71*Grad_par(Ve) + 0.71*(Vi-Ve)*Grad_par(logNe));
    }
    
    if(!mms) {
      // Source term. Note: Power source, so divide by Ne
      ddt(Te) += Sp*where(Sp, 1.0, Te*Ne) / Ne;
    
      // Source of particles shouldn't include energy, so Ne*Te=const
      // hence ddt(Te) = -(Te/Ne)*ddt(Ne)
      ddt(Te) -= (Te/Ne)*Sn*where(Sn, 1.0, 0.0);
    }
  }
  
  if(evolve_Vort) {
    // Vorticity
    ddt(Vort) = 
      - vE_Grad(Vort, phi)     // ExB term
      + 2.*coords->Bxy*C(Pe)/Ne + coords->Bxy*C(Gi)/(3.*Ne)
      + D(Vort, Dvort)
      + H(Vort, Hvort)
      ;

    if(parallel) {
      Field3D delV = Vi-Ve;
      mesh->communicate(delV);
      ddt(Vort) -= Vpar_Grad_par(Vi, Vort); // Parallel advection
      ddt(Vort) += SQ(coords->Bxy)*( Grad_par(delV) + (Vi - Ve)*Grad_par(logNe) );
    }
  }
  
  if(evolve_Ve) {
    // Electron velocity
    
    ddt(VePsi) = 
      - vE_Grad(Ve, phi)
      - Vpar_Grad_par(Ve, Ve)
      - mi_me*(2./3.)*Grad_par(Ge)
      - mi_me*nu*(Ve - Vi)
      + mi_me*Grad_par(phi)
      - mi_me*(  Te*Grad_par(logNe) + 1.71*Grad_par(Te) )
      + D(Ve, Dve)
      + H(Ve, Hve)
      ;
  }
  
  if(evolve_Vi) {
    // Ion velocity
    
    ddt(Vi) = 
      - vE_Grad(Vi, phi)
      - Vpar_Grad_par(Vi, Vi)
      - (2./3.)*Grad_par(Gi)
      - (Grad_par(Te) + Te*Grad_par(logNe)) // Parallel pressure
      + D(Vi, Dvi)
      + H(Vi, Hvi)
      ;
  }
  
  return 0;
}

const Field3D GBS::C(const Field3D &f) { // Curvature operator
  Field3D g; //Temporary in case we need to communicate
  switch(curv_method) {
  case 0:
    g = f ; 
    mesh->communicate(g);
    return V_dot_Grad(bxcv, g);
  case 1:
    g = f ; 
    mesh->communicate(g);
    return bxcv*Grad(g);
  }
  return coords->Bxy*bracket(logB, f, BRACKET_ARAKAWA);
}

const Field3D GBS::D(const Field3D &f, BoutReal d) { // Diffusion operator
  if(d < 0.0)
    return 0.0;
  return d * Delp2(f);
}

const Field3D GBS::H(const Field3D &f, BoutReal h) { // Numerical hyper-diffusion operator
  if(h < 0.0)
    return 0.0;
  return -h*(dx4*D4DX4(f) + dz4*D4DZ4(f));// + dy4*D4DY4(f)
}


// Standard main() function
BOUTMAIN(GBS);

