/****************************************************************
 * DALF3 model
 * 
 * Four-field model for electron pressure, vorticity, A|| and
 * parallel velocity
 *
 * References:
 *
 *   B.Scott, Plasma Phys. Contr. Fusion 39 (1997) 1635
 *
 *   B.Scott, "Drift Wave versus Interchange Turbulence in
 *             Tokamak Geometry: Linear versus Nonlinear
 *             Mode Structure"
 *             arXiv:physics/0207126  Feb 2001
 *
 * NOTE: The normalisation used here is different to in the above
 *       papers. See manual in doc/ subdirectory for details
 *
 ****************************************************************/

#include <bout.hxx>
#include <boutmain.hxx>

#include <utils.hxx>
#include <invert_laplace.hxx>
#include <math.h>

// Constants
const BoutReal MU0 = 4.0e-7*PI;
const BoutReal Charge = 1.60217646e-19; // electron charge e (C)
const BoutReal Mi = 2.0*1.67262158e-27; // Ion mass
const BoutReal Me = 9.1093816e-31;  // Electron mass
const BoutReal Me_Mi = Me / Mi; // Electron mass / Ion mass

// Normalisation factors
BoutReal Tenorm, Nenorm, Bnorm;
BoutReal Cs, rho_s, wci;

// Poisson brackets: b0 x Grad(f) dot Grad(g) / B = [f, g]
// Method to use: BRACKET_ARAKAWA, BRACKET_STD or BRACKET_SIMPLE
BRACKET_METHOD bm; // Bracket method for advection terms

// Evolving quantities
Field3D Vort, Ajpar, Pe, Vpar;

Field3D phi, apar, jpar;

Field2D B0, Pe0, Jpar0, P0;
Vector2D b0xcv;

Field2D eta; // Collisional damping (resistivity)
BoutReal beta_hat, mu_hat;
BoutReal viscosity_par;

int phi_flags, apar_flags;
bool ZeroElMass, estatic; 
bool curv_kappa;
bool flat_resist;
BoutReal mul_resist;
bool parallel_lc;
bool nonlinear;
bool jpar_noderiv; // Don't take Delp2(apar) to get jpar
bool warm_ion;
bool vpar_advect;  // Include Vpar advection terms

bool filter_z;

BoutReal viscosity, hyper_viscosity;

bool smooth_separatrix;
int jpar_boundary;

FieldGroup comms;

int physics_init(bool restarting) {
  
  /////////////////////////////////////////////////////
  // Load data from the grid
  
  GRID_LOAD(Jpar0);
  SAVE_ONCE(Jpar0);

  Field2D Ni0, Te0;
  GRID_LOAD2(Ni0, Te0);
  Ni0 *= 1e20; // To m^-3
  Pe0 = Charge * Ni0 * Te0; // Electron pressure in Pascals
  P0 = 2.*Pe0;
  if(!warm_ion)
    Pe0 = P0; // No ion pressure
  SAVE_ONCE2(Pe0, P0);
  
  // Load curvature term
  b0xcv.covariant = false; // Read contravariant components
  mesh->get(b0xcv, "bxcv"); // mixed units x: T y: m^-2 z: m^-2

  // Metric coefficients
  Field2D Rxy, Bpxy, Btxy, hthe;
  Field2D I; // Shear factor
  
  if(mesh->get(Rxy,  "Rxy")) { // m
    output_error.write("ERROR: Cannot read Rxy from grid\n");
    return 1;
  }
  if(mesh->get(Bpxy, "Bpxy")) { // T
    output_error.write("ERROR: Cannot read Bpxy from grid\n");
    return 1;
  }
  mesh->get(Btxy, "Btxy"); // T
  mesh->get(B0,   "Bxy");  // T
  mesh->get(hthe, "hthe"); // m
  mesh->get(I,    "sinty");// m^-2 T^-1

  //////////////////////////////////////////////////////////////
  // Options

  Options *globalOptions = Options::getRoot();
  Options *options = globalOptions->getSection("dalf3");
  
  OPTION(options, phi_flags, 0);
  OPTION(options, apar_flags, 0);
  OPTION(options, estatic, false);
  OPTION(options, ZeroElMass, false);
  OPTION(options, jpar_noderiv, true);
  OPTION(options, curv_kappa, false);
  OPTION(options, flat_resist, false);
  OPTION(options, mul_resist, 1.0);
  OPTION(options, viscosity, -1.0);
  OPTION(options, hyper_viscosity, -1.0);
  OPTION(options, viscosity_par, -1.0);
  OPTION(options, smooth_separatrix, false);
  OPTION(options, warm_ion, false);
  OPTION(options, vpar_advect, false);
  OPTION(options, jpar_boundary, 5); 
  OPTION(options, filter_z, false);

  OPTION(options, parallel_lc, true);
  OPTION(options, nonlinear, true);
  
  int bracket_method;
  OPTION(options, bracket_method, 0);
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
 
  // SHIFTED RADIAL COORDINATES

  if(mesh->ShiftXderivs) {
    if(mesh->IncIntShear) {
      // BOUT-06 style, using d/dx = d/dpsi + I * d/dz
      mesh->IntShiftTorsion = I;
      
    }else {
      // Dimits style, using local coordinate system
      b0xcv.z += I*b0xcv.x;
      I = 0.0;  // I disappears from metric
    }
  }
  
  ///////////////////////////////////////////////////
  // Normalisation
  
  Tenorm = max(Te0, true);
  Nenorm = max(Ni0, true);
  Bnorm  = max(B0, true);

  // Sound speed in m/s
  Cs = sqrt(Charge*Tenorm / Mi);
  
  // drift scale
  rho_s = Cs * Mi / (Charge * Bnorm);
 
  // Ion cyclotron frequency
  wci = Charge * Bnorm / Mi;
 
  beta_hat = 4.e-7*PI * Charge*Tenorm * Nenorm / (Bnorm*Bnorm);
  
  if(ZeroElMass) {
    mu_hat = 0.;
  }else
    mu_hat = Me / Mi;
  
  SAVE_ONCE3(Tenorm, Nenorm, Bnorm);
  SAVE_ONCE3(Cs, rho_s, wci);
  SAVE_ONCE2(beta_hat, mu_hat);

  // Spitzer resistivity 
  if(flat_resist) {
    // eta in Ohm-m. NOTE: ln(Lambda) = 20
    eta = 0.51*1.03e-4*20.*pow(Tenorm, -1.5);
  }else {
    eta = 0.51*1.03e-4*20.*(Te0^(-1.5)); 
  }
  if(mul_resist < 0.0)
    mul_resist = 0.0;
  eta *= mul_resist;

  // Plasma quantities
  Jpar0 /= Nenorm*Charge*Cs;
  Pe0 /= Nenorm*Charge*Tenorm;

  // Coefficients
  eta *= Charge * Nenorm / Bnorm;

  viscosity /= wci*SQ(rho_s);
  hyper_viscosity /= wci*SQ(SQ(rho_s));
  viscosity_par /= wci*SQ(rho_s);

  b0xcv.x /= Bnorm;
  b0xcv.y *= rho_s*rho_s;
  b0xcv.z *= rho_s*rho_s;

  // Metrics
  Rxy /= rho_s;
  hthe /= rho_s;
  I *= rho_s*rho_s*Bnorm;
  Bpxy /= Bnorm;
  Btxy /= Bnorm;
  B0 /= Bnorm;

  mesh->dx /= rho_s*rho_s*Bnorm;
  
  ///////////////////////////////////////////////////
  // CALCULATE METRICS
  
  mesh->g11 = (Rxy*Bpxy)^2;
  mesh->g22 = 1.0 / (hthe^2);
  mesh->g33 = (I^2)*mesh->g11 + (B0^2)/mesh->g11;
  mesh->g12 = 0.0;
  mesh->g13 = -I*mesh->g11;
  mesh->g23 = -Btxy/(hthe*Bpxy*Rxy);
  
  mesh->J = hthe / Bpxy;
  mesh->Bxy = B0;
  
  mesh->g_11 = 1.0/mesh->g11 + ((I*Rxy)^2);
  mesh->g_22 = (B0*hthe/Bpxy)^2;
  mesh->g_33 = Rxy*Rxy;
  mesh->g_12 = Btxy*hthe*I*Rxy/Bpxy;
  mesh->g_13 = I*Rxy*Rxy;
  mesh->g_23 = Btxy*hthe*Rxy/Bpxy;
  
  mesh->geometry(); // Calculate quantities from metric tensor

  SOLVE_FOR3(Vort, Pe, Vpar);
  comms.add(Vort, Pe, Vpar);
  if(!(estatic && ZeroElMass)) {
    SOLVE_FOR(Ajpar);
    // Never differentiate Ajpar -> don't communicate
  }
  if(estatic) {
    comms.add(jpar);
  }else {
    // Need to communicate apar first then jpar
    comms.add(apar);
  }
  
  comms.add(phi);

  phi.setBoundary("phi");
  apar.setBoundary("apar");
  jpar.setBoundary("jpar");

  SAVE_REPEAT3(jpar, apar, phi);

  return 0;
}

// Curvature operator
const Field3D Kappa(const Field3D &f) {
  if(curv_kappa) {
    // Use the b0xcv vector from grid file
    return -2.*b0xcv*Grad(f) / B0;
  }
  
  return 2.*bracket(log(B0), f, bm);
}

const Field3D Grad_parP_LtoC(const Field3D &f) {
  Field3D result;
  if(parallel_lc) {
    result = Grad_par_LtoC(f);
    if(nonlinear)
      result -= beta_hat * bracket(apar, f, BRACKET_ARAKAWA);
  }else {
    if(nonlinear) {
      result = Grad_parP(apar*beta_hat, f);
    }else {
      result = Grad_par(f);
    }
  }
  return result;
}

const Field3D Grad_parP_CtoL(const Field3D &f) {
  Field3D result;
  if(parallel_lc) {
    result = Grad_par_CtoL(f);
    if(nonlinear)
      result -= beta_hat * bracket(apar, f, BRACKET_ARAKAWA);
  }else {
    if(nonlinear) {
      result = Grad_parP(apar*beta_hat, f);
    }else {
      result = Grad_par(f);
    }
  }
  return result;
}

int physics_run(BoutReal time) {

  // Invert vorticity to get electrostatic potential
  phi = invert_laplace(Vort*B0, phi_flags);
  if(warm_ion) {
    phi -= 0.5*Pe / mesh->Bxy; // Pi = Pe
  }
  
  phi.applyBoundary();

  // Calculate apar and jpar
  if(estatic) {
    // Electrostatic
    apar = 0.;
    if(ZeroElMass) {
      // Not evolving Ajpar
      jpar = Grad_par_CtoL(Pe - phi) / eta;
      jpar.applyBoundary();
    }else {
      jpar = Ajpar / mu_hat;
    }
    mesh->communicate(comms);
  }else {
    // Electromagnetic
    if(ZeroElMass) {
      // Ignore electron inertia term
      apar = Ajpar / beta_hat;
      
      mesh->communicate(comms);
      jpar = -Delp2(apar);
      jpar.applyBoundary();
      mesh->communicate(jpar);
    }else {
      // All terms - solve Helmholtz equation
      // ajpar = beta_hat*apar + mu_hat*jpar
      Field2D a = beta_hat;
      Field2D d = -mu_hat;
      apar = invert_laplace(Ajpar, apar_flags, &a, NULL, &d);
      apar.applyBoundary();
      
      mesh->communicate(comms);
      if(jpar_noderiv) {
        // Already applied boundaries on Ajpar and apar
        jpar = (Ajpar - beta_hat*apar) / mu_hat;
      }else {
        jpar = -Delp2(apar);
        jpar.applyBoundary();
        mesh->communicate(jpar);
      }
    }
  }

  Field3D Pet = Pe0;
  if(nonlinear) {
    Pet += Pe;
  }
  
  Field3D P = Pe;
  if(warm_ion)
    P *= 2.; // Ti = Te
  
  Field3D Ptot = P0;
  if(nonlinear)
    Ptot += P;
  
  if(jpar_boundary > 0) {
    // Boundary in jpar
    if(mesh->firstX()) {
      for(int i=jpar_boundary-1;i>=0;i--)
        for(int j=0;j<mesh->LocalNy;j++)
  	  for(int k=0;k<mesh->LocalNz;k++) {
            jpar[i][j][k] = 0.0; //0.5*jpar[i+1][j][k];
	  }
    }
    if(mesh->lastX()) {
      for(int i=mesh->LocalNx-jpar_boundary;i<mesh->LocalNx;i++)
        for(int j=0;j<mesh->LocalNy;j++)
  	for(int k=0;k<mesh->LocalNz;k++) {
            jpar[i][j][k] = 0.0; //0.5*jpar[i-1][j][k];
  	}
    }
  }
  
  // Vorticity equation
  ddt(Vort) = 
    B0*B0*Grad_parP_LtoC(jpar/B0)
    - B0*Kappa(P) // Total perturbed pressure
    ;
 
  if(nonlinear) {
    ddt(Vort) -= bracket(phi, Vort, bm);    // ExB advection
    
    if(vpar_advect)
      ddt(Vort) -= Vpar_Grad_par(Vpar, Vort);
  }

  if(viscosity > 0.0) {
    ddt(Vort) +=  viscosity * Delp2(Vort);
  }
  if(hyper_viscosity > 0.0) {
    Field3D delp2_vort = Delp2(Vort);
    delp2_vort.applyBoundary("neumann");
    mesh->communicate(delp2_vort);
    
    ddt(Vort) += hyper_viscosity*Delp2(delp2_vort);
  }

  if(filter_z)
    ddt(Vort) = filter(ddt(Vort), 1);
 
  // Parallel Ohm's law
  if(!(estatic && ZeroElMass)) {
    // beta_hat*apar + mu_hat*jpar
    ddt(Ajpar) =
      Grad_parP_CtoL(Pe - phi) // Electron pressure only
      - beta_hat * bracket(apar, Pe0, BRACKET_ARAKAWA)
      - eta*jpar
      ;

    if(nonlinear) {
      ddt(Ajpar) -= mu_hat*bracket(phi, jpar, bm);
    }

    if(filter_z)
      ddt(Ajpar) = filter(ddt(Ajpar), 1);
  }
  
  // Parallel velocity
  ddt(Vpar) = 
    - Grad_parP_CtoL(P) // Total pressure
    + beta_hat * bracket(apar, P0, BRACKET_ARAKAWA)
    ;
  
  if(nonlinear) {
    ddt(Vpar) -= bracket(phi, Vpar, bm);
    
    if(vpar_advect)
      ddt(Vpar) -= Vpar_Grad_par(Vpar, Vpar);
  }

  if(viscosity_par > 0.) {
    ddt(Vpar) += viscosity_par * Grad2_par2(Vpar);
  }

  if(filter_z)
    ddt(Vpar) = filter(ddt(Vpar), 1);

  // Electron pressure
  ddt(Pe) =
    - bracket(phi, Pet, bm)
    + Pet * (
             Kappa(phi - Pe)
             + B0*Grad_parP_LtoC( (jpar - Vpar)/B0 )
             )
    ;
  
  if(nonlinear) {
    if(vpar_advect)
      ddt(Vort) -= Vpar_Grad_par(Vpar, Pe);
  }

  if(smooth_separatrix) {
    // Experimental smoothing across separatrix
    ddt(Vort) += mesh->smoothSeparatrix(Vort);
  }

  if(filter_z)
    ddt(Pe) = filter(ddt(Pe), 1);

  // Boundary in Vpar and vorticity
  
  if(mesh->firstX()) {
    for(int i=3;i>=0;i--)
      for(int j=0;j<mesh->LocalNy;j++)
	for(int k=0;k<mesh->LocalNz;k++) {
          ddt(Vpar)[i][j][k] = ddt(Vpar)[i+1][j][k];
          ddt(Vort)[i][j][k] = ddt(Vort)[i+1][j][k];
	}
    
    // Subtract DC component
    for(int i=0;i<10;i++)
      for(int j=0;j<mesh->LocalNy;j++) {
        BoutReal avg = 0.;
        for(int k=0;k<mesh->LocalNz;k++)
          avg += ddt(Vort)[i][j][k];
        avg /= (BoutReal) mesh->LocalNz;
        for(int k=0;k<mesh->LocalNz;k++)
          ddt(Vort)[i][j][k] -= avg;
      }
  }
  if(mesh->lastX()) {
    for(int i=mesh->LocalNx-3;i<mesh->LocalNx;i++)
      for(int j=0;j<mesh->LocalNy;j++)
	for(int k=0;k<mesh->LocalNz;k++) {
          ddt(Vpar)[i][j][k] = ddt(Vpar)[i-1][j][k];
          ddt(Vort)[i][j][k] = ddt(Vort)[i-1][j][k];
	}
  }
  
  return 0;
}

