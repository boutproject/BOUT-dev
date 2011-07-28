
#include <bout.hxx>
#include <boutmain.hxx>

#include <invert_laplace.hxx>
#include <math.h>

// Constants
const BoutReal MU0 = 4.0e-7*PI;
const BoutReal Charge = 1.60217646e-19; // electron charge e (C)
const BoutReal Mi = 2.0*1.67262158e-27; // Ion mass
const BoutReal Me = 9.1093816e-31;  // Electron mass
const BoutReal Me_Mi = Me / Mi; // Electron mass / Ion mass

// Evolving quantities
Field3D Vort, Ajpar, Pe, Upar;

Field3D phi, apar, jpar;

Field2D B0, Pe0, Jpar0;
Vector2D b0xcv;

BoutReal C; // Collisional damping (resistivity)
BoutReal beta_hat, mu_hat, eps_hat;
BoutReal mu_par;

int phi_flags, apar_flags;
bool ZeroElMass, estatic; 

FieldGroup comms;

int physics_init(bool restarting) {
  
  //////////////////////////////////////////////////////////////
  // Load data from the grid
  
  GRID_LOAD(Jpar0);
  
  Field2D Ni0, Te0, Ti0;
  GRID_LOAD3(Ni0, Te0, Ti0);
  Pe0 = Charge * 1e20*Ni0 * (Te0 + Ti0); // Pressure in Pascals
  SAVE_ONCE(Pe0);
  
  // Load curvature term
  b0xcv.covariant = false; // Read contravariant components
  mesh->get(b0xcv, "bxcv"); // mixed units x: T y: m^-2 z: m^-2

  // Metric coefficients
  Field2D Rxy, Bpxy, Btxy, hthe;
  Field2D I; // Shear factor
  
  if(mesh->get(Rxy,  "Rxy")) { // m
    output.write("Error: Cannot read Rxy from grid\n");
    return 1;
  }
  if(mesh->get(Bpxy, "Bpxy")) { // T
    output.write("Error: Cannot read Bpxy from grid\n");
    return 1;
  }
  mesh->get(Btxy, "Btxy"); // T
  mesh->get(B0,   "Bxy");  // T
  mesh->get(hthe, "hthe"); // m
  mesh->get(I,    "sinty");// m^-2 T^-1

  //////////////////////////////////////////////////////////////
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

  //////////////////////////////////////////////////////////////
  // Options

  Options *globalOptions = Options::getRoot();
  Options *options = globalOptions->getSection("dalf3");
  
  OPTION(options, phi_flags, 0);
  OPTION(options, estatic, false);
  OPTION(options, ZeroElMass, false);
  
  if(ZeroElMass) {
    mu_hat = 0.;
  }

  SOLVE_FOR3(Vort, Pe, Upar);
  comms.add(Vort, Pe, Upar);
  if(!(estatic && ZeroElMass)) {
    SOLVE_FOR(Ajpar);
    comms.add(Ajpar);
  }
  comms.add(phi);

  phi.setBoundary("phi");
  apar.setBoundary("apar");
  jpar.setBoundary("jpar");
}

// Curvature operator
const Field3D Kappa(const Field3D &x) {
  
}

const Field3D Grad_parP_LtoC(const Field3D &f) {
  return Grad_par_LtoC(f);
}

const Field3D Grad_parP_CtoL(const Field3D &f) {
  return Grad_par_CtoL(f);
}

int physics_run(BoutReal time) {

  // Invert vorticity to get electrostatic potential
  phi = invert_laplace(Vort, phi_flags);
  phi.applyBoundary();
  
  // Communicate evolving variables and phi
  mesh->communicate(comms);

  // Calculate apar and jpar
  if(estatic) {
    // Electrostatic
    apar = 0.;
    if(ZeroElMass) {
      // Not evolving Ajpar
      jpar = Grad_par_CtoL(Pe - phi) / C;
      jpar.applyBoundary();
    }else {
      jpar = Ajpar / mu_hat;
    }
  }else {
    // Electromagnetic
    if(ZeroElMass) {
      // Ignore electron inertia term
      apar = Ajpar / beta_hat;
    }else {
      // All terms - solve Helmholtz equation
      Field2D a = beta_hat;
      Field2D d = -mu_hat;
      apar = invert_laplace(Ajpar, apar_flags, &a, NULL, &d);
      apar.applyBoundary();
    }
    jpar = -Delp2(apar);
    jpar.applyBoundary();
  }
  
  mesh->communicate(apar, jpar);

  // Vorticity equation
  ddt(Vort) = 
    - b0xGrad_dot_Grad(phi, Vort)    // ExB advection
    + (B0^3)*Grad_parP_LtoC(jpar/B0)
    - B0*B0*Kappa(Pe)
    ;
  
  if(!(estatic && ZeroElMass)) {
    // beta_hat*apar + mu_hat*jpar
    ddt(Ajpar) =
      - mu_hat*b0xGrad_dot_Grad(phi, jpar)
      + Grad_parP_CtoL(Pe0 + Pe - phi)
      - C*jpar
      ;
  }
  
  ddt(Pe) =
    - b0xGrad_dot_Grad(phi, Pe + Pe0)
    + B0*Grad_parP_LtoC( (jpar - Upar)/B0 )
    - Kappa(Pe - phi)
    ;
  
  ddt(Upar) = 
    - b0xGrad_dot_Grad(phi, Upar)
    - Grad_parP_CtoL(Pe + Pe0) / eps_hat
    + (mu_par/eps_hat) * Grad2_par2(Upar)
    ;
}

