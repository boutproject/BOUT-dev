
#include <bout.hxx>
#include <boutmain.hxx>

#include <invert_laplace.hxx>
#include <math.h>

// Evolving quantities
Field3D Vort, Ajpar, Pe, Upar;


Field3D phi, apar, jpar;

Field2D B0, Pe0;
BoutReal C; // Collisional damping (resistivity)
BoutReal beta_hat, mu_hat, eps_hat;
BoutReal mu_par;

int phi_flags, apar_flags;
bool ZeroElMass, estatic; 

FieldGroup comms;

int physics_init(bool restarting) {
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

