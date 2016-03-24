
#include <bout/physicsmodel.hxx>
#include <bout/invert/laplacexy.hxx>
#include <bout/invert/laplacexz.hxx>
#include <invert_laplace.hxx>
#include <field_factory.hxx>

/// Fundamental constants
const BoutReal PI = 3.14159265;
const BoutReal e0  = 8.854e-12;      // Permittivity of free space
const BoutReal mu0 = 4.e-7*PI;       // Permeability of free space
const BoutReal qe  = 1.602e-19;      // Electron charge
const BoutReal Me  = 9.109e-31;      // Electron mass
const BoutReal Mp  = 1.67262158e-27; // Proton mass

class Alfven : public PhysicsModel {
private:
  Field3D Vort, Apar; // Evolving fields
  
  Field3D phi;  // Electrostatic potential
  Field3D jpar; // Parallel current
  
  Field2D phi2D; // Axisymmetric phi

  BoutReal Tnorm, Nnorm, Bnorm, AA; // Normalisation options
  BoutReal Cs0, Omega_ci, rho_s0, mi_me, beta_e;

  BoutReal mu_epar; // Electron parallel viscosity
  BoutReal resistivity;

  bool laplace_perp;    // Use Laplace_perp or Delp2?
  bool split_n0;        // Split solve into n=0 and n~=0?
  LaplaceXY *laplacexy; // Laplacian solver in X-Y (n=0)
  
  bool newXZsolver; 
  Laplacian *phiSolver; // Old Laplacian in X-Z
  LaplaceXZ *newSolver; // New Laplacian in X-Z
protected:
  
  int init(bool restarting) {
    
    // Normalisation
    Options *opt = Options::getRoot()->getSection("alfven");
    OPTION(opt, Tnorm, 100);  // Reference temperature [eV]
    OPTION(opt, Nnorm, 1e19); // Reference density [m^-3]
    OPTION(opt, Bnorm, 1.0);  // Reference magnetic field [T]
    OPTION(opt, AA, 2.0);     // Ion mass
  
    output.write("Normalisation Te=%e, Ne=%e, B=%e\n", Tnorm, Nnorm, Bnorm);
    SAVE_ONCE4(Tnorm, Nnorm, Bnorm, AA); // Save
    
    Cs0      = sqrt(qe*Tnorm / (AA*Mp)); // Reference sound speed [m/s]
    Omega_ci = qe*Bnorm / (AA*Mp);       // Ion cyclotron frequency [1/s]
    rho_s0   = Cs0 / Omega_ci;

    mi_me  = AA*Mp/Me;
    beta_e = qe*Tnorm*Nnorm / (SQ(Bnorm)/mu0);

    output.write("\tmi_me=%e, beta_e=%e\n", mi_me, beta_e);
    SAVE_ONCE2(mi_me, beta_e);
    
    output.write("\t Cs=%e, rho_s=%e, Omega_ci=%e\n", Cs0, rho_s0, Omega_ci);
    SAVE_ONCE3(Cs0, rho_s0, Omega_ci);
    
    OPTION(opt, mu_epar, -1e7); // Electron parallel viscosity [m^2/s]
    mu_epar /= rho_s0*rho_s0*Omega_ci * mi_me; // Normalise

    OPTION(opt, resistivity, 1e-7);

    // Load metric tensor from the mesh, passing length and B field normalisations
    LoadMetric(rho_s0, Bnorm);

    // Specify evolving variables
    SOLVE_FOR2(Vort, Apar);
   
    OPTION(opt, laplace_perp, true);  // Use Laplace_perp rather than Delp2
    OPTION(opt, split_n0, true); // Split into n=0 and n~=0

    if(split_n0) {
      // Create an XY solver for n=0 component
      laplacexy = new LaplaceXY(mesh);
      phi2D = 0.0; // Starting guess
    }
    
    // Create an XZ solver
    OPTION(opt, newXZsolver, false);
    if(newXZsolver) {
      // Test new LaplaceXZ solver
      newSolver = LaplaceXZ::create(mesh);
    }else {
      // Use older Laplacian solver
      phiSolver  = Laplacian::create();
    }
    phi = 0.0;

    // Set boundary condition on jpar from input file
    jpar.setBoundary("jpar");
    phi.setBoundary("phi");
    
    SAVE_REPEAT2(jpar,phi);
    
    //SAVE_REPEAT(phi2D);
    return 0;
  }

  int rhs(BoutReal time) {
    // Communicate evolving variables
    mesh->communicate(Vort, Apar);
    
    // Calculate parallel current from Apar
    if(laplace_perp) {
      jpar = Laplace_perp(Apar / (0.5*beta_e));
    }else {
      jpar = Delp2(Apar / (0.5*beta_e));
    }    

    // Apply boundary condition on jpar and phi
    jpar.applyBoundary(time);
    
    //Field2D Vort2D = Vort.DC(); // n=0 component
    //phi2D = laplacexy->solve(Vort2D, phi2D);

    // Calculate phi from potential
    if(split_n0) {
      // Split into axisymmetric and non-axisymmetric components
      Field2D Vort2D = Vort.DC(); // n=0 component
      phi2D = laplacexy->solve(Vort2D, phi2D);
      
      // Solve non-axisymmetric part using X-Z solver
      if(newXZsolver) {
        phi = newSolver->solve(Vort-Vort2D, phi);
      }else {
        phi = phiSolver->solve(Vort-Vort2D, phi);
      }
      phi.applyBoundary(time); // Apply Y boundaries
      phi += phi2D; // Add axisymmetric part
    }else {
      // Solve all components using X-Z solver
      if(newXZsolver) {
        // Use the new LaplaceXY solver
        phi = newSolver->solve(Vort, phi);
      }else {
        // Use older Laplacian solver
        phi = phiSolver->solve(Vort, phi);
      }
      phi.applyBoundary(time); // Apply Y boundaries
    }
    
    // Communicate jpar and phi
    mesh->communicate(jpar, phi);
    
    // Vorticity equation
    ddt(Vort) = Div_par(jpar);
    
    // Parallel electric field
    ddt(Apar) = Grad_par(phi);

    if(resistivity > 0.0) {
      ddt(Apar) += resistivity * jpar; 
    }   

    if(mu_epar > 0.0) {
      ddt(Apar) -= mu_epar * Laplace_par(jpar); // Parallel electron viscosity
    }
    
    return 0;
  }
  
  void LoadMetric(BoutReal Lnorm, BoutReal Bnorm) {
    // Load metric coefficients from the mesh
    Field2D Rxy, Bpxy, Btxy, hthe, sinty;
    GRID_LOAD5(Rxy, Bpxy, Btxy, hthe, sinty); // Load metrics
    
    // Checking for dpsi and qinty used in BOUT grids
    Field2D dx; 
    if(!mesh->get(dx,   "dpsi")) {
      output << "\tUsing dpsi as the x grid spacing\n";
      mesh->dx = dx; // Only use dpsi if found
    }else {
      // dx will have been read already from the grid
      output << "\tUsing dx as the x grid spacing\n";
    }
    Field2D qinty;
    if(!mesh->get(qinty, "qinty")) {
      output << "\tUsing qinty as the Z shift\n";
      mesh->zShift = qinty;
    }else {
      // Keep zShift
      output << "\tUsing zShift as the Z shift\n";
    }
    
    Rxy      /= Lnorm;
    hthe     /= Lnorm;
    sinty    *= SQ(Lnorm)*Bnorm;
    mesh->dx /= SQ(Lnorm)*Bnorm;
    
    Bpxy /= Bnorm;
    Btxy /= Bnorm;
    mesh->Bxy  /= Bnorm;
    
    // Calculate metric components
    if(mesh->ShiftXderivs) {
      sinty = 0.0;  // I disappears from metric
    }
    
    BoutReal sbp = 1.0; // Sign of Bp
    if(min(Bpxy, true) < 0.0)
      sbp = -1.0;
    
    mesh->g11 = (Rxy*Bpxy)^2;
    mesh->g22 = 1.0 / (hthe^2);
    mesh->g33 = (sinty^2)*mesh->g11 + (mesh->Bxy^2)/mesh->g11;
    mesh->g12 = 0.0;
    mesh->g13 = -sinty*mesh->g11;
    mesh->g23 = -sbp*Btxy/(hthe*Bpxy*Rxy);
    
    mesh->J = hthe / Bpxy;
    
    mesh->g_11 = 1.0/mesh->g11 + ((sinty*Rxy)^2);
    mesh->g_22 = (mesh->Bxy*hthe/Bpxy)^2;
    mesh->g_33 = Rxy*Rxy;
    mesh->g_12 = sbp*Btxy*hthe*sinty*Rxy/Bpxy;
    mesh->g_13 = sinty*Rxy*Rxy;
    mesh->g_23 = sbp*Btxy*hthe*Rxy/Bpxy;
    
    mesh->geometry();
  }
};

BOUTMAIN(Alfven);


