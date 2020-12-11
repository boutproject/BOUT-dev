
#include <bout/physicsmodel.hxx>
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
  
  bool newXZsolver; 
  std::unique_ptr<Laplacian> phiSolver{nullptr}; // Old Laplacian in X-Z
  std::unique_ptr<LaplaceXZ> newSolver{nullptr}; // New Laplacian in X-Z
protected:
  int init(bool) {
    // Normalisation
    auto opt = Options::root()["alfven"];
    Tnorm = opt["Tnorm"].withDefault(100);  // Reference temperature [eV]
    Nnorm = opt["Nnorm"].withDefault(1e19); // Reference density [m^-3]
    Bnorm = opt["Bnorm"].withDefault(1.0);  // Reference magnetic field [T]
    AA = opt["AA"].withDefault(2.0);        // Ion mass

    output.write("Normalisation Te={:e}, Ne={:e}, B={:e}\n", Tnorm, Nnorm, Bnorm);
    SAVE_ONCE4(Tnorm, Nnorm, Bnorm, AA); // Save
    
    Cs0      = sqrt(qe*Tnorm / (AA*Mp)); // Reference sound speed [m/s]
    Omega_ci = qe*Bnorm / (AA*Mp);       // Ion cyclotron frequency [1/s]
    rho_s0   = Cs0 / Omega_ci;

    mi_me  = AA*Mp/Me;
    beta_e = qe*Tnorm*Nnorm / (SQ(Bnorm)/mu0);

    output.write("\tmi_me={:e}, beta_e={:e}\n", mi_me, beta_e);
    SAVE_ONCE2(mi_me, beta_e);
    
    output.write("\t Cs={:e}, rho_s={:e}, Omega_ci={:e}\n", Cs0, rho_s0, Omega_ci);
    SAVE_ONCE3(Cs0, rho_s0, Omega_ci);

    mu_epar = opt["mu_epar"].withDefault(-1e7); // Electron parallel viscosity [m^2/s]
    mu_epar /= rho_s0*rho_s0*Omega_ci * mi_me; // Normalise

    resistivity = opt["resistivity"].withDefault(1e-7);

    // Load metric tensor from the mesh, passing length and B field normalisations
    LoadMetric(rho_s0, Bnorm);

    // Specify evolving variables
    SOLVE_FOR2(Vort, Apar);
   
    //////////////////////////////////////////
    // Solve potential as a constraint
    solver->constraint(phi, ddt(phi), "phi");
    phi = 0.0;
    
    // Specify the preconditioner function
    setPrecon( (preconfunc) &Alfven::precon );
    
    // Create an XZ solver
    newXZsolver = opt["newXZsolver"].withDefault(false);
    if(newXZsolver) {
      // Test new LaplaceXZ solver
      newSolver = LaplaceXZ::create(mesh);
    }else {
      // Use older Laplacian solver
      phiSolver  = Laplacian::create();
    }

    // Set boundary condition on jpar from input file
    jpar.setBoundary("jpar");
    phi.setBoundary("phi");
    
    SAVE_REPEAT(jpar);
    
    return 0;
  }

  int rhs(BoutReal time) {
    // Communicate evolving variables
    mesh->communicate(Vort, Apar);
    
    // Calculate parallel current from Apar
    jpar = Delp2(Apar / (0.5*beta_e));

    // Apply boundary condition on jpar and phi
    jpar.applyBoundary(time);
    
    // Communicate jpar and phi
    mesh->communicate(jpar, phi);

    // Calculate phi from potential as constraint
    phi.applyBoundary(time); // Apply boundaries
    //ddt(phi) = Laplace_perp(phi) - Vort;  // This constrained to be zero
    ddt(phi) = Delp2(phi) - Vort;

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
  
  /*!
   * Preconditioner. This inverts the operator (1 - gamma*J) 
   * where J is the Jacobian of the system
   *
   * The system state at time t is stored as usual
   * whilst the vector to be inverted is in ddt(f)
   * 
   * Inputs
   * ------
   *
   * t      = Current simulation time
   * gamma  = Coefficient proportional to timestep
   * delta  = Coefficient used in contrained problems
   * Vort,Apar,phi   = System state at current time
   * ddt(...) = Variables to be inverted
   *
   * Output
   * ------
   * 
   * ddt(f) = Result of the inversion
   */
  int precon(BoutReal, BoutReal, BoutReal) {
    if(newXZsolver) {
      ddt(phi) = newSolver->solve(ddt(phi) - ddt(Vort), 0.0);
    }else {
      ddt(phi) = phiSolver->solve(ddt(phi) - ddt(Vort), 0.0);
    }
    return 0;
  }

  void LoadMetric(BoutReal Lnorm, BoutReal Bnorm) {
    // Load metric coefficients from the mesh
    Field2D Rxy, Bpxy, Btxy, hthe, sinty;
    GRID_LOAD5(Rxy, Bpxy, Btxy, hthe, sinty); // Load metrics
    
    // Get the coordinates object
    Coordinates *coord = mesh->getCoordinates();

    // Checking for dpsi and qinty used in BOUT grids
    Field2D dx; 
    if(!mesh->get(dx,   "dpsi")) {
      output << "\tUsing dpsi as the x grid spacing\n";
      coord->dx = dx; // Only use dpsi if found
    }else {
      // dx will have been read already from the grid
      output << "\tUsing dx as the x grid spacing\n";
    }
    
    Rxy      /= Lnorm;
    hthe     /= Lnorm;
    sinty    *= SQ(Lnorm)*Bnorm;
    coord->dx /= SQ(Lnorm)*Bnorm;
    
    Bpxy /= Bnorm;
    Btxy /= Bnorm;
    coord->Bxy  /= Bnorm;
    
    // Check type of parallel transform
    std::string ptstr = Options::root()["mesh"]["paralleltransform"]["type"]
                                       .withDefault("identity");

    if(lowercase(ptstr) == "shifted") {
      // Using shifted metric method
      sinty = 0.0;  // I disappears from metric
    }
    
    BoutReal sbp = 1.0; // Sign of Bp
    if(min(Bpxy, true) < 0.0)
      sbp = -1.0;
    
    // Calculate metric components
    
    coord->g11 = SQ(Rxy*Bpxy);
    coord->g22 = 1.0 / SQ(hthe);
    coord->g33 = SQ(sinty)*coord->g11 + SQ(coord->Bxy)/coord->g11;
    coord->g12 = 0.0;
    coord->g13 = -sinty*coord->g11;
    coord->g23 = -sbp*Btxy/(hthe*Bpxy*Rxy);
    
    coord->J = hthe / Bpxy;
    
    coord->g_11 = 1.0/coord->g11 + SQ(sinty*Rxy);
    coord->g_22 = SQ(coord->Bxy*hthe/Bpxy);
    coord->g_33 = Rxy*Rxy;
    coord->g_12 = sbp*Btxy*hthe*sinty*Rxy/Bpxy;
    coord->g_13 = sinty*Rxy*Rxy;
    coord->g_23 = sbp*Btxy*hthe*Rxy/Bpxy;
    
    coord->geometry();
  }
};

BOUTMAIN(Alfven);


