/*******************************************************************************
 * Flute-Reduced MHD - including gravity term instead of curvature
 * Basically the same as Hazeltine-Meiss but different normalisations and have gravity intead of curvature.
 * Evolving Vorticity U, Parallel electric field Psi, Parallel velocity Vpar, Pressure p, and density rho.
 * Have included compressional terms in Vpar and in pressure and density evolution equations.
 *******************************************************************************/


#include <bout/physicsmodel.hxx>

#include <initialprofiles.hxx>
#include <invert_laplace.hxx>
#include <derivs.hxx>

const BoutReal PI = 3.14159265;

class GravityReduced : public PhysicsModel {
private:
  // 2D initial profiles
  
  Field2D rho0, p0;
  Field2D Jpar0; //calculated from equilibrium B field used in bbmhd Jpar0=b.curlB0
  Vector2D B0_vec;
  Field2D B0;
  Field2D G; //grad G will give us the gravity paramater.

  // Initial perturbations
  // Field3D U0; //calculated from intial velocity perturbation used in bbmhd.
  Field3D Vpar0; //parallel component of intial velocity perturbation.
  Field3D phi0;
  
  //3D evolving fields
  Field3D U, rho, p, Vpar, Psi;
  
  //Derived variables
  Field3D Jpar, phi;
  
  // Group of fields for communication
  FieldGroup comms;
  
  bool nonlinear;
  
  // metric coeffictients
  Coordinates *coord;
  
  // parameters
  BoutReal mu_0, Gamma;
  
  BoutReal viscos_par;  // Parallel viscosity
  BoutReal viscos_perp; // Perpendicular viscosity
  BoutReal hyperviscos; // Hyper-viscosity (radial)
  
  BRACKET_METHOD bm = BRACKET_ARAKAWA;

  /// Solver for inverting Laplacian
  std::unique_ptr<Laplacian> phiSolver{nullptr};
  
  int init(bool restarting) override {

    output << "Solving flute reduced MHD in a slab with gravity\n";
    
    //*************** LOAD DATE FROM GRID FILE ********************
  

    //   GRID_LOAD(U0);
    //   output << "Loaded U0\n";
    GRID_LOAD(Vpar0);
    output << "Loaded Vpar0\n";
    GRID_LOAD(rho0);
    output << "Loaded rho0\n";
    GRID_LOAD(p0);
    output << "Loaded p0\n";
    GRID_LOAD(Jpar0);
    output << "Loaded Jpar0\n";
    GRID_LOAD(G);
    G *= 1000.;
    output << "Loaded Gravity\n";
    GRID_LOAD(B0_vec);
    output << "Loaded B0_vec\n";
    GRID_LOAD(B0);
    output << "Loaded B0\n";
    
    GRID_LOAD(phi0);
    output << "Loaded phi0\n";
    
    // options stuff
    
    auto globalOptions = Options::root();
    auto options = globalOptions["gravity"];

    nonlinear = options["nonlinear"].withDefault(false);

    if (nonlinear) {
      output <<"Solving WITH nonlinear terms\n";
    } else {
      output <<"Solving WITHOUT nonlinear terms\n";
    }

    phi.setBoundary("phi");

    viscos_par = options["viscos_par"].withDefault(0.);
    viscos_perp = options["viscos_perp"].withDefault(0.);
    mu_0 = options["mu_0"].withDefault(1.);
    Gamma = options["Gamma"].withDefault(5. / 3.);

    // load metric tensor components

    coord = mesh->getCoordinates();
    
    BoutReal Lz; // Size of the Z box

    Lz = options["Lz"].withDefault(1.);

    // Set the metric tensor components to get Lz
    coord->g33 = SQ(2.*PI/Lz);
    coord->g_33 = 1. / coord->g33;
    
    /**************** SET EVOLVING VARIABLES *************/
    
    // Tell BOUT++ which variables to evolve
    // add evolving variables to the communication object
    
    SOLVE_FOR(rho, p, U, Psi, Vpar);

    if (!restarting) {
      // Set initial perturbation
      //     U = U0;
      //     U = Delp2(phi0);
      U = coord->g11*D2DX2(phi0) + coord->g33*D2DZ2(phi0);
      Vpar = Vpar0;
    }
    
    //******************Set up comms***************
    
    comms.add(rho, p, U, Psi, Vpar);
    
    // extra variables
    comms.add(phi);

    Jpar.setBoundary("jpar");
    
    // Add variables to output file
    SAVE_REPEAT(phi, Jpar);  // Save every output
    SAVE_ONCE(G, p0, rho0);

    // Save time derivatives
    SAVE_REPEAT(ddt(Psi));
    SAVE_REPEAT(ddt(U));
    SAVE_REPEAT(ddt(rho));

    // Create a solver for the Laplacian
    phiSolver = Laplacian::create();
    
    return 0;
  }
  
  int rhs(BoutReal UNUSED(t)) override {
    //   U = Delp2(phi);
    phi = phiSolver->solve(U); // Invert Laplacian
    phi.applyBoundary(); // Apply boundary condition in Y
    
    mesh->communicate(comms);
    
    Jpar = -(B0/mu_0)*Delp2(Psi);
    Jpar.applyBoundary();
    
    mesh->communicate(Jpar);
    
    //Parallel electric field
    ddt(Psi) = -(1/B0)*Grad_par_CtoL(B0*phi);// + 1e-2*Jpar; 
    
    if (nonlinear) {
      ddt(Psi) += (1/B0)*bracket(Psi, B0*phi, bm)*coord->Bxy;
    }
    
    //Parallel vorticity
    
    ddt(U) = (SQ(B0)/rho0)*(Grad_par_LtoC(Jpar/B0) );
    
    ddt(U) -= (1/rho0)*bracket(G,rho, bm)*coord->Bxy;
    
    ddt(U) -= (SQ(B0)/rho0)*bracket(Psi,Jpar0/B0, bm)*coord->Bxy;

    if (nonlinear) {
      ddt(U) -= bracket(phi,U, bm)*coord->Bxy;
      
      ddt(U) -= (SQ(B0)/rho0)*bracket(Psi,Jpar/B0, bm)*coord->Bxy;
    }
    
    // Viscosity terms 
    if (viscos_par > 0.0) {
      ddt(U) += viscos_par * Grad2_par2(U); // Parallel viscosity
    }
    
    if (viscos_perp > 0.0) {
      ddt(U) += viscos_perp * Delp2(U);     // Perpendicular viscosity
    }
    
    // Parallel velocity
    ddt(Vpar) = bracket(Psi,p0, bm)*coord->Bxy / rho0;
    
    ddt(Vpar) += -(Grad_par_CtoL(p))/rho0;
    
    ddt(Vpar) += bracket(G,Psi, bm)*coord->Bxy;
    
    if (nonlinear) {
      ddt(Vpar) -= bracket(phi,Vpar,bm)*coord->Bxy;
      
      ddt(Vpar) += bracket(Psi,p,bm)*coord->Bxy / rho0;
    }
    
    //Pressure
    ddt(p) = -bracket(phi,p0,bm);
    
    ddt(p) += -((Gamma*p0)/(1 + Gamma*p0*mu_0/SQ(B0)))*( (rho0*mu_0/SQ(B0))*bracket(G,phi,bm)*coord->Bxy + Grad_par_LtoC(Vpar)  - (Vpar/B0)*Grad_par(B0) );
    
    if (nonlinear) {
      ddt(p) -= bracket(phi,p,bm)*coord->Bxy;
      ddt(p) += ((Gamma*p0) / (1 + Gamma*p0*mu_0/SQ(B0))) * bracket(Psi, Vpar, bm)*coord->Bxy;
    }
    
    //Density
    ddt(rho) = -bracket(phi, rho0, bm)*coord->Bxy;
    
    ddt(rho) -= (rho0/(1 + Gamma*p0*mu_0/SQ(B0)))*( (rho0*mu_0/SQ(B0))*bracket(G,phi,bm)*coord->Bxy + Grad_par_LtoC(Vpar) - bracket(Psi,Vpar,bm)*coord->Bxy - (Vpar/B0)*Grad_par(B0) );
    
    if (nonlinear) {
      ddt(rho) -= bracket(phi, rho, bm)*coord->Bxy;
      ddt(rho) += ((rho0)/(1 + Gamma*p0*mu_0/SQ(B0)))*bracket(Psi, Vpar, bm)*coord->Bxy;
    }
    
    // Iterate over the lower Y boundary
    RangeIterator rlow = mesh->iterateBndryLowerY();
    for(rlow.first(); !rlow.isDone(); rlow.next()) {
      int x = rlow.ind;
      for(int y=2;y>=0;y--) 
        for(int z=0;z<mesh->LocalNz;z++) {
          ddt(rho)(x,y,z) = ddt(rho)(x,y+1,z);
          ddt(p)(x,y,z) = ddt(p)(x,y+1,z);
          ddt(Psi)(x,y,z) = ddt(Psi)(x,y+1,z);
        }
    }
    
    return 0;
  }
};

BOUTMAIN(GravityReduced);
