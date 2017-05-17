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
  //2D initial profiles
  
  Field2D rho0, p0;
  Field2D Jpar0; //calculated from equilibrium B field used in bbmhd Jpar0=b.curlB0
  Vector2D B0_vec;
  Field2D B0;
  Field2D G; //grad G will give us the gravity paramater.

  //Initial perturbations
  // Field3D U0; //calculated from intial velocity perturbation used in bbmhd.
  Field3D Vpar0; //parallel component of intial velocity perturbation.
  // Field3D psi1;
  // Field3D rho1;
  // Field3D p1;
  Field3D phi0;
  
  //testing variables
  Field3D testa;
  Field3D testb;
  Field3D testc;
  Field3D testd;
  Field3D teste;
  Field3D testf;
  Field3D testg;
  Field3D testh;
  
  
  //3D evolving fields
  Field3D U, rho, p, Vpar, Psi;
  
  //Derived variables
  Field3D Jpar, phi;
  // Field3D gam_U, gam_Vpar, gam_p, gam_rho, gam_Psi;
  
  // Group of fields for communication
  FieldGroup comms;
  
  
  bool nonlinear;
  
  //metric coeffictients
  Coordinates *coord;
  
  
  //parameters
  BoutReal mu_0, Gamma;
  
  BoutReal viscos_par;  // Parallel viscosity
  BoutReal viscos_perp; // Perpendicular viscosity
  BoutReal hyperviscos; // Hyper-viscosity (radial)
  
  // Number which specifies the boundary condition on phi in the inversion
  int phi_flags; 
  
  BRACKET_METHOD bm = BRACKET_ARAKAWA;
  
  int init(bool restarting) override {

    output <<"Solving flute reduced MHD in a slab with gravity\n";
    
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
    //   GRID_LOAD(p1);
    //   output << "Loaded p1\n";
    //   GRID_LOAD(rho1);
    //   output << "Loaded rho1\n";
    //   GRID_LOAD(psi1);
    //   output << "Loaded psi1\n";
    
    //options stuff
    
    Options *globalOptions = Options::getRoot();
    Options *options = globalOptions->getSection("gravity");
    
    OPTION(options, nonlinear, false);
    
    if (nonlinear) {
      output <<"Solving WITH nonlinear terms\n";
    } else {
      output <<"Solving WITHOUT nonlinear terms\n";
    }
    
    OPTION(options, phi_flags, 0);
    phi.setBoundary("phi");
    
    options->get("viscos_par", viscos_par, 0.);
    options->get("viscos_perp", viscos_par, 0.);
    
    options->get("mu_0", mu_0, 1.);
    options->get("gamma", Gamma, 5./3.);
    
    //load metric tensor components

    coord = mesh->coordinates();
    
    BoutReal Lz; // Size of the Z box
    
    OPTION(options, Lz, 1.);
    
    // Set the metric tensor components to get Lz
    coord->g33 = SQ(2.*PI/Lz);
    coord->g_33 = 1. / coord->g33;
    
    /**************** SET EVOLVING VARIABLES *************/
    
    // Tell BOUT++ which variables to evolve
    // add evolving variables to the communication object
    
    SOLVE_FOR5(rho, p, U, Psi, Vpar);

    if (!restarting) {
      // Set initial perturbation
      //     U = U0;
      //     U = Delp2(phi0);
      U = coord->g11*D2DX2(phi0) + coord->g33*D2DZ2(phi0);
      testa = coord->g11*D2DX2(phi0);
      testb = coord->g33*D2DZ2(phi0);
      testc = coord->g33*D2DZ2(Jpar0);
      
      testd = coord->g33*D2DZ2(rho0);
    
      Vpar = Vpar0;
      //     p = p1;
      //     rho = rho1;
      //     Psi = psi1;
    }
    
    //******************Set up comms***************
    
    comms.add(rho, p, U, Psi, Vpar);
    
    //extra variables
    comms.add(phi);

    Jpar.setBoundary("jpar");
  
    
    //add variables to output file
    
    dump.add(phi, "phi", 1);
    dump.add(Jpar, "Jpar", 1);
    dump.add(G, "G", 0);
    // dump.add(U0, "U0", 0);
    dump.add(p0, "p0", 0);
    dump.add(rho0, "rho0", 0);
    dump.add(testa, "testa", 0);
    dump.add(testb, "testb", 0);
    dump.add(testa, "testc", 0);
    dump.add(testb, "testd", 0);
    dump.add(teste, "teste", 1);
    dump.add(testf, "testf", 1);
    dump.add(testg, "testg", 1);
    dump.add(testh, "testh", 1);
    //   dump.add(p1, "p1", 0);
    //   dump.add(rho1, "rho1", 0);
    //   dump.add(psi1, "psi1", 0);
    //   dump.add(gam_Psi, "gam_psi", 1);
    //   dump.add(gam_p, "gam_p", 1);
    //   dump.add(gam_rho, "gam_rho", 1);
    //   dump.add(gam_U, "gam_U", 1);
    //   dump.add(gam_Vpar, "gam_Vpar", 1);
    //dump.add(Vpar0, "Vpar0", 0);
    
    SAVE_REPEAT(ddt(Psi));
    SAVE_REPEAT(ddt(U));
    SAVE_REPEAT(ddt(rho));
    
    return(0);
  }
  
  int rhs(BoutReal t) override {
    //   U = Delp2(phi);
    phi = invert_laplace(U, phi_flags); // Invert Laplacian, setting boundary condition in phi_flags
    phi.applyBoundary();
    
    
    mesh->communicate(comms);
    
    Jpar = -(B0/mu_0)*Delp2(Psi);
    Jpar.applyBoundary();
    
    mesh->communicate(Jpar);
    
    //Parallel electric field
    ddt(Psi) = -(1/B0)*Grad_par_CtoL(B0*phi);// + 1e-2*Jpar; 
    
    if (nonlinear) {
      //ddt(Psi) +=-(1/B0)*(- b0xGrad_dot_Grad(Psi,B0*phi));
      ddt(Psi) += (1/B0)*bracket(Psi, B0*phi, bm)*coord->Bxy;
    }
    
    //Parallel vorticity
    
    ddt(U) = (SQ(B0)/rho0)*(Grad_par_LtoC(Jpar/B0) );
    teste = (SQ(B0)/rho0)*(Grad_par_LtoC(Jpar/B0) );
    
    //ddt(U) += -(1/rho0)*b0xGrad_dot_Grad(G,rho);
    ddt(U) -= (1/rho0)*bracket(G,rho, bm)*coord->Bxy;
    testf = (1/rho0)*bracket(G,rho, bm)*coord->Bxy;
    
    ddt(U) -= (SQ(B0)/rho0)*bracket(Psi,Jpar0/B0, bm)*coord->Bxy; ////added 02/03/2011 WAS MISSING BEFORE - check effect - see if stable. - still stable...
    testg = (SQ(B0)/rho0)*bracket(Psi,Jpar0/B0, bm)*coord->Bxy;

    if (nonlinear) {
      ddt(U) -= bracket(phi,U, bm)*coord->Bxy;
      
      ddt(U) -= (SQ(B0)/rho0)*bracket(Psi,Jpar/B0, bm)*coord->Bxy;
    }
    
    // Viscosity terms 
    if (viscos_par > 0.0)
      ddt(U) += viscos_par * Grad2_par2(U); // Parallel viscosity
    
    if (viscos_perp > 0.0)
      ddt(U) += viscos_perp * Delp2(U);     // Perpendicular viscosity
    
    //Parallel velocity
    //ddt(Vpar) =  + (b0xGrad_dot_Grad(Psi,p0))/rho0;
    ddt(Vpar) = bracket(Psi,p0, bm)*coord->Bxy / rho0;
    
    ddt(Vpar) += -(Grad_par_CtoL(p))/rho0;
    
    //ddt(Vpar) += b0xGrad_dot_Grad(G,Psi);
    ddt(Vpar) += bracket(G,Psi, bm)*coord->Bxy;
    
    if (nonlinear) {
      //ddt(Vpar) += -b0xGrad_dot_Grad(phi,Vpar);
      ddt(Vpar) -= bracket(phi,Vpar,bm)*coord->Bxy;
      
      //ddt(Vpar) += (b0xGrad_dot_Grad(Psi,p))/rho0;
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
    
    
    //gam_U = 1/U0; //ddt(U)/U0;
    //   gam_Vpar = 1/Vpar0; //ddt(Vpar)/Vpar0;
    //   gam_p = 1/p1; //ddt(p)/p1;
    //   gam_rho = 1/rho1; //ddt(rho)/rho1;
    //   gam_Psi = 1/psi1; //ddt(Psi)/psi1;
    
    return 0;
  }
};

BOUTMAIN(GravityReduced);
