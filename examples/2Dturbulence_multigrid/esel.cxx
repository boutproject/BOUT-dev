#include <bout/physicsmodel.hxx>
#include <smoothing.hxx>
#include <invert_laplace.hxx>
#include <initialprofiles.hxx>
#include <derivs.hxx>

class ESEL : public PhysicsModel {
private:
  Field3D n/*, T*/, vort;                   // Evolving density, temp and vorticity
  Field3D N; // ln(n)
  Field3D phi ; 
  Field2D B ;                           // Magnetic field
  BoutReal D/*, chi*/, mu ;                 // Diffusion coefficients
  Field2D sigma_n, sigma_T, sigma_vort; // dissipation terms 
  BoutReal zeta ;                       // rho/R0   
  BRACKET_METHOD bm;               // Bracket method for advection terms
  class Laplacian* phiSolver;      // Laplacian solver for vort -> phi
  bool test_laplacian; // If true, compute the error on the Laplacian inversion and abort
  Field3D vort_error;

protected:
  int init(bool restart) {
    
    Options *options = Options::getRoot()->getSection("esel");
    
    OPTION(options, zeta, 2.15e-3) ;
    OPTION(options, D, 1.97e-3) ;
    //   OPTION(options, chi, 4.61e-3) ;
    OPTION(options, mu, 3.88e-2) ; 
    int bracket; 
    OPTION(options, bracket, 2);
    OPTION(options, test_laplacian, false);
    
    // Set sources and sinks from input profile
    initial_profile("sigma_n", sigma_n);
    initial_profile("sigma_T", sigma_T);
    initial_profile("sigma_vort", sigma_vort);
    initial_profile("B", B) ;
    
    SAVE_ONCE(sigma_n) ; 
    
    // Poisson brackets: b_hat x Grad(f) dot Grad(g) / B = [f, g]
    // Method to use: BRACKET_ARAKAWA, BRACKET_STD or BRACKET_SIMPLE
    // Choose method to use for Poisson bracket advection terms
    
    switch(bracket) {
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
      output << "ERROR: Invalid choice of bracket method. Must be 0-3\n";
      return 1;
    }

    Coordinates *coord = mesh->getCoordinates();
    
    // generate coordinate system 
    // coord->J = R0/B0 ;
    coord->Bxy = 1 ;
    
    coord->g11 = 1.0 ;
    coord->g22 = 1.0 ;
    coord->g33 = 1.0 ;
    coord->g12 = 0.0 ;
    coord->g13 = 0.0 ;
    coord->g23 = 0.0 ;
    
    coord->g_11 = 1.0 ;
    coord->g_22 = 1.0 ;
    coord->g_33 = 1.0 ;
    coord->g_12 = 0.0;
    coord->g_13 = 0.0;
    coord->g_23 = 0.0;
    
    coord->geometry();
    
    //   SOLVE_FOR3(n, T, vort);
    //   SOLVE_FOR(n);
    SOLVE_FOR(N);
    SOLVE_FOR(vort);
    SAVE_REPEAT(phi);
    if (test_laplacian) {
      SAVE_REPEAT(vort_error);
    }
    phiSolver = Laplacian::create();
    phi = 0.0 ; // Starting phi
    
    return 0;
  }
  
  Field3D C(const Field3D &f) {
    return zeta * DDZ(f) ;
  }
  
  int rhs(BoutReal time) {
    //   output<<"\r"<<time<<std::flush;
    //   output<<"\r"<<time-1.e5<<std::flush;
    //  output<<time-1.e4<<endl;
    
    //   mesh->communicate(n/*, T*/, vort);
    mesh->communicate(N/*, T*/, vort);
    
    // Solve for potential
    //   phiSolver->setCoefC(n);
    phiSolver->setCoefC2(N);
    phi = phiSolver->solve(vort, phi);
    
    // Communicate variables
    mesh->communicate(phi);
    
    if (test_laplacian) {
      
      Field3D vort2 = D2DX2(phi) + D2DZ2(phi) + DDX(N)*DDX(phi) + DDZ(N)*DDZ(phi);
      vort_error = (vort-vort2);
      
      dump.write();
      
      MPI_Barrier(BoutComm::get());
      
      return 1; // Abort execution
    }
    
    // Continuity equation:  
    //   ddt(n) = bracket(phi,n,bm) + n*C(phi) - C(n/**T*/) + D*Delp2(n) - sigma_n*n ; 
    ddt(N) = bracket(phi,N,bm) + C(phi) - C(N/**T*/) + D*Delp2(N) - sigma_n ; 
    
    // Energy equation:
    //   ddt(T) = bracket(phi,T,bm) + (2/3)*T*C(phi) - (7/3)*T*C(T) - (2/3)*(T*T/n)*C(n) + chi*Delp2(T) - sigma_T*T ;
    
    // Vorticity equation:
    //   ddt(vort) = bracket(phi, vort, bm) - C(n/**T*/) + mu*Delp2(vort) - sigma_vort*vort ; 
    ddt(vort) = bracket(phi, vort, bm) - C(exp(N)/**T*/) + mu*Delp2(vort) - sigma_vort*vort ; 
    // ddt(vort) += phi/1.e4; // Sheath dissipation term to try to stop potential getting too large...
    
    // n.b bracket terms do not have minus sign before them because 
    // B is pointing in -ve y direction in BOUT coordinates.  
    //// This may be wrong, but it is probably consistently wrong
    
    return 0;
  }
};

BOUTMAIN(ESEL);
