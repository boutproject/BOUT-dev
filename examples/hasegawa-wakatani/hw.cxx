
#include <bout/physicsmodel.hxx>
#include <bout/smoothing.hxx>
#include <bout/invert_laplace.hxx>
#include <bout/derivs.hxx>

class HW : public PhysicsModel {
private:
  Field3D n, vort;  // Evolving density and vorticity
  Field3D phi;      // Electrostatic potential

  // Model parameters
  BoutReal alpha;      // Adiabaticity (~conductivity)
  BoutReal kappa;      // Density gradient drive
  BoutReal Dvort, Dn;  // Diffusion 
  bool modified; // Modified H-W equations?
  
  // Poisson brackets: b0 x Grad(f) dot Grad(g) / B = [f, g]
  // Method to use: BRACKET_ARAKAWA, BRACKET_STD or BRACKET_SIMPLE
  BRACKET_METHOD bm; // Bracket method for advection terms
  
  class Laplacian* phiSolver; // Laplacian solver for vort -> phi

  // Simple implementation of 4th order perpendicular Laplacian
  Field3D Delp4(const Field3D &var) {
    Field3D tmp;
    tmp = Delp2(var, 0.0);
    mesh->communicate(tmp);
    tmp.applyBoundary("neumann");
    return Delp2(tmp, 0.0);
    
    //return Delp2(var);
  }
  
protected:
  int init(bool restart) {
  
    Options *options = Options::getRoot()->getSection("hw");
    OPTION(options, alpha, 1.0);
    OPTION(options, kappa, 0.1);
    OPTION(options, Dvort, 1e-2);
    OPTION(options, Dn,    1e-2);
  
    OPTION(options, modified, false);

    SOLVE_FOR2(n, vort);
    SAVE_REPEAT(phi);

    // Split into convective and diffusive parts
    setSplitOperator();
    
    phiSolver = Laplacian::create();
    phi = 0.; // Starting phi
    
    // Use default flags 
    
    // Choose method to use for Poisson bracket advection terms
    int bracket;
    OPTION(options, bracket, 0);
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
      output << "ERROR: Invalid choice of bracket method. Must be 0 - 3\n";
      return 1;
    }
    
    return 0;
  }

  int convective(BoutReal time) {
    // Non-stiff, convective part of the problem
    
    // Solve for potential
    phi = phiSolver->solve(vort, phi);
    
    // Communicate variables
    mesh->communicate(n, vort, phi);
    
    // Modified H-W equations, with zonal component subtracted from resistive coupling term
    Field3D nonzonal_n = n;
    Field3D nonzonal_phi = phi;
    if(modified) {
      // Subtract average in Y and Z
      nonzonal_n -= averageY(DC(n));
      nonzonal_phi -= averageY(DC(phi));
    }
    
    ddt(n) = -bracket(phi, n, bm) + alpha*(nonzonal_phi - nonzonal_n) - kappa*DDZ(phi);
    
    ddt(vort) = -bracket(phi, vort, bm) + alpha*(nonzonal_phi - nonzonal_n);
  
    return 0;
  }
  
  int diffusive(BoutReal time) {
    // Diffusive terms
    mesh->communicate(n, vort);
    ddt(n) = -Dn*Delp4(n);
    ddt(vort) = -Dvort*Delp4(vort);
    return 0;
  }
};

// Define a main() function
BOUTMAIN(HW);
