
#include <bout/physicsmodel.hxx>

#include <invert_laplace.hxx>

class DriftWave : public PhysicsModel { 
protected:
  int init(bool) {
    // Specify evolving variables
    solver->add(Vort, "Vort"); // Vorticity
    solver->add(Ne, "Ne");     // Electron density

    SAVE_REPEAT(phi);
    
    // Get the normalised resistivity
    nu = Options::root()["drift"]["nu"].doc("Normalised resistivity").withDefault(1.0);

    // Read background profile
    mesh->get(Ne0, "Ne0");

    setSplitOperator(); // Split into convective and diffusive (stiff)
    
    // Laplacian solver for potential
    phiSolver  = Laplacian::create();

    return 0;
  }
  
  int convective(BoutReal) {
    // Non-stiff parts of the problem here
    // Here just the nonlinear advection
    
    // Solve for potential
    phi = phiSolver->solve(Vort);
    mesh->communicate(phi, Vort, Ne);
    
    // Linear advection
    ddt(Ne) = -bracket(phi, Ne0, BRACKET_ARAKAWA);
    
    // Non-linear advection of density and vorticity
    ddt(Ne) -= bracket(phi, Ne, BRACKET_ARAKAWA);
    
    ddt(Vort) = -bracket(phi, Vort, BRACKET_ARAKAWA);
    
    return 0;
  }
  
  int diffusive(BoutReal) {
    // Parallel dynamics treated implicitly
    
    // Solve for potential
    phi = phiSolver->solve(Vort);
    mesh->communicate(phi, Vort, Ne);
    
    Ve = ( Grad_par(phi) - Grad_par(Ne) ) / nu;
    mesh->communicate(Ve);
    
    ddt(Ne) = -Div_par(Ve);
    ddt(Vort) = -Div_par(Ve);
    
    return 0;
  }
  
private:
  Field2D Ne0; // Background density profile
  
  Field3D Vort, Ne; // Vorticity and density

  Field3D phi; // Electrostatic potential
  Field3D Ve;  // parallel electron velocity
  
  Laplacian *phiSolver;

  BoutReal nu; // Resistivity parameter
};


BOUTMAIN(DriftWave);
