
#include <bout/physicsmodel.hxx>

#include <derivs.hxx>

#include <invert_laplace.hxx>
using bout::globals::mesh;

class DriftWave : public PhysicsModel { 
protected:
  int init(bool UNUSED(restart)) {
    // Specify evolving variables
    solver->add(Vort, "Vort"); // Vorticity
    solver->add(Ne, "Ne");     // Electron density
    
    // Get the normalised resistivity
    Options::getRoot()->getSection("drift")->get("nu", nu, 1.0);

    // Read background profile
    mesh->get(Ne0, "Ne0");
    
    // Split into convective and diffusive (stiff)
    setSplitOperator();

    // Solve potential as a constraint
    solver->constraint(phi, ddt(phi), "phi");
    phi = 0.0;
    
    // Coordinate system
    coord = mesh->getCoordinates();

    return 0;
  }
  
  int convective(BoutReal UNUSED(time)) {
    // Non-stiff parts of the problem here
    // Here just the nonlinear advection
    
    mesh->communicate(phi, Vort, Ne);
    
    // Linear advection
    ddt(Ne) = -bracket(phi, Ne0, BRACKET_ARAKAWA);
    
    // Non-linear advection of density and vorticity
    ddt(Ne) -= bracket(phi, Ne, BRACKET_ARAKAWA);
    
    ddt(Vort) = -bracket(phi, Vort, BRACKET_ARAKAWA);

    // potential is not evolved,
    // but solved as a constraint in the implicit stage
    ddt(phi) = 0.0;
    
    return 0;
  }
  
  int diffusive(BoutReal UNUSED(time)) {
    // Parallel dynamics treated implicitly
    
    mesh->communicate(phi, Vort, Ne);

    /*
      // This code results in ddt(Ne) depending on y+2, y-2
      // which are not (currently) included in the coloring
      // The result is that IMEX-BDF2 with coloring doesn't converge
      
    Ve = ( Grad_par(phi) - Grad_par(Ne) ) / nu;
    mesh->communicate(Ve);
    
    ddt(Ne) = -Div_par(Ve);
    ddt(Vort) = -Div_par(Ve);
    */
    
    // This version uses only y+1, y-1
    ddt(Ne) = Grad2_par2(Ne - phi)/nu;
    ddt(Vort) = ddt(Ne);
    
    // Calculate the constraint on potential
    // which is that ddt(phi) = 0
    
    phi.applyBoundary("dirichlet_o2");
    
    // This version uses FFTs, introducing dependencies in Z
    // ddt(phi) = Delp2(phi) - Vort;

    // This version uses central differencing for Delp2
    ddt(phi) = (coord->g11*D2DX2(phi) + coord->g33*D2DZ2(phi)) - Vort;
    
    return 0;
  }
  
private:
  Field2D Ne0; // Background density profile
  
  Field3D Vort, Ne; // Vorticity and density

  Field3D phi; // Electrostatic potential
  Field3D Ve;  // parallel electron velocity
  
  BoutReal nu; // Resistivity parameter
  
  Coordinates *coord; // Coordinate system metrics
};


BOUTMAIN(DriftWave);
