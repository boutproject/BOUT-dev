/*******************************************************************
 * Advection-Diffusion Example
 *
 * MVU 19-july-2011
 *******************************************************************/

#include <bout/physicsmodel.hxx>
#include <derivs.hxx>

class AdvDiff : public PhysicsModel {
private:
  // Evolving variables 
  Field3D V;

protected:
  int init(bool restarting) {
    // 2D initial profiles
    Field2D V0;
    
    // Read initial conditions
    
    Coordinates *coord = mesh->getCoordinates();
    
    mesh->get(V0, "V0");
    mesh->get(coord->dx,   "dx");
    mesh->get(coord->dy,   "dy");
    
    // read options
    
    // Set evolving variables
    SOLVE_FOR(V);
    
    if(!restarting) {
      // Set variables to these values (+ the initial perturbation)
      // NOTE: This must be after the calls to bout_solve
      V += V0;
    }
    return 0;
  }
  
  int rhs(BoutReal UNUSED(t)) {
    // Run communications
    mesh->communicate(V);
    
    //ddt(V) = D2DX2(V) + 0.5*DDX(V) + D2DY2(V);
    ddt(V) = DDX(V);
    
    return 0;
  }
};

BOUTMAIN(AdvDiff);
