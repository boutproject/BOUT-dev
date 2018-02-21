/*
 * Test of eigenvalue solver in box
 * 
 * Solves wave equation
 *  d^2f/dt^2 = d^2f/dx^2
 */

#include <bout/physicsmodel.hxx>
#include <bout/derivs.hxx>

class EigenBox : public PhysicsModel {
protected:
  int init(bool restarting) {
    solver->add(f, "f");
    solver->add(g, "g");
    return 0;
  }
  int rhs(BoutReal t) {
    mesh->communicate(f);
    
    ddt(g) = D2DX2(f);
    ddt(f) = g;
    
    return 0;
  }
private:
  Field3D f, g;
};

BOUTMAIN(EigenBox);
