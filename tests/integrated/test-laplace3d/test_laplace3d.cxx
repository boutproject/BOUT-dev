
#include <bout/physicsmodel.hxx>

#include <bout/invert/laplace3d.hxx>

class Laplace3DTest : public PhysicsModel {
protected:
  int init(bool restarting);
};


int Laplace3DTest::init(bool restarting) {
  
  // Create a 3D Laplacian solver
  Laplace3D *lap = Laplace3D::create();
  
  // Set coefficients
  
  // solve

  return 1;
}

// Generate a main() function to run test
BOUTMAIN(Laplace3DTest); 
