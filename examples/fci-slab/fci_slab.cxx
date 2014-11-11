#include <bout/physicsmodel.hxx>
#include <fci_derivs.hxx>
#include <utils.hxx>
#include <bout/mesh.hxx>

class FCISlab : public PhysicsModel {
public:

  // We need to initialise the FCI object with the mesh
  FCISlab() : fci(*mesh) {}

  int init(bool restarting) {

    D = 10;
    //OPTION(Options::getRoot(), D, 1.);

    mesh->get(mesh->g_22, "g_22");
    
    mesh->geometry();
    
    solver->add(f, "f");
    solver->add(g, "g");
    
    f.applyBoundary("dirichlet");
    g.applyBoundary("dirichlet");

    return 0;
  }

  int rhs(BoutReal time);

private:
  Field3D f, g;

  BoutReal D;

  FCI fci;
};

BOUTMAIN(FCISlab);

int FCISlab::rhs(BoutReal time) {
  mesh->communicate(f,g);
  ddt(f) = fci.Grad_par(g) + D*SQ(mesh->dy)*fci.Grad2_par2(f);

  //for(int i=0;i<mesh->ngx;i++)
  //  output.write("%i: %e\n", i, ddt(f)(i,16,0));

  ddt(g) = fci.Grad_par(f) + D*SQ(mesh->dy)*fci.Grad2_par2(g);

  return 0;
}
