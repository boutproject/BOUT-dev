#include <bout/physicsmodel.hxx>
#include <utils.hxx>
#include <bout/mesh.hxx>

class FCISlab : public PhysicsModel {
public:

  // We need to initialise the FCI object with the mesh
  FCISlab() {}

  int init(bool UNUSED(restarting)) {

    D = 10;

    Coordinates *coord = mesh->getCoordinates();

    mesh->get(coord->g_22, "g_22");

    coord->geometry();

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
};

BOUTMAIN(FCISlab);

int FCISlab::rhs(BoutReal time) {
  mesh->communicate(f,g);

  Coordinates *coord = mesh->getCoordinates();

  f.applyParallelBoundary(time);
  g.applyParallelBoundary(time);

  ddt(f) = Grad_par(g) + D*SQ(coord->dy)*Grad2_par2(f);

  ddt(g) = Grad_par(f) + D*SQ(coord->dy)*Grad2_par2(g);

  return 0;
}
