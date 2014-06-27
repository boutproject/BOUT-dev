#include <bout/physicsmodel.hxx>
#include <fci_derivs.hxx>
#include <utils.hxx>
#include <bout/mesh.hxx>

class FCISlab : public PhysicsModel {
public:

  // We need to initialise the FCI object with the mesh
  FCISlab() : fci(*mesh) {}

  int init(bool restarting) {

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

  FCI fci;
};

BOUTMAIN(FCISlab);

int FCISlab::rhs(BoutReal time) {
  mesh->communicate(f,g);
  ddt(f) = fci.Grad_par(g);
  ddt(g) = fci.Grad_par(f);

  return 0;
}
