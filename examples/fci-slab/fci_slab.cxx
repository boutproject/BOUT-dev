#include <bout/physicsmodel.hxx>
#include <fci_derivs.hxx>
#include <utils.hxx>
#include <bout/mesh.hxx>

#include <field_factory.hxx>

class FCISlab : public PhysicsModel {
public:

  // We need to initialise the FCI object with the mesh
  FCISlab() : fci(*mesh, false, false) {}

  int init(bool restarting) {

	D = 10;

    Coordinates *coord = mesh->coordinates();

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

  FCI fci;

  FieldGenerator* f_gen;
  FieldGenerator* g_gen;
};

BOUTMAIN(FCISlab);

int FCISlab::rhs(BoutReal time) {
  mesh->communicate(f,g);

  Coordinates *coord = mesh->coordinates();

  fci.calcYUpDown(f);
  f.applyParallelBoundary(time);
  fci.calcYUpDown(g);
  g.applyParallelBoundary(time);

  ddt(f) = fci.Grad_par(g) + D*SQ(coord->dy)*fci.Grad2_par2(f);

  ddt(g) = fci.Grad_par(f) + D*SQ(coord->dy)*fci.Grad2_par2(g);

  return 0;
}
