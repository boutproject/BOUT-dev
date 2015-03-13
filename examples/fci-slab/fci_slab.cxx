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

	Options *globalOptions = Options::getRoot();

	Options *f_options = globalOptions->getSection("f");

	string f_BC;
	OPTION(f_options, f_BC, "0.0");

	if(!f_BC.empty()) {
	  // First argument should be an expression
	  f_gen = FieldFactory::get()->parse(f_BC);
	}

	Options *g_options = globalOptions->getSection("g");

	string g_BC;
	OPTION(g_options, g_BC, "0.0");

	if(!g_BC.empty()) {
	  // First argument should be an expression
	  g_gen = FieldFactory::get()->parse(g_BC);
	}

	D = 10;
	//OPTION(Options::getRoot(), D, 1.);

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
  fci.applyBoundary(f, FCI::DIRICHLET, f_gen, time);
  fci.calcYUpDown(g);
  fci.applyBoundary(g, FCI::DIRICHLET, g_gen, time);

  ddt(f) = fci.Grad_par(g) + D*SQ(coord->dy)*fci.Grad2_par2(f);

  ddt(g) = fci.Grad_par(f) + D*SQ(coord->dy)*fci.Grad2_par2(g);

  return 0;
}
