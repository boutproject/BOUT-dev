#include <bout.hxx>

#include <invert_laplace.hxx>
#include <field_factory.hxx>
#include <bout/constants.hxx>

using bout::globals::mesh;

int main(int argc, char **argv) {
  int init_err = BoutInitialise(argc, argv);
  if (init_err < 0) {
    return 0;
  } else if (init_err > 0) {
    return init_err;
  }

  ////// Set mesh spacing
  Options *meshoptions = Options::getRoot()->getSection("mesh");

  BoutReal Lx;
  meshoptions->get("Lx",Lx,1.0);

  /*this assumes equidistant grid*/
  int nguard = mesh->xstart;
  mesh->getCoordinates()->dx = Lx/(mesh->GlobalNx - 2*nguard);
  mesh->getCoordinates()->dz = TWOPI*Lx/(mesh->LocalNz);
  /////

  // Create a Laplacian inversion solver
  auto lap = Laplacian::create();
  
  FieldFactory fact(mesh);

  std::shared_ptr<FieldGenerator> gen = fact.parse("input");
  output << "GEN = " << gen->str() << endl;

  Field3D input = fact.create3D("input");
  
  Field3D result = lap->solve(input);
  
  Field3D solution = fact.create3D("solution");

  Field3D error = result - solution;

  SAVE_ONCE4(input, result, solution, error);
  dump.write();
  
  BoutFinalise();
  return 0;
}
