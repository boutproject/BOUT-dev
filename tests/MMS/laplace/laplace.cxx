#include <bout/bout.hxx>
#include <bout/constants.hxx>
#include <bout/field3d.hxx>
#include <bout/field_factory.hxx>
#include <bout/invert_laplace.hxx>
#include <bout/output.hxx>

#include <string>

using bout::globals::mesh;
using namespace std::string_literals;

int main(int argc, char** argv) {
  int init_err = BoutInitialise(argc, argv);
  if (init_err < 0) {
    return 0;
  }
  if (init_err > 0) {
    return init_err;
  }

  ////// Set mesh spacing
  Options* meshoptions = Options::getRoot()->getSection("mesh");

  BoutReal Lx;
  meshoptions->get("Lx", Lx, 1.0);

  /*this assumes equidistant grid*/
  int nguard = mesh->xstart;
  mesh->getCoordinates()->dx = Lx / (mesh->GlobalNx - 2 * nguard);
  mesh->getCoordinates()->dz = TWOPI * Lx / (mesh->LocalNz);
  /////

  // Create a Laplacian inversion solver
  auto lap = Laplacian::create();

  FieldFactory fact(mesh);

  const auto input_name = "input_field"s;
  const auto gen = fact.parse(input_name);
  output.write("GEN = {}\n", gen->str());

  const Field3D input = fact.create3D(input_name);
  const Field3D result = lap->solve(input);
  const Field3D solution = fact.create3D("solution");
  const Field3D error = result - solution;

  Options dump;
  dump["input"] = input;
  dump["result"] = result;
  dump["solution"] = solution;
  dump["error"] = error;
  bout::writeDefaultOutputFile(dump);

  BoutFinalise();
  return 0;
}
