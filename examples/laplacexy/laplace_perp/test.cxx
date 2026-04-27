#include <bout/bout.hxx>
#include <bout/derivs.hxx>
#include <bout/field2d.hxx>
#include <bout/field_factory.hxx>
#include <bout/invert/laplacexy.hxx>
#include <bout/tokamak_coordinates.hxx>

using bout::globals::mesh;

int main(int argc, char** argv) {
  BoutInitialise(argc, argv);

  ///////////////////////////////////////
  const bool calc_metric = Options::root()["calc_metric"].withDefault(false);
  if (calc_metric) {
    bout::set_tokamak_coordinates(*mesh);
  }
  ///////////////////////////////////////

  // Read an analytic input
  const Field2D input =
      FieldFactory::get()->create2D("input_field", Options::getRoot(), mesh);

  // Create a LaplaceXY solver
  auto laplacexy = LaplaceXY::create(mesh);

  // Solve, using 0.0 as starting guess
  Field2D solved = laplacexy->solve(input, 0.0);

  // Need to communicate guard cells
  mesh->communicate(solved);

  // Now differentiate using Laplace_perp
  Options::root()["result"] = Laplace_perp(solved);

  // Write fields to output
  Options::root()["input"] = input;
  Options::root()["solved"] = solved;

  bout::writeDefaultOutputFile(Options::root());

  BoutFinalise();
  return 0;
}
