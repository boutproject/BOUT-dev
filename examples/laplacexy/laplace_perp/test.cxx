#include "../../common.hxx"

#include <bout/derivs.hxx>
#include <bout/field_factory.hxx>
#include <bout/invert/laplacexy.hxx>

using bout::globals::mesh;

int main(int argc, char** argv) {
  BoutInitialise(argc, argv);

  ///////////////////////////////////////
  bool calc_metric;
  calc_metric = Options::root()["calc_metric"].withDefault(false);
  if (calc_metric) {
    // Read metric tensor
    Field2D Rxy, Btxy, Bpxy, B0, hthe, I;
    mesh->get(Rxy, "Rxy");   // m
    mesh->get(Btxy, "Btxy"); // T
    mesh->get(Bpxy, "Bpxy"); // T
    mesh->get(B0, "Bxy");    // T
    mesh->get(hthe, "hthe"); // m
    mesh->get(I, "sinty");   // m^-2 T^-1

    Coordinates* coord = mesh->getCoordinates();

    tokamak_coordinates(coord, Rxy, Bpxy, hthe, I, B0, Btxy);
  }
  ///////////////////////////////////////

  // Read an analytic input
  Field2D input = FieldFactory::get()->create2D("input", Options::getRoot(), mesh);

  // Create a LaplaceXY solver
  LaplaceXY* laplacexy = new LaplaceXY(mesh);

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
