#include <bout/bout.hxx>

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

    // Calculate metrics
    const auto g11 = SQ(Rxy * Bpxy);
    const auto g22 = 1.0 / SQ(hthe);
    const auto g33 = SQ(I) * g11 + SQ(B0) / g11;
    const auto g12 = 0.0;
    const auto g13 = -I * g11;
    const auto g23 = -Btxy / (hthe * Bpxy * Rxy);

    const auto g_11 = 1.0 / g11 + SQ(I * Rxy);
    const auto g_22 = SQ(B0 * hthe / Bpxy);
    const auto g_33 = Rxy * Rxy;
    const auto g_12 = Btxy * hthe * I * Rxy / Bpxy;
    const auto g_13 = I * Rxy * Rxy;
    const auto g_23 = Btxy * hthe * Rxy / Bpxy;

    coord->setMetricTensor(ContravariantMetricTensor(g11, g22, g33, g12, g13, g23),
                           CovariantMetricTensor(g_11, g_22, g_33, g_12, g_13, g_23));

    coord->setJ(hthe / Bpxy);
    coord->setBxy(B0);
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
