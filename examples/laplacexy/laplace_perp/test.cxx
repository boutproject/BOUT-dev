#include <bout.hxx>

#include <bout/invert/laplacexy.hxx>
#include <derivs.hxx>
#include <field_factory.hxx>

using bout::globals::mesh;

int main(int argc, char** argv) {
  BoutInitialise(argc, argv);
  
  ///////////////////////////////////////
  bool calc_metric;
  calc_metric = Options::root()["calc_metric"].withDefault(false);
  if(calc_metric) {
    // Read metric tensor
    Field2D Rxy, Btxy, Bpxy, B0, hthe, I;
    mesh->get(Rxy,  "Rxy");  // m
    mesh->get(Btxy, "Btxy"); // T
    mesh->get(Bpxy, "Bpxy"); // T
    mesh->get(B0,   "Bxy");  // T
    mesh->get(hthe, "hthe"); // m
    mesh->get(I,    "sinty");// m^-2 T^-1

    Coordinates *coord = mesh->getCoordinates();

    // Calculate metrics
    coord->g11 = SQ(Rxy * Bpxy);
    coord->g22 = 1.0 / SQ(hthe);
    coord->g33 = SQ(I) * coord->g11 + SQ(B0) / coord->g11;
    coord->g12 = 0.0;
    coord->g13 = -I * coord->g11;
    coord->g23 = -Btxy / (hthe * Bpxy * Rxy);

    coord->J = hthe / Bpxy;
    coord->Bxy = B0;

    coord->g_11 = 1.0 / coord->g11 + SQ(I * Rxy);
    coord->g_22 = SQ(B0 * hthe / Bpxy);
    coord->g_33 = Rxy * Rxy;
    coord->g_12 = Btxy * hthe * I * Rxy / Bpxy;
    coord->g_13 = I * Rxy * Rxy;
    coord->g_23 = Btxy * hthe * Rxy / Bpxy;

    coord->geometry();
  }
  ///////////////////////////////////////
  
  // Read an analytic input
  Field2D input = FieldFactory::get()->create2D("input", Options::getRoot(), mesh);
  
  // Create a LaplaceXY solver
  LaplaceXY *laplacexy = new LaplaceXY(mesh);
  
  // Solve, using 0.0 as starting guess
  Field2D solved = laplacexy->solve(input, 0.0);
  
  // Need to communicate guard cells
  mesh->communicate(solved);
  
  // Now differentiate using Laplace_perp
  Options::root()["result"] = Laplace_perp(solved);

  // Write fields to output
  Options::root()["input"] = input;
  Options::root()["solved"] = solved;

  bout::writeDefaultOutputFile();

  BoutFinalise();
  return 0;
}

