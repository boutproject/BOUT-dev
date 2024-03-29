/*
 * Test LaplaceXZ solver
 *
 * Matrix assembly information:
 *    -mat_view ::ascii_info
 *
 * Useful diagnostics for SuperLU_dist solver
 * (pctype=lu, factor_package=superlu_dist)
 * -mat_superlu_dist_statprint
 */
#include <bout/bout.hxx>

#include <bout/derivs.hxx>
#include <bout/field_factory.hxx>
#include <bout/invert/laplacexz.hxx>

int main(int argc, char** argv) {
  BoutInitialise(argc, argv);

  auto inv = LaplaceXZ::create(bout::globals::mesh);

  auto coord = bout::globals::mesh->getCoordinates();
  coord->g13 = 1.8; // test off-diagonal components with nonzero value

  // create some input field
  Field3D f = FieldFactory::get()->create3D("f", Options::getRoot(), bout::globals::mesh);

  // Calculate the Laplacian with non-zero g13
  Field3D g = coord->g11 * D2DX2(f) + coord->g13 * D2DXDZ(f) + coord->g33 * D2DZ2(f);

  inv->setCoefs(Field2D(1.0), Field2D(0.0));

  Field3D f2 = inv->solve(g, 0.0); // Invert the Laplacian.

  coord->g13 = 0.0; // reset to 0.0 for original laplacexz test

  // Now the normal test.
  output.write("Setting coefficients\n");

  inv->setCoefs(Field3D(1.0), Field3D(0.0));

  output.write("First solve\n");

  Field3D rhs =
      FieldFactory::get()->create3D("rhs", Options::getRoot(), bout::globals::mesh);
  Field3D x = inv->solve(rhs, 0.0);

  output.write("Second solve\n");

  inv->setCoefs(Field3D(2.0), Field3D(0.1));

  Field3D rhs2 =
      FieldFactory::get()->create3D("rhs", Options::getRoot(), bout::globals::mesh);
  Field3D x2 = inv->solve(rhs2, 0.0);

  Options dump;

  dump["f"] = f;
  dump["f2"] = f2;
  dump["g"] = g;
  dump["rhs"] = rhs;
  dump["x"] = x;
  dump["rhs2"] = rhs2;
  dump["x2"] = x2;

  bout::writeDefaultOutputFile(dump);

  BoutFinalise();
  return 0;
}
