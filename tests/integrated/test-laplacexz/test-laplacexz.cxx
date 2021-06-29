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
#include <bout.hxx>

#include <bout/invert/laplacexz.hxx>
#include <field_factory.hxx>

using bout::globals::mesh;

int main(int argc, char** argv) {
  BoutInitialise(argc, argv);

  auto inv = LaplaceXZ::create(mesh);

  output.write("Setting coefficients\n");

  inv->setCoefs(Field3D(1.0),Field3D(0.0));

  output.write("First solve\n");

  Field3D rhs = FieldFactory::get()->create3D("rhs", Options::getRoot(), mesh);
  Field3D x = inv->solve(rhs, 0.0);

  SAVE_ONCE2(rhs, x);

  output.write("Second solve\n");

  inv->setCoefs(Field3D(2.0),Field3D(0.1));

  Field3D rhs2 = FieldFactory::get()->create3D("rhs", Options::getRoot(), mesh);
  Field3D x2 = inv->solve(rhs2, 0.0);
  SAVE_ONCE2(rhs2, x2);

  bout::globals::dump.write();

  BoutFinalise();
  return 0;
}
