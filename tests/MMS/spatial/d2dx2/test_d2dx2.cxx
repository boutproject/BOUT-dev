/*
 * Test D2DX2 operator without time integration
 */

#include <bout/bout.hxx>
#include <bout/derivs.hxx>
#include <bout/field_factory.hxx>

using bout::globals::mesh;

int main(int argc, char** argv) {

  BoutInitialise(argc, argv);

  Field3D input = FieldFactory::get()->create3D("input", Options::getRoot(), mesh);
  Field3D solution = FieldFactory::get()->create3D("solution", Options::getRoot(), mesh);
  // At this point the boundary cells are set to the analytic solution

  input.setBoundary("bndry");
  input.applyBoundary(0.0);

  // Boundaries of input now set using extrapolation around mid-point boundary

  Field3D result = D2DX2(input);
  // At this point result is not set in the boundary cells

  Options dump;
  dump["input"] = input;
  dump["solution"] = solution;
  dump["result"] = result;
  bout::writeDefaultOutputFile(dump);

  BoutFinalise();

  return 0;
}
