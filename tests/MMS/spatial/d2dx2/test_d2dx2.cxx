/*
 * Test D2DX2 operator without time integration
 */

#include <bout.hxx>
#include <bout/field_factory.hxx>
#include <bout/derivs.hxx>

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

  SAVE_ONCE3(input, solution, result);
  dump.write();
  
  BoutFinalise();

  return 0;
}
