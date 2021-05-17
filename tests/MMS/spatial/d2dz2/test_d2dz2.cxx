/*
 * Test D2DZ2 operator without time integration
 */

#include <bout.hxx>
#include <field_factory.hxx>
#include <derivs.hxx>

using bout::globals::mesh;
using bout::globals::dump;

int main(int argc, char** argv) {

  BoutInitialise(argc, argv);
  
  Field3D input = FieldFactory::get()->create3D("input", Options::getRoot(), mesh);
  Field3D solution = FieldFactory::get()->create3D("solution", Options::getRoot(), mesh);
  
  Field3D result = D2DZ2(input);

  SAVE_ONCE3(input, solution, result);
  dump.write();
  
  BoutFinalise();

  return 0;
}
