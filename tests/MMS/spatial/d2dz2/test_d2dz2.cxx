/*
 * Test D2DZ2 operator without time integration
 */

#include <bout.hxx>
#include <field_factory.hxx>
#include <derivs.hxx>

using bout::globals::mesh;

int main(int argc, char** argv) {

  BoutInitialise(argc, argv);

  Field3D input = FieldFactory::get()->create3D("input", Options::getRoot(), mesh);
  Field3D solution = FieldFactory::get()->create3D("solution", Options::getRoot(), mesh);

  Field3D result = D2DZ2(input);

  Options dump;
  dump["input"] = input;
  dump["solution"] = solution;
  dump["result"] = result;
  bout::writeDefaultOutputFile(dump);

  BoutFinalise();

  return 0;
}
