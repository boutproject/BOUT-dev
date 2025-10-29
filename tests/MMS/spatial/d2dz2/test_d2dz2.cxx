/*
 * Test D2DZ2 operator without time integration
 */

#include <bout/bout.hxx>
#include <bout/derivs.hxx>
#include <bout/field3d.hxx>
#include <bout/field_factory.hxx>
#include <bout/options.hxx>

using bout::globals::mesh;

int main(int argc, char** argv) {

  BoutInitialise(argc, argv);

  const Field3D input =
      FieldFactory::get()->create3D("input_field", Options::getRoot(), mesh);
  const Field3D solution =
      FieldFactory::get()->create3D("solution", Options::getRoot(), mesh);

  Field3D result = D2DZ2(input);

  Options dump;
  dump["input"] = input;
  dump["solution"] = solution;
  dump["result"] = result;
  bout::writeDefaultOutputFile(dump);

  BoutFinalise();

  return 0;
}
