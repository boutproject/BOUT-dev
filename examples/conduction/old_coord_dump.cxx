#include <bout/bout.hxx>
#include "bout/field_factory.hxx"

int main(int argc, char** argv) {

  BoutInitialise(argc, argv);
  using bout::globals::mesh;


  Options dump;
  FieldFactory factory {mesh, &Options::root()};
  dump["x"] = factory.create3D("x");
  dump["y"] = factory.create3D("y");
  bout::writeDefaultOutputFile(dump);

  BoutFinalise();

  return 0;
}