#include <bout/bout.hxx>
#include "bout/field_factory.hxx"

int main(int argc, char** argv) {

  BoutInitialise(argc, argv);
  using bout::globals::mesh;


  Options dump;
  FieldFactory factory {mesh, &Options::root()};
  dump["x"] = factory.create3D("x");
  dump["y"] = factory.create3D("y");

  // Write to BOUT_coord.dmp.* instead of the default BOUT.dmp.*
  Options coord_out_opts;
  coord_out_opts["prefix"] = "BOUT_coord.dmp";
  bout::OptionsIOFactory::getInstance().createOutput(&coord_out_opts)->write(dump);

  BoutFinalise();

  return 0;
}