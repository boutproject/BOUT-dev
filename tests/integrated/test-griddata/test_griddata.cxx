#include <bout/version.hxx>
#include <bout.hxx>

int main(int argc, char** argv) {
  BoutInitialise(argc, argv);

  Field2D Rxy, Bpxy;
  bout::globals::mesh->get(Rxy, "Rxy");
  bout::globals::mesh->get(Bpxy, "Bpxy");

  Options dump;
  dump["Rxy"] = Rxy;
  dump["Bpxy"] = Bpxy;
  bout::experimental::addBuildFlagsToOptions(dump);
  bout::globals::mesh->outputVars(dump);
  bout::OptionsNetCDF("data.nc").write(dump);

  BoutFinalise();
  return 0;
}
