#include <bout.hxx>
#include <bout/version.hxx>

int main(int argc, char** argv) {
  BoutInitialise(argc, argv);

  Datafile df(Options::getRoot()->getSection("output"));
  df.add(const_cast<BoutReal&>(bout::version::as_double), "BOUT_VERSION", false);

  bout::globals::mesh->outputVars(df);

  Field2D Rxy, Bpxy;
  bout::globals::mesh->get(Rxy, "Rxy");
  bout::globals::mesh->get(Bpxy, "Bpxy");

  df.add(Rxy, "Rxy");
  df.add(Bpxy, "Bpxy");

  df.write("data.nc");
  
  BoutFinalise();
  return 0;
}
