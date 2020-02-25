#include <bout.hxx>

int main(int argc, char** argv) {
  BoutInitialise(argc, argv);

  Datafile df(Options::getRoot()->getSection("output"));
  df.add(const_cast<BoutReal&>(bout::version::as_double), "BOUT_VERSION", false);

  mesh->outputVars(df);

  Field2D Rxy, Bpxy;
  mesh->get(Rxy, "Rxy");
  mesh->get(Bpxy, "Bpxy");

  df.add(Rxy, "Rxy");
  df.add(Bpxy, "Bpxy");

  df.write("data.nc");
  
  BoutFinalise();
  return 0;
}
