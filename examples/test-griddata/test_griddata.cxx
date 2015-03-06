#include <bout.hxx>

int main(int argc, char** argv) {
  BoutInitialise(argc, argv);

  Datafile df(Options::getRoot()->getSection("output"));

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
