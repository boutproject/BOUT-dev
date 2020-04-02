#include <bout.hxx>
#include <bout/version.hxx>

int main(int argc, char** argv) {
  BoutInitialise(argc, argv);

  auto options = Options::getRoot();//options["dump_format"].
  Datafile df(options->getSection("output"));
  df.add(const_cast<BoutReal&>(bout::version::as_double), "BOUT_VERSION", false);

  mesh->outputVars(df);

  Field2D Rxy, Bpxy;
  mesh->get(Rxy, "Rxy");
  mesh->get(Bpxy, "Bpxy");

  df.add(Rxy, "Rxy");
  df.add(Bpxy, "Bpxy");

  std::string ext = (*options)["dump_format"];
  df.write("data.{:s}",ext);
  
  BoutFinalise();
  return 0;
}
