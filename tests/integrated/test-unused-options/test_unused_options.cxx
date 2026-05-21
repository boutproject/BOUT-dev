#include <bout/bout.hxx>
#include <bout/field3d.hxx>
#include <bout/options.hxx>
#include <bout/output.hxx>

int main(int argc, char** argv) {
  BoutInitialise(argc, argv);

  Field3D f;
  f.setBoundary("f");

  bout::checkForUnusedOptions();

  BoutFinalise();

  // Note: Print message because sometimes MPI ranks error on tidyup.
  output.write("SUCCESS");

  return 0;
}
