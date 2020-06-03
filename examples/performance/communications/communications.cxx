#include <bout.hxx>
#include <bout/sys/timer.hxx>
#include <initialprofiles.hxx>

int main(int argc, char** argv) {

  BoutInitialise(argc, argv);

  Field3D f;
  initial_profile("f", f);
  int iterations = Options::root()["iterations"].withDefault(10000);

  Timer timer("comms");
  for (int i = 0; i < iterations; ++i) {
    f.getMesh()->communicate(f);
  }
  BoutReal run_length = timer.getTime();

  output << iterations << " iterations took " << run_length << "s";

  BoutFinalise();

  return 0;
}
