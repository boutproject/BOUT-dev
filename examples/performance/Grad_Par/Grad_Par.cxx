/*
 * Timing of serial tri solver
 *
 */

#include <bout/physicsmodel.hxx>
#include <bout/scorepwrapper.hxx>
#include <chrono>
#include <derivs.hxx>              // To use DDZ()

int main(int argc, char **argv) {

  SCOREP0();
  BoutInitialise(argc, argv);
  typedef std::chrono::time_point<std::chrono::steady_clock> SteadyClock;
  typedef std::chrono::duration<double> Duration;
  using namespace std::chrono;

  Field3D omega;
  Field3D phi;

  Options *globalOptions = Options::getRoot();
  Options *options = globalOptions->getSection("model");

  phi = 0.0; // Starting guess for first solve (if iterative)
  omega = 1.0;
  SteadyClock start1 = steady_clock::now();
  for(int i=0; i<1000; i++){
    phi = DDZ(omega);
  }
  Duration elapsed1 = steady_clock::now() - start1;
  
  output << "TIMING\n======\n";
  output << "DDZ: " << elapsed1.count() << endl;
  
  BoutFinalise();
  return 0;
};

