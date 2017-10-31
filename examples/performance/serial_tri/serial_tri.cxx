/*
 * Timing of serial tri solver
 *
 */

#include <bout/physicsmodel.hxx>
#include <bout/scorepwrapper.hxx>
#include <chrono>
#include <invert_laplace.hxx>      // Laplacian inversion

int main(int argc, char **argv) {

  SCOREP0();
  BoutInitialise(argc, argv);
  typedef std::chrono::time_point<std::chrono::steady_clock> SteadyClock;
  typedef std::chrono::duration<double> Duration;
  using namespace std::chrono;

  Field3D omega;
  Field3D phi;

  Laplacian *phiSolver;  ///< Performs Laplacian inversions to calculate phi

  Options *globalOptions = Options::getRoot();
  Options *options = globalOptions->getSection("model");

  phiSolver = Laplacian::create(Options::getRoot()->getSection("phiBoussinesq")); // BOUT.inp section "phiBoussinesq"
  phi = 0.0; // Starting guess for first solve (if iterative)

  omega = 1.0;
  SteadyClock start1 = steady_clock::now();
  for(int i=0; i<100; i++){
    phi = phiSolver->solve(omega);
  }
  Duration elapsed1 = steady_clock::now() - start1;
  
  output << "TIMING\n======\n";
  output << "Serial Tri: " << elapsed1.count() << endl;
  
  BoutFinalise();
  return 0;
};

