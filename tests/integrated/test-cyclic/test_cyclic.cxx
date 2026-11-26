/*
 * Test Cyclic Reduction parallel solver
 *
 */

#include <bout/bout.hxx>
#include <bout/bout_types.hxx>
#include <bout/cyclic_reduction.hxx>
#include <bout/options.hxx>
#include <bout/output.hxx>
#include <bout/utils.hxx>

#include <cstdlib>

#include <mpi.h>

// Change this to dcomplex to test complex matrix inversion
using T = BoutReal;

int main(int argc, char** argv) {

  // Initialise BOUT++, setting up mesh
  BoutInitialise(argc, argv);

  Options& options = Options::root();
  const int n = options["n"].withDefault(5);
  const int nsys = options["nsys"].withDefault(1);
  const BoutReal tol = options["tol"].withDefault(1e-10);
  const bool periodic = options["periodic"].withDefault(false);

  // Create a cyclic reduction object, operating on Ts
  auto cr = CyclicReduce<T>(BoutComm::get(), n);

  const int mype = BoutComm::rank();
  const int npe = BoutComm::size();

  Matrix<T> a(nsys, n);
  Matrix<T> b(nsys, n);
  Matrix<T> c(nsys, n);
  Matrix<T> rhs(nsys, n);
  Matrix<T> x(nsys, n);

  // Set coefficients to some random numbers
  for (int s = 0; s < nsys; s++) {
    for (int i = 0; i < n; i++) {
      a(s, i) = randomu();
      b(s, i) = 2. + randomu(); // So always diagonally dominant
      c(s, i) = randomu();

      // deterministic, so don't need to communicate
      x(s, i) = (((mype * n) + i) % 4) - 2.;
    }

    // Calculate RHS

    int i = 0;
    if ((mype == 0) && (!periodic)) {
      rhs(s, i) = (b(s, i) * x(s, i)) + (c(s, i) * x(s, i + 1));
    } else {
      const int pe = (mype - 1 + npe) % npe;
      const T xm = (((pe * n) + (n - 1)) % 4) - 2.;
      rhs(s, i) = (a(s, i) * xm) + (b(s, i) * x(s, i)) + (c(s, i) * x(s, i + 1));
    }

    for (i = 1; i < n - 1; i++) {
      rhs(s, i) = (a(s, i) * x(s, i - 1)) + (b(s, i) * x(s, i)) + (c(s, i) * x(s, i + 1));
    }

    if ((mype == (npe - 1)) && !periodic) {
      rhs(s, i) = (a(s, i) * x(s, i - 1)) + (b(s, i) * x(s, i));
    } else {
      const int pe = (mype + 1) % npe;
      const T xp = ((pe * n) % 4) - 2.;
      rhs(s, i) = (a(s, i) * x(s, i - 1)) + (b(s, i) * x(s, i)) + (c(s, i) * xp);
    }

    // Zero x, so that solve has to do something
    for (int i = 0; i < n; i++) {
      x(s, i) = 0.0;
    }
  }

  // Solve system
  cr.setPeriodic(periodic);
  cr.setCoefs(a, b, c);
  cr.solve(rhs, x);

  // Check result
  bool failed = false;
  for (int s = 0; s < nsys; s++) {
    output.write("System {}\n", s);
    for (int i = 0; i < n; i++) {
      const T val = (((mype * n) + i) % 4) - 2.;
      output.write("\t{} : {} ?= {}\n", i, val, x(s, i));
      failed |= (std::abs(val - x(s, i)) > tol);
    }
  }

  bool any_failed = false;
  MPI_Allreduce(&failed, &any_failed, 1, MPI_CXX_BOOL, MPI_LOR, BoutComm::get());

  output.write("******* Cyclic test case: ");
  if (any_failed) {
    output.write("FAILED\n");
  } else {
    output.write("PASSED\n");
  }

  MPI_Barrier(BoutComm::get());

  bout::checkForUnusedOptions();
  BoutFinalise();
  return static_cast<int>(any_failed);
}
