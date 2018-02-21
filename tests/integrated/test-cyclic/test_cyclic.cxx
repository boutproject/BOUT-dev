/*
 * Test Cyclic Reduction parallel solver
 *
 */

#include <bout.hxx>

#include <bout/cyclic_reduction.hxx>
#include <bout/dcomplex.hxx>

// Change this to dcomplex to test complex matrix inversion
typedef BoutReal T;

int main(int argc, char **argv) {

  // Initialise BOUT++, setting up mesh
  BoutInitialise(argc, argv);

  T **a, **b, **c, **rhs, **x;

  int nsys;
  int n;
  Options *options = Options::getRoot();
  OPTION(options, n, 5);
  OPTION(options, nsys, 1);
  BoutReal tol;
  OPTION(options, tol, 1e-10);
  bool periodic;
  OPTION(options, periodic, false);

  // Create a cyclic reduction object, operating on Ts
  CyclicReduce<T> *cr =
    new CyclicReduce<T>(BoutComm::get(), n);

  int mype, npe;
  MPI_Comm_rank(BoutComm::get(), &mype);
  MPI_Comm_size(BoutComm::get(), &npe);

  a   = matrix<T>(nsys, n);
  b   = matrix<T>(nsys, n);
  c   = matrix<T>(nsys, n);
  rhs = matrix<T>(nsys, n);
  x   = matrix<T>(nsys, n);

  // Set coefficients to some random numbers
  for(int s=0;s<nsys;s++) {
    for(int i=0;i<n;i++) {
      a[s][i] = randomu();
      b[s][i] = 2. + randomu(); // So always diagonally dominant
      c[s][i] = randomu();

      x[s][i] = ((mype*n + i) % 4) - 2.; // deterministic, so don't need to communicate
    }

    // Calculate RHS

    int i = 0;
    if((mype == 0) && (!periodic)) {
      rhs[s][i] = b[s][i]*x[s][i] + c[s][i]*x[s][i+1];
    }else {
      int pe = (mype - 1 + npe) % npe;
      T xm = ((pe*n + (n-1)) % 4) - 2.;
      rhs[s][i] = a[s][i]*xm + b[s][i]*x[s][i] + c[s][i]*x[s][i+1];
    }

    for(i=1;i<n-1;i++)
      rhs[s][i] = a[s][i]*x[s][i-1] + b[s][i]*x[s][i] + c[s][i]*x[s][i+1];

    if((mype == (npe-1)) && !periodic) {
      rhs[s][i] = a[s][i]*x[s][i-1] + b[s][i]*x[s][i];
    }else {
      int pe = (mype + 1) % npe;
      T xp = ((pe*n) % 4) - 2.;
      rhs[s][i] = a[s][i]*x[s][i-1] + b[s][i]*x[s][i] + c[s][i]*xp;
    }
  }

  // Solve system

  cr->setPeriodic(periodic);
  cr->setCoefs(nsys, a, b, c);
  cr->solve(nsys, rhs, x);

  // Check result

  int passed = 1;
  for(int s=0;s<nsys;s++) {
    output << "System " << s << endl;
    for(int i=0;i<n;i++) {
      T val = ((mype*n + i) % 4) - 2.;
      output << "\t" << i << " : " << val << " ?= " << x[s][i] << endl;
      if(abs(val - x[s][i]) > tol) {
        passed = 0;
      }
    }
  }

  int allpassed;
  MPI_Allreduce(&passed, &allpassed, 1, MPI_INT, MPI_MIN, BoutComm::get());

  // Saving state to file
  SAVE_ONCE(allpassed);

  output << "******* Cyclic test case: ";
  if(allpassed) {
    output << "PASSED" << endl;
  }else
    output << "FAILED" << endl;

  // Write data to file
  dump.write();
  dump.close();

  MPI_Barrier(BoutComm::get());

  BoutFinalise();
  return 0;
}
