#include <cstdio>

#include "gtest/gtest.h"
#include "bout/array.hxx"
#include "fft.hxx"
#include "output.hxx"
#include "bout/globalindexer.hxx"
// Note: petsclib included after globalindexer, or MPI_Waitall
// in mpi_wrapper.hxx is expanded as a macro
#include "bout/petsclib.hxx" 

GTEST_API_ int main(int argc, char** argv) {

  // Make sure fft functions are both quiet and deterministic by
  // setting fft_measure to false
  bout::fft::fft_init(false);

  printf("Running main() from bout_test_main.cxx\n");
  testing::InitGoogleTest(&argc, argv);

  // Explicitly setup and teardown PETSc to avoid reentry problems
  // with certain MPI implementations (see #1916 for details)
  output.disable();
  PetscLib petsclib{};
  output.enable();

  int result = RUN_ALL_TESTS();

  // Explicit cleanup of PetscLib because it might get destroyed
  // _after_ `output`
  output.disable();
  PetscLib::cleanup();
  output.enable();

  // Clean up the array store, so valgrind doesn't report false
  // positives
  Array<double>::cleanup();
  Array<int>::cleanup();
  Array<bool>::cleanup();

  return result;
}
