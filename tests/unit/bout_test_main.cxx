#include <stdio.h>

#include "gtest/gtest.h"
#include "bout/array.hxx"
#include "fft.hxx"

GTEST_API_ int main(int argc, char** argv) {

  // Make sure fft functions are both quiet and deterministic by
  // setting fft_measure to false
  bout::fft::fft_init(false);

  printf("Running main() from bout_test_main.cxx\n");
  testing::InitGoogleTest(&argc, argv);
  int result = RUN_ALL_TESTS();

  // Clean up the array store, so valgrind doesn't report false
  // positives
  Array<double>::cleanup();
  Array<int>::cleanup();
  Array<bool>::cleanup();
  return result;
}
