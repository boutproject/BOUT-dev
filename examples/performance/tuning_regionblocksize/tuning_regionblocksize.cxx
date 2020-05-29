/*
 * Testing performance of block region iterators over the mesh with
 * varying block sizes
 *
 */

#include <bout.hxx>

#include <algorithm>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <time.h>
#include <vector>

#include "bout/openmpwrap.hxx"
#include "bout/region.hxx"

using SteadyClock = std::chrono::time_point<std::chrono::steady_clock>;
using Duration = std::chrono::duration<double>;
using namespace std::chrono;
using namespace bout::globals;

#define ITERATOR_TEST_BLOCK(NAME, ...)                                                   \
  {                                                                                      \
    __VA_ARGS__                                                                          \
    names.push_back(NAME);                                                               \
    SteadyClock start = steady_clock::now();                                             \
    for (int repetitionIndex = 0; repetitionIndex < NUM_LOOPS; repetitionIndex++) {      \
      __VA_ARGS__;                                                                       \
    }                                                                                    \
    times.push_back(steady_clock::now() - start);                                        \
  }

int main(int argc, char **argv) {
  BoutInitialise(argc, argv);
  std::vector<std::string> names;
  std::vector<Duration> times;

  // Get options root
  auto globalOptions = Options::root();
  auto modelOpts = globalOptions["tuningRegionBlockSize"];
  const int NUM_LOOPS = modelOpts["NUM_LOOPS"].withDefault(100);
  const int numSteps = modelOpts["numSteps"].withDefault(16);

  ConditionalOutput time_output(Output::getInstance());
  time_output.enable(true);

  Field3D a = 1.0;
  Field3D b = 2.0;

  Field3D result;
  result.allocate();

  // Time simple task with different blocksizes
  for (int i = 0, blocksize = 1; i < numSteps; ++i, blocksize *= 2) {
    std::string name = "block size : " + std::to_string(blocksize);
    auto region =
        Region<Ind3D>(0, mesh->LocalNx - 1, 0, mesh->LocalNy - 1, 0, mesh->LocalNz - 1,
                      mesh->LocalNy, mesh->LocalNz, blocksize);

    ITERATOR_TEST_BLOCK(
        name, BOUT_FOR(i, region) { result[i] = a[i] + b[i]; });
  }

  // Report
  constexpr auto min_width = 5;
  const auto width =
      min_width
      + std::max_element(begin(names), end(names), [](const auto& a, const auto& b) {
          return a.size() < b.size();
        })->size();
  time_output << std::setw(width) << "Case name"
              << "\t"
              << "Time per iteration (s)"
              << "\n";
  for (std::size_t i = 0; i < names.size(); i++) {
    time_output << std::setw(width) << names[i] << "\t" << times[i].count() / NUM_LOOPS
                << "\n";
  }

  BoutFinalise();
  return 0;
}
