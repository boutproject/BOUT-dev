/*
 * Testing performance of iterators over the mesh
 *
 */

#include <bout/bout.hxx>

#include <chrono>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <time.h>
#include <vector>

#include <bout/field_factory.hxx>
#include <bout/initialprofiles.hxx>

#include "bout/openmpwrap.hxx"
#include "bout/region.hxx"

using SteadyClock = std::chrono::time_point<std::chrono::steady_clock>;
using Duration = std::chrono::duration<double>;
using namespace std::chrono;
using bout::globals::mesh;

#define ITERATOR_TEST_BLOCK(NAME, ...)                                              \
  {                                                                                 \
    __VA_ARGS__                                                                     \
    names.push_back(NAME);                                                          \
    SteadyClock start = steady_clock::now();                                        \
    for (int repetitionIndex = 0; repetitionIndex < NUM_LOOPS; repetitionIndex++) { \
      __VA_ARGS__;                                                                  \
    }                                                                               \
    times.push_back(steady_clock::now() - start);                                   \
  }

int main(int argc, char** argv) {
  BoutInitialise(argc, argv);
  std::vector<std::string> names;
  std::vector<Duration> times;

  // Get options root
  auto& globalOptions = Options::root();
  auto& modelOpts = globalOptions["performance"];
  int NUM_LOOPS;
  NUM_LOOPS = modelOpts["NUM_LOOPS"].withDefault(100);
  bool profileMode, includeHeader, do2D3D, do3D3D;
  profileMode = modelOpts["profileMode"].withDefault(false);
  includeHeader = modelOpts["includeHeader"].withDefault(false);
  do2D3D = modelOpts["do2D3D"].withDefault(false);
  do3D3D = modelOpts["do3D3D"].withDefault(false);

  ConditionalOutput time_output(Output::getInstance());
  time_output.enable(true);
  const int len = mesh->LocalNx * mesh->LocalNy * mesh->LocalNz;

  Field3D a;
  initial_profile("a", a);
  Field3D b;
  initial_profile("b", b);
  Field2D c;
  initial_profile("c", c);

  Field3D result;
  result.allocate();

  if (do2D3D) {
    ITERATOR_TEST_BLOCK("Bracket [2D,3D] ARAKAWA",
                        result = bracket(a, c, BRACKET_ARAKAWA););

    ITERATOR_TEST_BLOCK("Bracket [2D,3D] ARAKAWA_OLD",
                        result = bracket(a, c, BRACKET_ARAKAWA_OLD););

    ITERATOR_TEST_BLOCK("Bracket [2D,3D] SIMPLE",
                        result = bracket(a, c, BRACKET_SIMPLE););

    ITERATOR_TEST_BLOCK("Bracket [2D,3D] DEFAULT", result = bracket(a, c, BRACKET_STD););
  }

  if (do3D3D) {
    ITERATOR_TEST_BLOCK("Bracket [3D,3D] ARAKAWA",
                        result = bracket(a, b, BRACKET_ARAKAWA););

    ITERATOR_TEST_BLOCK("Bracket [3D,3D] ARAKAWA_OLD",
                        result = bracket(a, b, BRACKET_ARAKAWA_OLD););

    ITERATOR_TEST_BLOCK("Bracket [3D,3D] SIMPLE",
                        result = bracket(a, b, BRACKET_SIMPLE););

    ITERATOR_TEST_BLOCK("Bracket [3D,3D] DEFAULT", result = bracket(a, b, BRACKET_STD););
  }

  // Uncomment below for a "correctness" check
  // Field3D resNew = bracket(a, b, BRACKET_ARAKAWA); mesh->communicate(resNew);
  // Field3D resOld = bracket(a, b, BRACKET_ARAKAWA_OLD); mesh->communicate(resOld);
  // time_output << "Max abs diff is
  // "<<max(abs(resNew-resOld),true)/max(abs(resOld),true)<<std::endl;

  if (profileMode) {
    int nthreads = 0;
#if BOUT_USE_OPENMP
    nthreads = omp_get_max_threads();
#endif

    int width = 12;
    if (includeHeader) {
      time_output << "\n------------------------------------------------\n";
      time_output << "Case legend";
      time_output << "\n------------------------------------------------\n";

      for (std::size_t i = 0; i < names.size(); i++) {
        time_output << std::setw(width) << "Case " << i << ".\t" << names[i] << "\n";
      }
      time_output << "\n";
      time_output << std::setw(width) << "Nprocs"
                  << "\t";
      time_output << std::setw(width) << "Nthreads"
                  << "\t";
      time_output << std::setw(width) << "Num_loops"
                  << "\t";
      time_output << std::setw(width) << "Local grid"
                  << "\t";
      time_output << std::setw(width) << "Nx (global)"
                  << "\t";
      time_output << std::setw(width) << "Ny (global)"
                  << "\t";
      time_output << std::setw(width) << "Nz (global)"
                  << "\t";
      for (std::size_t i = 0; i < names.size(); i++) {
        time_output << std::setw(width) << "Case " << i << "\t";
      }
      time_output << "\n";
    }

    time_output << std::setw(width) << BoutComm::size() << "\t";
    time_output << std::setw(width) << nthreads << "\t";
    time_output << std::setw(width) << NUM_LOOPS << "\t";
    time_output << std::setw(width) << len << "\t";
    time_output << std::setw(width) << mesh->GlobalNx << "\t";
    time_output << std::setw(width) << mesh->GlobalNy << "\t";
    time_output << std::setw(width) << mesh->GlobalNz << "\t";
    for (std::size_t i = 0; i < names.size(); i++) {
      time_output << std::setw(width) << times[i].count() / NUM_LOOPS << "\t";
    }
    time_output << "\n";
  } else {
    std::size_t width = 0;
    for (const auto i : names) {
      width = i.size() > width ? i.size() : width;
    };
    width = width + 5;
    time_output << std::setw(width) << "Case name"
                << "\t"
                << "Time per iteration (s)"
                << "\n";
    for (std::size_t i = 0; i < names.size(); i++) {
      time_output << std::setw(width) << names[i] << "\t" << times[i].count() / NUM_LOOPS
                  << "\n";
    }
  };

  BoutFinalise();
  return 0;
}
