/*
 * Test performance of Laplacian inversion
 *
 */

#include <bout/bout.hxx>
#include <bout/field_factory.hxx>
#include <bout/invert_laplace.hxx>

#include <algorithm>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <numeric>
#include <time.h>
#include <vector>

using SteadyClock = std::chrono::time_point<std::chrono::steady_clock>;
using Duration = std::chrono::duration<double>;
using namespace std::chrono;

#define TEST_BLOCK(NAME, ...)                                                       \
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

  // Initialise BOUT++, setting up mesh
  BoutInitialise(argc, argv);
  std::vector<std::string> names;
  std::vector<Duration> times;

  // Get options root
  auto& globalOptions = Options::root();
  auto& modelOpts = globalOptions["LaplaceTest"];
  int NUM_LOOPS = modelOpts["NUM_LOOPS"].withDefault(1000);

  ConditionalOutput time_output(Output::getInstance());
  time_output.enable(true);

  FieldFactory f{bout::globals::mesh};

  Field3D input = f.create3D("(1-gauss(x-0.5,0.2))*gauss(y-pi)*gauss(z-pi)");
  Field2D a = f.create2D("gauss(x) * sin(y)");
  Field2D c = f.create2D("sin(x) * gauss(x-0.5) * gauss(y-pi)");
  Field2D d = f.create2D("y - pi/2");

  auto lap = std::unique_ptr<Laplacian>{Laplacian::create()};

  lap->setCoefA(0.0);
  lap->setCoefC(1.0);
  lap->setCoefD(1.0);
  Field3D flag0;
  TEST_BLOCK("flag0", flag0 = lap->solve(input););
  lap->setInnerBoundaryFlags(INVERT_DC_GRAD + INVERT_AC_GRAD);
  Field3D flag3;
  TEST_BLOCK("flag3", flag3 = lap->solve(input););

  lap->setCoefA(a);
  lap->setInnerBoundaryFlags(0);
  Field3D flag0a;
  TEST_BLOCK("flag0a", flag0a = lap->solve(input););
  lap->setInnerBoundaryFlags(INVERT_DC_GRAD + INVERT_AC_GRAD);
  Field3D flag3a = lap->solve(input);
  TEST_BLOCK("flag3a", flag3a = lap->solve(input););

  lap->setCoefC(c);
  lap->setInnerBoundaryFlags(0);
  Field3D flag0ac;
  TEST_BLOCK("flag0ac", flag0ac = lap->solve(input););
  lap->setInnerBoundaryFlags(INVERT_DC_GRAD + INVERT_AC_GRAD);
  Field3D flag3ac;
  TEST_BLOCK("flag3ac", flag3ac = lap->solve(input););

  lap->setCoefC(1.0);
  lap->setCoefD(d);
  lap->setInnerBoundaryFlags(0);
  Field3D flag0ad;
  TEST_BLOCK("flag0ad", flag0ad = lap->solve(input););
  lap->setInnerBoundaryFlags(INVERT_DC_GRAD + INVERT_AC_GRAD);
  Field3D flag3ad;
  TEST_BLOCK("flag3ad", flag3ad = lap->solve(input););

  /// Test new interface and INVERT_IN/OUT_SET flags

  Field2D set_to = f.create2D("cos(2*y)*(x - 0.5)");

  lap->setCoefA(0.0);
  lap->setCoefC(1.0);
  lap->setCoefD(1.0);

  lap->setInnerBoundaryFlags(INVERT_SET);
  Field3D flagis;
  TEST_BLOCK("flagis", flagis = lap->solve(input, set_to););
  lap->setInnerBoundaryFlags(0);
  lap->setOuterBoundaryFlags(INVERT_SET);
  Field3D flagos;
  TEST_BLOCK("flagos", flagos = lap->solve(input, set_to););

  lap->setCoefA(a);
  lap->setInnerBoundaryFlags(INVERT_SET);
  lap->setOuterBoundaryFlags(0);
  lap->setOuterBoundaryFlags(0);
  Field3D flagisa;
  TEST_BLOCK("flagisa", flagisa = lap->solve(input, set_to););
  lap->setInnerBoundaryFlags(0);
  lap->setOuterBoundaryFlags(INVERT_SET);
  Field3D flagosa;
  TEST_BLOCK("flagosa", flagosa = lap->solve(input, set_to););

  lap->setCoefC(c);
  lap->setInnerBoundaryFlags(INVERT_SET);
  lap->setOuterBoundaryFlags(0);
  Field3D flagisac;
  TEST_BLOCK("flagisac", flagisac = lap->solve(input, set_to););
  lap->setInnerBoundaryFlags(0);
  lap->setOuterBoundaryFlags(INVERT_SET);
  Field3D flagosac;
  TEST_BLOCK("flagosac", flagosac = lap->solve(input, set_to););

  lap->setCoefC(1.0);
  lap->setCoefD(d);
  lap->setInnerBoundaryFlags(INVERT_SET);
  lap->setOuterBoundaryFlags(0);
  Field3D flagisad;
  TEST_BLOCK("flagisad", flagisad = lap->solve(input, set_to););
  lap->setInnerBoundaryFlags(0);
  lap->setOuterBoundaryFlags(INVERT_SET);
  Field3D flagosad;
  TEST_BLOCK("flagosad", flagosad = lap->solve(input, set_to););

  MPI_Barrier(BoutComm::get()); // Wait for all processors to write data

  // Report
  constexpr auto min_width = 5;
  const auto width =
      min_width
      + std::max_element(begin(names), end(names), [](const auto& a, const auto& b) {
          return a.size() < b.size();
        })->size();
  time_output << std::setw(width) << "Case name"
              << "\t"
              << "Time (s)"
              << "\t"
              << "Time per iteration (s)"
              << "\n";
  for (std::size_t i = 0; i < names.size(); i++) {
    time_output << std::setw(width) << names[i] << "\t" << times[i].count() << "\t\t"
                << times[i].count() / NUM_LOOPS << "\n";
  }
  const auto sum_of_times =
      std::accumulate(begin(times), end(times), 0.0,
                      [](const auto& sum, const auto& t) { return sum + t.count(); });
  time_output << std::setw(width) << "Total"
              << "\t" << sum_of_times << "\n";
  time_output << std::setw(width) << "Average"
              << "\t" << sum_of_times / times.size() << "\t\t"
              << sum_of_times / NUM_LOOPS / times.size() << "\n";

  BoutFinalise();
  return 0;
}
