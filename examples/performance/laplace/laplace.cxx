/*
 * Test performance of Laplacian inversion
 *
 */

#include <bout.hxx>
#include <invert_laplace.hxx>
#include <field_factory.hxx>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <time.h>
#include <vector>

using SteadyClock = std::chrono::time_point<std::chrono::steady_clock>;
using Duration = std::chrono::duration<double>;
using namespace std::chrono;

#define TEST_BLOCK(NAME, ...)                                                            \
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

  // Initialise BOUT++, setting up mesh
  BoutInitialise(argc, argv);
  std::vector<std::string> names;
  std::vector<Duration> times;

  // Get options root
  auto globalOptions = Options::root();
  auto modelOpts = globalOptions["LaplaceTest"];
  int NUM_LOOPS = modelOpts["NUM_LOOPS"].withDefault(1000);

  ConditionalOutput time_output(Output::getInstance());
  time_output.enable(true);

  FieldFactory f(mesh);

  Field3D input = f.create3D("(1-gauss(x-0.5,0.2))*gauss(y-pi)*gauss(z-pi)");
  Field2D a = f.create2D("gauss(x) * sin(y)");
  Field2D c = f.create2D("sin(x) * gauss(x-0.5) * gauss(y-pi)");
  Field2D d = f.create2D("y - pi/2");

  Field3D flag0;
  TEST_BLOCK("flag0",
     flag0 = invert_laplace(input, 0);
  );

  Field3D flag3;
  TEST_BLOCK("flag3",
     flag3 = invert_laplace(input, 3);
  );

  Field3D flag0a;
  TEST_BLOCK("flag0a",
     flag0a = invert_laplace(input, 0, &a);
  );

  Field3D flag3a;
  TEST_BLOCK("flag3a",
     flag3a = invert_laplace(input, 3, &a);
  );

  Field3D flag0ac;
  TEST_BLOCK("flag0ac",
     flag0ac = invert_laplace(input, 0, &a, &c);
  );

  Field3D flag3ac;
  TEST_BLOCK("flag3ac",
     flag3ac = invert_laplace(input, 3, &a, &c);
  );

  Field3D flag0ad;
  TEST_BLOCK("flag0ad",
     flag0ad = invert_laplace(input, 0, &a, nullptr, &d);
  );

  Field3D flag3ad;
  TEST_BLOCK("flag3ad",
     flag3ad = invert_laplace(input, 3, &a, nullptr, &d);
  );

  /// Test new interface and INVERT_IN/OUT_SET flags

  Field2D set_to = f.create2D("cos(2*y)*(x - 0.5)");

  Laplacian *lap = Laplacian::create();
  lap->setFlags(4096);
  Field3D flagis;
  TEST_BLOCK("flagis",
     flagis = lap->solve(input, set_to);
  );

  lap->setFlags(8192);
  Field3D flagos;
  TEST_BLOCK("flagos",
     flagos = lap->solve(input, set_to);
  );

  lap->setCoefA(a);
  lap->setFlags(4096);
  Field3D flagisa;
  TEST_BLOCK("flagisa",
     flagisa = lap->solve(input, set_to);
  );

  lap->setFlags(8192);
  Field3D flagosa;
  TEST_BLOCK("flagosa",
     flagosa = lap->solve(input, set_to);
  );

  lap->setCoefC(c);
  lap->setFlags(4096);
  Field3D flagisac;
  TEST_BLOCK("flagisac",
     flagisac = lap->solve(input, set_to);
  );

  lap->setFlags(8192);
  Field3D flagosac;
  TEST_BLOCK("flagosac",
     flagosac = lap->solve(input, set_to);
  );

  lap->setCoefC(1.0);
  lap->setCoefD(d);
  lap->setFlags(4096);
  Field3D flagisad;
  TEST_BLOCK("flagisad",
     flagisad = lap->solve(input, set_to);
  );

  lap->setFlags(8192);
  Field3D flagosad;
  TEST_BLOCK("flagosad",
     flagosad = lap->solve(input, set_to);
  );

  // Delete Laplacian when done
  delete lap;

  // Write and close the output file
  dump.write();
  dump.close();

  MPI_Barrier(BoutComm::get()); // Wait for all processors to write data

  // Report
  int width = 0;
  for (const auto i : names) {
    width = i.size() > width ? i.size() : width;
  };
  width = width + 5;
  time_output << std::setw(width) << "Case name"
              << "\t"
              << "Time (s)"
              << "\t"
              << "Time per iteration (s)"
              << "\n";
  for (int i = 0; i < names.size(); i++) {
    time_output << std::setw(width) << names[i] << "\t" << times[i].count() << "\t\t" << times[i].count() / NUM_LOOPS 
                << "\n";
  }
  double sum_of_times = 0.0;
  for (auto t : times) {
        sum_of_times += t.count();
  }
  time_output << std::setw(width) << "Total" << "\t" << sum_of_times << "\n";
  time_output << std::setw(width) << "Average" << "\t" << sum_of_times / times.size() << "\t\t" << sum_of_times / NUM_LOOPS / times.size()
                << "\n";

  BoutFinalise();
  return 0;
}
