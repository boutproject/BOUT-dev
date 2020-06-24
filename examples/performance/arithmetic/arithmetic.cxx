/*
 * Timing of arithmetic operations
 *
 */

#include <benchmark/benchmark.h>

#include <bout/physicsmodel.hxx>

#include <bout/expr.hxx>

//#include <field3d.hxx>

static void BM_Inline(benchmark::State& state) {
  Field3D a = 1.0;
  Field3D b = 2.0;
  Field3D c = 3.0;
  Field3D result;
  for (auto _ : state) {
    //benchmark::DoNotOptimize(result = 2. * a + b * c);
    result = 2. * a + b * c;
  }
}

BENCHMARK(BM_Inline);

/*
static void BM_C_loop(benchmark::State& state) {
  Field3D a = 1.0;
  Field3D b = 2.0;
  Field3D c = 3.0;
  Field3D result;

  result.allocate();
  BoutReal *rd = &result2(0, 0, 0);
  BoutReal *ad = &a(0, 0, 0);
  BoutReal *bd = &b(0, 0, 0);
  BoutReal *cd = &c(0, 0, 0);

  for (auto _ : state) {
    for (int i = 0, iend = (mesh->LocalNx * mesh->LocalNy * mesh->LocalNz) - 1;
                  i != iend; i++) {
      *rd = 2. * (*ad) + (*bd) * (*cd);
      rd++;
      ad++;
      bd++;
      cd++;
    }
  }
}

BENCHMARK(BM_C_loop);
*/

static void BM_Templates(benchmark::State& state) {
  Field3D a = 1.0;
  Field3D b = 2.0;
  Field3D c = 3.0;
  Field3D result;
  for (auto _ : state) {
    result = eval3D(add(mul(2, a), mul(b, c)));
  }
}

BENCHMARK(BM_Templates);

static void BM_Range_For(benchmark::State& state) {
  Field3D a = 1.0;
  Field3D b = 2.0;
  Field3D c = 3.0;
  Field3D result;
  result.allocate();
  for (auto _ : state) {
    for (auto i : result) {
      result[i] = 2. * a[i] + b[i] * c[i];
    }
  }
}

BENCHMARK(BM_Range_For);

// Manual version replacing 
//BENCHMARK_MAIN();
// to allow call of BoutInitialise and BoutFinalise. This needs improving,
// as it ignores command line arguments to Google benchmarks.
int main(int argc, char** argv) {                                     
  BoutInitialise(argc, argv);
  ::benchmark::Initialize(&argc, argv);                               
  //if (::benchmark::ReportUnrecognizedArguments(argc, argv)) return 1; 
  ::benchmark::RunSpecifiedBenchmarks();
  BoutFinalise();
}
