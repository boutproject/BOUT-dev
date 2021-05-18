/*
 * Timing of arithmetic operations
 *
 */

#include <bout/physicsmodel.hxx>

#include <bout/expr.hxx>

#include <chrono>

using SteadyClock = std::chrono::time_point<std::chrono::steady_clock>;
using Duration = std::chrono::duration<double>;
using namespace std::chrono;

#define TIMEIT(elapsed, ...)                                                     \
  {                                                                              \
    SteadyClock start = steady_clock::now();                                     \
    { __VA_ARGS__; }                                                             \
    Duration diff = steady_clock::now() - start;                                 \
    diff *= 1000 * 1000;                                                         \
    elapsed.min = diff > elapsed.min ? elapsed.min : diff;                       \
    elapsed.max = diff < elapsed.max ? elapsed.max : diff;                       \
    elapsed.count++;                                                             \
    elapsed.avg = elapsed.avg * (1 - 1. / elapsed.count) + diff / elapsed.count; \
  }

struct Durations {
  Duration max;
  Duration min;
  Duration avg;
  int count;
};

class Arithmetic : public PhysicsModel {
protected:
  int init(bool) {

    Field3D a = 1.0;
    Field3D b = 2.0;
    Field3D c = 3.0;
    a.setRegion("RGN_ALL");
    b.setRegion("RGN_NOBNDRY");

    Field3D result1, result2, result3, result4;

    // Using Field methods (classic operator overloading)

    result1 = 2. * a + b * c;
#define dur_init {Duration::min(), Duration::max(), Duration::zero(), 0}
    Durations elapsed1 = dur_init, elapsed2 = dur_init, elapsed3 = dur_init,
              elapsed4 = dur_init;

    for (int ik = 0; ik < 1e3; ++ik) {
      TIMEIT(elapsed1, result1 = 2. * a + b * c;);

      // Using C loops
      result2.allocate();
      BoutReal *rd = &result2(0, 0, 0);
      BoutReal *ad = &a(0, 0, 0);
      BoutReal *bd = &b(0, 0, 0);
      BoutReal *cd = &c(0, 0, 0);
      TIMEIT(elapsed2,
             for (int i = 0, iend = (mesh->LocalNx * mesh->LocalNy * mesh->LocalNz) - 1;
                  i != iend; i++) {
               *rd = 2. * (*ad) + (*bd) * (*cd);
               rd++;
               ad++;
               bd++;
               cd++;
             });

      // Template expressions
      TIMEIT(elapsed3, result3 = eval3D(add(mul(2, a), mul(b, c))););

      // Range iterator
      result4.allocate();
      TIMEIT(elapsed4, for (auto i : result4) result4[i] = 2. * a[i] + b[i] * c[i];);
    }

    output.enable();
    output << "TIMING      |    minimum |       mean |    maximum\n"
           << "----------- | ---------- | ---------- | ----------\n";
    //#define PRINT(str,elapsed)   output << str << elapsed.min.count()<<
    //elapsed.avg.count()<< elapsed.max.count() << endl;
#define PRINT(str, elapsed)                                          \
  output.write("{:s} | {:7.3f} us | {:7.3f} us | {:7.3f} us\n", str, \
               elapsed.min.count(), elapsed.avg.count(), elapsed.max.count())
    PRINT("Fields:    ", elapsed1);
    PRINT("C loop:    ", elapsed2);
    PRINT("Templates: ", elapsed3);
    PRINT("Range For: ", elapsed4);
    output.disable();
    SOLVE_FOR(n);
    return 0;
  }

  int rhs(BoutReal) {
    ddt(n) = 0;
    return 0;
  }
  Field3D n;
};

BOUTMAIN(Arithmetic);
