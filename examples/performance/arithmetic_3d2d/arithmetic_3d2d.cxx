/*
 * Timing of arithmetic operations (Field3D/Field2D mixed)
 *
 */

#include <bout/physicsmodel.hxx>

#include <bout/expr.hxx>

#include <chrono>
#include <iomanip>

using SteadyClock = std::chrono::time_point<std::chrono::steady_clock>;
using Duration = std::chrono::duration<double>;
using namespace std::chrono;

#define TIMEIT(NAME, ...)                                                             \
  {                                                                                      \
    SteadyClock start = steady_clock::now();                                             \
    __VA_ARGS__                                                                          \
    Duration diff = steady_clock::now() - start;                                         \
    auto elapsed = elapsedMap[NAME];\
    elapsed.min = diff > elapsed.min ? elapsed.min : diff;                               \
    elapsed.max = diff < elapsed.max ? elapsed.max : diff;                               \
    elapsed.count++;                                                                     \
    elapsed.avg = elapsed.avg * (1 - 1. / elapsed.count) + diff / elapsed.count;         \
    elapsedMap[NAME] = elapsed;						\
  }

struct Durations {
  Duration max;
  Duration min;
  Duration avg;
  int count;
  Durations()
      : max(Duration::min()), min(Duration::max()), avg(Duration::zero()), count(0){};
};

class Arithmetic : public PhysicsModel {
protected:
  std::map<std::string, Durations> elapsedMap;

  int init(bool) {
    Field3D a = 1.0;
    Field3D b = 2.0;
    Field2D c = 3.0;

    Field3D result1, result2, result3, result4;

    // Using Field methods (classic operator overloading)
    result1 = 2. * a + b * c;

    for (int ik = 0; ik < 1e2; ++ik) {
      result1.allocate();
      TIMEIT("Fields", result1 = 2. * a + b * c;);

      // Using C loops
      result2.allocate();
      BoutReal *rd = &result2(0, 0, 0);
      BoutReal *ad = &a(0, 0, 0);
      BoutReal *bd = &b(0, 0, 0);
      BoutReal *cd = &c(0, 0, 0);
      TIMEIT("C loop",
             for (int i = 0, iend = (mesh->LocalNx * mesh->LocalNy) -1;
                  i != iend; i++) {
	       for (int j = 0, jend =  mesh->LocalNz - 1;  j != jend; j++){
		 *rd = 2. * (*ad) + (*bd) * (*cd);
		 rd++;
		 ad++;
		 bd++;
	       }
               cd++;
             });

      // Template expressions
      result3.allocate();
      TIMEIT("Templates", result3 = eval3D(add(mul(2, a), mul(b, c))););

      // Range iterator
      result4.allocate();
      TIMEIT("Range For", for (auto i : result4) result4[i] = 2. * a[i] + b[i] * c[i];);
    }

    output.enable();
    int width = 15; 
    output<<std::setw(width)<<"TIMING";
    output<<std::setw(width)<<"min";
    output<<std::setw(width)<<"avg";
    output<<std::setw(width)<<"max";
    output <<"\n======";
    for(int i=0;i<4*width;++i){output<<"=";};
    output<<"\n";
    
    for(const auto &approach: elapsedMap){
      output<<std::setw(width)<<approach.first;
      output<<std::setw(width)<<approach.second.min.count();
      output<<std::setw(width)<<approach.second.avg.count();
      output<<std::setw(width)<<approach.second.max.count();
      output<<"\n";
    }
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
