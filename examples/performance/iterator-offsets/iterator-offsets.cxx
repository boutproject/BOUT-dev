/*
 * Testing performance of iterators over the mesh
 *
 */

#include <bout.hxx>

#include <algorithm>
#include <chrono>
#include <iostream>
#include <iterator>
#include <vector>
#include <iomanip>

#include "bout/openmpwrap.hxx"
#include "bout/region.hxx"

using SteadyClock = std::chrono::time_point<std::chrono::steady_clock>;
using Duration = std::chrono::duration<double>;
using namespace std::chrono;

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

BoutReal test_ddy(stencil &s) {
  return (s.p - s.m);
}

using deriv_func = BoutReal (*)(stencil &);
deriv_func func_ptr = &test_ddy;

int main(int argc, char **argv) {
  BoutInitialise(argc, argv);
  std::vector<std::string> names;
  std::vector<Duration> times;

  //Get options root
  auto globalOptions = Options::root();
  auto modelOpts = globalOptions["performanceIterator"];
  int NUM_LOOPS;
  NUM_LOOPS = modelOpts["NUM_LOOPS"].withDefault(100);
  bool profileMode, includeHeader;
  profileMode = modelOpts["profileMode"].withDefault(false);
  includeHeader = modelOpts["includeHeader"].withDefault(false);

  ConditionalOutput time_output{Output::getInstance()};
  time_output.enable(true);

  const Field3D a{1.0};
  const Field3D b{2.0};

  Field3D result;
  result.allocate();

  const int len = mesh->LocalNx*mesh->LocalNy*mesh->LocalNz;

  // Nested loops over block data
  ITERATOR_TEST_BLOCK(
    "Nested loop",
    for(int i=0;i<mesh->LocalNx;++i) {
      for(int j=mesh->ystart;j<mesh->yend;++j) {
        for(int k=0;k<mesh->LocalNz;++k) {
          result(i,j,k) = (a(i,j+1,k) - a(i,j-1,k));
        }
      }
    }
    );

#ifdef _OPENMP
  ITERATOR_TEST_BLOCK(
    "Nested loop (omp)",
    BOUT_OMP(parallel for)
    for(int i=0;i<mesh->LocalNx;++i) {
      for(int j=mesh->ystart;j<mesh->yend;++j) {
        for(int k=0;k<mesh->LocalNz;++k) {
          result(i,j,k) = (a(i,j+1,k) - a(i,j-1,k));
        }
      }
    }
    );

#endif

  {
    auto deriv = DerivativeStore<Field3D>::getInstance().getStandardDerivative("C2",DIRECTION::Y, STAGGER::None);
    ITERATOR_TEST_BLOCK(
			"DerivativeStore without fetching",
			deriv(a, result, "RGN_NOY");
			);
  };
  
  ITERATOR_TEST_BLOCK(
		      "DerivativeStore with fetch",
		      auto deriv = DerivativeStore<Field3D>::getInstance().getStandardDerivative("C2",DIRECTION::Y, STAGGER::None);    
		      deriv(a, result, "RGN_NOY");
		      );

  ITERATOR_TEST_BLOCK(
		      "Region with stencil",
		      BOUT_OMP(parallel)
		      {
			stencil s;
			BOUT_FOR_INNER(i, mesh->getRegion3D("RGN_NOY")) {
			  s.m = a[i.ym()];
			  s.c = a[i];
			  s.p = a[i.yp()];
	
			  result[i] = (s.p - s.m);
			}
		      }
		      );
  
  ITERATOR_TEST_BLOCK(
		      "Region with stencil and function pointer",
		      BOUT_OMP(parallel)
		      {
			stencil s;
			BOUT_FOR_INNER(i, mesh->getRegion3D("RGN_NOY")) {
			  s.m = a[i.ym()];
			  s.c = a[i];
			  s.p = a[i.yp()];

			  result[i] = func_ptr(s);
			}
		      }
		      );

   if (profileMode) {
     int nthreads = 0;
#ifdef _OPENMP
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
       time_output << std::setw(width) << "Nprocs" << "\t";
       time_output << std::setw(width) << "Nthreads" << "\t";
       time_output << std::setw(width) << "Num_loops" << "\t";
       time_output << std::setw(width) << "Local grid" << "\t";
       time_output << std::setw(width) << "Nx (global)" << "\t";
       time_output << std::setw(width) << "Ny (global)" << "\t";
       time_output << std::setw(width) << "Nz (global)" << "\t";
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
     for (const auto &time : times) {
       time_output << std::setw(width) << time.count() / NUM_LOOPS << "\t";
     }
     time_output << "\n";
   } else {
     std::vector<int> sizes(names.size());
     std::transform(names.begin(), names.end(), sizes.begin(),
                    [](const std::string &name) -> int { return static_cast<int>(name.size()); });
     int width = *std::max_element(sizes.begin(), sizes.end());
     width += 5;
     time_output << std::setw(width) << "Case name" << "\t" << "Time per iteration (s)" << "\n";
     for (std::size_t i = 0; i < names.size(); i++) {
       time_output << std::setw(width) << names[i] << "\t" << times[i].count() / NUM_LOOPS
                   << "\n";
     }
   };

   BoutFinalise();
   return 0;
}
