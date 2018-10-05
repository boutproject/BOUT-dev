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

Mesh::deriv_func func_ptr = &test_ddy;

int main(int argc, char **argv) {
  BoutInitialise(argc, argv);
  std::vector<std::string> names;
  std::vector<Duration> times;

  //Get options root
  Options *globalOptions = Options::getRoot();
  Options *modelOpts = globalOptions->getSection("performanceIterator");
  int NUM_LOOPS;
  OPTION(modelOpts, NUM_LOOPS, 100);
  bool profileMode, includeHeader;
  OPTION(modelOpts, profileMode, false);
  OPTION(modelOpts, includeHeader, false);

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

  // Range based for DataIterator with indices
  ITERATOR_TEST_BLOCK(
    "C++11 range-based for (omp)",
    BOUT_OMP(parallel)
    for(auto i : result.region(RGN_NOY)){
      result(i.x,i.y,i.z) = (a(i.x,i.y+1,i.z) - a(i.x,i.y-1,i.z));
    }
    );

  // Range based DataIterator
  ITERATOR_TEST_BLOCK(
    "C++11 range-based for [i] (omp)",
    BOUT_OMP(parallel)
    for(const auto &i : result.region(RGN_NOY)){
      result[i] = (a[i.yp()] - a[i.ym()]);
    }
    );

  ITERATOR_TEST_BLOCK(
    "C++11 range-based for over Region [i] (omp)",
    BOUT_OMP(parallel)
    {
      for(const auto &i : mesh->getRegion3D("RGN_NOY")){
        result[i] = (a[i.yp()] - a[i.ym()]);
      }
    }
    );

  ITERATOR_TEST_BLOCK(
    "C++11 range-based for [i] with stencil (omp)",
    BOUT_OMP(parallel)
    {
      stencil s;
      for(const auto &i : result.region(RGN_NOY)){
        s.mm = nan("");
        s.m = a[i.ym()];
        s.c = a[i];
        s.p = a[i.yp()];
        s.pp = nan("");
        result[i] = (s.p - s.m);
      }
    }
    );
#endif

  ITERATOR_TEST_BLOCK(
    "C++11 range-based for [i] with stencil (serial)",
    stencil s;
    for(const auto &i : result.region(RGN_NOY)){
      s.mm = nan("");
      s.m = a[i.ym()];
      s.c = a[i];
      s.p = a[i.yp()];
      s.pp = nan("");
      result[i] = (s.p - s.m);
    }
    );

  // Region macro
  ITERATOR_TEST_BLOCK(
      "Region with stencil (serial)",
      stencil s;
      BOUT_FOR_SERIAL(i, mesh->getRegion3D("RGN_NOY")) {
        s.mm = nan("");
        s.m = a[i.ym()];
        s.c = a[i];
        s.p = a[i.yp()];
        s.pp = nan("");

        result[i] = (s.p - s.m);
      }
    );

#ifdef _OPENMP

  ITERATOR_TEST_BLOCK(
    "Region with stencil (parallel section nowait omp)",
    BOUT_OMP(parallel)
    {
      stencil s;
      BOUT_FOR_INNER(i, mesh->getRegion3D("RGN_NOY")) {
        s.mm = nan("");
        s.m = a[i.ym()];
        s.c = a[i];
        s.p = a[i.yp()];
        s.pp = nan("");

        result[i] = (s.p - s.m);
      }
    }
    );

   ITERATOR_TEST_BLOCK(
    "Region with stencil (parallel section wait omp)",
    BOUT_OMP(parallel)
    {
      stencil s;
      BOUT_FOR_OMP(i, mesh->getRegion3D("RGN_NOY"), for schedule(guided)) {
        s.mm = nan("");
        s.m = a[i.ym()];
        s.c = a[i];
        s.p = a[i.yp()];
        s.pp = nan("");

        result[i] = (s.p - s.m);
      }
    }
    );

  ITERATOR_TEST_BLOCK(
    "Region with stencil & raw ints (parallel section omp)",
    const auto nz = mesh->LocalNz;
    BOUT_OMP(parallel)
    {
      stencil s;
      BOUT_FOR_INNER(i, mesh->getRegion3D("RGN_NOY")) {
        s.mm = nan("");
        s.m = a[i.ind - (nz)];
        s.c = a[i.ind];
        s.p = a[i.ind + (nz)];
        s.pp = nan("");

        result[i] = (s.p - s.m);
      }
    }
    );

  ITERATOR_TEST_BLOCK(
    "Region with stencil (single loop omp)",
    BOUT_OMP(parallel)
    {
      stencil s;
      const auto &region = mesh->getRegion3D("RGN_NOY").getIndices();
      BOUT_OMP(for schedule(guided))
        for (auto i = region.cbegin(); i < region.cend(); ++i) {
          s.mm = nan("");
          s.m = a[i->ym()];
          s.c = a[*i];
          s.p = a[i->yp()];
          s.pp = nan("");

          result[*i] = (s.p - s.m);
        }
    }
    );

  ITERATOR_TEST_BLOCK(
    "Region with stencil and function pointer (parallel section nowait omp)",
    BOUT_OMP(parallel)
    {
      stencil s;
      BOUT_FOR_INNER(i, mesh->getRegion3D("RGN_NOY")) {
        s.mm = nan("");
        s.m = a[i.ym()];
        s.c = a[i];
        s.p = a[i.yp()];
        s.pp = nan("");

        result[i] = func_ptr(s);
      }
    }
    );

#endif

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
