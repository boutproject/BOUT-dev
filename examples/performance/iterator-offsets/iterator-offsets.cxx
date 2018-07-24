/*
 * Testing performance of iterators over the mesh
 *
 */

#include <bout.hxx>

#include <chrono>
#include <iostream>
#include <iterator>
#include <time.h>
#include <vector>
#include <iomanip>

#include "bout/openmpwrap.hxx"
#include "bout/region.hxx"
#include <bout/indexoffset.hxx>

typedef std::chrono::time_point<std::chrono::steady_clock> SteadyClock;
typedef std::chrono::duration<double> Duration;
using namespace std::chrono;

#define ITERATOR_TEST_BLOCK(NAME, ...)		\
  {__VA_ARGS__								\
      names.push_back(NAME);						\
    SteadyClock start = steady_clock::now();				\
    for (int repetitionIndex = 0 ; repetitionIndex < NUM_LOOPS ; repetitionIndex++){ \
      __VA_ARGS__;							\
    }									\
    times.push_back(steady_clock::now()-start);				\
  }

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

  ConditionalOutput time_output(Output::getInstance());
  time_output.enable(true);
  
  Field3D a = 1.0;
  Field3D b = 2.0;

  Field3D result;
  result.allocate();

  //Some small setup for C loop case
  BoutReal *ad = &a(0,0,0);
  BoutReal *bd = &b(0,0,0);
  BoutReal *rd = &result(0,0,0);
  
  const int len = mesh->LocalNx*mesh->LocalNy*mesh->LocalNz;

  //Raw C loop
///  ITERATOR_TEST_BLOCK("C loop",
///		      for(int j=0;j<len;++j) {
///			rd[j] = ad[j] + bd[j];
///		      };
///		      );
#ifdef _OPENMP  
///  ITERATOR_TEST_BLOCK("C loop (omp)",
///		      BOUT_OMP(parallel for)
///		      for(int j=0;j<len;++j) {
///			rd[j] = ad[j] + bd[j];
///		      };
///		    );
#endif
  
  // Nested loops over block data
  ITERATOR_TEST_BLOCK("Nested loop",
		    for(int i=0;i<mesh->LocalNx;++i) {
		      for(int j=mesh->ystart;j<mesh->yend;++j) {
			for(int k=0;k<mesh->LocalNz;++k) {
			  result(i,j,k) = (a(i,j+1,k) - a(i,j-1,k))/(2.*mesh->coordinates()->dy(i,j));
			}
		      }
		    }
		    );

#ifdef _OPENMP  
  ITERATOR_TEST_BLOCK("Nested loop (omp)",
		      BOUT_OMP(parallel for)
		      for(int i=0;i<mesh->LocalNx;++i) {
			for(int j=mesh->ystart;j<mesh->yend;++j) {
			  for(int k=0;k<mesh->LocalNz;++k) {
			    result(i,j,k) = (a(i,j+1,k) - a(i,j-1,k))/(2.*mesh->coordinates()->dy(i,j));
			  }
		      }
		      }
		      );
#endif
  
  // DataIterator using begin(), end()
///  IndexRange r{0, mesh->LocalNx, mesh->ystart, mesh->yend, 0, mesh->LocalNz};
///  ITERATOR_TEST_BLOCK("DI begin/end",
///		    for(DataIterator i = r.begin(), rend=r.end(); i != rend; ++i){
///		      result[i] = (a[i.yp()] - a[i.ym()])/(2.*mesh->coordinates()->dy[i]);
///		    }
///		    );

  // DataIterator with done()
///  IndexRange r{0, mesh->LocalNx, mesh->ystart, mesh->yend, 0, mesh->LocalNz};
///  ITERATOR_TEST_BLOCK("DI begin/done",
///		    for(DataIterator i = r.begin(); !i.done() ; ++i){
///		      result(i.x,i.y,i.z);// = a(i.x,i.y,i.z) + b(i.x,i.y,i.z);
///		    }
///		    );
  
  // Range based for DataIterator with indices
  ITERATOR_TEST_BLOCK("C++11 range-based for",
                    BOUT_OMP(parallel)
		    for(auto i : result.region(RGN_NOY)){
		      result(i.x,i.y,i.z) = (a(i.x,i.y+1,i.z) - a(i.x,i.y-1,i.z))/(2.*mesh->coordinates()->dy(i.x,i.y));
		    }
		    );

  // Range based DataIterator 
  ITERATOR_TEST_BLOCK("C++11 range-based for [i]", 
                    BOUT_OMP(parallel)
		    for(const auto &i : result.region(RGN_NOY)){
		      result[i] = (a[i.yp()] - a[i.ym()])/(2.*mesh->coordinates()->dy[i]);
		    }
		    );
  
  // DataIterator over fields
///  ITERATOR_TEST_BLOCK("DI (done) [i]",
///                    BOUT_OMP(parallel)
///		    for(DataIterator d = result.iterator(); !d.done(); d++)
///		      result[d] = a[d] + b[d];
///		    );

  //Raw C loop
///  ITERATOR_TEST_BLOCK("C loop repeat",
///		    for(int j=0;j<len;++j) {
///		      rd[j] = ad[j] + bd[j];
///		    };
///		    );

  // Region macro
  stencil s;
  IndexOffset<Ind3D> offset(*mesh);
  ITERATOR_TEST_BLOCK(
      "Region (serial)",
      BLOCK_REGION_LOOP_SERIAL(mesh->getRegion3D("RGN_NOY"), i,
        s.mm = nan("");
        s.m = a[offset.ym(i)];
        s.c = a[i];
        s.p = a[offset.yp(i)];
        s.pp = nan("");

        result[i] = (s.p - s.m)/(2.*mesh->coordinates()->dy[i]);
);
		      );
#ifdef _OPENMP
  ITERATOR_TEST_BLOCK("Region (omp)",
BOUT_OMP(parallel)
{
      BLOCK_REGION_LOOP_PARALLEL_SECTION(mesh->getRegion3D("RGN_NOY"), i,
        s.mm = nan("");
        s.m = a[offset.ym(i)];
        s.c = a[i];
        s.p = a[offset.yp(i)];
        s.pp = nan("");

        result[i] = (s.p - s.m)/(2.*mesh->coordinates()->dy[i]);
);
}
		      );
#endif
  
  if(profileMode){
    int nthreads=0;
#ifdef _OPENMP
    nthreads = omp_get_max_threads();
#endif

    int width = 12;
    if(includeHeader){
      time_output << "\n------------------------------------------------\n";
      time_output << "Case legend";
      time_output <<"\n------------------------------------------------\n";
      
      for (int i = 0 ; i < names.size(); i++){	
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
      for (int i = 0 ; i < names.size(); i++){	
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
    for (int i = 0 ; i < names.size(); i++){	
      time_output << std::setw(width) << times[i].count()/NUM_LOOPS << "\t";
    }
    time_output << "\n";
  }else{
    int width = 0;
    for (const auto i: names){ width = i.size() > width ? i.size() : width;};
    width = width + 5;
    time_output << std::setw(width) << "Case name" << "\t" << "Time per iteration (s)" << "\n";
    for(int i = 0 ; i < names.size(); i++){
      time_output <<  std::setw(width) << names[i] << "\t" << times[i].count()/NUM_LOOPS << "\n";
    }
  };

  BoutFinalise();
  return 0;
}
