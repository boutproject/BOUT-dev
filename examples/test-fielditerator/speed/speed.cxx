/*
 * Test Cyclic Reduction parallel solver
 *
 */

#include <bout.hxx>
#include <boutmain.hxx>

#include <fielditerator.hxx>
#include <assert.h>

#include <unistd.h>
#include <stdlib.h>

Field3D n;


int physics_init(bool restarting) {

  
  Field3D d3=1;
  

  const int iter=1e3;
  double t1=omp_get_wtime();
  for (int i=0;i<iter;++i){
    double myrand=rand();
#pragma omp parallel
    for (FieldIteratorCIndex cxit( *mesh);cxit;cxit.next3()){
      d3[cxit]= myrand;
    }
  }
  double t2=omp_get_wtime();
  for (int i=0;i<iter;++i){
    double myrand=rand();
#pragma omp parallel
    for (int jx=0;jx < mesh->ngx;++jx){
      for (int jy=0;jy < mesh->ngy;++jy){
	for (int jz=0;jz < mesh->ngz;++jz){
	  d3[jx][jy][jz]=myrand;
	}
      }
    }
  }
  double t3=omp_get_wtime();

  output.write("\n\n\n\t***** RESULTS *****\n\n  ***********************************\n");
  output.write("    FieldIteratorCIndex: %f\n",t2-t1);
  output.write("    triple for loop    : %f\n",t3-t2);
  output.write("  ***********************************\n\n\n\n\n\n");

  
  n=1;
  SOLVE_FOR(n);
  return 0;
}

int physics_run(BoutReal t) {
  // Doesn't do anything
  ddt(n)=0;
  return 0;
}
