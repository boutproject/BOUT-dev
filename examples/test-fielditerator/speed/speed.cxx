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

#ifndef _OPENMP
#include <time.h>
#include <sys/time.h>
double omp_get_wtime(){
  struct timeval time;
  if (gettimeofday(&time,NULL)){
    //  Handle error
    return 0;
  }
  return (double)time.tv_sec + (double)time.tv_usec * .000001;
}
#else
#include <omp.h>
#endif

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
#pragma omp parallel for collapse(3)
    for (int jx=0;jx < mesh->ngx;++jx){
      for (int jy=0;jy < mesh->ngy;++jy){
	for (int jz=0;jz < mesh->ngz;++jz){
	  d3[jx][jy][jz]=myrand;
	}
      }
    }
  }
  double t3=omp_get_wtime();
  {
    const int mymax=mesh->ngx*mesh->ngy*mesh->ngz;
    double * data=d3[0][0];
    for (int i=0;i<iter;++i){
      double myrand=rand();
#pragma omp parallel
      for (int j=0;j < mymax;++j){
	data[j]=myrand;
      }
    }
  }
  double t4=omp_get_wtime();

  output.write("\n\n\n\t***** RESULTS *****\n\n  ***********************************\n");
  output.write("    FieldIteratorCIndex: %f\n",t2-t1);
  output.write("    triple for loop    : %f\n",t3-t2);
  output.write("    single for loop    : %f\n",t4-t3);
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
