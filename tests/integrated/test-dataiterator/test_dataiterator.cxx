/*
 * Test Cyclic Reduction parallel solver
 *
 */

#include <bout.hxx>
#include <bout/boutmain.hxx>

//#include <fielditerator.hxx>
#include <assert.h>

#include <unistd.h>

Field3D n;

#ifdef _OPENMP
#define PRINT_DEBUG {\
    usleep(2e5*omp_get_thread_num());					\
    output.write("FieldIterator failed at ???");			\
    output.write("- expected %.1f but got %.1f. ",exp,d3[i]);		\
    output.write("I am %d !\n",omp_get_thread_num());			\
    usleep(2e5*omp_get_num_threads());					\
    abort(); }
#else
#define PRINT_DEBUG {							\
    output.write("FieldIterator failed at ");				\
    output.write("- expected %.1f but got %.1f. ",exp,d3[i]);		\
    abort(); }
#endif

int spread_work(int num_work, int thread, int max_thread);

int physics_init(bool restarting) {

  // for (int w=2;w<8;w++){
  //   int we=pow(2,w)+1.1;
  //   output.write("spreading %d\n",we);
  //   for (int t=1;t<5;++t){
  //     for (int ct=0;ct<=t;++ct){
  // output.write("%d ",spread_work(we,ct,t));
  //     }
  //     output.write("\n");
  //   }
  // }


  /**********************************************
   * * *                                    * * *
   * * *          Field3D test              * * *
   * * *                                    * * *
   **********************************************/

  Field3D d3=1;

  BoutReal exp=1;
  // for (FieldIteratorCIndex cxit( *mesh);cxit;cxit.next3()){
  //   if (d3[cxit] != exp) PRINT_DEBUG;
  //   d3[cxit]=2;
  // }
  // exp=2;

#pragma omp parallel
  for (auto i: d3){
    //output.print("%d %d %d\n",i.x,i.y,i.z);
    if (d3[i] != exp) PRINT_DEBUG;
    //if (i.x != 8)
    d3[i]=3;
  }
  exp=3;

  for (int jx=0;jx<mesh->LocalNx;++jx)
    for (int jy=0;jy<mesh->LocalNy;++jy)
      for (int jz=0;jz<mesh->LocalNz;++jz){
	if (d3(jx,jy,jz) != exp) {
	  output.print("%d %d %d: expected %g got %g",jx,jy,jz,exp,d3(jx,jy,jz)) ;
	};
	assert(d3(jx,jy,jz)==exp);
      }

  /*
#pragma omp parallel
  for (auto i: d3){
    if (d3[i] != exp) PRINT_DEBUG;
    d3[i]=4;
  }

  for (int jx=0;jx<mesh->LocalNx;++jx){
    if (jx < mesh->xstart || jx > mesh->xend)
      exp=3;
    else
      exp=4;
    for (int jy=0;jy<mesh->LocalNy;++jy)
      for (int jz=0;jz<mesh->LocalNz;++jz)
	assert(d3(jx,jy,jz)==exp);
  }
  d3=3;
  #pragma omp parallel
  for (auto i : d3){
    if (d3[i] != exp) PRINT_DEBUG;
    d3[i]=4;
  }
  for (int jx=0;jx<mesh->LocalNx;++jx)
    for (int jy=0;jy<mesh->LocalNy;++jy){
      if (jy < mesh->ystart || jy > mesh->yend)
	exp=3;
      else
	exp=4;
      for (int jz=0;jz<mesh->LocalNz;++jz)
	assert(d3(jx,jy,jz)==exp);
    }
  d3=3;exp=3;
  #pragma omp parallel
  for (auto i: d3){
    if (d3[i] != exp) PRINT_DEBUG;
    d3[i]=4;
  }
  for (int jx=0;jx<mesh->LocalNx;++jx)
    for (int jy=0;jy<mesh->LocalNy;++jy)
      for (int jz=0;jz<mesh->LocalNz;++jz){
	if (jz == mesh->LocalNz)
	  exp=3;
	else
	  exp=4;
	assert(d3(jx,jy,jz)==exp);
      }
  d3=3;exp=3;
  #pragma omp parallel
  for (auto i:d3){
    if (d3[i] != exp) PRINT_DEBUG;
    d3[i]=4;
  }
  for (int jx=0;jx<mesh->LocalNx;++jx)
    for (int jy=0;jy<mesh->LocalNy;++jy)
      for (int jz=0;jz<mesh->LocalNz;++jz){
	if (jz == mesh->LocalNz ||
	    jy < mesh->ystart || jy > mesh->yend ||
	    jx < mesh->xstart || jx > mesh->xend)
	  exp=3;
	else
	  exp=4;
	assert(d3(jx,jy,jz)==exp);
      }
*/

  /**********************************************
   * * *                                    * * *
   * * *          Field2D test              * * *
   * * *                                    * * *
   **********************************************/

  Field2D d2=1;
  exp=1;
  #pragma omp parallel
  for (auto i: d2){
    if (d2[i] != exp) PRINT_DEBUG;
    d2[i]=3;
  }
  exp=3;
  for (int jx=0;jx<mesh->LocalNx;++jx)
    for (int jy=0;jy<mesh->LocalNy;++jy)
      assert(d2(jx,jy)==exp);
  /*
  d2=1;exp=1;
  for (auto i : d2){
    if (d2[i] != exp) PRINT_DEBUG;
    d2[ i]=2;
  }
  for (int jx=0;jx<mesh->LocalNx;++jx){
    if (jx < mesh->xstart || jx > mesh->xend)
      exp=1;
    else
      exp=2;
    for (int jy=0;jy<mesh->LocalNy;++jy)
      assert(d2(jx,jy)==exp);
  }
  d2=1;exp=1;
  for (auto i:d2){
    if (d2[i] != exp) PRINT_DEBUG;
    d2[i]=2;
  }
  for (int jx=0;jx<mesh->LocalNx;++jx)
    for (int jy=0;jy<mesh->LocalNy;++jy){
      if (jy < mesh->ystart || jy > mesh->yend)
	exp=1;
      else
	exp=2;
      assert(d2(jx,jy)==exp);
    }
  d2=1;exp=1;
  for (auto i:d2){
    if (d2[i] != exp) PRINT_DEBUG;
    d2[i]=2;
  }
  for (int jx=0;jx<mesh->LocalNx;++jx)
    for (int jy=0;jy<mesh->LocalNy;++jy){
      if (jx < mesh->xstart || jx > mesh->xend ||
	  jy < mesh->ystart || jy > mesh->yend)
	exp=1;
      else
	exp=2;
      assert(d2(jx,jy)==exp);
    }

  */

  bool fail;
  OPTION(Options::getRoot(), fail,                  false) ;
  n=1;
  SOLVE_FOR(n);
  if (fail){
    return 1;
  }
  return 0;
}

int physics_run(BoutReal t) {
  // Doesn't do anything
  ddt(n)=0;
  return 0;
}
