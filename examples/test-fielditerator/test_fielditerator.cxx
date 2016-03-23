/*
 * Test Cyclic Reduction parallel solver
 *
 */

#include <bout.hxx>
#include <boutmain.hxx>

#include <fielditerator.hxx>
#include <assert.h>

#include <unistd.h>

Field3D n;

#define PRINT_DEBUG {						\
    usleep(2e5*omp_get_thread_num());				\
    output.write("FieldIterator failed at ");			\
    cxit.printState();						\
    output.write("- expected %.1f but got %.1f. ",exp,d3[cxit]);	\
    output.write("I am %d !\n",omp_get_thread_num());		\
    usleep(2e5*omp_get_num_threads());				\
    abort(); }

int spread_work(int num_work, int thread, int max_thread);

int physics_init(bool restarting) {

  // for (int w=2;w<8;w++){
  //   int we=pow(2,w)+1.1;
  //   output.write("spreading %d\n",we);
  //   for (int t=1;t<5;++t){
  //     for (int ct=0;ct<=t;++ct){
  // 	output.write("%d ",spread_work(we,ct,t));
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
  for (FieldIteratorCIndex cxit( *mesh);cxit;cxit.next3()){
    if (d3[cxit] != exp) PRINT_DEBUG;
    d3[cxit]=3;
  }
  exp=3;

  for (int jx=0;jx<mesh->ngx;++jx)
    for (int jy=0;jy<mesh->ngy;++jy)
      for (int jz=0;jz<mesh->ngz;++jz)
	assert(d3(jx,jy,jz)==exp);
  
#pragma omp parallel
  for (FieldIteratorCIndex cxit( *mesh,NO_X);cxit;cxit.next3()){
    if (d3[cxit] != exp) PRINT_DEBUG;
    d3[cxit]=4;
  }

  for (int jx=0;jx<mesh->ngx;++jx){
    if (jx < mesh->xstart || jx > mesh->xend)
      exp=3;
    else
      exp=4;
    for (int jy=0;jy<mesh->ngy;++jy)
      for (int jz=0;jz<mesh->ngz;++jz)
	assert(d3(jx,jy,jz)==exp);
  }
  d3=3;
#pragma omp parallel
  for (FieldIteratorCIndex cxit( *mesh,NO_Y);cxit;cxit.next3()){
    if (d3[cxit] != exp) PRINT_DEBUG;
    d3[cxit]=4;
  }
  for (int jx=0;jx<mesh->ngx;++jx)
    for (int jy=0;jy<mesh->ngy;++jy){
      if (jy < mesh->ystart || jy > mesh->yend)
	exp=3;
      else
	exp=4;
      for (int jz=0;jz<mesh->ngz;++jz)
	assert(d3(jx,jy,jz)==exp);
  }
  d3=3;exp=3;
#pragma omp parallel
  for (FieldIteratorCIndex cxit( *mesh,NO_Z);cxit;cxit.next3()){
    if (d3[cxit] != exp) PRINT_DEBUG;
    d3[cxit]=4;
  }
  for (int jx=0;jx<mesh->ngx;++jx)
    for (int jy=0;jy<mesh->ngy;++jy)
      for (int jz=0;jz<mesh->ngz;++jz){
	if (jz == mesh->ngz-1)
	  exp=3;
	else
	  exp=4;
	assert(d3(jx,jy,jz)==exp);
      }
  d3=3;exp=3;
#pragma omp parallel
  for (FieldIteratorCIndex cxit( *mesh,NO_X|NO_Y|NO_Z);cxit;cxit.next3()){
    if (d3[cxit] != exp) PRINT_DEBUG;
    d3[cxit]=4;
  }
  for (int jx=0;jx<mesh->ngx;++jx)
    for (int jy=0;jy<mesh->ngy;++jy)
      for (int jz=0;jz<mesh->ngz;++jz){
	if (jz == mesh->ngz-1 ||
	    jy < mesh->ystart || jy > mesh->yend ||
	    jx < mesh->xstart || jx > mesh->xend) 
	  exp=3;
	else
	  exp=4;
	assert(d3(jx,jy,jz)==exp);
      }
  
  /**********************************************
   * * *                                    * * *
   * * *          Field2D test              * * *
   * * *                                    * * *
   **********************************************/
   
  Field2D d2=1;
  exp=1;
#pragma omp parallel
  for (FieldIteratorCIndex cxit(*mesh,FIELD2D);cxit;cxit.next2()){
    if (d2[cxit] != exp) PRINT_DEBUG;
    d2[cxit]=3;
  }
  exp=3;
  for (int jx=0;jx<mesh->ngx;++jx)
    for (int jy=0;jy<mesh->ngy;++jy)
      assert(d2(jx,jy)==exp);
  d2=1;exp=1;
  for (FieldIteratorCIndex cxit(*mesh,FIELD2D|NO_X);cxit;cxit.next2()){
    if (d2[cxit] != exp) PRINT_DEBUG;
    d2[cxit]=2;
  }
  for (int jx=0;jx<mesh->ngx;++jx){
    if (jx < mesh->xstart || jx > mesh->xend)
      exp=1;
    else
      exp=2;
    for (int jy=0;jy<mesh->ngy;++jy)
      assert(d2(jx,jy)==exp);
  }
  d2=1;exp=1;
  for (FieldIteratorCIndex cxit(*mesh,FIELD2D|NO_Y);cxit;cxit.next2()){
    if (d2[cxit] != exp) PRINT_DEBUG;
    d2[cxit]=2;
  }
  for (int jx=0;jx<mesh->ngx;++jx)
    for (int jy=0;jy<mesh->ngy;++jy){
      if (jy < mesh->ystart || jy > mesh->yend)
	exp=1;
      else
	exp=2;
      assert(d2(jx,jy)==exp);
  }
  d2=1;exp=1;
  for (FieldIteratorCIndex cxit(*mesh,FIELD2D|NO_X|NO_Y);cxit;cxit.next2()){
    if (d2[cxit] != exp) PRINT_DEBUG;
    d2[cxit]=2;
  }
  for (int jx=0;jx<mesh->ngx;++jx)
    for (int jy=0;jy<mesh->ngy;++jy){
      if (jx < mesh->xstart || jx > mesh->xend ||
	  jy < mesh->ystart || jy > mesh->yend)
	exp=1;
      else
	exp=2;
      assert(d2(jx,jy)==exp);
  }

  
  n=1;
  SOLVE_FOR(n);
  return 0;
}

int physics_run(BoutReal t) {
  // Doesn't do anything
  ddt(n)=0;
  return 0;
}
