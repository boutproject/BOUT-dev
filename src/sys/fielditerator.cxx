// The iterator are not iterators in the C++ sense, as they are not
// lists. However, they behave similar, and should allow to use
// OpenMP.

#include "fielditerator.hxx"

// int calcStart(int);
// int calcEnd(int);

// FieldIteratorIndex::FieldIteratorIndex(int flags,Mesh & mesh): mesh(mesh){
//   int offset=0;
//   int nz=mesh.ngz;
//   int ny=mesh.ngy;
//   if (flags&NO_Y){
//     //ny-=2*mesh->ystart;
//     //end    = - mesh->ystart*nz;
//     offset =   mesh.ystart*nz;
//   }
//   int nx=mesh.ngx;
//   if (flags&NO_X){
//     //nx =mesh->xstart;
//     offset += mesh.xstart*ny*nz;
//   }
//   int total=nx*ny*nz-2*offset;
//   if (flags&PARALLEL){
//     index = calcStart(total);
//     end   = calcEnd(total);
//   } else {
//     index = 0;
//     end   = total;
//   }
//   index +=offset;
//   start  =index;
//   end   -=offset;
// }

// void FieldIteratorIndex::reset(){
//   index=start;
// }

// FieldIteratorCIndex::FieldIteratorCIndex(int flags,Mesh & mesh): flags(flags),mesh(mesh){
//   xend=mesh.ngx;
//   xstart=0;
//   if (flags&NO_X){
//     xend-=mesh.xstart;
//     xstart=mesh.xstart;
//   }
//   yend=mesh.ngy;
//   ystart=0;
//   if (flags&NO_Y){
//     yend-=mesh.ystart;
//     ystart=mesh.ystart;
//   }
//   zend=mesh.ngz;
//   if (flags&NO_Z){
//     zend--;
//   }
  
  
// }


// int calcEnd(int _max){
// #ifndef _OPENMP
//     return _max;
// #else
//     int ipp    = _max/(omp_get_num_threads());
//     int my_max = ipp*(omp_get_thread_num()+1);
//     int rest   = _max%omp_get_num_threads();
//     if (omp_get_thread_num() < rest){
//       my_max += omp_get_thread_num()+1;
//     } else{
//       my_max += rest;
//     }
//     return my_max;
// #endif
// }

// int calcStart(int _max){
// #ifndef _OPENMP
//     return 0;
// #else
//     int ipp = _max/omp_get_num_threads();
//     int a=ipp*omp_get_thread_num();
//     int rest=_max%omp_get_num_threads();
//     if (omp_get_thread_num() < rest){
//       a      += omp_get_thread_num();
//     } else{
//       a+=rest;
//     }
//     return a;
// #endif
// }


FieldIteratorCIndex::FieldIteratorCIndex(int flags, Mesh & mesh):flags(flags),mesh(mesh){
  init();
}

FieldIteratorCIndex::FieldIteratorCIndex(Mesh & mesh, int flags):flags(flags),mesh(mesh){
  init();
}

static int spread_work(int num_work, int thread, int max_thread);

void FieldIteratorCIndex::init(){
  if (flags & NO_X){
    start.jx=mesh.xstart;
    end.jx = mesh.xend+1;
  } else {
    start.jx=0;
    end.jx  = mesh.ngx;
  }
  if (flags & NO_Y){
    start.jy=mesh.ystart;
    end.jy  =mesh.yend+1;
  } else {
    start.jy=0;
    end.jy  = mesh.ngy;
  }
  if (flags & NO_Z){
    start.jz=0;
    end.jz  =mesh.ngz-1;
  } else {
    start.jz=0;
    end.jz  =mesh.ngz;
  }
  ymin=start.jy;
#ifdef _OPENMP
  if (flags & PARALLEL){
    CIndex work = end-start;
    int iwork   = work.jx*work.jy*work.jz;
    int np=omp_get_num_threads();
    int cp=omp_get_thread_num();
    int istart=spread_work(iwork,cp,np);
    int iend=spread_work(iwork,cp+1,np);
    end.jz=start.jz+iend%work.jz;
    iend/=work.jz;
    end.jy=start.jy+iend%work.jy;
    iend/=work.jy;
    end.jx=start.jx+iend;
    start.jz+=istart%work.jz;
    istart/=work.jz;
    start.jy+=istart%work.jy;
    istart/=work.jy;
    start.jx+=istart;
  }
#endif
  start.flags=flags;
  current=start;
}

void FieldIteratorCIndex::reset(){
  current=start;
}


int spread_work(int work,int cp,int np){
  int pp=work/np;
  int rest=work%np;
  int result=pp*cp;
  if (rest > cp){
    result +=cp;
  } else {
    result +=rest;
  }
  return result;
};
  
