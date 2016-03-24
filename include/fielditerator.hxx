/**************************************************************************
 * Several FieldIterators.
 * 
 * FieldIteratorIndex 
 * does include some boundaries, even if it is told to not include
 * them. Therefore it should only be used, if the calculation itself
 * is fast, or anyway all boundaries are included. It is mainly
 * intended for the multiplication etc. in Field3D 
 *
 **************************************************************************
 * This file is part of BOUT++.
 *
 * BOUT++ is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * BOUT++ is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with BOUT++.  If not, see <http://www.gnu.org/licenses/>.
 * 
 **************************************************************************/
#pragma once

#include "mpi.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#include "stencils.hxx"
#include "bout/mesh.hxx"
#include "bout/assert.hxx"


enum flags_fielditerator_t{
  NONE=0,
  NO_X=1,
  NO_Y=2,
  NO_Z=4,
  PARALLEL=0,
  NOT_PARALLEL=8,
  CALC_INDEX=16,
  NO_BNDRY=NO_X|NO_Y|NO_Z,
  FIELD2D=32
};


// class FieldIteratorIndex{
// public:
//   //FieldIteratorIndex & operator++(){++index;};
//   FieldIteratorIndex & operator++(){next();};
//   void next(){++index;};
//   //next(){++index;};
//   operator int() {return index;};
//   FieldIteratorIndex(int flags=0);
//   bool finished() {return !notFinished();};
//   bool notFinished() {return index < end;};
//   //FieldIteratorIndex(mesh);
//   FieldIteratorIndex(int flags,Mesh & mesh);
//   void reset();
// private:
//   // please use the prefix iterator
//   FieldIteratorIndex & operator++(int);
//   // no need to copy
//   FieldIteratorIndex ( FieldIteratorIndex &);
//   // actuall private stuff
//   int start;
//   int end;
//   int index;
//   Mesh & mesh;
// };


class CIndex{
public:
  int jx,jy,jz;
  int flags;
  //Mesh & mesh;
  //CIndex(int _jx,int _jy, int _jz,int _flags, Mesh &
  //_mesh):mesh(_mesh){
  CIndex(int _jx,int _jy, int _jz,int _flags){
    jx=_jx;
    jy=_jy;
    jz=_jz;
    flags=_flags;
  };
  CIndex(int _jx,int _jy, int _jz){
    jx=_jx;
    jy=_jy;
    jz=_jz;
  };
  CIndex(){};// no init
  operator bindex(){
    bindex bx;
    bx.jx=jx;
    bx.jy=jy;
    bx.jz=jz;
    if (flags&CALC_INDEX){
      calc_index(&bx);
    }
    return bx;
  };
  CIndex& operator-=(const CIndex& rhs){
    jx-=rhs.jx;
    jy-=rhs.jy;
    jz-=rhs.jz;
  }
  CIndex& operator+=(const CIndex& rhs){
    jx+=rhs.jx;
    jy+=rhs.jy;
    jz+=rhs.jz;
  }
  /*int getIndex(){
    int ny=mesh.ngy;
    int nz=mesh.ngz;
    return jz +(jx*ny+jy)*nz;
  };
  void set(int index){
    ASSERT1(index > 0);
    int ny=mesh.ngy;
    int nz=mesh.ngz;
    jz=index % nz;
    index /= nz;
    jy=index % ny;
    index /= ny;
    jx=index;
    ASSERT1(jx < mesh.ngx);
    };*/
};

inline bool operator< (const CIndex& lhs,const CIndex& rhs){
  if (lhs.jx < rhs.jx){
    return true;
  }
  if (lhs.jx > rhs.jx){
    return false;
  }
  if (lhs.jy < rhs.jy){
    return true;
  }
  if (lhs.jy > rhs.jy){
    return false;
  }
  return (lhs.jz < rhs.jz);
};
// bool operator > (const CIndex &rhs){
//   if (jx > rhs.jx){
//     return true;
//   }
//   if (jx < rhs.jx){
//     return false;
//   }
//   if (jy > rhs.jy){
//     return true;
//   }
//   if (jy < rhs.jy){
//     return false;
//   }
//   return (jz > rhs.jz);
// };
inline bool operator==(const CIndex&lhs, const CIndex& rhs){
  return (lhs.jx==rhs.jx && lhs.jy == rhs.jy && lhs.jz==rhs.jz);
};
inline bool operator!=(const CIndex& lhs, const CIndex& rhs){return !operator==(lhs,rhs);}
inline bool operator> (const CIndex& lhs, const CIndex& rhs){return  operator< (rhs,lhs);}
inline bool operator<=(const CIndex& lhs, const CIndex& rhs){return !operator> (lhs,rhs);}
inline bool operator>=(const CIndex& lhs, const CIndex& rhs){return !operator< (lhs,rhs);}
inline CIndex operator-(CIndex lhs, const CIndex& rhs){
  lhs-=rhs;
  return lhs;
};
inline CIndex operator+(CIndex lhs, const CIndex& rhs){
  lhs+=rhs;
  return lhs;
};

class FieldIteratorCIndex{
public:
  // needs to be first
  // otherwise it cannot be casted as CIndex
  CIndex current;
  operator bindex(){return current;};
  operator CIndex&() {return current;};
  void next3(){
    ASSERT2(!(flags&FIELD2D));
    if (++current.jz == max.jz){
      current.jz=0;
      if (++current.jy == max.jy){
	current.jy=ymin;
	if (++current.jx == max.jx){
#ifndef _OPENMP
	  notfinished=false;
#endif
	};
      }
    }
#ifdef _OPENMP
    notfinished= current <= end;
#endif
    //return notfinished;
  };
  void next2(){
    ASSERT2(flags&FIELD2D);
    if (++current.jy == max.jy){
      current.jy=ymin;
#ifdef _OPENMP
      ++current.jx;
#else
      if (++current.jx == max.jx){
	notfinished=false;
      }
#endif
    }
#ifdef _OPENMP
    notfinished = current < end;
#endif
    //return notfinished;
  };
  FieldIteratorCIndex & operator++(){
    this->next3();
  };
  FieldIteratorCIndex(int flags, Mesh & mesh);
  FieldIteratorCIndex(Mesh & mesh,int flags=0);
  bool finished() {return !notfinished;};
  bool notFinished() {return notfinished ;};
  operator bool() {return notfinished;};
  //FieldIteratorCIndex(int flags,Mesh & mesh);
  void reset();
  void printState();
private:
  // please use the prefix iterator
  FieldIteratorCIndex & operator++(int);
  // no need to copy
  FieldIteratorCIndex ( FieldIteratorCIndex &);
  // actuall private stuff
  //int zmax,ymax,xmax;
  CIndex max;
  CIndex start;
  CIndex end;
  int ymin;//zmin=0 -> not needed
  int flags;
  Mesh & mesh;
  bool notfinished;
  //int xend,yend,zend;
  //int xstart,ystart,zstart;
  void init();
};


