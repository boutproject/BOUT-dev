/**************************************************************************
 * Class for 2D X-Z slices
 *
 **************************************************************************
 * Copyright 2010 B.D.Dudson, S.Farley, M.V.Umansky, X.Q.Xu
 *
 * Contact: Ben Dudson, bd512@york.ac.uk
 * 
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

#include <globals.hxx>

#include <stdlib.h>
#include <math.h>

#include <fieldperp.hxx>
#include <utils.hxx>
#include <boutexception.hxx>

FieldPerp::FieldPerp() {
  // Get mesh size
  nx = mesh->ngx;
  nz = mesh->ngz;
  
  yindex = -1;
}


FieldPerp* FieldPerp::clone() const {
  return new FieldPerp(*this);
}

/***************************************************************
 *                         ASSIGNMENT 
 ***************************************************************/

FieldPerp & FieldPerp::operator=(const FieldPerp &rhs) {
  nx = rhs.nx;
  nz = rhs.nz;
  yindex = rhs.yindex;
  data = rhs.data;
}

FieldPerp & FieldPerp::operator=(const BoutReal rhs) {
  allocate();

  for(BoutReal *it = data.begin(); it != data.end(); it++)
    *it = rhs;
  
  return *this;
}

/***************************************************************
 *                         OPERATORS 
 ***************************************************************/

#define FPERP_OP_FIELD(op, bop, ftype)			\
  FieldPerp& FieldPerp::operator op(const ftype &rhs) { \
    if(data.unique()) {                                 \
      /* Only reference to the data */			\
      for(int i=0;i<nx;i++)                             \
        for(int k=0;k<nz;k++)                           \
          (*this)(i,k) op rhs(i, yindex, k);            \
    }else {  			                        \
      /* Shared with another FieldPerp */		\
      (*this) = (*this) bop rhs;                        \
    }                                                   \
    return *this;                                       \
  }

FPERP_OP_FIELD(+=, +, FieldPerp);
FPERP_OP_FIELD(+=, +, Field3D);
FPERP_OP_FIELD(+=, +, Field2D);

FPERP_OP_FIELD(-=, -, FieldPerp);
FPERP_OP_FIELD(-=, -, Field3D);
FPERP_OP_FIELD(-=, -, Field2D);

FPERP_OP_FIELD(*=, *, FieldPerp);
FPERP_OP_FIELD(*=, *, Field3D);
FPERP_OP_FIELD(*=, *, Field2D);

FPERP_OP_FIELD(/=, /, FieldPerp);
FPERP_OP_FIELD(/=, /, Field3D);
FPERP_OP_FIELD(/=, /, Field2D);

#define FPERP_OP_REAL(op, bop)  			\
  FieldPerp& FieldPerp::operator op(const BoutReal &rhs) { \
    if(data.unique()) {                                 \
      /* Only reference to the data */           	\
      for(int i=0;i<nx;i++)                             \
        for(int k=0;k<nz;k++)                           \
          (*this)(i,k) op rhs;				\
    }else {  			                        \
      /* Shared with another FieldPerp */		\
      (*this) = (*this) bop rhs;                        \
    }                                                   \
  }

FPERP_OP_REAL(+=, +);
FPERP_OP_REAL(-=, -);
FPERP_OP_REAL(*=, *);
FPERP_OP_REAL(/=, /);

////////////////////// STENCILS //////////////////////////

/*
void FieldPerp::setStencil(bstencil *fval, bindex *bx) const {
  fval->cc = (*this)(bx->jx,bx->jz);

  fval->xp = (*this)(bx->jxp,bx->jz);
  fval->xm = (*this)(bx->jxm,bx->jz);
  fval->x2p = (*this)(bx->jx2p,bx->jz);
  fval->x2m = (*this)(bx->jx2m,bx->jz);
  
  fval->yp = (*this)(bx->jx,bx->jz);
  fval->ym = (*this)(bx->jx,bx->jz);
  fval->zp = (*this)(bx->jx,bx->jzp);
  fval->zm = (*this)(bx->jx,bx->jzm);

  fval->y2p = (*this)(bx->jx,bx->jy);
  fval->y2m = (*this)(bx->jx,bx->jy);
  fval->z2p = (*this)(bx->jx,bx->jz2p);
  fval->z2m = (*this)(bx->jx,bx->jz2m);
}
*/

void FieldPerp::setXStencil(stencil &fval, const bindex &bx, CELL_LOC loc) const {
  fval.p = (*this)(bx.jxp,bx.jz);
  fval.m = (*this)(bx.jxm,bx.jz);
  fval.pp = (*this)(bx.jx2p,bx.jz);
  fval.mm = (*this)(bx.jx2m,bx.jz);
}

void FieldPerp::setYStencil(stencil &fval, const bindex &bx, CELL_LOC loc) const {
  fval = (*this)(bx.jx,bx.jz);
}

void FieldPerp::setZStencil(stencil &fval, const bindex &bx, CELL_LOC loc) const {
  fval.p = (*this)(bx.jx,bx.jzp);
  fval.m = (*this)(bx.jx,bx.jzm);
  fval.pp = (*this)(bx.jx,bx.jz2p);
  fval.mm = (*this)(bx.jx,bx.jz2m);
}

////////////// NON-MEMBER OVERLOADED OPERATORS //////////////

// Operator on FieldPerp and another field
#define FPERP_FPERP_OP_FIELD(op, ftype)                     	          \
  const FieldPerp operator op(const FieldPerp &lhs, const ftype &rhs) {   \
    FieldPerp result;                                                     \
    result.allocate();                                                    \
                                                                          \
    int y = lhs.getIndex();            					\
    result.setIndex(y);                                                   \
                                                                          \
    for(int i=0; i<mesh->ngx; i++)                                        \
      for(int j=0; j<mesh->ngz; j++)                                      \
        result(i,j) = lhs(i,j) op rhs(i,y,j);                             \
                                                                          \
    return result;                                                        \
  }

FPERP_FPERP_OP_FIELD(+, FieldPerp);
FPERP_FPERP_OP_FIELD(+, Field3D);
FPERP_FPERP_OP_FIELD(+, Field2D);

FPERP_FPERP_OP_FIELD(-, FieldPerp);
FPERP_FPERP_OP_FIELD(-, Field3D);
FPERP_FPERP_OP_FIELD(-, Field2D);

FPERP_FPERP_OP_FIELD(*, FieldPerp);
FPERP_FPERP_OP_FIELD(*, Field3D);
FPERP_FPERP_OP_FIELD(*, Field2D);

FPERP_FPERP_OP_FIELD(/, FieldPerp);
FPERP_FPERP_OP_FIELD(/, Field3D);
FPERP_FPERP_OP_FIELD(/, Field2D);

// Operator on FieldPerp and BoutReal
#define FPERP_FPERP_OP_REAL(op)                     	                   \
  const FieldPerp operator op(const FieldPerp &lhs, const BoutReal &rhs) { \
    FieldPerp result;                                                     \
    result.allocate();                                                    \
                                                                          \
    int y = lhs.getIndex();						\
    result.setIndex(y);                                                   \
                                                                          \
    for(int i=0; i<mesh->ngx; i++)                                        \
      for(int j=0; j<mesh->ngz; j++)                                      \
        result(i,j) = lhs(i,j) op rhs;                                    \
                                                                          \
    return result;                                                        \
  }

FPERP_FPERP_OP_REAL(+);
FPERP_FPERP_OP_REAL(-);
FPERP_FPERP_OP_REAL(*);
FPERP_FPERP_OP_REAL(/);

#define FPERP_REAL_OP_FPERP(op)                     	                   \
  const FieldPerp operator op(const BoutReal &lhs, const FieldPerp &rhs) { \
    FieldPerp result;                                                     \
    result.allocate();                                                    \
                                                                          \
    int y = rhs.getIndex();						\
    result.setIndex(y);                                                   \
                                                                          \
    for(int i=0; i<mesh->ngx; i++)                                        \
      for(int j=0; j<mesh->ngz; j++)                                      \
        result(i,j) = lhs op rhs(i,j);                                    \
                                                                          \
    return result;                                                        \
  }

// Only need the asymmetric operators
FPERP_REAL_OP_FPERP(-);
FPERP_REAL_OP_FPERP(/);

const FieldPerp copy(const FieldPerp &f) {
  FieldPerp fcopy = f;
  fcopy.allocate();
  return fcopy;
}

const FieldPerp SQ(const FieldPerp &f) {
  return f*f;
}
 
const FieldPerp sliceXZ(const Field3D& f, int y) {
  // Source field should be valid
  ASSERT1(f.isAllocated());
  
  FieldPerp result;

  // Allocate memory
  result.allocate();

  result.setIndex(y);

  for(int jx=0;jx<mesh->ngx;jx++)
    for(int jz=0;jz<mesh->ngz;jz++)
      result(jx,jz) = f(jx,y,jz);
  
  return result;
}

