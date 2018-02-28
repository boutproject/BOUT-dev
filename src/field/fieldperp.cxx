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
#include <msg_stack.hxx>
#include <bout/scorepwrapper.hxx>

FieldPerp::FieldPerp(Mesh *localmesh) {
  SCOREP0();
  // Get mesh size
  fieldmesh = localmesh;
  if (localmesh) {
    nx = localmesh->LocalNx;
    nz = localmesh->LocalNz;
  }
  
#if CHECK > 0
  else {
    nx=-1;
    nz=-1;
  }
#endif

  yindex = -1;
}

/***************************************************************
 *                         ASSIGNMENT 
 ***************************************************************/

FieldPerp & FieldPerp::operator=(const FieldPerp &rhs) {
  SCOREP0();
  nx = rhs.nx;
  nz = rhs.nz;
  yindex = rhs.yindex;
  data = rhs.data;
  return *this;
}

FieldPerp & FieldPerp::operator=(const BoutReal rhs) {
  SCOREP0();
  allocate();

  for(auto&& d : data) {
    d = rhs;
  }
  
  return *this;
}

/***************************************************************
 *                         ITERATORS
 ***************************************************************/

const DataIterator FieldPerp::begin() const {
  SCOREP0();
  return DataIterator( 0, nx-1,
                      yindex, yindex,
                      0, nz-1);
}

const DataIterator FieldPerp::end() const {
  SCOREP0();
  return DataIterator( 0, nx-1,
                      yindex, yindex,
		       0, nz-1,DI_GET_END);
}


/***************************************************************
 *                         OPERATORS 
 ***************************************************************/

#define FPERP_OP_FIELD(op, bop, ftype)			\
  FieldPerp& FieldPerp::operator op(const ftype &rhs) { \
    SCOREP0();                                          \
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
  FieldPerp& FieldPerp::operator op(BoutReal rhs) { \
    SCOREP0();                                          \
    if(data.unique()) {                                 \
      /* Only reference to the data */           	\
      for(int i=0;i<nx;i++)                             \
        for(int k=0;k<nz;k++)                           \
          (*this)(i,k) op rhs;				\
    }else {  			                        \
      /* Shared with another FieldPerp */		\
      (*this) = (*this) bop rhs;                        \
    }                                                   \
    return *this;                                       \
  }

FPERP_OP_REAL(+=, +);
FPERP_OP_REAL(-=, -);
FPERP_OP_REAL(*=, *);
FPERP_OP_REAL(/=, /);

////////////// NON-MEMBER OVERLOADED OPERATORS //////////////

// Operator on FieldPerp and another field
#define FPERP_FPERP_OP_FIELD(op, ftype)                                                  \
  const FieldPerp operator op(const FieldPerp &lhs, const ftype &rhs) {                  \
    SCOREP0();                                                                           \
    FieldPerp result(lhs.getMesh());                                                     \
    result.allocate();                                                                   \
                                                                                         \
    int y = lhs.getIndex();                                                              \
    result.setIndex(y);                                                                  \
                                                                                         \
    for (auto i : result)                                                                \
      result[i] = lhs[i] op rhs[i];                                                      \
                                                                                         \
    return result;                                                                       \
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
#define FPERP_FPERP_OP_REAL(op)                                                          \
  const FieldPerp operator op(const FieldPerp &lhs, BoutReal rhs) {                      \
    SCOREP0();                                                                           \
    FieldPerp result(lhs.getMesh());                                                     \
    result.allocate();                                                                   \
                                                                                         \
    int y = lhs.getIndex();                                                              \
    result.setIndex(y);                                                                  \
                                                                                         \
    for (auto i : result)                                                                \
      result[i] = lhs[i] op rhs;                                                         \
                                                                                         \
    return result;                                                                       \
  }

FPERP_FPERP_OP_REAL(+);
FPERP_FPERP_OP_REAL(-);
FPERP_FPERP_OP_REAL(*);
FPERP_FPERP_OP_REAL(/);

#define FPERP_REAL_OP_FPERP(op)                                                          \
  const FieldPerp operator op(BoutReal lhs, const FieldPerp &rhs) {                      \
    SCOREP0();                                                                           \
    FieldPerp result(rhs.getMesh());                                                     \
    result.allocate();                                                                   \
                                                                                         \
    int y = rhs.getIndex();                                                              \
    result.setIndex(y);                                                                  \
                                                                                         \
    for (auto i : result)                                                                \
      result[i] = lhs op rhs[i];                                                         \
                                                                                         \
    return result;                                                                       \
  }

// Only need the asymmetric operators
FPERP_REAL_OP_FPERP(-);
FPERP_REAL_OP_FPERP(/);

const FieldPerp copy(const FieldPerp &f) {
  SCOREP0();
  FieldPerp fcopy = f;
  fcopy.allocate();
  return fcopy;
}

const FieldPerp sliceXZ(const Field3D& f, int y) {
  SCOREP0();
  // Source field should be valid
  ASSERT1(f.isAllocated());

  FieldPerp result(f.getMesh());

  // Allocate memory
  result.allocate();
  result.setIndex(y);

  for(auto i : result)
    result[i] = f[i];
  
  return result;
}

