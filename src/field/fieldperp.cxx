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

#include <boutcomm.hxx>
#include <globals.hxx>

#include <stdlib.h>
#include <math.h>

#include <fieldperp.hxx>
#include <utils.hxx>
#include <boutexception.hxx>
#include <msg_stack.hxx>
#include <bout/scorepwrapper.hxx>

FieldPerp::FieldPerp(Mesh *localmesh) : Field(localmesh), yindex(-1) {
  SCOREP0();
  if (fieldmesh) {
    nx = fieldmesh->LocalNx;
    nz = fieldmesh->LocalNz;
  }
  
#if CHECK > 0
  else {
    nx=-1;
    nz=-1;
  }
#endif
}

FieldPerp::FieldPerp(BoutReal val, Mesh *localmesh) : Field(localmesh), yindex(-1) {
  nx = fieldmesh->LocalNx;
  nz = fieldmesh->LocalNz;
  *this = val;
}

void FieldPerp::allocate() {
  if (data.empty()) {
    if (!fieldmesh) {
      /// If no mesh, use the global
      fieldmesh = mesh;
      nx = fieldmesh->LocalNx;
      nz = fieldmesh->LocalNz;
    }
    data = Array<BoutReal>(nx * nz);
#if CHECK > 2
    invalidateGuards(*this);
#endif
  } else
    data.ensureUnique();
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
  TRACE("FieldPerp = BoutReal");

  allocate();

#if CHECK > 0
  if (!finite(rhs))
    throw BoutException("FieldPerp: Assignment from non-finite BoutReal\n");
#endif

  for (const auto &i : (*this))
    (*this)[i] = rhs;

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

const IndexRange FieldPerp::region(REGION rgn) const {
  switch (rgn) {
  case RGN_ALL:
  case RGN_NOZ:
    return IndexRange{0, nx - 1, 0, 0, 0, nz - 1};
    break;
  case RGN_NOX:
    return IndexRange{getMesh()->xstart, getMesh()->xend, 0, 0, 0, nz - 1};
    break;
  default:
    throw BoutException("FieldPerp::region() : Requested region not implemented");
    break;
  };
}

//////////////// NON-MEMBER FUNCTIONS //////////////////

////////////// NON-MEMBER OVERLOADED OPERATORS //////////////

// Unary minus
FieldPerp operator-(const FieldPerp &f) { return -1.0 * f; }

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

/////////////////////////////////////////////////
// functions

/*!
 * This macro takes a function \p func, which is
 * assumed to operate on a single BoutReal and return
 * a single BoutReal, and wraps it up into a function
 * of a FieldPerp called \p name.
 *
 * @param name  The name of the function to define
 * @param func  The function to apply to each value
 *
 * If CHECK >= 1, checks if the FieldPerp is allocated
 *
 * Loops over the entire domain, applies function,
 * and uses checkData() to, if CHECK >= 3, check
 * result for non-finite numbers
 *
 */
#define FPERP_FUNC(name, func)                                                           \
  const FieldPerp name(const FieldPerp &f, REGION rgn) {                                 \
    TRACE(#name "(FieldPerp)");                                                          \
    /* Check if the input is allocated */                                                \
    ASSERT1(f.isAllocated());                                                            \
    /* Define and allocate the output result */                                          \
    FieldPerp result(f.getMesh());                                                       \
    result.allocate();                                                                   \
    /* Loop over domain */                                                               \
    for (const auto &d : result.region(rgn)) {                                           \
      result[d] = func(f[d]);                                                            \
    }                                                                                    \
    checkData(result);                                                                   \
    return result;                                                                       \
  }

FPERP_FUNC(abs, ::fabs);

FPERP_FUNC(sqrt, ::sqrt);

FPERP_FUNC(exp, ::exp);
FPERP_FUNC(log, ::log);

FPERP_FUNC(sin, ::sin);
FPERP_FUNC(cos, ::cos);
FPERP_FUNC(tan, ::tan);

FPERP_FUNC(sinh, ::sinh);
FPERP_FUNC(cosh, ::cosh);
FPERP_FUNC(tanh, ::tanh);

const FieldPerp copy(const FieldPerp &f) {
  SCOREP0();
  FieldPerp fcopy = f;
  fcopy.allocate();
  return fcopy;
}

const FieldPerp floor(const FieldPerp &var, BoutReal f, REGION rgn) {
  FieldPerp result = copy(var);

  for (const auto &d : result.region(rgn))
    if (result[d] < f)
      result[d] = f;

  return result;
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

BoutReal min(const FieldPerp &f, bool allpe, REGION rgn) {
  TRACE("FieldPerp::Min() %s", allpe ? "over all PEs" : "");

  ASSERT2(f.isAllocated());

  BoutReal result = f[f.region(rgn).begin()];

  for (const auto &i : f.region(rgn))
    if (f[i] < result)
      result = f[i];

  if (allpe) {
    // MPI reduce
    BoutReal localresult = result;
    MPI_Allreduce(&localresult, &result, 1, MPI_DOUBLE, MPI_MIN, BoutComm::get());
  }

  return result;
}

BoutReal max(const FieldPerp &f, bool allpe, REGION rgn) {
  TRACE("FieldPerp::Max() %s", allpe ? "over all PEs" : "");

  ASSERT2(f.isAllocated());

  BoutReal result = f[f.region(rgn).begin()];

  for (const auto &i : f.region(rgn))
    if (f[i] > result)
      result = f[i];

  if (allpe) {
    // MPI reduce
    BoutReal localresult = result;
    MPI_Allreduce(&localresult, &result, 1, MPI_DOUBLE, MPI_MAX, BoutComm::get());
  }

  return result;
}

bool finite(const FieldPerp &f, REGION rgn) {
  TRACE("finite(FieldPerp)");

  if (!f.isAllocated()) {
    return false;
  }

  for (const auto &i : f.region(rgn)) {
    if (!::finite(f[i])) {
      return false;
    }
  }

  return true;
}

FieldPerp pow(const FieldPerp &lhs, const FieldPerp &rhs, REGION rgn) {
  TRACE("pow(FieldPerp, FieldPerp)");
  // Check if the inputs are allocated
  ASSERT1(lhs.isAllocated());
  ASSERT1(rhs.isAllocated());

  // Define and allocate the output result
  ASSERT1(lhs.getMesh() == rhs.getMesh());
  FieldPerp result(lhs.getMesh());
  result.allocate();

  // Loop over domain
  for (const auto &i : result.region(rgn)) {
    result[i] = ::pow(lhs[i], rhs[i]);
  }

  checkData(result);
  return result;
}

FieldPerp pow(const FieldPerp &lhs, BoutReal rhs, REGION rgn) {
  TRACE("pow(FieldPerp, BoutReal)");
  // Check if the inputs are allocated
  ASSERT1(lhs.isAllocated());

  // Define and allocate the output result
  FieldPerp result(lhs.getMesh());
  result.allocate();

  // Loop over domain
  for (const auto &i : result.region(rgn)) {
    result[i] = ::pow(lhs[i], rhs);
  }

  checkData(result);
  return result;
}

FieldPerp pow(BoutReal lhs, const FieldPerp &rhs, REGION rgn) {
  TRACE("pow(lhs, FieldPerp)");
  // Check if the inputs are allocated
  ASSERT1(rhs.isAllocated());

  // Define and allocate the output result
  FieldPerp result(rhs.getMesh());
  result.allocate();

  // Loop over domain
  for (const auto &i : result.region(rgn)) {
    result[i] = ::pow(lhs, rhs[i]);
  }

  checkData(result);
  return result;
}

#if CHECK > 0
/// Check if the data is valid
void checkData(const FieldPerp &f, REGION region) {
  if (!f.isAllocated()) {
    throw BoutException("FieldPerp: Operation on empty data\n");
  }

#if CHECK > 2
  // Do full checks
  for (const auto &i : f.region(region)) {
    if (!::finite(f[i])) {
      throw BoutException("FieldPerp: Operation on non-finite data at [%d][%d]\n", i.x,
                          i.z);
    }
  }
#endif
}
#endif

void invalidateGuards(FieldPerp &var) {
#if CHECK > 2
  Mesh *localmesh = var.getMesh();

  // Inner x -- all z
  for (int ix = 0; ix < localmesh->xstart; ix++) {
    for (int iz = 0; iz < localmesh->LocalNz; iz++) {
      var(ix, iz) = std::nan("");
    }
  }

  // Outer x -- all z
  for (int ix = localmesh->xend + 1; ix < localmesh->LocalNx; ix++) {
    for (int iz = 0; iz < localmesh->LocalNz; iz++) {
      var(ix, iz) = std::nan("");
    }
  }
#endif
  return;
}
