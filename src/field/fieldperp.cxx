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

#include <bout/boutcomm.hxx>
#include <bout/boutexception.hxx>
#include <bout/fieldperp.hxx>
#include <bout/globals.hxx>
#include <bout/msg_stack.hxx>
#include <bout/utils.hxx>

#include <stdlib.h>
#include <math.h>

FieldPerp::FieldPerp(Mesh *localmesh) : Field(localmesh) {
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

FieldPerp::FieldPerp(BoutReal val, Mesh *localmesh) : Field(localmesh) {
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

FieldPerp &FieldPerp::operator=(const FieldPerp &rhs) {
  /// Check for self-assignment
  if (this == &rhs) {
    return (*this); // skip this assignment
  }

  checkData(rhs);

  nx = rhs.nx;
  nz = rhs.nz;
  yindex = rhs.yindex;
  data = rhs.data;
  return *this;
}

FieldPerp & FieldPerp::operator=(const BoutReal rhs) {
  TRACE("FieldPerp = BoutReal");

  allocate();

  checkData(rhs);

  const Region<IndPerp> &region_all = fieldmesh->getRegionPerp("RGN_ALL");

  BOUT_FOR(i, region_all) {
    (*this)[i] = rhs;
  }

  return *this;
}

/***************************************************************
 *                         OPERATORS 
 ***************************************************************/

#define FPERP_OP_FPERP(op, bop, ftype)                                                   \
  FieldPerp &FieldPerp::operator op(const ftype &rhs) {                                  \
    if (data.unique()) {                                                                 \
      checkData(rhs);                                                                    \
      /* Only reference to the data */                                                   \
      const Region<IndPerp> &region = fieldmesh->getRegionPerp("RGN_ALL");               \
      BOUT_FOR(i, region) {                                                              \
        (*this)[i] op rhs[i];                                                            \
      }                                                                                  \
      checkData(*this);                                                                  \
    } else {                                                                             \
      /* Shared with another FieldPerp */                                                \
      (*this) = (*this)bop rhs;                                                          \
    }                                                                                    \
    return *this;                                                                        \
  }

FPERP_OP_FPERP(+=, +, FieldPerp);
FPERP_OP_FPERP(-=, -, FieldPerp);
FPERP_OP_FPERP(*=, *, FieldPerp);
FPERP_OP_FPERP(/=, /, FieldPerp);

#define FPERP_OP_FIELD(op, bop, ftype)                                                   \
  FieldPerp &FieldPerp::operator op(const ftype &rhs) {                                  \
    if (data.unique()) {                                                                 \
      checkData(*this);                                                                  \
      checkData(rhs);                                                                    \
      /* Only reference to the data */                                                   \
      const Region<IndPerp> &region = fieldmesh->getRegionPerp("RGN_ALL");               \
      BOUT_FOR(i, region) {                                                              \
        (*this)[i] op rhs(i.x(), yindex, i.z());                                         \
      }                                                                                  \
      checkData(*this);                                                                  \
    } else {                                                                             \
      /* Shared with another FieldPerp */                                                \
      (*this) = (*this)bop rhs;                                                          \
    }                                                                                    \
    return *this;                                                                        \
  }

FPERP_OP_FIELD(+=, +, Field3D);
FPERP_OP_FIELD(+=, +, Field2D);

FPERP_OP_FIELD(-=, -, Field3D);
FPERP_OP_FIELD(-=, -, Field2D);

FPERP_OP_FIELD(*=, *, Field3D);
FPERP_OP_FIELD(*=, *, Field2D);

FPERP_OP_FIELD(/=, /, Field3D);
FPERP_OP_FIELD(/=, /, Field2D);

#define FPERP_OP_REAL(op, bop)                                                           \
  FieldPerp &FieldPerp::operator op(BoutReal rhs) {                                      \
    if (data.unique()) {                                                                 \
      checkData(rhs);                                                                    \
      /* Only reference to the data */                                                   \
      const Region<IndPerp> &region = fieldmesh->getRegionPerp("RGN_ALL");               \
      BOUT_FOR(i, region) {                                                              \
        (*this)[i] op rhs;                                                               \
      }                                                                                  \
      checkData(*this);                                                                  \
    } else {                                                                             \
      /* Shared with another FieldPerp */                                                \
      (*this) = (*this)bop rhs;                                                          \
    }                                                                                    \
    return *this;                                                                        \
  }

FPERP_OP_REAL(+=, +);
FPERP_OP_REAL(-=, -);
FPERP_OP_REAL(*=, *);
FPERP_OP_REAL(/=, /);

const Region<IndPerp> &FieldPerp::getRegion(REGION region) const {
  return fieldmesh->getRegionPerp(REGION_STRING(region));
};
const Region<IndPerp> &FieldPerp::getRegion(const std::string &region_name) const {
  return fieldmesh->getRegionPerp(region_name);
};

//////////////// NON-MEMBER FUNCTIONS //////////////////

////////////// NON-MEMBER OVERLOADED OPERATORS //////////////

// Unary minus
FieldPerp operator-(const FieldPerp &f) { return -1.0 * f; }

// Operator on FieldPerp and another field
#define FPERP_FPERP_OP_FPERP(op, ftype)                                                  \
  const FieldPerp operator op(const FieldPerp &lhs, const ftype &rhs) {                  \
    checkData(lhs);                                                                      \
    checkData(rhs);                                                                      \
    FieldPerp result(lhs.getMesh());                                                     \
    result.allocate();                                                                   \
    result.setIndex(lhs.getIndex());                                                     \
                                                                                         \
    const Region<IndPerp> &region = lhs.getMesh()->getRegionPerp("RGN_ALL");             \
    BOUT_FOR(i, region) {                                                                \
      result[i] = lhs[i] op rhs[i];                                                      \
    }                                                                                    \
    checkData(result);                                                                   \
    return result;                                                                       \
  }

FPERP_FPERP_OP_FPERP(+, FieldPerp);
FPERP_FPERP_OP_FPERP(-, FieldPerp);
FPERP_FPERP_OP_FPERP(*, FieldPerp);
FPERP_FPERP_OP_FPERP(/, FieldPerp);

// Operator on FieldPerp and another field
#define FPERP_FPERP_OP_FIELD(op, ftype)                                                  \
  const FieldPerp operator op(const FieldPerp &lhs, const ftype &rhs) {                  \
    checkData(lhs);                                                                      \
    checkData(rhs);                                                                      \
    FieldPerp result(lhs.getMesh());                                                     \
    result.allocate();                                                                   \
    result.setIndex(lhs.getIndex());                                                     \
                                                                                         \
    const Region<IndPerp> &region = lhs.getMesh()->getRegionPerp("RGN_ALL");             \
    BOUT_FOR(i, region) {                                                                \
      result[i] = lhs[i] op rhs(i.x(), lhs.getIndex(), i.z());                           \
    }                                                                                    \
    checkData(result);                                                                   \
    return result;                                                                       \
  }

FPERP_FPERP_OP_FIELD(+, Field3D);
FPERP_FPERP_OP_FIELD(+, Field2D);

FPERP_FPERP_OP_FIELD(-, Field3D);
FPERP_FPERP_OP_FIELD(-, Field2D);

FPERP_FPERP_OP_FIELD(*, Field3D);
FPERP_FPERP_OP_FIELD(*, Field2D);

FPERP_FPERP_OP_FIELD(/, Field3D);
FPERP_FPERP_OP_FIELD(/, Field2D);

// Operator on FieldPerp and BoutReal
#define FPERP_FPERP_OP_REAL(op)                                                          \
  const FieldPerp operator op(const FieldPerp &lhs, BoutReal rhs) {                      \
    checkData(lhs);                                                                      \
    checkData(rhs);                                                                      \
    FieldPerp result(lhs.getMesh());                                                     \
    result.allocate();                                                                   \
    result.setIndex(lhs.getIndex());                                                     \
                                                                                         \
    const Region<IndPerp> &region = result.getMesh()->getRegionPerp("RGN_ALL");          \
    BOUT_FOR (i, region) {                                                               \
      result[i] = lhs[i] op rhs;                                                         \
    }                                                                                    \
                                                                                         \
    checkData(result);                                                                   \
    return result;                                                                       \
  }

FPERP_FPERP_OP_REAL(+);
FPERP_FPERP_OP_REAL(-);
FPERP_FPERP_OP_REAL(*);
FPERP_FPERP_OP_REAL(/);

#define FPERP_REAL_OP_FPERP(op)                                                          \
  const FieldPerp operator op(BoutReal lhs, const FieldPerp &rhs) {                      \
    checkData(lhs);                                                                      \
    checkData(rhs);                                                                      \
    FieldPerp result(rhs.getMesh());                                                     \
    result.allocate();                                                                   \
    result.setIndex(rhs.getIndex());                                                     \
                                                                                         \
    const Region<IndPerp> &region = result.getMesh()->getRegionPerp("RGN_ALL");          \
    BOUT_FOR (i, region) {                                                               \
      result[i] = lhs op rhs[i];                                                         \
    }                                                                                    \
                                                                                         \
    checkData(result);                                                                   \
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
    checkData(f);                                                                        \
    TRACE(#name "(FieldPerp)");                                                          \
    /* Check if the input is allocated */                                                \
    ASSERT1(f.isAllocated());                                                            \
    /* Define and allocate the output result */                                          \
    FieldPerp result(f.getMesh());                                                       \
    result.allocate();                                                                   \
    result.setIndex(f.getIndex());                                                       \
    const Region<IndPerp> &region = f.getMesh()->getRegionPerp(REGION_STRING(rgn));      \
    BOUT_FOR (d, region) {                                                               \
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
  FieldPerp fcopy = f;
  fcopy.allocate();
  return fcopy;
}

const FieldPerp floor(const FieldPerp &var, BoutReal f, REGION rgn) {
  checkData(var);
  FieldPerp result = copy(var);

  const Region<IndPerp> &region = var.getMesh()->getRegionPerp(REGION_STRING(rgn));
  BOUT_FOR(d, region) {
    if (result[d] < f) {
      result[d] = f;
    }
  }

  checkData(result);
  return result;
}

const FieldPerp sliceXZ(const Field3D& f, int y) {
  // Source field should be valid
  checkData(f);

  FieldPerp result(f.getMesh());

  // Allocate memory
  result.allocate();
  result.setIndex(y);

  const Region<IndPerp> &region_all = f.getMesh()->getRegionPerp("RGN_ALL");
  BOUT_FOR(i, region_all) {
    result[i] = f(i, y);
  }
  
  checkData(result);
  return result;
}

BoutReal min(const FieldPerp &f, bool allpe, REGION rgn) {
  TRACE("FieldPerp::Min() %s", allpe ? "over all PEs" : "");

  checkData(f);

  const Region<IndPerp> &region = f.getMesh()->getRegionPerp(REGION_STRING(rgn));

  BoutReal result = f[*region.cbegin()];

  BOUT_FOR_OMP(i, region, parallel for reduction(min:result)) {
    if (f[i] < result) {
      result = f[i];
    }
  }

  if (allpe) {
    // MPI reduce
    BoutReal localresult = result;
    MPI_Allreduce(&localresult, &result, 1, MPI_DOUBLE, MPI_MIN, BoutComm::get());
  }

  return result;
}

BoutReal max(const FieldPerp &f, bool allpe, REGION rgn) {
  TRACE("FieldPerp::Max() %s", allpe ? "over all PEs" : "");

  checkData(f);

  const Region<IndPerp> &region = f.getMesh()->getRegionPerp(REGION_STRING(rgn));

  BoutReal result = f[*region.cbegin()];

  BOUT_FOR_OMP(i, region, parallel for reduction(max:result)) {
    if (f[i] > result) {
      result = f[i];
    }
  }

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

  const Region<IndPerp> &region = f.getMesh()->getRegionPerp(REGION_STRING(rgn));

  BOUT_FOR_SERIAL(i, region) {
    if (!::finite(f[i])) {
      return false;
    }
  }

  return true;
}

FieldPerp pow(const FieldPerp &lhs, const FieldPerp &rhs, REGION rgn) {
  TRACE("pow(FieldPerp, FieldPerp)");
  // Check if the inputs are allocated
  checkData(lhs);
  checkData(rhs);

  // Define and allocate the output result
  ASSERT1(lhs.getMesh() == rhs.getMesh());
  ASSERT1(lhs.getIndex() == rhs.getIndex());
  FieldPerp result(lhs.getMesh());
  result.allocate();
  result.setIndex(lhs.getIndex());

  const Region<IndPerp> &region = result.getMesh()->getRegionPerp(REGION_STRING(rgn));

  BOUT_FOR(i, region) {
    result[i] = ::pow(lhs[i], rhs[i]);
  }

  checkData(result);
  return result;
}

FieldPerp pow(const FieldPerp &lhs, BoutReal rhs, REGION rgn) {
  TRACE("pow(FieldPerp, BoutReal)");
  // Check if the inputs are allocated
  checkData(lhs);
  checkData(rhs);

  // Define and allocate the output result
  FieldPerp result(lhs.getMesh());
  result.allocate();
  result.setIndex(lhs.getIndex());

  const Region<IndPerp> &region = result.getMesh()->getRegionPerp(REGION_STRING(rgn));

  BOUT_FOR(i, region) {
    result[i] = ::pow(lhs[i], rhs);
  }

  checkData(result);
  return result;
}

FieldPerp pow(BoutReal lhs, const FieldPerp &rhs, REGION rgn) {
  TRACE("pow(lhs, FieldPerp)");
  // Check if the inputs are allocated
  checkData(lhs);
  checkData(rhs);

  // Define and allocate the output result
  FieldPerp result(rhs.getMesh());
  result.allocate();
  result.setIndex(rhs.getIndex());

  const Region<IndPerp> &region = result.getMesh()->getRegionPerp(REGION_STRING(rgn));

  BOUT_FOR(i, region) {
    result[i] = ::pow(lhs, rhs[i]);
  }

  checkData(result);
  return result;
}

#if CHECK > 2
void checkDataIsFiniteOnRegion(const FieldPerp &f, REGION region) {
  const Region<IndPerp> &new_region = f.getMesh()->getRegionPerp(REGION_STRING(region));
  
  // Do full checks
  BOUT_FOR_SERIAL(i, new_region) {
    if (!::finite(f[i])) {
      throw BoutException("FieldPerp: Operation on non-finite data at [%d][%d]\n", i.x(),
                          i.z());
    }
  }
}
#else
void checkDataIsFiniteOnRegion(const FieldPerp &UNUSED(f), REGION UNUSED(region)) {}
#endif


#if CHECK > 0
/// Check if the data is valid
void checkData(const FieldPerp &f, REGION region) {
  if (!f.isAllocated()) {
    throw BoutException("FieldPerp: Operation on empty data\n");
  }

  ASSERT3(f.getIndex() >= 0 && f.getIndex() < f.getMesh()->LocalNy);

  checkDataIsFiniteOnRegion(f, region);
}
#endif

#if CHECK > 2
void invalidateGuards(FieldPerp &var) {
  Mesh *localmesh = var.getMesh();

  const Region<IndPerp> &region_guards = localmesh->getRegionPerp("RGN_GUARDS");

  BOUT_FOR(i, region_guards) {
    var[i] = BoutNaN;
  }
}
#endif
