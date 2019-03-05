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

#include <cmath>

#include <bout/mesh.hxx>
#include <fieldperp.hxx>
#include <utils.hxx>
#include <boutexception.hxx>
#include <msg_stack.hxx>

FieldPerp::FieldPerp(Mesh *localmesh) : Field(localmesh) {
  if (fieldmesh) {
    nx = fieldmesh->LocalNx;
    nz = fieldmesh->LocalNz;
  }
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
      fieldmesh = bout::globals::mesh;
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

void FieldPerp::setLocation(CELL_LOC new_location) {
  AUTO_TRACE();
  if (getMesh()->StaggerGrids) {
    if (new_location == CELL_VSHIFT) {
      throw BoutException(
          "FieldPerp: CELL_VSHIFT cell location only makes sense for vectors");
    }
    if (new_location == CELL_DEFAULT) {
      new_location = CELL_CENTRE;
    }
    
    location = new_location;
  } else {
#if CHECK > 0
    if (new_location != CELL_CENTRE && new_location != CELL_DEFAULT) {
      throw BoutException("FieldPerp: Trying to set off-centre location on "
                          "non-staggered grid\n"
                          "         Did you mean to enable staggered grids?");
    }
#endif
    location = CELL_CENTRE;
  }

  // Ensures Coordinates object is initialized for this Field's location
  getCoordinates();
}

CELL_LOC FieldPerp::getLocation() const {
  AUTO_TRACE();
  return location;
}

/***************************************************************
 *                         ASSIGNMENT 
 ***************************************************************/

FieldPerp &FieldPerp::operator=(const FieldPerp &rhs) {
  /// Check for self-assignment
  if (this == &rhs) {
    return (*this); // skip this assignment
  }

  nx = rhs.nx;
  nz = rhs.nz;
  yindex = rhs.yindex;
  data = rhs.data;

  setLocation(rhs.location);
  return *this;
}

FieldPerp & FieldPerp::operator=(const BoutReal rhs) {
  TRACE("FieldPerp = BoutReal");

  allocate();

  BOUT_FOR(i, getRegion("RGN_ALL")) { (*this)[i] = rhs; }

  return *this;
}

/***************************************************************
 *                         OPERATORS 
 ***************************************************************/

#define FPERP_OP_FPERP(op, bop)                                   \
  FieldPerp& FieldPerp::operator op(const FieldPerp& rhs) {       \
    ASSERT1(getMesh() == rhs.getMesh());                          \
    ASSERT1(location == rhs.getLocation());                       \
    if (data.unique()) {                                          \
      checkData(rhs);                                             \
      /* Only reference to the data */                            \
      BOUT_FOR(i, getRegion("RGN_ALL")) { (*this)[i] op rhs[i]; } \
      checkData(*this);                                           \
    } else {                                                      \
      /* Shared with another FieldPerp */                         \
      (*this) = (*this)bop rhs;                                   \
    }                                                             \
    return *this;                                                 \
  }

FPERP_OP_FPERP(+=, +);
FPERP_OP_FPERP(-=, -);
FPERP_OP_FPERP(*=, *);
FPERP_OP_FPERP(/=, /);

#define FPERP_OP_FIELD(op, bop, ftype)                                               \
  FieldPerp& FieldPerp::operator op(const ftype& rhs) {                              \
    ASSERT1(getMesh() == rhs.getMesh());                                             \
    ASSERT1(location == rhs.getLocation());                                          \
    if (data.unique()) {                                                             \
      checkData(*this);                                                              \
      checkData(rhs);                                                                \
      /* Only reference to the data */                                               \
      BOUT_FOR(i, getRegion("RGN_ALL")) { (*this)[i] op rhs(i.x(), yindex, i.z()); } \
      checkData(*this);                                                              \
    } else {                                                                         \
      /* Shared with another FieldPerp */                                            \
      (*this) = (*this)bop rhs;                                                      \
    }                                                                                \
    return *this;                                                                    \
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
      BOUT_FOR(i, getRegion("RGN_ALL")) { (*this)[i] op rhs; }                           \
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
#define FPERP_FPERP_OP_FPERP(op)                                               \
  const FieldPerp operator op(const FieldPerp& lhs, const FieldPerp& rhs) {    \
    ASSERT1(lhs.getMesh() == rhs.getMesh());                                   \
    ASSERT1(rhs.getLocation() == rhs.getLocation());                           \
    checkData(lhs);                                                            \
    checkData(rhs);                                                            \
    FieldPerp result(lhs.getMesh());                                           \
    result.allocate();                                                         \
    result.setIndex(lhs.getIndex());                                           \
    result.setLocation(rhs.getLocation());                                     \
                                                                               \
    BOUT_FOR(i, result.getRegion("RGN_ALL")) { result[i] = lhs[i] op rhs[i]; } \
    checkData(result);                                                         \
    return result;                                                             \
  }

FPERP_FPERP_OP_FPERP(+);
FPERP_FPERP_OP_FPERP(-);
FPERP_FPERP_OP_FPERP(*);
FPERP_FPERP_OP_FPERP(/);

// Operator on FieldPerp and another field
#define FPERP_FPERP_OP_FIELD(op, ftype)                                 \
  const FieldPerp operator op(const FieldPerp& lhs, const ftype& rhs) { \
    ASSERT1(lhs.getMesh() == rhs.getMesh());                            \
    ASSERT1(rhs.getLocation() == rhs.getLocation());                    \
    checkData(lhs);                                                     \
    checkData(rhs);                                                     \
    FieldPerp result(lhs.getMesh());                                    \
    result.allocate();                                                  \
    result.setIndex(lhs.getIndex());                                    \
    result.setLocation(rhs.getLocation());                              \
                                                                        \
    BOUT_FOR(i, result.getRegion("RGN_ALL")) {                          \
      result[i] = lhs[i] op rhs(i.x(), lhs.getIndex(), i.z());          \
    }                                                                   \
    checkData(result);                                                  \
    return result;                                                      \
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
#define FPERP_FPERP_OP_REAL(op)                                             \
  const FieldPerp operator op(const FieldPerp& lhs, BoutReal rhs) {         \
    checkData(lhs);                                                         \
    checkData(rhs);                                                         \
    FieldPerp result(lhs.getMesh());                                        \
    result.allocate();                                                      \
    result.setIndex(lhs.getIndex());                                        \
    result.setLocation(lhs.getLocation());                                  \
                                                                            \
    BOUT_FOR(i, result.getRegion("RGN_ALL")) { result[i] = lhs[i] op rhs; } \
                                                                            \
    checkData(result);                                                      \
    return result;                                                          \
  }

FPERP_FPERP_OP_REAL(+);
FPERP_FPERP_OP_REAL(-);
FPERP_FPERP_OP_REAL(*);
FPERP_FPERP_OP_REAL(/);

#define FPERP_REAL_OP_FPERP(op)                                             \
  const FieldPerp operator op(BoutReal lhs, const FieldPerp& rhs) {         \
    checkData(lhs);                                                         \
    checkData(rhs);                                                         \
    FieldPerp result(rhs.getMesh());                                        \
    result.allocate();                                                      \
    result.setIndex(rhs.getIndex());                                        \
    result.setLocation(rhs.getLocation());                                  \
                                                                            \
    BOUT_FOR(i, result.getRegion("RGN_ALL")) { result[i] = lhs op rhs[i]; } \
                                                                            \
    checkData(result);                                                      \
    return result;                                                          \
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
#define FPERP_FUNC(name, func)                                     \
  const FieldPerp name(const FieldPerp& f, REGION rgn) {           \
    checkData(f);                                                  \
    TRACE(#name "(FieldPerp)");                                    \
    /* Check if the input is allocated */                          \
    ASSERT1(f.isAllocated());                                      \
    /* Define and allocate the output result */                    \
    FieldPerp result(f.getMesh());                                 \
    result.allocate();                                             \
    result.setIndex(f.getIndex());                                 \
    result.setLocation(f.getLocation());                           \
    BOUT_FOR(d, result.getRegion(rgn)) { result[d] = func(f[d]); } \
    checkData(result);                                             \
    return result;                                                 \
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

  BOUT_FOR(d, var.getRegion(rgn)) {
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
  result.setLocation(f.getLocation());
  BOUT_FOR(i, result.getRegion("RGN_ALL")) { result[i] = f(i, y); }

  checkData(result);
  return result;
}

BoutReal min(const FieldPerp &f, bool allpe, REGION rgn) {
  TRACE("FieldPerp::Min() %s", allpe ? "over all PEs" : "");

  checkData(f);

  const auto region = f.getRegion(rgn);
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

  const auto region = f.getRegion(rgn);
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

  BOUT_FOR_SERIAL(i, f.getRegion(rgn)) {
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
  ASSERT1(lhs.getLocation() == rhs.getLocation());
  FieldPerp result(lhs.getMesh());
  result.allocate();
  result.setIndex(lhs.getIndex());
  result.setLocation(lhs.getLocation());
  BOUT_FOR(i, result.getRegion(rgn)) { result[i] = ::pow(lhs[i], rhs[i]); }

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
  result.setLocation(lhs.getLocation());
  BOUT_FOR(i, result.getRegion(rgn)) { result[i] = ::pow(lhs[i], rhs); }

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
  result.setLocation(rhs.getLocation());
  BOUT_FOR(i, result.getRegion(rgn)) { result[i] = ::pow(lhs, rhs[i]); }

  checkData(result);
  return result;
}

#if CHECK > 2
void checkDataIsFiniteOnRegion(const FieldPerp &f, REGION region) {
  // Do full checks
  BOUT_FOR_SERIAL(i, f.getRegion(region)) {
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
  BOUT_FOR(i, var.getRegion("RGN_GUARDS")) { var[i] = BoutNaN; }
}
#endif
