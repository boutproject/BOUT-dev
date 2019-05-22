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

FieldPerp::FieldPerp(Mesh *localmesh, CELL_LOC location_in, int yindex_in,
      DirectionTypes directions)
    : Field(localmesh, location_in, directions),
      yindex(yindex_in) {
  if (fieldmesh) {
    nx = fieldmesh->LocalNx;
    nz = fieldmesh->LocalNz;
  }
}

FieldPerp::FieldPerp(BoutReal val, Mesh *localmesh) : FieldPerp(localmesh) {
  *this = val;
}

FieldPerp& FieldPerp::allocate() {
  if (data.empty()) {
    if (!fieldmesh) {
      // fieldmesh was not initialized when this field was initialized, so use
      // the global mesh and set some members to default values
      fieldmesh = bout::globals::mesh;
      nx = fieldmesh->LocalNx;
      nz = fieldmesh->LocalNz;
    }
    data.reallocate(nx * nz);
#if CHECK > 2
    invalidateGuards(*this);
#endif
  } else
    data.ensureUnique();

  return *this;
}

/***************************************************************
 *                         ASSIGNMENT 
 ***************************************************************/

FieldPerp &FieldPerp::operator=(const FieldPerp &rhs) {
  /// Check for self-assignment
  if (this == &rhs) {
    return (*this); // skip this assignment
  }

  copyFieldMembers(rhs);

  nx = rhs.nx;
  nz = rhs.nz;
  yindex = rhs.yindex;
  data = rhs.data;

  return *this;
}

FieldPerp & FieldPerp::operator=(const BoutReal rhs) {
  TRACE("FieldPerp = BoutReal");

  allocate();

  BOUT_FOR(i, getRegion("RGN_ALL")) { (*this)[i] = rhs; }

  return *this;
}

const Region<IndPerp> &FieldPerp::getRegion(REGION region) const {
  return fieldmesh->getRegionPerp(toString(region));
};
const Region<IndPerp> &FieldPerp::getRegion(const std::string &region_name) const {
  return fieldmesh->getRegionPerp(region_name);
};

//////////////// NON-MEMBER FUNCTIONS //////////////////

FieldPerp toFieldAligned(const FieldPerp& f, const REGION region) {
  return f.getCoordinates()->getParallelTransform().toFieldAligned(f, region);
}

FieldPerp fromFieldAligned(const FieldPerp& f, const REGION region) {
  return f.getCoordinates()->getParallelTransform().fromFieldAligned(f, region);
}

////////////// NON-MEMBER OVERLOADED OPERATORS //////////////

// Unary minus
FieldPerp operator-(const FieldPerp &f) { return -1.0 * f; }

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
    FieldPerp result{emptyFrom(f)};                                \
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

  FieldPerp result(f.getMesh(), f.getLocation(), y,
                   {f.getDirectionY(), f.getDirectionZ()});

  // Allocate memory
  result.allocate();
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
  ASSERT1(areFieldsCompatible(lhs, rhs));
  FieldPerp result{emptyFrom(lhs)};
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
  FieldPerp result{emptyFrom(lhs)};
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
  FieldPerp result{emptyFrom(rhs)};
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
