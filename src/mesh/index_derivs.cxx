/**************************************************************************
 * Basic derivative methods in mesh index space
 *
 *
 * Four kinds of differencing methods:
 *
 * 1. First derivative DD*
 *    Central differencing e.g. Div(f)
 *
 * 2. Second derivatives D2D*2
 *    Central differencing e.g. Delp2(f)
 *
 * 3. Upwinding VDD*
 *    Terms like v*Grad(f)
 *
 * 4. Flux methods FDD* (e.g. flux conserving, limiting)
 *    Div(v*f)
 *
 * Changelog
 * =========
 *
 * 2014-11-22   Ben Dudson  <benjamin.dudson@york.ac.uk>
 *    o Moved here from sys/derivs, made part of Mesh
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

#include <bout/constants.hxx>
#include <derivs.hxx>
#include <fft.hxx>
#include <globals.hxx>
#include <interpolation.hxx>
#include <bout/constants.hxx>
#include <bout/derivs.hxx>
#include <bout/openmpwrap.hxx>

#include <msg_stack.hxx>
#include <stencils.hxx>
#include <utils.hxx>
#include <unused.hxx>

#include <cmath>
#include <stdlib.h>
#include <string.h>

#include <output.hxx>

#include <bout/mesh.hxx>

/*******************************************************************************
 * Limiters
 *******************************************************************************/

/// Van Leer limiter. Used in TVD code
BoutReal VANLEER(BoutReal r) { return r + fabs(r) / (1.0 + fabs(r)); }

// Superbee limiter
BoutReal SUPERBEE(BoutReal r) {
  return BOUTMAX(0.0, BOUTMIN(2. * r, 1.0), BOUTMIN(r, 2.));
}

/*******************************************************************************
 * Basic derivative methods.
 * All expect to have an input grid cell at the same location as the output
 * Hence convert cell centred values -> centred values, or left -> left
 *******************************************************************************/

//////////////////////// MUSCL scheme ///////////////////////

void DDX_KT_LR(const stencil &f, BoutReal &fLp, BoutReal &fRp, BoutReal &fLm,
               BoutReal &fRm) {
  // Limiter functions
  BoutReal phi = SUPERBEE((f.c - f.m) / (f.p - f.c));
  BoutReal phi_m = SUPERBEE((f.m - f.mm) / (f.c - f.m));
  BoutReal phi_p = SUPERBEE((f.p - f.c) / (f.pp - f.p));

  fLp = f.c + 0.5 * phi * (f.p - f.c);
  fRp = f.p - 0.5 * phi_p * (f.pp - f.p);

  fLm = f.m + 0.5 * phi_m * (f.c - f.m);
  fRm = f.c - 0.5 * phi * (f.p - f.c);
}

// du/dt = d/dx(f)  with maximum local velocity Vmax
BoutReal DDX_KT(const stencil &f, const stencil &u, const BoutReal Vmax) {
  BoutReal uLp, uRp, uLm, uRm;
  BoutReal fLp, fRp, fLm, fRm;

  DDX_KT_LR(u, uLp, uRp, uLm, uRm);
  DDX_KT_LR(f, fLp, fRp, fLm, fRm);

  BoutReal Fm = 0.5 * (fRm + fLm - Vmax * (uRm - uLm));
  BoutReal Fp = 0.5 * (fRp + fLp - Vmax * (uRp - uLp));

  return Fm - Fp;
}

/*******************************************************************************
 * Lookup tables of functions. Map between names, codes and functions
 *******************************************************************************/

/// Translate between DIFF_METHOD codes, and functions
struct DiffLookup {
  DIFF_METHOD method;
  Mesh::deriv_func func;     // Single-argument differencing function
  Mesh::upwind_func up_func; // Upwinding function
  Mesh::flux_func fl_func;   // Flux function
  operator Mesh::deriv_func (){
    return func;
  }
  operator Mesh::upwind_func (){
    return up_func;
  }
  operator Mesh::flux_func (){
    return fl_func;
  }
};

/// Translate between short names, long names and DIFF_METHOD codes
struct DiffNameLookup {
  DIFF_METHOD method;
  const char *label; // Short name
  const char *name;  // Long name
};

/// Differential function name/code lookup
static DiffNameLookup DiffNameTable[] = {
    {DIFF_U1, "U1", "First order upwinding"},
    {DIFF_U2, "U2", "Second order upwinding"},
    {DIFF_C2, "C2", "Second order central"},
    {DIFF_W2, "W2", "Second order WENO"},
    {DIFF_W3, "W3", "Third order WENO"},
    {DIFF_C4, "C4", "Fourth order central"},
    {DIFF_U3, "U3", "Third order upwinding"},
    {DIFF_U3, "U4", "Third order upwinding (Can't do 4th order yet)."},
    {DIFF_S2, "S2", "Smoothing 2nd order"},
    {DIFF_FFT, "FFT", "FFT"},
    {DIFF_SPLIT, "SPLIT", "Split into upwind and central"},
    {DIFF_DEFAULT, nullptr, nullptr}}; // Use to terminate the list

/// First derivative lookup table
static DiffLookup FirstDerivTable[] = {};

/// Second derivative lookup table
static DiffLookup SecondDerivTable[] = {};

/// Upwinding functions lookup table
static DiffLookup UpwindTable[] = {};

/// Flux functions lookup table
static DiffLookup FluxTable[] = {};

/// First staggered derivative lookup
static DiffLookup FirstStagDerivTable[] = {};

/// Second staggered derivative lookup
static DiffLookup SecondStagDerivTable[] = {};

/// Upwinding staggered lookup
static DiffLookup UpwindStagTable[] = {};

/// Flux staggered lookup
static DiffLookup FluxStagTable[] = {};

/*******************************************************************************
 * Routines to use the above tables to map between function codes, names
 * and pointers
 *******************************************************************************/


/// Test if a given DIFF_METHOD exists in a table
bool isImplemented(DiffLookup *table, DIFF_METHOD method) {
  int i = 0;
  do {
    if (table[i].method == method)
      return true;
    i++;
  } while (table[i].method != DIFF_DEFAULT);

  return false;
}


DiffLookup lookupFunc(DiffLookup * table, DIFF_METHOD method) {
  for (int i=0; ; ++i){
    if (table[i].method == method) {
      return table[i];
    }
    if (table[i].method == DIFF_DEFAULT){
      return table[i];
    }
  }
}

void printFuncName(DIFF_METHOD method) {
  // Find this entry
  int i = 0;
  do {
    if (DiffNameTable[i].method == method) {
      output_info.write(" %s (%s)\n", DiffNameTable[i].name, DiffNameTable[i].label);
      return;
    }
    i++;
  } while (DiffNameTable[i].method != DIFF_DEFAULT);

  // None
  output_error.write(" == INVALID DIFFERENTIAL METHOD ==\n");
}

/// This function is used during initialisation only (i.e. doesn't need to be particularly
/// fast) Returns DIFF_METHOD, rather than function so can be applied to central and
/// upwind tables
DiffLookup lookupFunc(DiffLookup *table, const std::string & label){

  // Loop through the name lookup table
  for (int i = 0; DiffNameTable[i].method != DIFF_DEFAULT; ++i) {
    if (strcasecmp(label.c_str(), DiffNameTable[i].label) == 0) { // Whole match
      auto method=DiffNameTable[i].method;
      if (isImplemented(table, method)) {
        printFuncName(method);
        for (int j=0;;++j){
          if (table[j].method == method){
            return table[j];
          }
        }
      }
    }
  }

  // No exact match, so throw
  std::string avail{};
  for (int i = 0; DiffNameTable[i].method != DIFF_DEFAULT; ++i) {
    avail += DiffNameTable[i].label;
    avail += "\n";
  }
  throw BoutException("Unknown option %s.\nAvailable options are:\n%s", label.c_str(),
                      avail.c_str());
}


/*******************************************************************************
 * Default functions
 *
 *
 *******************************************************************************/

// Central -> Central (or Left -> Left) functions
Mesh::deriv_func fDDX, fDDY, fDDZ;       ///< Differencing methods for each dimension
Mesh::deriv_func fD2DX2, fD2DY2, fD2DZ2; ///< second differential operators
Mesh::upwind_func fVDDX, fVDDY, fVDDZ;   ///< Upwind functions in the three directions
Mesh::flux_func fFDDX, fFDDY, fFDDZ;     ///< Default flux functions

// Central -> Left (or Left -> Central) functions
Mesh::deriv_func sfDDX, sfDDY, sfDDZ;
Mesh::deriv_func sfD2DX2, sfD2DY2, sfD2DZ2;
Mesh::flux_func sfVDDX, sfVDDY, sfVDDZ;
Mesh::flux_func sfFDDX, sfFDDY, sfFDDZ;

/*******************************************************************************
 * Initialisation
 *******************************************************************************/

/// Set the derivative method, given a table and option name
template <typename T, typename Ts>
void derivs_set(std::vector<Options *> options, DiffLookup *table, DiffLookup *stable,
                const std::string &name, const std::string &def, T &f, Ts &sf,
                bool staggerGrids) {
  TRACE("derivs_set()");
  output_info.write("\t%-12s: ", name.c_str());
  string label = def;
  for (auto &opts : options) {
    if (opts->isSet(name)) {
      opts->get(name, label, "");
      break;
    }
  }

  f = lookupFunc(table, label); // Find the function

  label = def;
  if (staggerGrids) {
    output_info.write("\tStag. %-6s: ", name.c_str());
    for (auto &_name : {name + "stag", name}) {
      for (auto &opts : options) {
        if (opts->isSet(_name)) {
          opts->get(_name, label, "");
          sf = lookupFunc(stable, label); // Find the function
          return;
        }
      }
    }
  }
  sf = lookupFunc(stable, label); // Find the function
}

/// Initialise derivatives from options
void derivs_initialise(Options *optionbase, std::string sec, bool staggerGrids,
                       Mesh::deriv_func &fdd, Mesh::deriv_func &sfdd,
                       Mesh::deriv_func &fd2d, Mesh::deriv_func &sfd2d,
                       Mesh::upwind_func &fu, Mesh::flux_func &sfu, Mesh::flux_func &ff,
                       Mesh::flux_func &sff) {
  std::vector<Options *> options = {optionbase->getSection(sec),
                                    optionbase->getSection("diff")};
  derivs_set(options, FirstDerivTable, FirstStagDerivTable, "First", "C2", fdd, sfdd,
             staggerGrids);

  derivs_set(options, SecondDerivTable, SecondStagDerivTable, "Second", "C2", fd2d, sfd2d,
             staggerGrids);

  derivs_set(options, UpwindTable, UpwindStagTable, "Upwind", "U1", fu, sfu,
             staggerGrids);

  derivs_set(options, FluxTable, FluxStagTable, "Flux", "U1", ff, sff, staggerGrids);
}

/// Initialise the derivative methods. Must be called before any derivatives are used
void Mesh::derivs_init(Options *options) {
  TRACE("Initialising derivatives");

  output_info.write("Setting X differencing methods\n");
  derivs_initialise(options, "ddx", StaggerGrids, fDDX, sfDDX, fD2DX2, sfD2DX2, fVDDX,
                    sfVDDX, fFDDX, sfFDDX);

  if ((fDDX == nullptr) || (fD2DX2 == nullptr))
    throw BoutException("FFT cannot be used in X\n");

  output_info.write("Setting Y differencing methods\n");
  derivs_initialise(options, "ddy", StaggerGrids, fDDY, sfDDY, fD2DY2, sfD2DY2, fVDDY,
                    sfVDDY, fFDDY, sfFDDY);

  if ((fDDY == nullptr) || (fD2DY2 == nullptr))
    throw BoutException("FFT cannot be used in Y\n");

  output_info.write("Setting Z differencing methods\n");
  derivs_initialise(options, "ddz", StaggerGrids, fDDZ, sfDDZ, fD2DZ2, sfD2DZ2, fVDDZ,
                    sfVDDZ, fFDDZ, sfFDDZ);

  // Get the fraction of modes filtered out in FFT derivatives
  options->getSection("ddz")->get("fft_filter", fft_derivs_filter, 0.0);
}

/*******************************************************************************
 * Apply differential operators. These are fairly brain-dead functions
 * which apply a derivative function to a field (sort of like map). Decisions
 * of what to apply are made in the DDX,DDY and DDZ functions lower down.
 *
 * loc  is the cell location of the result
 *******************************************************************************/

// X derivative
template<typename T>
const T Mesh::applyXdiff(const T &var, Mesh::deriv_func func,
                               CELL_LOC outloc, REGION region) {

  static_assert(std::is_base_of<Field2D, T>::value || std::is_base_of<Field3D, T>::value,
                "applyXdiff only works on Field2D or Field3D input");

  // Check that the mesh is correct
  ASSERT1(this == var.getMesh());
  // Check that the input variable has data
  ASSERT1(var.isAllocated());

  CELL_LOC inloc = var.getLocation();
  if (outloc == CELL_DEFAULT)
    outloc = inloc;

  // Determine what the non-centre allowed location is for outloc
  CELL_LOC allowedStaggerLoc;
  allowedStaggerLoc = CELL_XLOW;

  // Allowed staggers:
  ASSERT1(outloc == inloc || (outloc == CELL_CENTRE && inloc == allowedStaggerLoc) ||
          (outloc == allowedStaggerLoc && inloc == CELL_CENTRE));

  int nPoint, nGuard;
  nPoint = var.getNx();
  nGuard = xstart;

  if (nPoint == 1) {
    auto tmp = T(0., this);
    tmp.setLocation(outloc);
    return tmp;
  }

  T result(this);
  result.allocate(); // Make sure data allocated
  result.setLocation(outloc);

  if (StaggerGrids && (outloc != inloc)) {
    // Staggered differencing

    if (nGuard > 1) {
      // More than one guard cell, so set pp and mm values
      // This allows higher-order methods to be used
      if (outloc == allowedStaggerLoc) {
        applyDiffKernel<DIRECTION::X, STAGGER::C2L, 2>(var, func, result, region);
      } else {
        applyDiffKernel<DIRECTION::X, STAGGER::L2C, 2>(var, func, result, region);
      }
    } else {
      // Only one guard cell, so no pp or mm values
      if (outloc == allowedStaggerLoc) {
        applyDiffKernel<DIRECTION::X, STAGGER::C2L, 1>(var, func, result, region);
      } else {
        applyDiffKernel<DIRECTION::X, STAGGER::L2C, 1>(var, func, result, region);
      }
    }
  } else {
    // Non-staggered differencing
    if (nGuard > 1) {
      // More than one guard cell, so set pp and mm values
      // This allows higher-order methods to be used
      applyDiffKernel<DIRECTION::X, STAGGER::None, 2>(var, func, result, region);
    } else {
      // Only one guard cell, so no pp or mm values
      applyDiffKernel<DIRECTION::X, STAGGER::None, 1>(var, func, result, region);
    }
  }

  return result;
}

// Y derivative
template <typename T>
const T Mesh::applyYdiff(const T &var, Mesh::deriv_func func, CELL_LOC outloc,
                         REGION region) {

  static_assert(std::is_base_of<Field2D, T>::value || std::is_base_of<Field3D, T>::value,
                "applyYdiff only works on Field2D or Field3D input");

  ASSERT1(this == var.getMesh());
  // Check that the input variable has data
  ASSERT1(var.isAllocated());

  // Cell location of the input field
  CELL_LOC inloc = var.getLocation();
  if (outloc == CELL_DEFAULT)
    outloc = inloc;

  // Determine what the non-centre allowed location is for outloc
  CELL_LOC allowedStaggerLoc;
  allowedStaggerLoc = CELL_YLOW;

  // Allowed staggers:
  ASSERT1(outloc == inloc || (outloc == CELL_CENTRE && inloc == allowedStaggerLoc) ||
          (outloc == allowedStaggerLoc && inloc == CELL_CENTRE));

  int nPoint, nGuard;
  nPoint = var.getNy();
  nGuard = ystart;

  if (nPoint == 1) {
    auto tmp = T(0., this);
    tmp.setLocation(outloc);
    return tmp;
  }

  T result(this);
  result.allocate(); // Make sure data allocated
  result.setLocation(outloc);

  if (std::is_base_of<Field3D, T>::value && var.hasYupYdown() &&
      ((&var.yup() != &var) || (&var.ydown() != &var))) {
    // Field "var" has distinct yup and ydown fields which
    // will be used to calculate a derivative along
    // the magnetic field

    if (StaggerGrids && (outloc != inloc)) {
      // Staggered differencing
      if (outloc == allowedStaggerLoc) {
        applyDiffKernel<DIRECTION::YOrthogonal, STAGGER::C2L, 1>(var, func, result,
                                                                 region);
      } else {
        applyDiffKernel<DIRECTION::YOrthogonal, STAGGER::L2C, 1>(var, func, result,
                                                                 region);
      }
    } else {
      // Non-staggered
      applyDiffKernel<DIRECTION::YOrthogonal, STAGGER::None, 1>(var, func, result,
                                                                region);
    }
  } else {
    // var has no yup/ydown fields, so we need to shift into field-aligned coordinates
    T var_fa = toFieldAligned(var);

    if (StaggerGrids && (outloc != inloc)) {
      // Staggered differencing
      if (nGuard > 1) {
        // More than one guard cell, so set pp and mm values
        // This allows higher-order methods to be used
        if (outloc == allowedStaggerLoc) {
          applyDiffKernel<DIRECTION::Y, STAGGER::C2L, 2>(var_fa, func, result, region);
        } else {
          applyDiffKernel<DIRECTION::Y, STAGGER::L2C, 2>(var_fa, func, result, region);
        }
      } else {
        // Only one guard cell, so no pp or mm values
        if (outloc == allowedStaggerLoc) {
          applyDiffKernel<DIRECTION::Y, STAGGER::C2L, 1>(var_fa, func, result, region);
        } else {
          applyDiffKernel<DIRECTION::Y, STAGGER::L2C, 1>(var_fa, func, result, region);
        }
      }
    } else {
      // Non-staggered differencing

      if (nGuard > 1) {
        // More than one guard cell, so set pp and mm values
        // This allows higher-order methods to be used
        applyDiffKernel<DIRECTION::Y, STAGGER::None, 2>(var_fa, func, result, region);
      } else {
        // Only one guard cell, so no pp or mm values
        applyDiffKernel<DIRECTION::Y, STAGGER::None, 1>(var_fa, func, result, region);
      }
    }

    // Shift result back
    result = fromFieldAligned(result);
  }

  return result;
}

// Z derivative
template <typename T>
const T Mesh::applyZdiff(const T &var, Mesh::deriv_func func, CELL_LOC outloc,
                         REGION region) {

  static_assert(std::is_base_of<Field2D, T>::value || std::is_base_of<Field3D, T>::value,
                "applyZdiff only works on Field2D or Field3D input");

  ASSERT1(this == var.getMesh());
  // Check that the input variable has data
  ASSERT1(var.isAllocated());

  CELL_LOC inloc = var.getLocation();
  if (outloc == CELL_DEFAULT)
    outloc = inloc;

  // Determine what the non-centre allowed location is for outloc
  CELL_LOC allowedStaggerLoc;
  allowedStaggerLoc = CELL_ZLOW;

  // Allowed staggers:
  ASSERT1(outloc == inloc || (outloc == CELL_CENTRE && inloc == allowedStaggerLoc) ||
          (outloc == allowedStaggerLoc && inloc == CELL_CENTRE));

  int nPoint, nGuard;
  nPoint = var.getNz();
  nGuard = 2;

  if (nPoint == 1) {
    auto tmp = T(0., this);
    tmp.setLocation(outloc);
    return tmp;
  }

  T result(this);
  result.allocate(); // Make sure data allocated
  result.setLocation(outloc);

  if (StaggerGrids && (outloc != inloc)) {
    // Staggered differencing
    if (nGuard > 1) {
      // More than one guard cell, so set pp and mm values
      // This allows higher-order methods to be used
      if (outloc == allowedStaggerLoc) {
        applyDiffKernel<DIRECTION::Z, STAGGER::C2L, 2>(var, func, result, region);
      } else {
        applyDiffKernel<DIRECTION::Z, STAGGER::L2C, 2>(var, func, result, region);
      }
    } else {
      // Only one guard cell, so no pp or mm values
      if (outloc == allowedStaggerLoc) {
        applyDiffKernel<DIRECTION::Z, STAGGER::C2L, 1>(var, func, result, region);
      } else {
        applyDiffKernel<DIRECTION::Z, STAGGER::L2C, 1>(var, func, result, region);
      }
    }
  } else {
    // Non-staggered differencing
    if (nGuard > 1) {
      // More than one guard cell, so set pp and mm values
      // This allows higher-order methods to be used
      applyDiffKernel<DIRECTION::Z, STAGGER::None, 2>(var, func, result, region);
    } else {
      // Only one guard cell, so no pp or mm values
      applyDiffKernel<DIRECTION::Z, STAGGER::None, 1>(var, func, result, region);
    }
  }

  return result;
}

template <DIRECTION direction, STAGGER stagger, int nGuard, typename T>
void Mesh::applyDiffKernel(const T &var, Mesh::deriv_func func, T &result,
                           REGION region) {
  BOUT_OMP(parallel) {
    stencil s;
    BOUT_FOR_INNER(i, result.getRegion(region)) {
      populateStencil<direction, stagger, nGuard, T>(s, var, i);
      result[i] = func(s);
    }
  }
  return;
}
/*******************************************************************************
 * First central derivatives
 *******************************************************************************/

const STAGGER Mesh::getStagger(const CELL_LOC inloc, const CELL_LOC outloc, const CELL_LOC allowedStaggerLoc) {
  ASSERT1(outloc == inloc || (outloc == CELL_CENTRE && inloc == allowedStaggerLoc) ||
          (outloc == allowedStaggerLoc && inloc == CELL_CENTRE));

  if ( (!StaggerGrids) || outloc == inloc) return STAGGER::None;
  if (outloc == allowedStaggerLoc) {
    return STAGGER::C2L;
  } else {
    return STAGGER::L2C;
  }
}

const STAGGER Mesh::getStagger(const CELL_LOC vloc, const CELL_LOC inloc, const CELL_LOC outloc, const CELL_LOC allowedStaggerLoc) {
  ASSERT1(vloc == inloc);
  return getStagger(inloc, outloc, allowedStaggerLoc);
}

template<DIRECTION direction>
const CELL_LOC Mesh::getAllowedStaggerLoc() {
  switch(direction) {
  case(DIRECTION::X):
    return CELL_XLOW;
  case(DIRECTION::Y):
  case(DIRECTION::YOrthogonal):
  case(DIRECTION::YAligned):    
    return CELL_YLOW;
  case(DIRECTION::Z):
    return CELL_ZLOW;
  }
};


template<DIRECTION direction>
const int Mesh::getNpoints() {
  switch(direction) {
  case(DIRECTION::X):
    return LocalNx;
  case(DIRECTION::Y):
  case(DIRECTION::YOrthogonal):
  case(DIRECTION::YAligned):    
    return LocalNy;
  case(DIRECTION::Z):
    return LocalNz;
  }
};

template<typename T, DIRECTION direction, int order>
const T Mesh::indexStandardDerivative(const T &f, CELL_LOC outloc, DIFF_METHOD method, REGION region) {
  // Checks
  static_assert(std::is_base_of<Field2D, T>::value || std::is_base_of<Field3D, T>::value,
                "indexDDX only works on Field2D or Field3D input");
  // Check that the mesh is correct
  ASSERT1(this == f.getMesh());
  // Check that the input variable has data
  ASSERT1(f.isAllocated());

  // Define properties of this approach
  const CELL_LOC allowedStaggerLoc = getAllowedStaggerLoc<direction>();

  // Handle the staggering
  const CELL_LOC inloc = f.getLocation(); // Input location
  if (outloc == CELL_DEFAULT)
    outloc = inloc;
  const STAGGER stagger = getStagger(inloc, outloc, allowedStaggerLoc);

  // Check for early exit
  const int nPoint = getNpoints<direction>();

  if (nPoint == 1) {
    auto tmp = T(0., this);
    tmp.setLocation(outloc);
    return tmp;
  }
  
  // Lookup the method
  auto derivativeStore = DerivativeStore<T>{}.getInstance();
  typename DerivativeStore<T>::standardFunc derivativeMethod;
  
  if (order == 1) {
    derivativeMethod = derivativeStore.getStandardDerivative(DIFF_METHOD_STRING(method), direction, stagger);
  } else if (order == 2) {
    derivativeMethod = derivativeStore.getStandard2ndDerivative(DIFF_METHOD_STRING(method), direction, stagger);
  } else if (order == 4) {
    derivativeMethod = derivativeStore.getStandard4thDerivative(DIFF_METHOD_STRING(method), direction, stagger);
  } else {
    throw BoutException("Invalid order used in indexStandardDerivative.");
  }
  
  // Create the result field
  T result(this);
  result.allocate(); // Make sure data allocated
  result.setLocation(outloc);

  // Apply method
  derivativeMethod(f, result, region);

  return result;
}

////////////// X DERIVATIVE /////////////////

template<typename T>
const T Mesh::indexDDX(const T &f, CELL_LOC outloc, DIFF_METHOD method, REGION region) {
  return indexStandardDerivative<T, DIRECTION::X, 1>(f, outloc, method, region);
}

template<typename T>
const T Mesh::indexD2DX2(const T &f, CELL_LOC outloc, DIFF_METHOD method, REGION region) {
  return indexStandardDerivative<T, DIRECTION::X, 2>(f, outloc, method, region);
}

template<typename T>
const T Mesh::indexD4DX4(const T &f, CELL_LOC outloc, DIFF_METHOD method, REGION region) {
  return indexStandardDerivative<T, DIRECTION::X, 4>(f, outloc, method, region);
}

////////////// Y DERIVATIVE /////////////////

template<typename T>
const T Mesh::indexDDY(const T &f, CELL_LOC outloc, DIFF_METHOD method, REGION region) {
  return indexStandardDerivative<T, DIRECTION::Y, 1>(f, outloc, method, region);
}

template<typename T>
const T Mesh::indexD2DY2(const T &f, CELL_LOC outloc, DIFF_METHOD method, REGION region) {
  return indexStandardDerivative<T, DIRECTION::Y, 2>(f, outloc, method, region);
}

template<typename T>
const T Mesh::indexD4DY4(const T &f, CELL_LOC outloc, DIFF_METHOD method, REGION region) {
  return indexStandardDerivative<T, DIRECTION::Y, 4>(f, outloc, method, region);
}

////////////// Z DERIVATIVE /////////////////
template<typename T>
const T Mesh::indexDDZ(const T &f, CELL_LOC outloc, DIFF_METHOD method, REGION region) {
  return indexStandardDerivative<T, DIRECTION::Z, 1>(f, outloc, method, region);
}

template<typename T>
const T Mesh::indexD2DZ2(const T &f, CELL_LOC outloc, DIFF_METHOD method, REGION region) {
  return indexStandardDerivative<T, DIRECTION::Z, 2>(f, outloc, method, region);
}

template<typename T>
const T Mesh::indexD4DZ4(const T &f, CELL_LOC outloc, DIFF_METHOD method, REGION region) {
  return indexStandardDerivative<T, DIRECTION::Z, 4>(f, outloc, method, region);
}

/*******************************************************************************
 * Advection schemes
 *
 * Jan 2018  - Re-written to use iterators and handle staggering as different cases
 * Jan 2009  - Re-written to use Set*Stencil routines
 *******************************************************************************/

/*******************************************************************************
 * Flux conserving schemes
 *******************************************************************************/


template<typename T, DIRECTION direction, DERIV derivType>
const T Mesh::indexFlowDerivative(const T &vel, const T &f, CELL_LOC outloc, DIFF_METHOD method, REGION region) {
  // Checks
  static_assert(std::is_base_of<Field2D, T>::value || std::is_base_of<Field3D, T>::value,
                "indexDDX only works on Field2D or Field3D input");
  // Check that the mesh is correct
  ASSERT1(this == f.getMesh());
  ASSERT1(this == v.getMesh());  
  // Check that the input variable has data
  ASSERT1(f.isAllocated());
  ASSERT1(v.isAllocated());  

  // Define properties of this approach
  const CELL_LOC allowedStaggerLoc = getAllowedStaggerLoc<direction>();

  // Handle the staggering
  const CELL_LOC inloc = f.getLocation(); // Input locations
  const CELL_LOC vloc = vel.getLocation();
  if (outloc == CELL_DEFAULT)
    outloc = inloc;
  const STAGGER stagger = getStagger(vloc, inloc, outloc, allowedStaggerLoc);

  // Check for early exit
  const int nPoint = getNpoints<direction>();

  if (nPoint == 1) {
    auto tmp = T(0., this);
    tmp.setLocation(outloc);
    return tmp;
  }
  
  // Lookup the method
  auto derivativeStore = DerivativeStore<T>{}.getInstance();
  typename DerivativeStore<T>::upwindFunc derivativeMethod;
  if (derivType == DERIV::Upwind) {
    derivativeMethod = derivativeStore.getUpwindDerivative(DIFF_METHOD_STRING(method), direction, stagger);
  } else if (derivType == DERIV::Flux) {
    derivativeMethod = derivativeStore.getFluxDerivative(DIFF_METHOD_STRING(method), direction, stagger);
  } else {
    throw BoutException("Invalid derivative type in call to indexFlowDerivative.");
  }
  
  // Create the result field
  T result(this);
  result.allocate(); // Make sure data allocated
  result.setLocation(outloc);

  // Apply method
  derivativeMethod(vel, f, result, region);

  return result;
}

////////////// X DERIVATIVE /////////////////

template<typename T>
const T Mesh::indexVDDX(const T& vel, const T &f, CELL_LOC outloc, DIFF_METHOD method, REGION region) {
  return indexFlowDerivative<T, DIRECTION::X, DERIV::Upwind>(vel, f, outloc, method, region);
}

template<typename T>
const T Mesh::indexFDDX(const T& vel, const T &f, CELL_LOC outloc, DIFF_METHOD method, REGION region) {
  return indexFlowDerivative<T, DIRECTION::X, DERIV::Flux>(vel, f, outloc, method, region);
}

////////////// Y DERIVATIVE /////////////////

template<typename T>
const T Mesh::indexVDDY(const T& vel, const T &f, CELL_LOC outloc, DIFF_METHOD method, REGION region) {
  return indexFlowDerivative<T, DIRECTION::Y, DERIV::Upwind>(vel, f, outloc, method, region);
}

template<typename T>
const T Mesh::indexFDDY(const T& vel, const T &f, CELL_LOC outloc, DIFF_METHOD method, REGION region) {
  return indexFlowDerivative<T, DIRECTION::Y, DERIV::Flux>(vel, f, outloc, method, region);
}

////////////// Z DERIVATIVE /////////////////

template<typename T>
const T Mesh::indexVDDZ(const T& vel, const T &f, CELL_LOC outloc, DIFF_METHOD method, REGION region) {
  return indexFlowDerivative<T, DIRECTION::Z, DERIV::Upwind>(vel, f, outloc, method, region);
}

template<typename T>
const T Mesh::indexFDDZ(const T& vel, const T &f, CELL_LOC outloc, DIFF_METHOD method, REGION region) {
  return indexFlowDerivative<T, DIRECTION::Z, DERIV::Flux>(vel, f, outloc, method, region);
}
