/**************************************************************************
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

#include "bout/build_config.hxx"

#include "bout/traits.hxx"
#include <bout/index_derivs.hxx>
#include <bout/mesh.hxx>
#include <msg_stack.hxx>
#include <unused.hxx>

/*******************************************************************************
 * Helper routines
 *******************************************************************************/

/// Initialise the derivative methods. Must be called before any derivatives are used
void Mesh::derivs_init(Options* options) {
  TRACE("Initialising derivatives");
  // For each direction need to set what the default method is for each type
  // of derivative.
  DerivativeStore<Field3D>::getInstance().initialise(options);
  DerivativeStore<Field2D>::getInstance().initialise(options);
  // Get the fraction of modes filtered out in FFT derivatives
  options->getSection("ddz")->get("fft_filter", fft_derivs_filter, 0.0);
}

STAGGER Mesh::getStagger(const CELL_LOC inloc, const CELL_LOC outloc,
                         const CELL_LOC allowedStaggerLoc) const {
  TRACE("Mesh::getStagger -- three arguments");
  ASSERT1(outloc == inloc || (outloc == CELL_CENTRE && inloc == allowedStaggerLoc)
          || (outloc == allowedStaggerLoc && inloc == CELL_CENTRE));

  if ((!StaggerGrids) || outloc == inloc)
    return STAGGER::None;
  if (outloc == allowedStaggerLoc) {
    return STAGGER::C2L;
  } else {
    return STAGGER::L2C;
  }
}

STAGGER Mesh::getStagger(const CELL_LOC vloc, MAYBE_UNUSED(const CELL_LOC inloc),
                         const CELL_LOC outloc, const CELL_LOC allowedStaggerLoc) const {
  TRACE("Mesh::getStagger -- four arguments");
  ASSERT1(inloc == outloc);
  ASSERT1(vloc == inloc || (vloc == CELL_CENTRE && inloc == allowedStaggerLoc)
          || (vloc == allowedStaggerLoc && inloc == CELL_CENTRE));
  return getStagger(vloc, outloc, allowedStaggerLoc);
}

////////////////////// FIRST DERIVATIVES /////////////////////

/// central, 2nd order
REGISTER_STANDARD_DERIVATIVE(DDX_C2, "C2", 1, DERIV::Standard) {
  return 0.5 * (f.p - f.m);
}

/// central, 4th order
REGISTER_STANDARD_DERIVATIVE(DDX_C4, "C4", 2, DERIV::Standard) {
  return (8. * f.p - 8. * f.m + f.mm - f.pp) / 12.;
}

/// Central WENO method, 2nd order (reverts to 1st order near shocks)
REGISTER_STANDARD_DERIVATIVE(DDX_CWENO2, "W2", 1, DERIV::Standard) {
  BoutReal isl, isr, isc;  // Smoothness indicators
  BoutReal al, ar, ac, sa; // Un-normalised weights
  BoutReal dl, dr, dc;     // Derivatives using different stencils

  dc = 0.5 * (f.p - f.m);
  dl = f.c - f.m;
  dr = f.p - f.c;

  isl = SQ(dl);
  isr = SQ(dr);
  isc = (13. / 3.) * SQ(f.p - 2. * f.c + f.m) + 0.25 * SQ(f.p - f.m);

  al = 0.25 / SQ(WENO_SMALL + isl);
  ar = 0.25 / SQ(WENO_SMALL + isr);
  ac = 0.5 / SQ(WENO_SMALL + isc);
  sa = al + ar + ac;

  return (al * dl + ar * dr + ac * dc) / sa;
}

// Smoothing 2nd order derivative
REGISTER_STANDARD_DERIVATIVE(DDX_S2, "S2", 2, DERIV::Standard) {

  // 4th-order differencing
  BoutReal result = (8. * f.p - 8. * f.m + f.mm - f.pp) / 12.;

  result += SIGN(f.c) * (f.pp - 4. * f.p + 6. * f.c - 4. * f.m + f.mm) / 12.;

  return result;
}

/// Also CWENO3 but needs an upwind op so define later.

//////////////////////////////
//--- Second order derivatives
//////////////////////////////

/// Second derivative: Central, 2nd order
REGISTER_STANDARD_DERIVATIVE(D2DX2_C2, "C2", 1, DERIV::StandardSecond) {
  return f.p + f.m - 2. * f.c;
}

/// Second derivative: Central, 4th order
REGISTER_STANDARD_DERIVATIVE(D2DX2_C4, "C4", 2, DERIV::StandardSecond) {
  return (-f.pp + 16. * f.p - 30. * f.c + 16. * f.m - f.mm) / 12.;
}

//////////////////////////////
//--- Fourth order derivatives
//////////////////////////////
REGISTER_STANDARD_DERIVATIVE(D4DX4_C2, "C2", 2, DERIV::StandardFourth) {
  return (f.pp - 4. * f.p + 6. * f.c - 4. * f.m + f.mm);
}

////////////////////////////////////////////////////////////////////////////////
/// Staggered methods
///
/// Map Centre -> Low or Low -> Centre
///
/// These expect the output grid cell to be at a different location to the input
///
/// The stencil no longer has a value in 'C' (centre)
/// instead, points are shifted as follows:
///
///  mm  -> -3/2 h
///  m   -> -1/2 h
///  p   -> +1/2 h
///  pp  -? +3/2 h
///
/// NOTE: Cell widths (dx, dy, dz) are currently defined as centre->centre
/// for the methods above. This is currently not taken account of, so large
/// variations in cell size will cause issues.
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
/// Standard methods -- first order
////////////////////////////////////////////////////////////////////////////////
REGISTER_STANDARD_DERIVATIVE_STAGGERED(DDX_C2_stag, "C2", 1, DERIV::Standard) {
  return f.p - f.m;
}

REGISTER_STANDARD_DERIVATIVE_STAGGERED(DDX_C4_stag, "C4", 2, DERIV::Standard) {
  return (27. * (f.p - f.m) - (f.pp - f.mm)) / 24.;
}

////////////////////////////////////////////////////////////////////////////////
/// Standard methods -- second order
////////////////////////////////////////////////////////////////////////////////
REGISTER_STANDARD_DERIVATIVE_STAGGERED(D2DX2_C2_stag, "C2", 2, DERIV::StandardSecond) {
  return (f.pp + f.mm - f.p - f.m) / 2.;
}
