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

/// Initialise the derivative methods. Must be called before any derivatives are used
void Mesh::derivs_init(Options *options) {
  TRACE("Initialising derivatives");
  // Get the fraction of modes filtered out in FFT derivatives
  options->getSection("ddz")->get("fft_filter", fft_derivs_filter, 0.0);
}


/*******************************************************************************
 * Actual derivative operators
 *******************************************************************************/

/*******************************************************************************
 * Advection schemes
 *
 * Jan 2018  - Re-written to use iterators and handle staggering as different cases
 * Jan 2009  - Re-written to use Set*Stencil routines
 *******************************************************************************/

/*******************************************************************************
 * Flux conserving schemes
 *******************************************************************************/

/*******************************************************************************
 * Helper routines
 *******************************************************************************/

STAGGER Mesh::getStagger(const CELL_LOC inloc, const CELL_LOC outloc, const CELL_LOC allowedStaggerLoc) const {
  ASSERT1(outloc == inloc || (outloc == CELL_CENTRE && inloc == allowedStaggerLoc) ||
          (outloc == allowedStaggerLoc && inloc == CELL_CENTRE));

  if ( (!StaggerGrids) || outloc == inloc) return STAGGER::None;
  if (outloc == allowedStaggerLoc) {
    return STAGGER::C2L;
  } else {
    return STAGGER::L2C;
  }
}

STAGGER Mesh::getStagger(const CELL_LOC vloc, const CELL_LOC inloc, const CELL_LOC outloc, const CELL_LOC allowedStaggerLoc) const {
  ASSERT1(vloc == inloc);
  return getStagger(inloc, outloc, allowedStaggerLoc);
}
