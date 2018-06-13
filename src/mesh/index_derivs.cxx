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
#include <bout/scorepwrapper.hxx>
#include <bout/indexoffset.hxx>

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

const BoutReal WENO_SMALL = 1.0e-8; // Small number for WENO schemes

////////////////////// FIRST DERIVATIVES /////////////////////

/// central, 2nd order
BoutReal DDX_C2(stencil &f) { return 0.5 * (f.p - f.m); }

/// central, 4th order
BoutReal DDX_C4(stencil &f) { return (8. * f.p - 8. * f.m + f.mm - f.pp) / 12.; }

/// Central WENO method, 2nd order (reverts to 1st order near shocks)
BoutReal DDX_CWENO2(stencil &f) {
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
BoutReal DDX_S2(stencil &f) {

  // 4th-order differencing
  BoutReal result = (8. * f.p - 8. * f.m + f.mm - f.pp) / 12.;

  result += SIGN(f.c) * (f.pp - 4. * f.p + 6. * f.c - 4. * f.m + f.mm) / 12.;

  return result;
}

///////////////////// SECOND DERIVATIVES ////////////////////

/// Second derivative: Central, 2nd order
BoutReal D2DX2_C2(stencil &f) { return f.p + f.m - 2. * f.c; }

/// Second derivative: Central, 4th order
BoutReal D2DX2_C4(stencil &f) {
  return (-f.pp + 16. * f.p - 30. * f.c + 16. * f.m - f.mm) / 12.;
}

//////////////////////// UPWIND METHODS ///////////////////////

/// Upwinding: Central, 2nd order
BoutReal VDDX_C2(BoutReal vc, stencil &f) { return vc * 0.5 * (f.p - f.m); }

/// Upwinding: Central, 4th order
BoutReal VDDX_C4(BoutReal vc, stencil &f) {
  return vc * (8. * f.p - 8. * f.m + f.mm - f.pp) / 12.;
}

/// upwind, 1st order
BoutReal VDDX_U1(BoutReal vc, stencil &f) {
  return vc >= 0.0 ? vc * (f.c - f.m) : vc * (f.p - f.c);
}

/// upwind, 2nd order
BoutReal VDDX_U2(BoutReal vc, stencil &f) {
  return vc >= 0.0 ? vc * (1.5 * f.c - 2.0 * f.m + 0.5 * f.mm)
                   : vc * (-0.5 * f.pp + 2.0 * f.p - 1.5 * f.c);
}

/// upwind, 3rd order
BoutReal VDDX_U3(BoutReal vc, stencil &f) {
  return vc >= 0.0 ? vc*(4.*f.p - 12.*f.m + 2.*f.mm + 6.*f.c)/12.
    : vc*(-4.*f.m + 12.*f.p - 2.*f.pp - 6.*f.c)/12.;
}

/// 3rd-order WENO scheme
BoutReal VDDX_WENO3(BoutReal vc, stencil &f) {
  BoutReal deriv, w, r;

  if (vc > 0.0) {
    // Left-biased stencil

    r = (WENO_SMALL + SQ(f.c - 2.0 * f.m + f.mm)) /
        (WENO_SMALL + SQ(f.p - 2.0 * f.c + f.m));
    w = 1.0 / (1.0 + 2.0 * r * r);

    deriv = 0.5 * (f.p - f.m) - 0.5 * w * (-f.mm + 3. * f.m - 3. * f.c + f.p);

  } else {
    // Right-biased

    r = (WENO_SMALL + SQ(f.pp - 2.0 * f.p + f.c)) /
        (WENO_SMALL + SQ(f.p - 2.0 * f.c + f.m));
    w = 1.0 / (1.0 + 2.0 * r * r);

    deriv = 0.5 * (f.p - f.m) - 0.5 * w * (-f.m + 3. * f.c - 3. * f.p + f.pp);
  }

  return vc * deriv;
}

/// 3rd-order CWENO. Uses the upwinding code and split flux
BoutReal DDX_CWENO3(stencil &f) {
  BoutReal a, ma = fabs(f.c);
  // Split flux
  a = fabs(f.m);
  if (a > ma)
    ma = a;
  a = fabs(f.p);
  if (a > ma)
    ma = a;
  a = fabs(f.mm);
  if (a > ma)
    ma = a;
  a = fabs(f.pp);
  if (a > ma)
    ma = a;

  stencil sp, vp, sm, vm;

  sp = f + ma;
  sm = ma - f;

  return VDDX_WENO3(0.5, sp) + VDDX_WENO3(-0.5, sm);
}

//////////////////////// FLUX METHODS ///////////////////////

BoutReal FDDX_U1(stencil &v, stencil &f) {
  // Velocity at lower end
  BoutReal vs = 0.5 * (v.m + v.c);
  BoutReal result = (vs >= 0.0) ? vs * f.m : vs * f.c;
  // and at upper
  vs = 0.5 * (v.c + v.p);
  result -= (vs >= 0.0) ? vs * f.c : vs * f.p;

  return - result;
}

BoutReal FDDX_C2(stencil &v, stencil &f) { return 0.5 * (v.p * f.p - v.m * f.m); }

BoutReal FDDX_C4(stencil &v, stencil &f) {
  return (8. * v.p * f.p - 8. * v.m * f.m + v.mm * f.mm - v.pp * f.pp) / 12.;
}

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
 * Staggered differencing methods
 * These expect the output grid cell to be at a different location to the input
 *
 * The stencil no longer has a value in 'C' (centre)
 * instead, points are shifted as follows:
 *
 * mm  -> -3/2 h
 * m   -> -1/2 h
 * p   -> +1/2 h
 * pp  -? +3/2 h
 *
 * NOTE: Cell widths (dx, dy, dz) are currently defined as centre->centre
 * for the methods above. This is currently not taken account of, so large
 * variations in cell size will cause issues.
 *******************************************************************************/

/////////////////////// FIRST DERIVATIVES //////////////////////
// Map Centre -> Low or Low -> Centre

// Second order differencing (staggered)
BoutReal DDX_C2_stag(stencil &f) { return f.p - f.m; }

BoutReal DDX_C4_stag(stencil &f) { return (27. * (f.p - f.m) - (f.pp - f.mm)) / 24.; }

BoutReal D2DX2_C2_stag(stencil &f) { return (f.pp + f.mm - f.p - f.m) / 2.; }
/////////////////////////// UPWINDING ///////////////////////////
// Map (Low, Centre) -> Centre  or (Centre, Low) -> Low
// Hence v contains only (mm, m, p, pp) fields whilst f has 'c' too
//
// v.p is v at +1/2, v.m is at -1/2

BoutReal VDDX_U1_stag(stencil &v, stencil &f) {
  // Lower cell boundary
  BoutReal result = (v.m >= 0) ? v.m * f.m : v.m * f.c;

  // Upper cell boundary
  result -= (v.p >= 0) ? v.p * f.c : v.p * f.p;

  result *= -1;

  // result is now d/dx(v*f), but want v*d/dx(f) so subtract f*d/dx(v)
  result -= f.c * (v.p - v.m);

  return result;
}

BoutReal VDDX_U2_stag(stencil &v, stencil &f) {
  BoutReal result;

  if (v.p > 0 && v.m > 0) {
    // Extrapolate v to centre from below, use 2nd order backward difference on f
    result = (1.5 * v.m - .5 * v.mm) * (.5 * f.mm - 2. * f.m + 1.5 * f.c);
  } else if (v.p < 0 && v.m < 0) {
    // Extrapolate v to centre from above, use 2nd order forward difference on f
    result = (1.5 * v.p - .5 * v.pp) * (-1.5 * f.c + 2. * f.p - .5 * f.pp);
  } else {
    // Velocity changes sign, hence is almost zero: use centred interpolation/differencing
    result = .25 * (v.p + v.m) * (f.p - f.m);
  }

  return result;
}

BoutReal VDDX_C2_stag(stencil &v, stencil &f) {
  // Result is needed at location of f: interpolate v to f's location and take an
  // unstaggered derivative of f
  return 0.5 * (v.p + v.m) * 0.5 * (f.p - f.m);
}

BoutReal VDDX_C4_stag(stencil &v, stencil &f) {
  // Result is needed at location of f: interpolate v to f's location and take an
  // unstaggered derivative of f
  return (9. * (v.m + v.p) - v.mm - v.pp) / 16. * (8. * f.p - 8. * f.m + f.mm - f.pp) /
         12.;
}

/////////////////////////// FLUX ///////////////////////////
// Map (Low, Centre) -> Centre  or (Centre, Low) -> Low
// Hence v contains only (mm, m, p, pp) fields whilst f has 'c' too
//
// v.p is v at +1/2, v.m is at -1/2

BoutReal FDDX_U1_stag(stencil &v, stencil &f) {
  // Lower cell boundary
  BoutReal result = (v.m >= 0) ? v.m * f.m : v.m * f.c;

  // Upper cell boundary
  result -= (v.p >= 0) ? v.p * f.c : v.p * f.p;

  return - result;
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
    {DIFF_DEFAULT, NULL, NULL}}; // Use to terminate the list

/// First derivative lookup table
static DiffLookup FirstDerivTable[] = {
    {DIFF_C2, DDX_C2, NULL, NULL},     {DIFF_W2, DDX_CWENO2, NULL, NULL},
    {DIFF_W3, DDX_CWENO3, NULL, NULL}, {DIFF_C4, DDX_C4, NULL, NULL},
    {DIFF_S2, DDX_S2, NULL, NULL},     {DIFF_FFT, NULL, NULL, NULL},
    {DIFF_DEFAULT, NULL, NULL, NULL}};

/// Second derivative lookup table
static DiffLookup SecondDerivTable[] = {{DIFF_C2, D2DX2_C2, NULL, NULL},
                                        {DIFF_C4, D2DX2_C4, NULL, NULL},
                                        {DIFF_FFT, NULL, NULL, NULL},
                                        {DIFF_DEFAULT, NULL, NULL, NULL}};

/// Upwinding functions lookup table
static DiffLookup UpwindTable[] = {
    {DIFF_U1, NULL, VDDX_U1, NULL},    {DIFF_U2, NULL, VDDX_U2, NULL},
    {DIFF_C2, NULL, VDDX_C2, NULL},    {DIFF_U3, NULL, VDDX_U3, NULL},
    {DIFF_W3, NULL, VDDX_WENO3, NULL}, {DIFF_C4, NULL, VDDX_C4, NULL},
    {DIFF_DEFAULT, NULL, NULL, NULL}};

/// Flux functions lookup table
static DiffLookup FluxTable[] = {
    {DIFF_SPLIT, NULL, NULL, NULL},   {DIFF_U1, NULL, NULL, FDDX_U1},
    {DIFF_C2, NULL, NULL, FDDX_C2},   {DIFF_C4, NULL, NULL, FDDX_C4},
    {DIFF_DEFAULT, NULL, NULL, NULL}};

/// First staggered derivative lookup
static DiffLookup FirstStagDerivTable[] = {{DIFF_C2, DDX_C2_stag, NULL, NULL},
                                           {DIFF_C4, DDX_C4_stag, NULL, NULL},
                                           {DIFF_DEFAULT, NULL, NULL, NULL}};

/// Second staggered derivative lookup
static DiffLookup SecondStagDerivTable[] = {{DIFF_C2, D2DX2_C2_stag, NULL, NULL},
                                            {DIFF_DEFAULT, NULL, NULL, NULL}};

/// Upwinding staggered lookup
static DiffLookup UpwindStagTable[] = {{DIFF_U1, NULL, NULL, VDDX_U1_stag},
                                       {DIFF_U2, NULL, NULL, VDDX_U2_stag},
                                       {DIFF_C2, NULL, NULL, VDDX_C2_stag},
                                       {DIFF_C4, NULL, NULL, VDDX_C4_stag},
                                       {DIFF_DEFAULT, NULL, NULL, NULL}};

/// Flux staggered lookup
static DiffLookup FluxStagTable[] = {{DIFF_SPLIT, NULL, NULL, NULL},
                                     {DIFF_U1, NULL, NULL, FDDX_U1_stag},
                                     {DIFF_DEFAULT, NULL, NULL, NULL}};

/*******************************************************************************
 * Routines to use the above tables to map between function codes, names
 * and pointers
 *******************************************************************************/

Mesh::deriv_func lookupFunc(DiffLookup *table, DIFF_METHOD method) {
  int i = 0;
  do {
    if (table[i].method == method)
      return table[i].func;
    i++;
  } while (table[i].method != DIFF_DEFAULT);
  // Not found in list. Return the first

  return table[0].func;
}

Mesh::upwind_func lookupUpwindFunc(DiffLookup *table, DIFF_METHOD method) {
  int i = 0;
  do {
    if (table[i].method == method)
      return table[i].up_func;
    i++;
  } while (table[i].method != DIFF_DEFAULT);
  // Not found in list. Return the first

  return table[0].up_func;
}

Mesh::flux_func lookupFluxFunc(DiffLookup *table, DIFF_METHOD method) {
  SCOREP0();
  int i = 0;
  do {
    if (table[i].method == method)
      return table[i].fl_func;
    i++;
  } while (table[i].method != DIFF_DEFAULT);
  // Not found in list. Return the first

  return table[0].fl_func;
}

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

/// This function is used during initialisation only (i.e. doesn't need to be particularly
/// fast) Returns DIFF_METHOD, rather than function so can be applied to central and
/// upwind tables
DIFF_METHOD lookupFunc(DiffLookup *table, const string &label) {

  if (label.empty())
    return table[0].method;

  // Loop through the name lookup table
  for (int i = 0; DiffNameTable[i].method != DIFF_DEFAULT; ++i) {
    if (strcasecmp(label.c_str(), DiffNameTable[i].label) == 0) { // Whole match
      if (isImplemented(table, DiffNameTable[i].method)) {
        return DiffNameTable[i].method;
      } else {
        std::string avail{};

        for (int i = 0; DiffNameTable[i].method != DIFF_DEFAULT; ++i) {
          if (isImplemented(table, DiffNameTable[i].method)) {
            avail += DiffNameTable[i].label;
            avail += "\n";
          }
        }
        throw BoutException("Option %s is known but not valid for this differencing "
                            "type.\nAvailable options are:\n%s",
                            label.c_str(), avail.c_str());
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
void derivs_set(Options *options, DiffLookup *table, const char *name,
                Mesh::deriv_func &f) {
  TRACE("derivs_set( deriv_func )");
  string label;
  options->get(name, label, "C2");

  DIFF_METHOD method = lookupFunc(table, label); // Find the function
  printFuncName(method);                         // Print differential function name
  f = lookupFunc(table, method);                 // Find the function pointer
}

void derivs_set(Options *options, DiffLookup *table, const char *name,
                Mesh::upwind_func &f) {
  TRACE("derivs_set( upwind_func )");
  string label;
  options->get(name, label, "U1");

  DIFF_METHOD method = lookupFunc(table, label); // Find the function
  printFuncName(method);                         // Print differential function name
  f = lookupUpwindFunc(table, method);
}

void derivs_set(Options *options, DiffLookup *table, const char *name,
                Mesh::flux_func &f) {
  TRACE("derivs_set( flux_func )");
  string label;
  options->get(name, label, "U1");

  DIFF_METHOD method = lookupFunc(table, label); // Find the function
  printFuncName(method);                         // Print differential function name
  f = lookupFluxFunc(table, method);
}

/// Initialise derivatives from options
void derivs_initialise(Options *options, bool StaggerGrids, Mesh::deriv_func &fdd,
                       Mesh::deriv_func &sfdd, Mesh::deriv_func &fd2d,
                       Mesh::deriv_func &sfd2d, Mesh::upwind_func &fu,
                       Mesh::flux_func &sfu, Mesh::flux_func &ff, Mesh::flux_func &sff) {
  output_info.write("\tFirst       : ");
  derivs_set(options, FirstDerivTable, "first", fdd);
  if (StaggerGrids) {
    output_info.write("\tStag. First : ");
    derivs_set(options, FirstStagDerivTable, "first", sfdd);
  }
  output_info.write("\tSecond      : ");
  derivs_set(options, SecondDerivTable, "second", fd2d);
  if (StaggerGrids) {
    output_info.write("\tStag. Second: ");
    derivs_set(options, SecondStagDerivTable, "second", sfd2d);
  }
  output_info.write("\tUpwind      : ");
  derivs_set(options, UpwindTable, "upwind", fu);
  if (StaggerGrids) {
    output_info.write("\tStag. Upwind: ");
    derivs_set(options, UpwindStagTable, "upwind", sfu);
  }
  output_info.write("\tFlux        : ");
  derivs_set(options, FluxTable, "flux", ff);
  if (StaggerGrids) {
    output_info.write("\tStag. Flux  : ");
    derivs_set(options, FluxStagTable, "flux", sff);
  }
}

/// Initialise the derivative methods. Must be called before any derivatives are used
void Mesh::derivs_init(Options *options) {
  TRACE("Initialising derivatives");

  output_info.write("Setting X differencing methods\n");
  derivs_initialise(options->getSection("ddx"), StaggerGrids, fDDX, sfDDX, fD2DX2,
                    sfD2DX2, fVDDX, sfVDDX, fFDDX, sfFDDX);

  if ((fDDX == NULL) || (fD2DX2 == NULL))
    throw BoutException("FFT cannot be used in X\n");

  output_info.write("Setting Y differencing methods\n");
  derivs_initialise(options->getSection("ddy"), StaggerGrids, fDDY, sfDDY, fD2DY2,
                    sfD2DY2, fVDDY, sfVDDY, fFDDY, sfFDDY);

  if ((fDDY == NULL) || (fD2DY2 == NULL))
    throw BoutException("FFT cannot be used in Y\n");

  output_info.write("Setting Z differencing methods\n");
  derivs_initialise(options->getSection("ddz"), StaggerGrids, fDDZ, sfDDZ, fD2DZ2,
                    sfD2DZ2, fVDDZ, sfVDDZ, fFDDZ, sfFDDZ);

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

const Field2D Mesh::applyXdiff(const Field2D &var, Mesh::deriv_func func,
                               CELL_LOC loc, REGION region) {
  ASSERT1(this == var.getMesh());
  ASSERT1(var.isAllocated());

  if (var.getNx() == 1) {
    return Field2D(0., this);
  }

  CELL_LOC diffloc = var.getLocation();

  Field2D result(this);
  result.allocate(); // Make sure data allocated

  if (this->StaggerGrids && (loc != CELL_DEFAULT) && (loc != var.getLocation())) {
    // Staggered differencing

    CELL_LOC location = var.getLocation();

    if (this->xstart > 1) {
      // More than one guard cell, so set pp and mm values
      // This allows higher-order methods to be used
BOUT_OMP(parallel)
{
      stencil s;
      IndexOffset<Ind3D> offset(*mesh);
      BLOCK_REGION_LOOP_PARALLEL_SECTION(mesh->getRegion3D(region), i,
        s.mm = var[offset.xmm(i)];
        s.m = var[offset.xm(i)];
        s.c = var[i];
        s.p = var[offset.xp(i)];
        s.pp = var[offset.xpp(i)];

        if ((location == CELL_CENTRE) && (loc == CELL_XLOW)) {
          // Producing a stencil centred around a lower X value
          s.pp = s.p;
          s.p = s.c;
        } else if (location == CELL_XLOW) {
          // Stencil centred around a cell centre
          s.mm = s.m;
          s.m = s.c;
        }

        result[i] = func(s);
      );
}
    } else {
      // Only one guard cell, so no pp or mm values
BOUT_OMP(parallel)
{
      stencil s;
      IndexOffset<Ind3D> offset(*mesh);
      BLOCK_REGION_LOOP_PARALLEL_SECTION(mesh->getRegion3D(region), i,
        s.mm = nan("");
        s.m = var[offset.xm(i)];
        s.c = var[i];
        s.p = var[offset.xp(i)];
        s.pp = nan("");

        if ((location == CELL_CENTRE) && (loc == CELL_XLOW)) {
          // Producing a stencil centred around a lower X value
          s.pp = s.p;
          s.p = s.c;
        } else if (location == CELL_XLOW) {
          // Stencil centred around a cell centre
          s.mm = s.m;
          s.m = s.c;
        }

        result[i] = func(s);
      );
}
    }

  } else {
    // Non-staggered differencing

    if (this->xstart > 1) {
      // More than one guard cell, so set pp and mm values
      // This allows higher-order methods to be used
BOUT_OMP(parallel)
{
      stencil s;
      IndexOffset<Ind3D> offset(*mesh);
      BLOCK_REGION_LOOP_PARALLEL_SECTION(mesh->getRegion3D(region), i,
        s.mm = var[offset.xmm(i)];
        s.m = var[offset.xm(i)];
        s.c = var[i];
        s.p = var[offset.xp(i)];
        s.pp = var[offset.xpp(i)];

        result[i] = func(s);
      );
}
    } else {
      // Only one guard cell, so no pp or mm values
BOUT_OMP(parallel)
{
      stencil s;
      IndexOffset<Ind3D> offset(*mesh);
      BLOCK_REGION_LOOP_PARALLEL_SECTION(mesh->getRegion3D(region), i,
        s.mm = nan("");
        s.m = var[offset.xm(i)];
        s.c = var[i];
        s.p = var[offset.xp(i)];
        s.pp = nan("");

        result[i] = func(s);
      );
}
    }
  }

  result.setLocation(diffloc);

#if CHECK > 0
  // Mark boundaries as invalid
  result.bndry_xin = result.bndry_xout = result.bndry_yup = result.bndry_ydown = false;
#endif

  return result;
}

const Field3D Mesh::applyXdiff(const Field3D &var, Mesh::deriv_func func,
                               CELL_LOC loc, REGION region) {
  // Check that the mesh is correct
  ASSERT1(this == var.getMesh());
  // Check that the input variable has data
  ASSERT1(var.isAllocated());

  if (var.getNx() == 1) {
    return Field3D(0., this);
  }

  CELL_LOC diffloc = var.getLocation();

  Field3D result(this);
  result.allocate(); // Make sure data allocated

  if (this->StaggerGrids && (loc != CELL_DEFAULT) && (loc != var.getLocation())) {
    // Staggered differencing

    CELL_LOC location = var.getLocation();

    if (this->xstart > 1) {
      // More than one guard cell, so set pp and mm values
      // This allows higher-order methods to be used
BOUT_OMP(parallel)
{
      stencil s;
      IndexOffset<Ind3D> offset(*mesh);
      //for (const auto &i : result.region(region)) {
      BLOCK_REGION_LOOP_PARALLEL_SECTION(mesh->getRegion3D(region), i,
        s.mm = var[offset.xmm(i)];
        s.m = var[offset.xm(i)];
        s.c = var[i];
        s.p = var[offset.xp(i)];
        s.pp = var[offset.xpp(i)];

        if ((location == CELL_CENTRE) && (loc == CELL_XLOW)) {
          // Producing a stencil centred around a lower X value
          s.pp = s.p;
          s.p = s.c;
        } else if (location == CELL_XLOW) {
          // Stencil centred around a cell centre
          s.mm = s.m;
          s.m = s.c;
        }

        result[i] = func(s);
      );
}
    } else {
      // Only one guard cell, so no pp or mm values
BOUT_OMP(parallel)
{
      stencil s;
      IndexOffset<Ind3D> offset(*mesh);
      //for (const auto &i : result.region(region)) {
      BLOCK_REGION_LOOP_PARALLEL_SECTION(mesh->getRegion3D(region), i,
        s.mm = nan("");
        s.m = var[offset.xm(i)];
        s.c = var[i];
        s.p = var[offset.xp(i)];
        s.pp = nan("");

        if ((location == CELL_CENTRE) && (loc == CELL_XLOW)) {
          // Producing a stencil centred around a lower X value
          s.pp = s.p;
          s.p = s.c;
        } else if (location == CELL_XLOW) {
          // Stencil centred around a cell centre
          s.mm = s.m;
          s.m = s.c;
        }

        result[i] = func(s);
      );
}
    }

  } else {
    // Non-staggered differencing

    if (this->xstart > 1) {
      // More than one guard cell, so set pp and mm values
      // This allows higher-order methods to be used
BOUT_OMP(parallel)
{
      stencil s;
      IndexOffset<Ind3D> offset(*mesh);
      //for (const auto &i : result.region(region)) {
      BLOCK_REGION_LOOP_PARALLEL_SECTION(mesh->getRegion3D(region), i,
        s.mm = var[offset.xmm(i)];
        s.m = var[offset.xm(i)];
        s.c = var[i];
        s.p = var[offset.xp(i)];
        s.pp = var[offset.xpp(i)];

        result[i] = func(s);
      );
}
    } else {
      // Only one guard cell, so no pp or mm values
BOUT_OMP(parallel)
{
      stencil s;
      IndexOffset<Ind3D> offset(*mesh);
      s.mm = nan("");
      s.pp = nan("");
      //for (const auto &i : result.region(region)) {
      BLOCK_REGION_LOOP_PARALLEL_SECTION(mesh->getRegion3D(region), i,
        s.m = var[offset.xm(i)];
        s.c = var[i];
        s.p = var[offset.xp(i)];

        result[i] = func(s);
      );
}
    }
  }

  result.setLocation(diffloc);

#if CHECK > 0
  // Mark boundaries as invalid
  result.bndry_xin = result.bndry_xout = result.bndry_yup = result.bndry_ydown = false;
#endif

  return result;
}

// Y derivative

const Field2D Mesh::applyYdiff(const Field2D &var, Mesh::deriv_func func, CELL_LOC UNUSED(loc),
                               REGION region) {
SCOREP0();
  ASSERT1(this == var.getMesh());
  // Check that the input variable has data
  ASSERT1(var.isAllocated());

  if (var.getNy() == 1) {
    return Field2D(0., this);
  }

  CELL_LOC diffloc = var.getLocation();

  Field2D result(this);
  result.allocate(); // Make sure data allocated
BOUT_OMP(parallel)
{

  if (this->ystart > 1) {
    // More than one guard cell, so set pp and mm values
    // This allows higher-order methods to be used

    stencil s;
    IndexOffset<Ind3D> offset(*mesh);
    BLOCK_REGION_LOOP_PARALLEL_SECTION(mesh->getRegion3D(region), i,
      // Set stencils
      s.mm = var[offset.ymm(i)];
      s.m = var[offset.ym(i)];
      s.c = var[i];
      s.p = var[offset.yp(i)];
      s.pp = var[offset.ypp(i)];

      result[i] = func(s);
    );
  } else {
    // Only one guard cell, so no pp or mm values
    stencil s;
    IndexOffset<Ind3D> offset(*mesh);
    s.mm = nan("");
    s.pp = nan("");
    BLOCK_REGION_LOOP_PARALLEL_SECTION(mesh->getRegion3D(region), i,
      // Set stencils
      s.m = var[offset.ym(i)];
      s.c = var[i];
      s.p = var[offset.yp(i)];

      result[i] = func(s);
    );
  }

  result.setLocation(diffloc);

#if CHECK > 0
  // Mark boundaries as invalid
  result.bndry_yup = result.bndry_ydown = false;
#endif

}
  return result;
}

const Field3D Mesh::applyYdiff(const Field3D &var, Mesh::deriv_func func, CELL_LOC loc,
                               REGION region) {
  SCOREP0();
  ASSERT1(this == var.getMesh());
  // Check that the input variable has data
  ASSERT1(var.isAllocated());

  if (var.getNy() == 1) {
    return Field3D(0., this);
  }

  CELL_LOC diffloc = var.getLocation();

  Field3D result(this);
  result.allocate(); // Make sure data allocated

  if (var.hasYupYdown() && ((&var.yup() != &var) || (&var.ydown() != &var))) {
    // Field "var" has distinct yup and ydown fields which
    // will be used to calculate a derivative along
    // the magnetic field

    if (this->StaggerGrids && (loc != CELL_DEFAULT) && (loc != var.getLocation())) {
      // Staggered differencing

      // Cell location of the input field
      CELL_LOC location = var.getLocation();

BOUT_OMP(parallel)
{
      stencil s;
      IndexOffset<Ind3D> offset(*mesh);
      //for (const auto &i : result.region(region)) {
      BLOCK_REGION_LOOP_PARALLEL_SECTION(mesh->getRegion3D(region), i,
        s.mm = nan("");
        s.m = var.ydown()[offset.ym(i)];
        s.c = var[i];
        s.p = var.yup()[offset.yp(i)];
        s.pp = nan("");

        if ((location == CELL_CENTRE) && (loc == CELL_YLOW)) {
          // Producing a stencil centred around a lower Y value
          s.pp = s.p;
          s.p = s.c;
        } else if (location == CELL_YLOW) {
          // Stencil centred around a cell centre
          s.mm = s.m;
          s.m = s.c;
        }

        result[i] = func(s);
);
}
    } else {
      // Non-staggered
BOUT_OMP(parallel)
{
      stencil s;
      IndexOffset<Ind3D> offset(*mesh);
      //for (const auto &i : result.region(region)) {
      BLOCK_REGION_LOOP_PARALLEL_SECTION(mesh->getRegion3D(region), i,
        s.mm = nan("");
        s.m = var.ydown()[offset.ym(i)];
        s.c = var[i];
        s.p = var.yup()[offset.yp(i)];
        s.pp = nan("");

        result[i] = func(s);
);
}
    }
  } else {
    // var has no yup/ydown fields, so we need to shift into field-aligned coordinates

    Field3D var_fa = this->toFieldAligned(var);

    if (this->StaggerGrids && (loc != CELL_DEFAULT) && (loc != var.getLocation())) {
      // Staggered differencing

      // Cell location of the input field
      CELL_LOC location = var.getLocation();

      if (this->ystart > 1) {
BOUT_OMP(parallel)
{
        stencil s;
        IndexOffset<Ind3D> offset(*mesh);
        // More than one guard cell, so set pp and mm values
        // This allows higher-order methods to be used
        //for (const auto &i : result.region(region)) {
        BLOCK_REGION_LOOP_PARALLEL_SECTION(mesh->getRegion3D(region), i,
          s.mm = var_fa[offset.ymm(i)];
          s.m = var_fa[offset.ym(i)];
          s.c = var_fa[i];
          s.p = var_fa[offset.yp(i)];
          s.pp = var_fa[offset.ypp(i)];

          if ((location == CELL_CENTRE) && (loc == CELL_YLOW)) {
            // Producing a stencil centred around a lower Y value
            s.pp = s.p;
            s.p = s.c;
          } else if (location == CELL_YLOW) {
            // Stencil centred around a cell centre
            s.mm = s.m;
            s.m = s.c;
          }

          result[i] = func(s);
);
}
      } else {
BOUT_OMP(parallel)
{
        stencil s;
        IndexOffset<Ind3D> offset(*mesh);
        //for (const auto &i : result.region(region)) {
        BLOCK_REGION_LOOP_PARALLEL_SECTION(mesh->getRegion3D(region), i,
          s.mm = nan("");
          s.m = var_fa[offset.ym(i)];
          s.c = var_fa[i];
          s.p = var_fa[offset.yp(i)];
          s.pp = nan("");

          if ((location == CELL_CENTRE) && (loc == CELL_YLOW)) {
            // Producing a stencil centred around a lower Y value
            s.pp = s.p;
            s.p = s.c;
          } else if (location == CELL_YLOW) {
            // Stencil centred around a cell centre
            s.mm = s.m;
            s.m = s.c;
          }

          result[i] = func(s);
);
}
      }

    } else {
      // Non-staggered differencing

      if (this->ystart > 1) {
BOUT_OMP(parallel)
{
        stencil s;
        IndexOffset<Ind3D> offset(*mesh);
        // More than one guard cell, so set pp and mm values
        // This allows higher-order methods to be used
        //for (const auto &i : result.region(region)) {
        BLOCK_REGION_LOOP_PARALLEL_SECTION(mesh->getRegion3D(region), i,
          s.mm = var_fa[offset.ymm(i)];
          s.m = var_fa[offset.ym(i)];
          s.c = var_fa[i];
          s.p = var_fa[offset.yp(i)];
          s.pp = var_fa[offset.ypp(i)];

          result[i] = func(s);
);
}
      } else {
BOUT_OMP(parallel)
{
        stencil s;
        IndexOffset<Ind3D> offset(*mesh);
        // Only one guard cell, so no pp or mm values
        //for (const auto &i : result.region(region)) {
        BLOCK_REGION_LOOP_PARALLEL_SECTION(mesh->getRegion3D(region), i,
          s.mm = nan("");
          s.m = var_fa[offset.ym(i)];
          s.c = var_fa[i];
          s.p = var_fa[offset.yp(i)];
          s.pp = nan("");

          result[i] = func(s);
);
}
      }
    }

    // Shift result back

    result = this->fromFieldAligned(result);
  }

  result.setLocation(diffloc);

#if CHECK > 0
  // Mark boundaries as invalid
  result.bndry_xin = result.bndry_xout = result.bndry_yup = result.bndry_ydown = false;
#endif

  return result;
}

// Z derivative

const Field3D Mesh::applyZdiff(const Field3D &var, Mesh::deriv_func func, CELL_LOC loc,
                               REGION region) {
  SCOREP0();
  ASSERT1(this == var.getMesh());
  // Check that the input variable has data
  ASSERT1(var.isAllocated());

  if (var.getNz() == 1) {
    return Field3D(0., this);
  }


  CELL_LOC diffloc = var.getLocation();

  if (this->StaggerGrids && (loc != CELL_DEFAULT) && (loc != var.getLocation())) {
    // Staggered differencing
    throw BoutException("No one used this before. And no one implemented it.");
  }

  Field3D result(this);
  result.allocate(); // Make sure data allocated

  // Check that the input variable has data
  ASSERT1(var.isAllocated());

///  BOUT_OMP(parallel)
///  {
///    stencil s;
///    IndexOffset<Ind3D> offset(*mesh);
///    //for (const auto &i : result.region(region)) {
///    BLOCK_REGION_LOOP_PARALLEL_SECTION( mesh->getRegion3D(region), i,
///      s.mm = var[offset.zmm(i)];
///      s.m = var[offset.zm(i)];
///      s.c = var[i];
///      s.p = var[offset.zp(i)];
///      s.pp = var[offset.zpp(i)];
///
///      result[i] = func(s);
///    );
///  }

  stencil s;
  for (const auto &i : result.region(region)) {
    s.c = var[i];
    s.p = var[i.zp()];
    s.m = var[i.zm()];
    s.pp = var[i.offset(0, 0, 2)];
    s.mm = var[i.offset(0, 0, -2)];

    result[i] = func(s);
  }

  result.setLocation(diffloc);

  return result;
}

/*******************************************************************************
 * First central derivatives
 *******************************************************************************/

////////////// X DERIVATIVE /////////////////

const Field3D Mesh::indexDDX(const Field3D &f, CELL_LOC outloc, DIFF_METHOD method, REGION region) {

  Mesh::deriv_func func = fDDX; // Set to default function
  DiffLookup *table = FirstDerivTable;

  CELL_LOC inloc = f.getLocation(); // Input location
  CELL_LOC diffloc = inloc;         // Location of differential result

  Field3D result(this);

  if (this->StaggerGrids && (outloc == CELL_DEFAULT)) {
    // Take care of CELL_DEFAULT case
    outloc = diffloc; // No shift (i.e. same as no stagger case)
  }

  if (this->StaggerGrids && (outloc != inloc)) {
    // Shifting to a new location

    if (((inloc == CELL_CENTRE) && (outloc == CELL_XLOW)) ||
        ((inloc == CELL_XLOW) && (outloc == CELL_CENTRE))) {
      // Shifting in X. Centre -> Xlow, or Xlow -> Centre

      func = sfDDX;                // Set default
      table = FirstStagDerivTable; // Set table for others
      diffloc = (inloc == CELL_CENTRE) ? CELL_XLOW : CELL_CENTRE;

    } else {
      // Derivative of interpolated field or interpolation of derivative field
      // cannot be taken without communicating and applying boundary
      // conditions, so throw an exception instead
      throw BoutException("Unsupported combination of {inloc =%s} and {outloc =%s} in Mesh:indexDDX(Field3D).", strLocation(inloc), strLocation(outloc));
    }
  }

  if (method != DIFF_DEFAULT) {
    // Lookup function
    func = lookupFunc(table, method);
    if (func == NULL)
      throw BoutException("Cannot use FFT for X derivatives");
  }

  result = applyXdiff(f, func, diffloc, region);

  result.setLocation(diffloc); // Set the result location

  return result;
}

const Field2D Mesh::indexDDX(const Field2D &f, CELL_LOC outloc,
                             DIFF_METHOD method, REGION region) {
  ASSERT1(outloc == CELL_DEFAULT || outloc == f.getLocation());
  ASSERT1(method == DIFF_DEFAULT);
  return applyXdiff(f, fDDX, f.getLocation(), region);
}

////////////// Y DERIVATIVE /////////////////

const Field3D Mesh::indexDDY(const Field3D &f, CELL_LOC outloc,
                             DIFF_METHOD method, REGION region) {
  Mesh::deriv_func func = fDDY; // Set to default function
  DiffLookup *table = FirstDerivTable;

  CELL_LOC inloc = f.getLocation(); // Input location
  CELL_LOC diffloc = inloc;         // Location of differential result

  Field3D result(this);

  if (this->StaggerGrids && (outloc == CELL_DEFAULT)) {
    // Take care of CELL_DEFAULT case
    outloc = diffloc; // No shift (i.e. same as no stagger case)
  }

  if (this->StaggerGrids && (outloc != inloc)) {
    // Shifting to a new location
    if (((inloc == CELL_CENTRE) && (outloc == CELL_YLOW)) ||
        ((inloc == CELL_YLOW) && (outloc == CELL_CENTRE))) {
      // Shifting in Y. Centre -> Ylow, or Ylow -> Centre

      func = sfDDY;                // Set default
      table = FirstStagDerivTable; // Set table for others
      diffloc = (inloc == CELL_CENTRE) ? CELL_YLOW : CELL_CENTRE;

    } else {
      // Derivative of interpolated field or interpolation of derivative field
      // cannot be taken without communicating and applying boundary
      // conditions, so throw an exception instead
      throw BoutException("Unsupported combination of {inloc =%s} and {outloc =%s} in Mesh:indexDDY(Field3D).", strLocation(inloc), strLocation(outloc));
    }
  }

  if (method != DIFF_DEFAULT) {
    // Lookup function
    func = lookupFunc(table, method);
    if (func == NULL)
      throw BoutException("Cannot use FFT for Y derivatives");
  }

  result = applyYdiff(f, func, diffloc, region);

  result.setLocation(diffloc); // Set the result location

  return result;
}

const Field2D Mesh::indexDDY(const Field2D &f, CELL_LOC outloc,
                             DIFF_METHOD method, REGION region) {
  ASSERT1(outloc == CELL_DEFAULT || outloc == f.getLocation());
  ASSERT1(method == DIFF_DEFAULT);
  return applyYdiff(f, fDDY, f.getLocation(), region);
}

////////////// Z DERIVATIVE /////////////////

const Field3D Mesh::indexDDZ(const Field3D &f, CELL_LOC outloc,
                             DIFF_METHOD method, REGION region) {
  SCOREP0();
  Mesh::deriv_func func = fDDZ; // Set to default function
  DiffLookup *table = FirstDerivTable;

  CELL_LOC inloc = f.getLocation(); // Input location
  CELL_LOC diffloc = inloc;         // Location of differential result

  Field3D result(this);

  if (this->StaggerGrids && (outloc == CELL_DEFAULT)) {
    // Take care of CELL_DEFAULT case
    outloc = diffloc; // No shift (i.e. same as no stagger case)
  }

  if (this->StaggerGrids && (outloc != inloc)) {
    // Shifting to a new location

    if (((inloc == CELL_CENTRE) && (outloc == CELL_ZLOW)) ||
        ((inloc == CELL_ZLOW) && (outloc == CELL_CENTRE))) {
      // Shifting in Z. Centre -> Zlow, or Zlow -> Centre

      func = sfDDZ;                // Set default
      table = FirstStagDerivTable; // Set table for others
      diffloc = (inloc == CELL_CENTRE) ? CELL_ZLOW : CELL_CENTRE;

    } else {
      // Derivative of interpolated field or interpolation of derivative field
      // cannot be taken without communicating and applying boundary
      // conditions, so throw an exception instead
      throw BoutException("Unsupported combination of {inloc =%s} and {outloc =%s} in Mesh:indexDDZ(Field3D).", strLocation(inloc), strLocation(outloc));
    }
  }

  if (method != DIFF_DEFAULT) {
    // Lookup function
    func = lookupFunc(table, method);
  }

  if (func == NULL) {
    // Use FFT

    BoutReal shift = 0.; // Shifting result in Z?
    if (this->StaggerGrids) {
      if ((inloc == CELL_CENTRE) && (diffloc == CELL_ZLOW)) {
        // Shifting down - multiply by exp(-0.5*i*k*dz)
        shift = -1.;
        throw BoutException("Not tested - probably broken");
      } else if ((inloc == CELL_ZLOW) && (diffloc == CELL_CENTRE)) {
        // Shifting up
        shift = 1.;
        throw BoutException("Not tested - probably broken");
      }
    }

    result.allocate(); // Make sure data allocated

    auto region_index = f.region(region);
    int xs = region_index.xstart;
    int xe = region_index.xend;
    int ys = region_index.ystart;
    int ye = region_index.yend;
    ASSERT2(region_index.zstart == 0);
    int ncz = region_index.zend + 1;

    BOUT_OMP(parallel)
    {
      Array<dcomplex> cv(ncz / 2 + 1);


      // Calculate how many Z wavenumbers will be removed
      int kfilter =
          static_cast<int>(fft_derivs_filter * ncz / 2); // truncates, rounding down
      if (kfilter < 0)
        kfilter = 0;
      if (kfilter > (ncz / 2))
        kfilter = ncz / 2;
      int kmax = ncz / 2 - kfilter; // Up to and including this wavenumber index

      BOUT_OMP(for)
      for (int jx = xs; jx <= xe; jx++) {
        for (int jy = ys; jy <= ye; jy++) {
          rfft(f(jx, jy), ncz, cv.begin()); // Forward FFT

          for (int jz = 0; jz <= kmax; jz++) {
            BoutReal kwave = jz * 2.0 * PI / ncz; // wave number is 1/[rad]

            cv[jz] *= dcomplex(0.0, kwave);
            if (shift)
              cv[jz] *= exp(Im * (shift * kwave));
          }
          for (int jz = kmax + 1; jz <= ncz / 2; jz++) {
            cv[jz] = 0.0;
          }

          irfft(cv.begin(), ncz, result(jx, jy)); // Reverse FFT
        }
      }
    }
      // End of parallel section

#if CHECK > 0
    // Mark boundaries as invalid
    result.bndry_xin = false;
    result.bndry_xout = false;
    result.bndry_yup = false;
    result.bndry_ydown = false;
#endif

  } else {
    // All other (non-FFT) functions
    result = applyZdiff(f, func, diffloc, region);
  }

  result.setLocation(diffloc);

  return result;
}

const Field2D Mesh::indexDDZ(const Field2D &f, CELL_LOC outloc,
                             DIFF_METHOD method, REGION region) {
  SCOREP0();
  ASSERT1(this == f.getMesh());
  return Field2D(0., this);
}

/*******************************************************************************
 * 2nd derivatives
 *******************************************************************************/

////////////// X DERIVATIVE /////////////////

/*!
 * @brief Calculates second X derivative on Mesh in index space
 *
 * @param[in] f        3D scalar field to be differentiated.
 *                     Must be allocated and finite
 *
 * @param[in] outloc   The cell location of the result
 *
 * @param[in] method   The numerical method to use
 *
 * @return  A 3D scalar field with invalid data in the
 *          guard cells
 *
 */
const Field3D Mesh::indexD2DX2(const Field3D &f, CELL_LOC outloc,
                               DIFF_METHOD method, REGION region) {
  Mesh::deriv_func func = fD2DX2; // Set to default function
  DiffLookup *table = SecondDerivTable;

  CELL_LOC inloc = f.getLocation(); // Input location
  CELL_LOC diffloc = inloc;         // Location of differential result

  ASSERT1(this == f.getMesh());

  Field3D result(this);

  if (StaggerGrids && (outloc == CELL_DEFAULT)) {
    // Take care of CELL_DEFAULT case
    outloc = diffloc; // No shift (i.e. same as no stagger case)
  }

  if (StaggerGrids && (outloc != inloc)) {
    // Shifting to a new location

    if (((inloc == CELL_CENTRE) && (outloc == CELL_XLOW)) ||
        ((inloc == CELL_XLOW) && (outloc == CELL_CENTRE))) {
      // Shifting in X. Centre -> Xlow, or Xlow -> Centre

      func = sfD2DX2;               // Set default
      table = SecondStagDerivTable; // Set table for others
      diffloc = (inloc == CELL_CENTRE) ? CELL_XLOW : CELL_CENTRE;

    } else {
      // Derivative of interpolated field or interpolation of derivative field
      // cannot be taken without communicating and applying boundary
      // conditions, so throw an exception instead
      throw BoutException("Unsupported combination of {inloc =%s} and {outloc =%s} in Mesh:indexD2DX2(Field3D).", strLocation(inloc), strLocation(outloc));
    }
  }

  if (method != DIFF_DEFAULT) {
    // Lookup function
    func = lookupFunc(table, method);
    if (func == NULL)
      throw BoutException("Cannot use FFT for X derivatives");
  }

  result = applyXdiff(f, func, diffloc, region);

  result.setLocation(diffloc);

  return result;
}

/*!
 * @brief Calculates second X derivative on Mesh in index space
 *
 * @param[in] f        2D scalar field to be differentiated.
 *                     Must be allocated and finite
 *
 * @return  A 2D scalar field with invalid data in the
 *          guard cells
 *
 */
const Field2D Mesh::indexD2DX2(const Field2D &f,  CELL_LOC outloc,
                               DIFF_METHOD method, REGION region) {
  ASSERT1(outloc == CELL_DEFAULT || outloc == f.getLocation());
  ASSERT1(method == DIFF_DEFAULT);
  return applyXdiff(f, fD2DX2, f.getLocation(), region);
}

////////////// Y DERIVATIVE /////////////////

/*!
 * @brief Calculates second Y derivative on Mesh in index space
 *
 * @param[in] f        3D scalar field to be differentiated.
 *                     Must be allocated and finite
 *
 * @return  A 3D scalar field with invalid data in the
 *          guard cells
 *
 */
const Field3D Mesh::indexD2DY2(const Field3D &f, CELL_LOC outloc,
                               DIFF_METHOD method, REGION region) {
  Mesh::deriv_func func = fD2DY2; // Set to default function
  DiffLookup *table = SecondDerivTable;

  CELL_LOC inloc = f.getLocation(); // Input location
  CELL_LOC diffloc = inloc;         // Location of differential result

  ASSERT1(this == f.getMesh());

  Field3D result(this);

  if (StaggerGrids && (outloc == CELL_DEFAULT)) {
    // Take care of CELL_DEFAULT case
    outloc = diffloc; // No shift (i.e. same as no stagger case)
  }

  if (StaggerGrids && (outloc != inloc)) {
    // Shifting to a new location

    if (((inloc == CELL_CENTRE) && (outloc == CELL_YLOW)) ||
        ((inloc == CELL_YLOW) && (outloc == CELL_CENTRE))) {
      // Shifting in Y. Centre -> Ylow, or Ylow -> Centre

      func = sfD2DY2;               // Set default
      table = SecondStagDerivTable; // Set table for others
      diffloc = (inloc == CELL_CENTRE) ? CELL_YLOW : CELL_CENTRE;

    } else {
      // Derivative of interpolated field or interpolation of derivative field
      // cannot be taken without communicating and applying boundary
      // conditions, so throw an exception instead
      throw BoutException("Unsupported combination of {inloc =%s} and {outloc =%s} in Mesh:indexD2DY2(Field3D).", strLocation(inloc), strLocation(outloc));
    }
  }

  if (method != DIFF_DEFAULT) {
    // Lookup function
    func = lookupFunc(table, method);
    if (func == NULL)
      throw BoutException("Cannot use FFT for Y derivatives");
  }

  result = applyYdiff(f, func, diffloc, region);

  result.setLocation(diffloc);

  return result;
}

/*!
 * @brief Calculates second Y derivative on Mesh in index space
 *
 * @param[in] f        2D scalar field to be differentiated.
 *                     Must be allocated and finite
 *
 * @return  A 2D scalar field with invalid data in the
 *          guard cells
 *
 */
const Field2D Mesh::indexD2DY2(const Field2D &f, CELL_LOC outloc,
                               DIFF_METHOD method, REGION region) {
  ASSERT1(outloc == CELL_DEFAULT || outloc == f.getLocation());
  ASSERT1(method == DIFF_DEFAULT);
  return applyYdiff(f, fD2DY2, f.getLocation(), region);
}

////////////// Z DERIVATIVE /////////////////

/*!
 * @brief Calculates second Z derivative on Mesh in index space
 *
 * @param[in] f        3D scalar field to be differentiated.
 *                     Must be allocated and finite
 *
 * @return  A 3D scalar field with invalid data in the
 *          guard cells
 *
 */
const Field3D Mesh::indexD2DZ2(const Field3D &f, CELL_LOC outloc,
                               DIFF_METHOD method, REGION region) {
  Mesh::deriv_func func = fD2DZ2; // Set to default function
  DiffLookup *table = SecondDerivTable;

  CELL_LOC inloc = f.getLocation(); // Input location
  CELL_LOC diffloc = inloc;         // Location of differential result

  ASSERT1(this == f.getMesh());

  Field3D result(this);

  if (StaggerGrids && (outloc == CELL_DEFAULT)) {
    // Take care of CELL_DEFAULT case
    outloc = diffloc; // No shift (i.e. same as no stagger case)
  }

  if (StaggerGrids && (outloc != inloc)) {
    // Shifting to a new location

    if (((inloc == CELL_CENTRE) && (outloc == CELL_ZLOW)) ||
        ((inloc == CELL_ZLOW) && (outloc == CELL_CENTRE))) {
      // Shifting in Z. Centre -> Zlow, or Zlow -> Centre

      func = sfD2DZ2;               // Set default
      table = SecondStagDerivTable; // Set table for others
      diffloc = (inloc == CELL_CENTRE) ? CELL_ZLOW : CELL_CENTRE;

    } else {
      // Derivative of interpolated field or interpolation of derivative field
      // cannot be taken without communicating and applying boundary
      // conditions, so throw an exception instead
      throw BoutException("Unsupported combination of {inloc =%s} and {outloc =%s} in Mesh:indexD2DZ2(Field3D).", strLocation(inloc), strLocation(outloc));
    }
  }

  if (method != DIFF_DEFAULT) {
    // Lookup function
    func = lookupFunc(table, method);
  }

  if (func == NULL) {
    // Use FFT

    BoutReal shift = 0.; // Shifting result in Z?
    if (StaggerGrids) {
      if ((inloc == CELL_CENTRE) && (diffloc == CELL_ZLOW)) {
	      // Shifting down - multiply by exp(-0.5*i*k*dz) 
        throw BoutException("Not tested - probably broken");
      } else if((inloc == CELL_ZLOW) && (diffloc == CELL_CENTRE)) {
	      // Shifting up
        throw BoutException("Not tested - probably broken");

      } else if (diffloc != CELL_DEFAULT && diffloc != inloc){
        throw BoutException("Not implemented!");
      }
    }

    result.allocate(); // Make sure data allocated

    auto region_index = f.region(region);
    int xs = region_index.xstart;
    int xe = region_index.xend;
    int ys = region_index.ystart;
    int ye = region_index.yend;
    ASSERT2(region_index.zstart == 0);
    int ncz = region_index.zend + 1;

    // TODO: The comment does not match the check
    ASSERT1(ncz % 2 == 0); // Must be a power of 2
    Array<dcomplex> cv(ncz / 2 + 1);
    
    for (int jx = xs; jx <= xe; jx++) {
      for (int jy = ys; jy <= ye; jy++) {

        rfft(f(jx, jy), ncz, cv.begin()); // Forward FFT

        for (int jz = 0; jz <= ncz / 2; jz++) {
          BoutReal kwave = jz * 2.0 * PI / ncz; // wave number is 1/[rad]

          cv[jz] *= -SQ(kwave);
          if (shift)
            cv[jz] *= exp(0.5 * Im * (shift * kwave));
        }

        irfft(cv.begin(), ncz, result(jx, jy)); // Reverse FFT
      }
    }

#if CHECK > 0
    // Mark boundaries as invalid
    result.bndry_xin = false;
    result.bndry_xout = false;
    result.bndry_yup = false;
    result.bndry_ydown = false;
#endif

  } else {
    // All other (non-FFT) functions
    result = applyZdiff(f, func, diffloc, region);
  }

  result.setLocation(diffloc);

  return result;
}

/*******************************************************************************
 * Fourth derivatives
 *******************************************************************************/

BoutReal D4DX4_C2(stencil &f) { return (f.pp - 4. * f.p + 6. * f.c - 4. * f.m + f.mm); }

const Field3D Mesh::indexD4DX4(const Field3D &f, CELL_LOC outloc,
                               DIFF_METHOD method, REGION region) {
  ASSERT1(outloc == CELL_DEFAULT || outloc == f.getLocation());
  ASSERT1(method == DIFF_DEFAULT);
  return applyXdiff(f, D4DX4_C2, f.getLocation(), region);
}

const Field2D Mesh::indexD4DX4(const Field2D &f, CELL_LOC outloc,
                               DIFF_METHOD method, REGION region) {
  ASSERT1(outloc == CELL_DEFAULT || outloc == f.getLocation());
  ASSERT1(method == DIFF_DEFAULT);
  return applyXdiff(f, D4DX4_C2, f.getLocation(), region);
}

const Field3D Mesh::indexD4DY4(const Field3D &f, CELL_LOC outloc,
                               DIFF_METHOD method, REGION region) {
  ASSERT1(outloc == CELL_DEFAULT || outloc == f.getLocation());
  ASSERT1(method == DIFF_DEFAULT);
  return applyYdiff(f, D4DX4_C2, f.getLocation(), region);
}

const Field2D Mesh::indexD4DY4(const Field2D &f, CELL_LOC outloc,
                               DIFF_METHOD method, REGION region) {
  ASSERT1(outloc == CELL_DEFAULT || outloc == f.getLocation());
  ASSERT1(method == DIFF_DEFAULT);
  return applyYdiff(f, D4DX4_C2, f.getLocation(), region);
}

const Field3D Mesh::indexD4DZ4(const Field3D &f, CELL_LOC outloc,
                               DIFF_METHOD method, REGION region){
  ASSERT1(outloc == CELL_DEFAULT || outloc == f.getLocation());
  ASSERT1(method == DIFF_DEFAULT);
  return applyZdiff(f, D4DX4_C2, f.getLocation(), region);
}

const Field2D Mesh::indexD4DZ4(const Field2D &f, CELL_LOC outloc,
                               DIFF_METHOD method, REGION region){
  ASSERT1(outloc == CELL_DEFAULT || outloc == f.getLocation());
  return Field2D(0.,this);
}


/*******************************************************************************
 * Mixed derivatives
 *******************************************************************************/

/*******************************************************************************
 * Advection schemes
 *
 * Jan 2018  - Re-written to use iterators and handle staggering as different cases
 * Jan 2009  - Re-written to use Set*Stencil routines
 *******************************************************************************/

////////////// X DERIVATIVE /////////////////

/// Special case where both arguments are 2D. Output location ignored for now
const Field2D Mesh::indexVDDX(const Field2D &v, const Field2D &f, CELL_LOC outloc,
                              DIFF_METHOD method, REGION region) {
  TRACE("Mesh::indexVDDX(Field2D, Field2D)");

  CELL_LOC diffloc = f.getLocation();

  Mesh::upwind_func func = fVDDX;

  if (method != DIFF_DEFAULT) {
    // Lookup function
    func = lookupUpwindFunc(UpwindTable, method);
  }

  ASSERT1(this == f.getMesh());
  ASSERT1(this == v.getMesh());
  ASSERT2((v.getLocation() == f.getLocation()) && ((outloc == CELL_DEFAULT) || (outloc == f.getLocation()))); // No staggering allowed for Field2D

  Field2D result(this);
  result.allocate(); // Make sure data allocated

  if (this->xstart > 1) {
    // Two or more guard cells

    stencil s;
    for (const auto &i : result.region(region)) {
      s.c = f[i];
      s.p = f[i.xp()];
      s.m = f[i.xm()];
      s.pp = f[i.offset(2, 0, 0)];
      s.mm = f[i.offset(-2, 0, 0)];

      result[i] = func(v[i], s);
    }

  } else if (this->xstart == 1) {
    // Only one guard cell

    stencil s;
    s.pp = nan("");
    s.mm = nan("");

    for (const auto &i : result.region(region)) {
      s.c = f[i];
      s.p = f[i.xp()];
      s.m = f[i.xm()];

      result[i] = func(v[i], s);
    }
  } else {
    // No guard cells
    throw BoutException("Error: Derivatives in X requires at least one guard cell");
  }

#if CHECK > 0
  // Mark boundaries as invalid
  result.bndry_xin = result.bndry_xout = false;
#endif

  result.setLocation(diffloc);

  return result;
}

/// General version for 3D objects.
/// 2D objects passed as input will result in copying
const Field3D Mesh::indexVDDX(const Field3D &v, const Field3D &f, CELL_LOC outloc,
                              DIFF_METHOD method, REGION region) {
  TRACE("Mesh::indexVDDX(Field3D, Field3D)");

  ASSERT1(this == v.getMesh());
  ASSERT1(this == f.getMesh());

  Field3D result(this);
  result.allocate(); // Make sure data allocated

  CELL_LOC vloc = v.getLocation();
  CELL_LOC inloc = f.getLocation(); // Input location
  CELL_LOC diffloc = inloc;         // Location of differential result

  if (StaggerGrids && (outloc == CELL_DEFAULT)) {
    // Take care of CELL_DEFAULT case
    outloc = diffloc; // No shift (i.e. same as no stagger case)
  }

  if (StaggerGrids && (vloc != inloc)) {
    // Staggered grids enabled, and velocity at different location to value

    Mesh::flux_func func = sfVDDX;
    DiffLookup *table = UpwindTable;

    if (vloc == CELL_XLOW) {
      // V staggered w.r.t. variable
      func = sfVDDX;
      table = UpwindStagTable;
      diffloc = CELL_CENTRE;
    } else if ((vloc == CELL_CENTRE) && (inloc == CELL_XLOW)) {
      // Shifted
      func = sfVDDX;
      table = UpwindStagTable;
      diffloc = CELL_XLOW;
    } else {
      // More complicated shifting. The user should probably
      // be explicit about what interpolation should be done

      throw BoutException("Unhandled shift in Mesh::indexVDDX");
    }

    if (method != DIFF_DEFAULT) {
      // Lookup function
      func = lookupFluxFunc(table, method);
    }

    // Note: The velocity stencil contains only (mm, m, p, pp)
    // v.p is v at +1/2, v.m is at -1/2 relative to the field f

    if (this->xstart > 1) {
      // Two or more guard cells

      if ((vloc == CELL_XLOW) && (diffloc == CELL_CENTRE)) {
        stencil fs, vs;
        vs.c = nan("");

        for (const auto &i : result.region(region)) {
          fs.c = f[i];
          fs.p = f[i.xp()];
          fs.m = f[i.xm()];
          fs.pp = f[i.offset(2, 0, 0)];
          fs.mm = f[i.offset(-2, 0, 0)];

          vs.mm = v[i.xm()];
          vs.m = v[i];
          vs.p = v[i.xp()];
          vs.pp = v[i.offset(2, 0, 0)];

          result[i] = func(vs, fs);
        }

      } else if ((vloc == CELL_CENTRE) && (diffloc == CELL_XLOW)) {
        stencil fs, vs;
        vs.c = nan("");

        for (const auto &i : result.region(region)) {
          fs.c = f[i];
          fs.p = f[i.xp()];
          fs.m = f[i.xm()];
          fs.pp = f[i.offset(2, 0, 0)];
          fs.mm = f[i.offset(-2, 0, 0)];

          vs.mm = v[i.offset(-2, 0, 0)];
          vs.m = v[i.xm()];
          vs.p = v[i];
          vs.pp = v[i.xp()];

          result[i] = func(vs, fs);
        }
      } else {
        throw BoutException("Unhandled shift in Mesh::indexVDDX");
      }
    } else if (this->xstart == 1) {
      // One guard cell

      if ((vloc == CELL_XLOW) && (diffloc == CELL_CENTRE)) {
        stencil fs, vs;
        vs.c = nan("");
        vs.pp = nan("");
        fs.pp = nan("");
        fs.mm = nan("");

        for (const auto &i : result.region(region)) {
          fs.c = f[i];
          fs.p = f[i.xp()];
          fs.m = f[i.xm()];

          vs.mm = v[i.xm()];
          vs.m = v[i];
          vs.p = v[i.xp()];

          result[i] = func(vs, fs);
        }

      } else if ((vloc == CELL_CENTRE) && (diffloc == CELL_XLOW)) {
        stencil fs, vs;

        fs.pp = nan("");
        fs.mm = nan("");
        vs.c = nan("");
        vs.mm = nan("");

        for (const auto &i : result.region(region)) {
          fs.c = f[i];
          fs.p = f[i.xp()];
          fs.m = f[i.xm()];

          vs.m = v[i.xm()];
          vs.p = v[i];
          vs.pp = v[i.xp()];

          result[i] = func(vs, fs);
        }
      } else {
        throw BoutException("Unhandled shift in Mesh::indexVDDX");
      }
    } else {
      // No guard cells
      throw BoutException("Error: Derivatives in X requires at least one guard cell");
    }

  } else {
    // Not staggered
    Mesh::upwind_func func = fVDDX;
    DiffLookup *table = UpwindTable;

    if (method != DIFF_DEFAULT) {
      // Lookup function
      func = lookupUpwindFunc(table, method);
    }

    if (this->xstart > 1) {
      // Two or more guard cells
      stencil fs;
      for (const auto &i : result.region(region)) {
        fs.c = f[i];
        fs.p = f[i.xp()];
        fs.m = f[i.xm()];
        fs.pp = f[i.offset(2, 0, 0)];
        fs.mm = f[i.offset(-2, 0, 0)];

        result[i] = func(v[i], fs);
      }
    } else if (this->xstart == 1) {
      // Only one guard cell
      stencil fs;
      fs.pp = nan("");
      fs.mm = nan("");
      for (const auto &i : result.region(region)) {
        fs.c = f[i];
        fs.p = f[i.xp()];
        fs.m = f[i.xm()];

        result[i] = func(v[i], fs);
      }
    } else {
      // No guard cells
      throw BoutException("Error: Derivatives in X requires at least one guard cell");
    }
  }

  result.setLocation(diffloc);

#if CHECK > 0
  // Mark boundaries as invalid
  result.bndry_xin = result.bndry_xout = result.bndry_yup = result.bndry_ydown = false;
#endif

  return result;
}

////////////// Y DERIVATIVE /////////////////

// special case where both are 2D
const Field2D Mesh::indexVDDY(const Field2D &v, const Field2D &f, CELL_LOC outloc,
                              DIFF_METHOD method, REGION region) {
  SCOREP0();
  TRACE("Mesh::indexVDDY");

  ASSERT1(this == v.getMesh());
  ASSERT1(this == f.getMesh());
  ASSERT2((v.getLocation() == f.getLocation()) && ((outloc == CELL_DEFAULT) || (outloc == f.getLocation()))); // No staggering allowed for Field2D

  ASSERT1(this->ystart > 0); // Must have at least one guard cell

  Field2D result(this);
  result.allocate(); // Make sure data allocated

  CELL_LOC vloc = v.getLocation();
  CELL_LOC inloc = f.getLocation(); // Input location
  CELL_LOC diffloc = inloc;         // Location of differential result

  if (StaggerGrids && (outloc == CELL_DEFAULT)) {
    // Take care of CELL_DEFAULT case
    outloc = diffloc; // No shift (i.e. same as no stagger case)
  }

  if (StaggerGrids && (vloc != inloc)) {
    // Staggered grids enabled, and velocity at different location to value

    Mesh::flux_func func = sfVDDY;
    DiffLookup *table = UpwindTable;
    if ((vloc == CELL_YLOW) && (diffloc == CELL_CENTRE)) {
      // V staggered w.r.t. variable
      func = sfVDDY;
      table = UpwindStagTable;
    } else if ((vloc == CELL_CENTRE) && (inloc == CELL_YLOW)) {
      // Shifted
      func = sfVDDY;
      table = UpwindStagTable;
      diffloc = CELL_YLOW;
    } else {
      // More complicated. Deciding what to do here isn't straightforward

      throw BoutException("Unhandled shift in indexVDDY(Field2D, Field2D)");
    }

    if (method != DIFF_DEFAULT) {
      // Lookup function
      func = lookupFluxFunc(table, method);
    }

    // Note: vs.c not used for staggered differencing
    // vs.m is at i-1/2, vs.p is as i+1/2
    if ((vloc == CELL_YLOW) && (diffloc == CELL_CENTRE)) {
      if (this->ystart > 1) {
        // Two or more guard cells
        stencil fs, vs;
        vs.c = nan("");
        for (const auto &i : result.region(region)) {

          fs.c = f[i];
          fs.p = f[i.yp()];
          fs.m = f[i.ym()];
          fs.pp = f[i.offset(0, 2, 0)];
          fs.mm = f[i.offset(0, -2, 0)];

          vs.pp = v[i.offset(0, 2, 0)];
          vs.p = v[i.yp()];
          vs.m = v[i];
          vs.mm = v[i.ym()];

          result[i] = func(vs, fs);
        }

      } else {
        // Only one guard cell

        stencil fs, vs;
        fs.pp = nan("");
        fs.mm = nan("");
        vs.c = nan("");
        vs.pp = nan("");

        for (const auto &i : result.region(region)) {
          fs.c = f[i];
          fs.p = f[i.yp()];
          fs.m = f[i.ym()];

          vs.p = v[i.yp()];
          vs.m = v[i];
          vs.mm = v[i.ym()];

          result[i] = func(vs, fs);
        }
      }
    } else if ((vloc == CELL_CENTRE) && (inloc == CELL_YLOW)) {
      if (this->ystart > 1) {
        // Two or more guard cells
        stencil fs, vs;
        vs.c = nan("");
        for (const auto &i : result.region(region)) {

          fs.c = f[i];
          fs.p = f[i.yp()];
          fs.m = f[i.ym()];
          fs.pp = f[i.offset(0, 2, 0)];
          fs.mm = f[i.offset(0, -2, 0)];

          vs.pp = v[i.yp()];
          vs.p = v[i];
          vs.m = v[i.ym()];
          vs.mm = v[i.offset(0, -2, 0)];

          result[i] = func(vs, fs);
        }

      } else {
        // Only one guard cell

        stencil fs, vs;
        fs.pp = nan("");
        fs.mm = nan("");
        vs.c = nan("");
        vs.mm = nan("");

        for (const auto &i : result.region(region)) {
          fs.c = f[i];
          fs.p = f[i.yp()];
          fs.m = f[i.ym()];

          vs.pp = v[i.yp()];
          vs.p = v[i];
          vs.m = v[i.ym()];

          result[i] = func(vs, fs);
        }
      }
    } else {
      throw BoutException("Unhandled shift in indexVDDY(Field2D, Field2D)");
    }

  } else {
    // Not staggered

    Mesh::upwind_func func = fVDDY;
    DiffLookup *table = UpwindTable;

    if (method != DIFF_DEFAULT) {
      // Lookup function
      func = lookupUpwindFunc(table, method);
    }

    if (this->ystart > 1) {
      // Two or more guard cells
      stencil fs;
      for (const auto &i : result.region(region)) {

        fs.c = f[i];
        fs.p = f[i.yp()];
        fs.m = f[i.ym()];
        fs.pp = f[i.offset(0, 2, 0)];
        fs.mm = f[i.offset(0, -2, 0)];

        result[i] = func(v[i], fs);
      }

    } else {
      // Only one guard cell

      stencil fs;
      fs.pp = nan("");
      fs.mm = nan("");

      for (const auto &i : result.region(region)) {
        fs.c = f[i];
        fs.p = f[i.yp()];
        fs.m = f[i.ym()];

        result[i] = func(v[i], fs);
      }
    }
  }

  result.setLocation(diffloc);

#if CHECK > 0
  // Mark boundaries as invalid
  result.bndry_xin = result.bndry_xout = result.bndry_yup = result.bndry_ydown = false;
#endif

  return result;
}

// general case
const Field3D Mesh::indexVDDY(const Field3D &v, const Field3D &f, CELL_LOC outloc,
                              DIFF_METHOD method, REGION region) {
  SCOREP0();
  TRACE("Mesh::indexVDDY(Field3D, Field3D, ..., REGION<Ind3D>)");

  ASSERT1(this == v.getMesh());
  ASSERT1(this == f.getMesh());

  ASSERT1(this->ystart > 0); // Need at least one guard cell

  Field3D result(this);
  result.allocate(); // Make sure data allocated

  CELL_LOC vloc = v.getLocation();
  CELL_LOC inloc = f.getLocation(); // Input location
  CELL_LOC diffloc = inloc;         // Location of differential result

  if (StaggerGrids && (outloc == CELL_DEFAULT)) {
    // Take care of CELL_DEFAULT case
    outloc = diffloc; // No shift (i.e. same as no stagger case)
  }

  if (StaggerGrids && (vloc != inloc)) {
    // Staggered grids enabled, and velocity at different location to value

    Mesh::flux_func func = sfVDDY;
    DiffLookup *table = UpwindTable;

    if (vloc == CELL_YLOW) {
      // V staggered w.r.t. variable
      func = sfVDDY;
      table = UpwindStagTable;
      diffloc = CELL_CENTRE;
    } else if ((vloc == CELL_CENTRE) && (inloc == CELL_YLOW)) {
      // Shifted
      func = sfVDDY;
      table = UpwindStagTable;
      diffloc = CELL_YLOW;
    } else {
      // More complicated. Deciding what to do here isn't straightforward

      throw BoutException("Unhandled shift in VDDY(Field, Field)");
    }

    if (method != DIFF_DEFAULT) {
      // Lookup function
      func = lookupFluxFunc(table, method);
    }

    // There are four cases, corresponding to whether or not f and v
    // have yup, ydown fields.

    // If vUseUpDown is true, field "v" has distinct yup and ydown fields which
    // will be used to calculate a derivative along
    // the magnetic field
    bool vUseUpDown = (v.hasYupYdown() && ((&v.yup() != &v) || (&v.ydown() != &v)));
    bool fUseUpDown = (f.hasYupYdown() && ((&f.yup() != &f) || (&f.ydown() != &f)));

    if (vUseUpDown && fUseUpDown) {
      // Both v and f have up/down fields

BOUT_OMP(parallel)
{
      stencil vval, fval;
      IndexOffset<Ind3D> offset(*mesh);
      //for (const auto &i : result.region(region)) {
      BLOCK_REGION_LOOP_PARALLEL_SECTION( result.getMesh()->getRegion3D(region), i,

        vval.mm = nan("");
        vval.m = v.ydown()[offset.ym(i)];
        vval.c = v[i];
        vval.p = v.yup()[offset.yp(i)];
        vval.pp = nan("");

        fval.mm = nan("");
        fval.m = f.ydown()[offset.ym(i)];
        fval.c = f[i];
        fval.p = f.yup()[offset.yp(i)];
        fval.pp = nan("");

        if (diffloc != CELL_DEFAULT) {
          // Non-centred stencil
          if ((vloc == CELL_CENTRE) && (diffloc == CELL_YLOW)) {
            // Producing a stencil centred around a lower Y value
            vval.pp = vval.p;
            vval.p = vval.c;
          } else if (vloc == CELL_YLOW) {
            // Stencil centred around a cell centre
            vval.mm = vval.m;
            vval.m = vval.c;
          }
          // Shifted in one direction -> shift in another
          // Could produce warning
        }
        result[i] = func(vval, fval);
      //}
      );
}
    } else if (vUseUpDown) {
      // Only v has up/down fields
      // f must shift to field aligned coordinates
      Field3D f_fa = this->toFieldAligned(f);

BOUT_OMP(parallel)
{
      stencil vval, fval;
      IndexOffset<Ind3D> offset(*mesh);
      BLOCK_REGION_LOOP_PARALLEL_SECTION( result.getMesh()->getRegion3D(region), i,
        vval.mm = nan("");
        vval.m = v.ydown()[offset.ym(i)];
        vval.c = v[i];
        vval.p = v.yup()[offset.yp(i)];
        vval.pp = nan("");

        fval.mm = f_fa[offset.ymm(i)];
        fval.m = f_fa[offset.ym(i)];
        fval.c = f_fa[i];
        fval.p = f_fa[offset.yp(i)];
        fval.pp = f_fa[offset.ypp(i)];

        if (diffloc != CELL_DEFAULT) {
          // Non-centred stencil
          if ((vloc == CELL_CENTRE) && (diffloc == CELL_YLOW)) {
            // Producing a stencil centred around a lower Y value
            vval.pp = vval.p;
            vval.p = vval.c;
          } else if (vloc == CELL_YLOW) {
            // Stencil centred around a cell centre
            vval.mm = vval.m;
            vval.m = vval.c;
          }
          // Shifted in one direction -> shift in another
          // Could produce warning
        }
        result[i] = func(vval, fval);
      );
}
    } else if (fUseUpDown) {
      // Only f has up/down fields
      // v must shift to field aligned coordinates
      Field3D v_fa = this->toFieldAligned(v);

BOUT_OMP(parallel)
{
      stencil vval, fval;
      IndexOffset<Ind3D> offset(*mesh);
      BLOCK_REGION_LOOP_PARALLEL_SECTION( result.getMesh()->getRegion3D(region), i,

        vval.mm = v_fa[offset.ymm(i)];
        vval.m = v_fa[offset.ym(i)];
        vval.c = v_fa[i];
        vval.p = v_fa[offset.yp(i)];
        vval.pp = v_fa[offset.ypp(i)];

        fval.mm = nan("");
        fval.m = f.ydown()[offset.ym(i)];
        fval.c = f[i];
        fval.p = f.yup()[offset.yp(i)];
        fval.pp = nan("");

        if (diffloc != CELL_DEFAULT) {
          // Non-centred stencil
          if ((vloc == CELL_CENTRE) && (diffloc == CELL_YLOW)) {
            // Producing a stencil centred around a lower Y value
            vval.pp = vval.p;
            vval.p = vval.c;
          } else if (vloc == CELL_YLOW) {
            // Stencil centred around a cell centre
            vval.mm = vval.m;
            vval.m = vval.c;
          }
          // Shifted in one direction -> shift in another
          // Could produce warning
        }
        result[i] = func(vval, fval);
      );
}
    } else {
      // Both must shift to field aligned
      Field3D v_fa = this->toFieldAligned(v);
      Field3D f_fa = this->toFieldAligned(f);

      BOUT_OMP(parallel)
      {
        stencil vval, fval;
        IndexOffset<Ind3D> offset(*mesh);

        //for (const auto &i : result.region(region)) {
        BLOCK_REGION_LOOP_PARALLEL_SECTION( result.getMesh()->getRegion3D(region), i,

          vval.mm = v_fa[offset.ymm(i)];
          vval.m = v_fa[offset.ym(i)];
          vval.c = v_fa[i];
          vval.p = v_fa[offset.yp(i)];
          vval.pp = v_fa[offset.ypp(i)];

          fval.mm = f_fa[offset.ymm(i)];
          fval.m = f_fa[offset.ym(i)];
          fval.c = f[i];
          fval.p = f_fa[offset.yp(i)];
          fval.pp = f_fa[offset.ypp(i)];

          if (diffloc != CELL_DEFAULT) {
            // Non-centred stencil
            if ((vloc == CELL_CENTRE) && (diffloc == CELL_YLOW)) {
              // Producing a stencil centred around a lower Y value
              vval.pp = vval.p;
              vval.p = vval.c;
            } else if (vloc == CELL_YLOW) {
              // Stencil centred around a cell centre
              vval.mm = vval.m;
              vval.m = vval.c;
            }
            // Shifted in one direction -> shift in another
            // Could produce warning
          }
          result[i] = func(vval, fval);
        //}
        );
      }
    }
  } else {
    // Non-staggered case

    Mesh::upwind_func func = fVDDY;
    DiffLookup *table = UpwindTable;

    if (method != DIFF_DEFAULT) {
      // Lookup function
      func = lookupUpwindFunc(table, method);
    }

    if (f.hasYupYdown() && ((&f.yup() != &f) || (&f.ydown() != &f))) {
      // f has yup and ydown fields which are distinct

      BOUT_OMP(parallel)
      {
        stencil fs;
        fs.pp = nan("");
        fs.mm = nan("");

        Field3D f_yup = f.yup();
        Field3D f_ydown = f.ydown();

        IndexOffset<Ind3D> offset(*mesh);

        //for (const auto &i : result.region(region)) {
        BLOCK_REGION_LOOP_PARALLEL_SECTION( result.getMesh()->getRegion3D(region), i,

          fs.m = f_ydown[offset.ym(i)];
          fs.c = f[i];
          fs.p = f_yup[offset.yp(i)];

          result[i] = func(v[i], fs);
        );
      }

    } else {
      // Not using yup/ydown fields, so first transform to field-aligned coordinates

      Field3D f_fa = this->toFieldAligned(f);
      Field3D v_fa = this->toFieldAligned(v);

      if (this->ystart > 1) {
        BOUT_OMP(parallel)
        {

          stencil fs;
          IndexOffset<Ind3D> offset(*mesh);

          //for (const auto &i : result.region(region)) {
          BLOCK_REGION_LOOP_PARALLEL_SECTION( result.getMesh()->getRegion3D(region), i,

            fs.mm = f_fa[offset.ymm(i)];
            fs.m = f_fa[offset.ym(i)];
            fs.c = f_fa[i];
            fs.p = f_fa[offset.yp(i)];
            fs.pp = f_fa[offset.ypp(i)];

            result[i] = func(v_fa[i], fs);
          );
        }
      } else {
        BOUT_OMP(parallel)
        {

          stencil fs;
          IndexOffset<Ind3D> offset(*mesh);
          fs.mm = nan("");
          fs.pp = nan("");

          //for (const auto &i : result.region(region)) {
          BLOCK_REGION_LOOP_PARALLEL_SECTION( result.getMesh()->getRegion3D(region), i,

            fs.m = f_fa[offset.ym(i)];
            fs.c = f_fa[i];
            fs.p = f_fa[offset.yp(i)];

            result[i] = func(v_fa[i], fs);
          );
        }
      }
      // Shift result back
      result = this->fromFieldAligned(result);
    }
  }

  result.setLocation(diffloc);

#if CHECK > 0
  // Mark boundaries as invalid
  result.bndry_xin = result.bndry_xout = result.bndry_yup = result.bndry_ydown = false;
#endif

  return result;
}

////////////// Z DERIVATIVE /////////////////

// general case
const Field3D Mesh::indexVDDZ(const Field3D &v, const Field3D &f, CELL_LOC outloc,
                              DIFF_METHOD method, REGION region) {
  SCOREP0();
  TRACE("Mesh::indexVDDZ");

  ASSERT1(this == v.getMesh());
  ASSERT1(this == f.getMesh());

  Field3D result(this);
  result.allocate(); // Make sure data allocated

  CELL_LOC vloc = v.getLocation();
  CELL_LOC inloc = f.getLocation(); // Input location
  CELL_LOC diffloc = inloc;         // Location of differential result

  if (StaggerGrids && (outloc == CELL_DEFAULT)) {
    // Take care of CELL_DEFAULT case
    outloc = diffloc; // No shift (i.e. same as no stagger case)
  }

  if (StaggerGrids && (vloc != inloc)) {
    // Staggered grids enabled, and velocity at different location to value

    Mesh::flux_func func = sfVDDZ;
    DiffLookup *table = UpwindTable;

    if (vloc == CELL_ZLOW) {
      // V staggered w.r.t. variable
      func = sfVDDZ;
      table = UpwindStagTable;
      diffloc = CELL_CENTRE;
    } else if ((vloc == CELL_CENTRE) && (inloc == CELL_ZLOW)) {
      // Shifted
      func = sfVDDZ;
      table = UpwindStagTable;
      diffloc = CELL_ZLOW;
    } else {
      // More complicated. Deciding what to do here isn't straightforward

      throw BoutException("Unhandled shift in indexVDDZ");
    }

    if (method != DIFF_DEFAULT) {
      // Lookup function
      func = lookupFluxFunc(table, method);
    }

    stencil vval, fval;
    for (const auto &i : result.region(region)) {
      fval.mm = f[i.offset(0,0,-2)];
      fval.m = f[i.zm()];
      fval.c = f[i];
      fval.p = f[i.zp()];
      fval.pp = f[i.offset(0,0,2)];

      vval.mm = v[i.offset(0,0,-2)];
      vval.m = v[i.zm()];
      vval.c = v[i];
      vval.p = v[i.zp()];
      vval.pp = v[i.offset(0,0,2)];

      if((diffloc != CELL_DEFAULT) && (diffloc != vloc)) {
        // Non-centred stencil

        if((vloc == CELL_CENTRE) && (diffloc == CELL_ZLOW)) {
          // Producing a stencil centred around a lower Z value
          vval.pp = vval.p;
          vval.p  = vval.c;

        }else if(vloc == CELL_ZLOW) {
          // Stencil centred around a cell centre

          vval.mm = vval.m;
          vval.m  = vval.c;
        }
        // Shifted in one direction -> shift in another
        // Could produce warning
      }
      result[i] = func(vval, fval);
    }
  } else {
    Mesh::upwind_func func = fVDDZ;
    DiffLookup *table = UpwindTable;

    if (method != DIFF_DEFAULT) {
      // Lookup function
      func = lookupUpwindFunc(table, method);
    }

    stencil fval;
    for (const auto &i : result.region(region)) {
      fval.mm = f[i.offset(0,0,-2)];
      fval.m = f[i.zm()];
      fval.c = f[i];
      fval.p = f[i.zp()];
      fval.pp = f[i.offset(0,0,2)];

      result[i] = func(v[i], fval);
    }
  }

  result.setLocation(diffloc);

#if CHECK > 0
  // Mark boundaries as invalid
  result.bndry_xin = result.bndry_xout = result.bndry_yup = result.bndry_ydown = false;
#endif

  return result;
}

/*******************************************************************************
 * Flux conserving schemes
 *******************************************************************************/

const Field2D Mesh::indexFDDX(const Field2D &v, const Field2D &f, CELL_LOC outloc,
                              DIFF_METHOD method, REGION region) {
  TRACE("Mesh::::indexFDDX(Field2D, Field2D)");

  CELL_LOC diffloc = f.getLocation();

  if ((method == DIFF_SPLIT) || ((method == DIFF_DEFAULT) && (fFDDX == NULL))) {
    // Split into an upwind and a central differencing part
    // d/dx(v*f) = v*d/dx(f) + f*d/dx(v)
    return indexVDDX(v, f, outloc, DIFF_DEFAULT) + f * indexDDX(v);
  }

  Mesh::flux_func func = fFDDX;
  if (method != DIFF_DEFAULT) {
    // Lookup function
    func = lookupFluxFunc(FluxTable, method);
  }

  Field2D result(this);
  result.allocate(); // Make sure data allocated

  if (StaggerGrids &&
      ((v.getLocation() != CELL_CENTRE) || (f.getLocation() != CELL_CENTRE))) {
    // Staggered differencing
    throw BoutException("Unhandled staggering");
  }
  else {
    result.setLocation(CELL_CENTRE);
  }

  ASSERT1(this == v.getMesh());
  ASSERT1(this == f.getMesh());

  if (this->xstart > 1) {
    // Two or more guard cells

    stencil fs;
    stencil vs;
    for (const auto &i : result.region(region)) {
      fs.c = f[i];
      fs.p = f[i.xp()];
      fs.m = f[i.xm()];
      fs.pp = f[i.offset(2, 0, 0)];
      fs.mm = f[i.offset(-2, 0, 0)];

      vs.c = v[i];
      vs.p = v[i.xp()];
      vs.m = v[i.xm()];
      vs.pp = v[i.offset(2, 0, 0)];
      vs.mm = v[i.offset(-2, 0, 0)];

      result[i] = func(vs, fs);
    }
  } else if (this->xstart == 1) {
    // Only one guard cell

    stencil fs;
    fs.pp = nan("");
    fs.mm = nan("");
    stencil vs;
    vs.pp = nan("");
    vs.mm = nan("");

    for (const auto &i : result.region(region)) {
      fs.c = f[i];
      fs.p = f[i.xp()];
      fs.m = f[i.xm()];

      vs.c = v[i];
      vs.p = v[i.xp()];
      vs.m = v[i.xm()];

      result[i] = func(vs, fs);
    }
  } else {
    // No guard cells
    throw BoutException("Error: Derivatives in X requires at least one guard cell");
  }

#if CHECK > 0
  // Mark boundaries as invalid
  result.bndry_xin = result.bndry_xout = false;
#endif

  return result;
}

const Field3D Mesh::indexFDDX(const Field3D &v, const Field3D &f, CELL_LOC outloc,
                              DIFF_METHOD method, REGION region) {
  TRACE("Mesh::indexFDDX(Field3D, Field3D)");

  if ((method == DIFF_SPLIT) || ((method == DIFF_DEFAULT) && (fFDDX == NULL))) {
    // Split into an upwind and a central differencing part
    // d/dx(v*f) = v*d/dx(f) + f*d/dx(v)
    return indexVDDX(v, f, outloc, DIFF_DEFAULT) + indexDDX(v, outloc, DIFF_DEFAULT) * f;
  }

  Mesh::flux_func func = fFDDX;
  DiffLookup *table = FluxTable;

  CELL_LOC vloc = v.getLocation();
  CELL_LOC inloc = f.getLocation(); // Input location
  CELL_LOC diffloc = inloc;         // Location of differential result

  if (StaggerGrids && (outloc == CELL_DEFAULT)) {
    // Take care of CELL_DEFAULT case
    outloc = diffloc; // No shift (i.e. same as no stagger case)
  }

  if (StaggerGrids && (vloc != inloc)) {
    // Staggered grids enabled, and velocity at different location to value
    if (vloc == CELL_XLOW) {
      // V staggered w.r.t. variable
      func = sfFDDX;
      table = FluxStagTable;
      diffloc = CELL_CENTRE;
    } else if ((vloc == CELL_CENTRE) && (inloc == CELL_XLOW)) {
      // Shifted
      func = sfFDDX;
      table = FluxStagTable;
      diffloc = CELL_XLOW;
    } else {
      // More complicated. Deciding what to do here isn't straightforward
      throw BoutException("Unhandled shift in indexFDDX");
    }
  }

  if (method != DIFF_DEFAULT) {
    // Lookup function
    func = lookupFluxFunc(table, method);
  }

  ASSERT1(this == f.getMesh());
  ASSERT1(this == v.getMesh());

  Field3D result(this);
  result.allocate(); // Make sure data allocated

  if (this->xstart > 1) {
    // Two or more guard cells
    if (StaggerGrids) {
      if ((vloc == CELL_CENTRE) && (diffloc == CELL_XLOW)) {
        // Producing a stencil centred around a lower X value

        stencil fs, vs;
        vs.c = nan("");
        for (const auto &i : result.region(region)) {
          // Location of f always the same as the output
          fs.c = f[i];
          fs.p = f[i.xp()];
          fs.m = f[i.xm()];
          fs.pp = f[i.offset(2, 0, 0)];
          fs.mm = f[i.offset(-2, 0, 0)];

          // Note: Location in diffloc

          vs.mm = v[i.offset(-2, 0, 0)];
          vs.m = v[i.xm()];
          vs.p = v[i];
          vs.pp = v[i.xp()];

          result[i] = func(vs, fs);
        }
      } else if ((vloc == CELL_XLOW) && (diffloc == CELL_CENTRE)) {
        // Stencil centred around a cell centre
        stencil fs, vs;
        vs.c = nan("");
        for (const auto &i : result.region(region)) {
          // Location of f always the same as the output
          fs.c = f[i];
          fs.p = f[i.xp()];
          fs.m = f[i.xm()];
          fs.pp = f[i.offset(2, 0, 0)];
          fs.mm = f[i.offset(-2, 0, 0)];

          vs.mm = v[i.xm()];
          vs.m = v[i];
          vs.p = v[i.xp()];
          vs.pp = v[i.offset(2, 0, 0)];

          result[i] = func(vs, fs);
        }
      } else {
        throw BoutException("Unhandled staggering");
      }
    } else {
      // Non-staggered, two or more guard cells
      stencil fs;
      stencil vs;
      for (const auto &i : result.region(region)) {
        // Location of f always the same as the output
        fs.c = f[i];
        fs.p = f[i.xp()];
        fs.m = f[i.xm()];
        fs.pp = f[i.offset(2, 0, 0)];
        fs.mm = f[i.offset(-2, 0, 0)];

        // Note: Location in diffloc
        vs.c = v[i];
        vs.p = v[i.xp()];
        vs.m = v[i.xm()];
        vs.pp = v[i.offset(2, 0, 0)];
        vs.mm = v[i.offset(-2, 0, 0)];

        result[i] = func(vs, fs);
      }
    }
  } else if (this->xstart == 1) {
    // One guard cell

    stencil fs;
    fs.pp = nan("");
    fs.mm = nan("");

    stencil vs;
    vs.pp = nan("");
    vs.mm = nan("");
    vs.c = nan("");

    if (StaggerGrids) {
      if ((vloc == CELL_CENTRE) && (diffloc == CELL_XLOW)) {
        // Producing a stencil centred around a lower X value

        for (const auto &i : result.region(region)) {
          // Location of f always the same as the output
          fs.c = f[i];
          fs.p = f[i.xp()];
          fs.m = f[i.xm()];

          // Note: Location in diffloc
          vs.m = v[i.xm()];
          vs.p = v[i];
          vs.pp = v[i.xp()];

          result[i] = func(vs, fs);
        }
      } else if ((vloc == CELL_XLOW) && (diffloc == CELL_CENTRE)) {
        // Stencil centred around a cell centre
        for (const auto &i : result.region(region)) {
          // Location of f always the same as the output
          fs.c = f[i];
          fs.p = f[i.xp()];
          fs.m = f[i.xm()];

          vs.mm = v[i.xm()];
          vs.m = v[i];
          vs.p = v[i.xp()];

          result[i] = func(vs, fs);
        }
      } else {
        throw BoutException("Unhandled staggering");
      }
    } else {
      // Non-staggered, one guard cell
      for (const auto &i : result.region(region)) {
        // Location of f always the same as the output
        fs.c = f[i];
        fs.p = f[i.xp()];
        fs.m = f[i.xm()];

        vs.c = v[i];
        vs.p = v[i.xp()];
        vs.m = v[i.xm()];

        result[i] = func(vs, fs);
      }
    }
  } else {
    // No guard cells
    throw BoutException("Error: Derivatives in X requires at least one guard cell");
  }

  result.setLocation(diffloc);

#if CHECK > 0
  // Mark boundaries as invalid
  result.bndry_xin = result.bndry_xout = result.bndry_yup = result.bndry_ydown = false;
#endif

  return result;
}

/////////////////////////////////////////////////////////////////////////

const Field2D Mesh::indexFDDY(const Field2D &v, const Field2D &f, CELL_LOC outloc,
                              DIFF_METHOD method, REGION region) {
  TRACE("Mesh::indexFDDY(Field2D, Field2D)");

  ASSERT1(this == v.getMesh());
  ASSERT1(this == f.getMesh());

  CELL_LOC diffloc = f.getLocation();

  if ((method == DIFF_SPLIT) || ((method == DIFF_DEFAULT) && (fFDDY == NULL))) {
    // Split into an upwind and a central differencing part
    // d/dx(v*f) = v*d/dx(f) + f*d/dx(v)
    return indexVDDY(v, f, outloc, DIFF_DEFAULT) + f * indexDDY(v);
  }

  Mesh::flux_func func = fFDDY;
  if (method != DIFF_DEFAULT) {
    // Lookup function
    func = lookupFluxFunc(FluxTable, method);
  }

  Field2D result(this);
  result.allocate(); // Make sure data allocated

  if (StaggerGrids &&
      ((v.getLocation() != CELL_CENTRE) || (f.getLocation() != CELL_CENTRE))) {
    // Staggered differencing
    throw BoutException("Unhandled staggering");
  }
  else {
    result.setLocation(CELL_CENTRE);
  }

  if (this->ystart > 1) {
    // Two or more guard cells
    stencil fs, vs;
    for (const auto &i : result.region(region)) {

      fs.c = f[i];
      fs.p = f[i.yp()];
      fs.m = f[i.ym()];
      fs.pp = f[i.offset(0, 2, 0)];
      fs.mm = f[i.offset(0, -2, 0)];

      vs.c = v[i];
      vs.p = v[i.yp()];
      vs.m = v[i.ym()];
      vs.pp = v[i.offset(0, 2, 0)];
      vs.mm = v[i.offset(0, -2, 0)];

      result[i] = func(vs, fs);
    }

  } else if (this->ystart == 1) {
    // Only one guard cell

    stencil fs;
    fs.pp = nan("");
    fs.mm = nan("");
    stencil vs;
    vs.pp = nan("");
    vs.mm = nan("");

    for (const auto &i : result.region(region)) {
      fs.c = f[i];
      fs.p = f[i.yp()];
      fs.m = f[i.ym()];

      vs.c = v[i];
      vs.p = v[i.yp()];
      vs.m = v[i.ym()];
      result[i] = func(vs, fs);
    }
  } else {
    // No guard cells
    throw BoutException("Error: Derivatives in Y requires at least one guard cell");
  }

  result.setLocation(diffloc);

#if CHECK > 0
  // Mark boundaries as invalid
  result.bndry_xin = result.bndry_xout = false;
#endif

  return result;
}

const Field3D Mesh::indexFDDY(const Field3D &v, const Field3D &f, CELL_LOC outloc,
                              DIFF_METHOD method, REGION region) {
  TRACE("Mesh::indexFDDY");

  if ((method == DIFF_SPLIT) || ((method == DIFF_DEFAULT) && (fFDDY == NULL))) {
    // Split into an upwind and a central differencing part
    // d/dx(v*f) = v*d/dx(f) + f*d/dx(v)
    return indexVDDY(v, f, outloc, DIFF_DEFAULT) + indexDDY(v, outloc, DIFF_DEFAULT) * f;
  }
  Mesh::flux_func func = fFDDY;
  DiffLookup *table = FluxTable;

  CELL_LOC vloc = v.getLocation();
  CELL_LOC inloc = f.getLocation(); // Input location
  CELL_LOC diffloc = inloc;         // Location of differential result

  if (StaggerGrids && (outloc == CELL_DEFAULT)) {
    // Take care of CELL_DEFAULT case
    outloc = diffloc; // No shift (i.e. same as no stagger case)
  }

  if (StaggerGrids && (vloc != inloc)) {
    // Staggered grids enabled, and velocity at different location to value
    if (vloc == CELL_YLOW) {
      // V staggered w.r.t. variable
      func = sfFDDY;
      table = FluxStagTable;
      diffloc = CELL_CENTRE;
    } else if ((vloc == CELL_CENTRE) && (inloc == CELL_YLOW)) {
      // Shifted
      func = sfFDDY;
      table = FluxStagTable;
      diffloc = CELL_YLOW;
    } else {
      // More complicated. Deciding what to do here isn't straightforward
      throw BoutException("Unhandled shift in indexFDDY");
    }
  }

  if (method != DIFF_DEFAULT) {
    // Lookup function
    func = lookupFluxFunc(table, method);
  }

  if (func == NULL) {
    // To catch when no function
    return indexVDDY(v, f, outloc, DIFF_DEFAULT) + indexDDY(v, outloc, DIFF_DEFAULT) * f;
  }

  ASSERT1(this == v.getMesh());
  ASSERT1(this == f.getMesh());

  Field3D result(this);
  result.allocate(); // Make sure data allocated

  // There are four cases, corresponding to whether or not f and v
  // have yup, ydown fields.

  // If vUseUpDown is true, field "v" has distinct yup and ydown fields which
  // will be used to calculate a derivative along
  // the magnetic field
  bool vUseUpDown = (v.hasYupYdown() && ((&v.yup() != &v) || (&v.ydown() != &v)));
  bool fUseUpDown = (f.hasYupYdown() && ((&f.yup() != &f) || (&f.ydown() != &f)));

  if (vUseUpDown && fUseUpDown) {
    // Both v and f have up/down fields
    stencil vval, fval;
    vval.mm = nan("");
    vval.pp = nan("");
    fval.mm = nan("");
    fval.pp = nan("");
    for (const auto &i : result.region(region)) {

      fval.m = f.ydown()[i.ym()];
      fval.c = f[i];
      fval.p = f.yup()[i.yp()];

      vval.m = v.ydown()[i.ym()];
      vval.c = v[i];
      vval.p = v.yup()[i.yp()];

      if(StaggerGrids && (diffloc != CELL_DEFAULT) && (diffloc != vloc)) {
        // Non-centred stencil
        if((vloc == CELL_CENTRE) && (diffloc == CELL_YLOW)) {
          // Producing a stencil centred around a lower Y value
          vval.pp = vval.p;
          vval.p  = vval.c;
        }else if(vloc == CELL_YLOW) {
          // Stencil centred around a cell centre
          vval.mm = vval.m;
          vval.m  = vval.c;
        }
        // Shifted in one direction -> shift in another
        // Could produce warning
      }
      result[i] = func(vval, fval);
    }
  }
  else if (vUseUpDown) {
    // Only v has up/down fields
    // f must shift to field aligned coordinates
    Field3D f_fa = this->toFieldAligned(f);

    stencil vval;
    vval.mm = nan("");
    vval.pp = nan("");

    stencil fval;
    for (const auto &i : result.region(region)) {

      fval.mm = f_fa[i.offset(0, -2, 0)];
      fval.m = f_fa[i.ym()];
      fval.c = f_fa[i];
      fval.p = f_fa[i.yp()];
      fval.pp = f_fa[i.offset(0, 2, 0)];

      vval.m = v.ydown()[i.ym()];
      vval.c = v[i];
      vval.p = v.yup()[i.yp()];

      if(StaggerGrids && (diffloc != CELL_DEFAULT) && (diffloc != vloc)) {
        // Non-centred stencil
        if((vloc == CELL_CENTRE) && (diffloc == CELL_YLOW)) {
          // Producing a stencil centred around a lower Y value
          vval.pp = vval.p;
          vval.p  = vval.c;
        }else if(vloc == CELL_YLOW) {
          // Stencil centred around a cell centre
          vval.mm = vval.m;
          vval.m  = vval.c;
        }
        // Shifted in one direction -> shift in another
        // Could produce warning
      }
      result[i] = func(vval, fval);
    }
  }
  else if (fUseUpDown) {
    // Only f has up/down fields
    // v must shift to field aligned coordinates
    Field3D v_fa = this->toFieldAligned(v);

    stencil vval;

    stencil fval;
    fval.pp = nan("");
    fval.mm = nan("");

    for (const auto &i : result.region(region)) {

      fval.m = f.ydown()[i.ym()];
      fval.c = f[i];
      fval.p = f.yup()[i.yp()];

      vval.mm = v_fa[i.offset(0,-2,0)];
      vval.m = v_fa[i.ym()];
      vval.c = v_fa[i];
      vval.p = v_fa[i.yp()];
      vval.pp = v_fa[i.offset(0,2,0)];

      if(StaggerGrids && (diffloc != CELL_DEFAULT) && (diffloc != vloc)) {
        // Non-centred stencil
        if((vloc == CELL_CENTRE) && (diffloc == CELL_YLOW)) {
          // Producing a stencil centred around a lower Y value
          vval.pp = vval.p;
          vval.p  = vval.c;
        }else if(vloc == CELL_YLOW) {
          // Stencil centred around a cell centre
          vval.mm = vval.m;
          vval.m  = vval.c;
        }
        // Shifted in one direction -> shift in another
        // Could produce warning
      }
      result[i] = func(vval, fval);
    }
  }
  else {
    // Both must shift to field aligned
    Field3D v_fa = this->toFieldAligned(v);
    Field3D f_fa = this->toFieldAligned(f);

    stencil vval, fval;

    for (const auto &i : result.region(region)) {

      fval.mm = f_fa[i.offset(0,-2,0)];
      fval.m = f_fa[i.ym()];
      fval.c = f_fa[i];
      fval.p = f_fa[i.yp()];
      fval.pp = f_fa[i.offset(0,2,0)];

      vval.mm = v_fa[i.offset(0,-2,0)];
      vval.m = v_fa[i.ym()];
      vval.c = v_fa[i];
      vval.p = v_fa[i.yp()];
      vval.pp = v_fa[i.offset(0,2,0)];

      if(StaggerGrids && (diffloc != CELL_DEFAULT) && (diffloc != vloc)) {
        // Non-centred stencil
        if((vloc == CELL_CENTRE) && (diffloc == CELL_YLOW)) {
          // Producing a stencil centred around a lower Y value
          vval.pp = vval.p;
          vval.p  = vval.c;
        }else if(vloc == CELL_YLOW) {
          // Stencil centred around a cell centre
          vval.mm = vval.m;
          vval.m  = vval.c;
        }
        // Shifted in one direction -> shift in another
        // Could produce warning
      }
      result[i] = func(vval, fval);
    }
  }

  result.setLocation(diffloc);

#if CHECK > 0
  // Mark boundaries as invalid
  result.bndry_xin = result.bndry_xout = result.bndry_yup = result.bndry_ydown = false;
#endif

  return result;
}

/////////////////////////////////////////////////////////////////////////

const Field3D Mesh::indexFDDZ(const Field3D &v, const Field3D &f, CELL_LOC outloc,
                              DIFF_METHOD method, REGION region) {
  TRACE("Mesh::indexFDDZ(Field3D, Field3D)");
  if ((method == DIFF_SPLIT) || ((method == DIFF_DEFAULT) && (fFDDZ == NULL))) {
    // Split into an upwind and a central differencing part
    // d/dx(v*f) = v*d/dx(f) + f*d/dx(v)
    return indexVDDZ(v, f, outloc, DIFF_DEFAULT) +
           indexDDZ(v, outloc, DIFF_DEFAULT, true) * f;
  }

  Mesh::flux_func func = fFDDZ;
  DiffLookup *table = FluxTable;

  CELL_LOC vloc = v.getLocation();
  CELL_LOC inloc = f.getLocation(); // Input location
  CELL_LOC diffloc = inloc;         // Location of differential result

  if (StaggerGrids && (outloc == CELL_DEFAULT)) {
    // Take care of CELL_DEFAULT case
    outloc = diffloc; // No shift (i.e. same as no stagger case)
  }

  if (StaggerGrids && (vloc != inloc)) {
    // Staggered grids enabled, and velocity at different location to value
    if (vloc == CELL_ZLOW) {
      // V staggered w.r.t. variable
      func = sfFDDZ;
      table = FluxStagTable;
      diffloc = CELL_CENTRE;
    } else if ((vloc == CELL_CENTRE) && (inloc == CELL_ZLOW)) {
      // Shifted
      func = sfFDDZ;
      table = FluxStagTable;
      diffloc = CELL_ZLOW;
    } else {
      // More complicated. Deciding what to do here isn't straightforward
      throw BoutException("Unhandled shift in indexFDDZ");
    }
  }

  if (method != DIFF_DEFAULT) {
    // Lookup function
    func = lookupFluxFunc(table, method);
  }

  ASSERT1(this == v.getMesh());
  ASSERT1(this == f.getMesh());

  Field3D result(this);
  result.allocate(); // Make sure data allocated

  stencil vval, fval;
  for (const auto &i : result.region(region)) {

    fval.mm = f[i.offset(0,0,-2)];
    fval.m = f[i.zm()];
    fval.c = f[i];
    fval.p = f[i.zp()];
    fval.pp = f[i.offset(0,0,2)];

    vval.mm = v[i.offset(0,0,-2)];
    vval.m = v[i.zm()];
    vval.c = v[i];
    vval.p = v[i.zp()];
    vval.pp = v[i.offset(0,0,2)];

    if(StaggerGrids && (diffloc != CELL_DEFAULT) && (diffloc != vloc)) {
      // Non-centred stencil

      if((vloc == CELL_CENTRE) && (diffloc == CELL_ZLOW)) {
      // Producing a stencil centred around a lower Z value
        vval.pp = vval.p;
        vval.p  = vval.c;

      }else if(vloc == CELL_ZLOW) {
        // Stencil centred around a cell centre

        vval.mm = vval.m;
        vval.m  = vval.c;
      }
      // Shifted in one direction -> shift in another
      // Could produce warning
    }
    result[i] = func(vval, fval);
  }

  result.setLocation(diffloc);

#if CHECK > 0
  // Mark boundaries as invalid
  result.bndry_xin = result.bndry_xout = result.bndry_yup = result.bndry_ydown = false;
#endif

  return result;
}
