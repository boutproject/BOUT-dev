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

  stencil sp, sm;

  sp.mm = f.mm + ma;
  sp.m = f.m + ma;
  sp.c = f.c + ma;
  sp.p = f.p + ma;
  sp.pp = f.pp + ma;

  sm.mm = ma - f.mm;
  sm.m = ma - f.m;
  sm.c = ma - f.c;
  sm.p = ma - f.p;
  sm.pp = ma - f.pp;

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
  // Calculate d(v*f)/dx = (v*f)[i+1/2] - (v*f)[i-1/2]

  // Upper cell boundary
  BoutReal result = (v.p >= 0.) ? v.p * (1.5*f.c - 0.5*f.m) : v.p * (1.5*f.p - 0.5*f.pp);

  // Lower cell boundary
  result -= (v.m >= 0.) ? v.m * (1.5*f.m - 0.5*f.mm) : v.m * (1.5*f.c - 0.5*f.p);

  // result is now d/dx(v*f), but want v*d/dx(f) so subtract f*d/dx(v)
  result -= f.c * (v.p - v.m);

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
static DiffLookup FirstDerivTable[] = {
    {DIFF_C2, DDX_C2, nullptr, nullptr},     {DIFF_W2, DDX_CWENO2, nullptr, nullptr},
    {DIFF_W3, DDX_CWENO3, nullptr, nullptr}, {DIFF_C4, DDX_C4, nullptr, nullptr},
    {DIFF_S2, DDX_S2, nullptr, nullptr},     {DIFF_FFT, nullptr, nullptr, nullptr},
    {DIFF_DEFAULT, nullptr, nullptr, nullptr}};

/// Second derivative lookup table
static DiffLookup SecondDerivTable[] = {{DIFF_C2, D2DX2_C2, nullptr, nullptr},
                                        {DIFF_C4, D2DX2_C4, nullptr, nullptr},
                                        {DIFF_FFT, nullptr, nullptr, nullptr},
                                        {DIFF_DEFAULT, nullptr, nullptr, nullptr}};

/// Upwinding functions lookup table
static DiffLookup UpwindTable[] = {
    {DIFF_U1, nullptr, VDDX_U1, nullptr},    {DIFF_U2, nullptr, VDDX_U2, nullptr},
    {DIFF_C2, nullptr, VDDX_C2, nullptr},    {DIFF_U3, nullptr, VDDX_U3, nullptr},
    {DIFF_W3, nullptr, VDDX_WENO3, nullptr}, {DIFF_C4, nullptr, VDDX_C4, nullptr},
    {DIFF_DEFAULT, nullptr, nullptr, nullptr}};

/// Flux functions lookup table
static DiffLookup FluxTable[] = {
    {DIFF_SPLIT, nullptr, nullptr, nullptr},   {DIFF_U1, nullptr, nullptr, FDDX_U1},
    {DIFF_C2, nullptr, nullptr, FDDX_C2},   {DIFF_C4, nullptr, nullptr, FDDX_C4},
    {DIFF_DEFAULT, nullptr, nullptr, nullptr}};

/// First staggered derivative lookup
static DiffLookup FirstStagDerivTable[] = {{DIFF_C2, DDX_C2_stag, nullptr, nullptr},
                                           {DIFF_C4, DDX_C4_stag, nullptr, nullptr},
                                           {DIFF_DEFAULT, nullptr, nullptr, nullptr}};

/// Second staggered derivative lookup
static DiffLookup SecondStagDerivTable[] = {{DIFF_C2, D2DX2_C2_stag, nullptr, nullptr},
                                            {DIFF_DEFAULT, nullptr, nullptr, nullptr}};

/// Upwinding staggered lookup
static DiffLookup UpwindStagTable[] = {{DIFF_U1, nullptr, nullptr, VDDX_U1_stag},
                                       {DIFF_U2, nullptr, nullptr, VDDX_U2_stag},
                                       {DIFF_C2, nullptr, nullptr, VDDX_C2_stag},
                                       {DIFF_C4, nullptr, nullptr, VDDX_C4_stag},
                                       {DIFF_DEFAULT, nullptr, nullptr, nullptr}};

/// Flux staggered lookup
static DiffLookup FluxStagTable[] = {{DIFF_SPLIT, nullptr, nullptr, nullptr},
                                     {DIFF_U1, nullptr, nullptr, FDDX_U1_stag},
                                     {DIFF_DEFAULT, nullptr, nullptr, nullptr}};

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

const Field2D Mesh::applyXdiff(const Field2D &var, Mesh::deriv_func func,
                               CELL_LOC outloc, REGION region) {
  ASSERT1(this == var.getMesh());
  ASSERT1(var.isAllocated());
  CELL_LOC inloc = var.getLocation();
  if (outloc == CELL_DEFAULT)
    outloc = inloc;
  // Allowed staggers:
  ASSERT1(outloc == inloc || (outloc == CELL_CENTRE && inloc == CELL_XLOW) ||
          (outloc == CELL_XLOW && inloc == CELL_CENTRE));

  if (var.getNx() == 1) {
    auto tmp = Field2D(0., this);
    tmp.setLocation(var.getLocation());
    return tmp;
  }

  /// Convert REGION enum to a Region string identifier
  const auto region_str = REGION_STRING(region);
  
  Field2D result(this);
  result.allocate(); // Make sure data allocated

  if (this->StaggerGrids && (outloc != inloc)) {
    // Staggered differencing

    if (this->xstart > 1) {
      // More than one guard cell, so set pp and mm values
      // This allows higher-order methods to be used
      BOUT_OMP(parallel)
      {
        stencil s;
        BOUT_FOR_INNER(i, this->getRegion2D(region_str)) {
          s.mm = var[i.xmm()];
          s.m = var[i.xm()];
          s.c = var[i];
          s.p = var[i.xp()];
          s.pp = var[i.xpp()];

          if (outloc == CELL_XLOW) {
            // Producing a stencil centred around a lower X value
            s.pp = s.p;
            s.p = s.c;
          } else {
            // Stencil centred around a cell centre
            s.mm = s.m;
            s.m = s.c;
          }

          result[i] = func(s);
        }
      }
    } else {
      // Only one guard cell, so no pp or mm values
      BOUT_OMP(parallel)
      {
        stencil s;
        BOUT_FOR_INNER(i, this->getRegion2D(region_str)) {
          s.m = var[i.xm()];
          s.c = var[i];
          s.p = var[i.xp()];

          if (outloc == CELL_XLOW) {
            // Producing a stencil centred around a lower X value
            s.pp = s.p;
            s.p = s.c;
          } else {
            // Stencil centred around a cell centre
            s.mm = s.m;
            s.m = s.c;
          }

          result[i] = func(s);
        }
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
        BOUT_FOR_INNER(i, this->getRegion2D(region_str)) {
          s.mm = var[i.xmm()];
          s.m = var[i.xm()];
          s.c = var[i];
          s.p = var[i.xp()];
          s.pp = var[i.xpp()];

          result[i] = func(s);
        }
      }
    } else {
      // Only one guard cell, so no pp or mm values
      BOUT_OMP(parallel)
      {
        stencil s;
        BOUT_FOR_INNER(i, this->getRegion2D(region_str)) {
          s.m = var[i.xm()];
          s.c = var[i];
          s.p = var[i.xp()];

          result[i] = func(s);
        }
      }
    }
  }

  result.setLocation(outloc);

#if CHECK > 0
  // Mark boundaries as invalid
  result.bndry_xin = result.bndry_xout = result.bndry_yup = result.bndry_ydown = false;
#endif

  return result;
}

const Field3D Mesh::applyXdiff(const Field3D &var, Mesh::deriv_func func,
                               CELL_LOC outloc, REGION region) {
  // Check that the mesh is correct
  ASSERT1(this == var.getMesh());
  // Check that the input variable has data
  ASSERT1(var.isAllocated());

  CELL_LOC inloc = var.getLocation();
  if (outloc == CELL_DEFAULT)
    outloc = inloc;
  // Allowed staggers:
  ASSERT1(outloc == inloc || (outloc == CELL_CENTRE && inloc == CELL_XLOW) ||
          (outloc == CELL_XLOW && inloc == CELL_CENTRE));

  if (var.getNx() == 1) {
    auto tmp = Field3D(0., this);
    tmp.setLocation(var.getLocation());
    return tmp;
  }

  /// Convert REGION enum to a Region string identifier
  const auto region_str = REGION_STRING(region);

  Field3D result(this);
  result.allocate(); // Make sure data allocated

  if (this->StaggerGrids && (outloc != inloc)) {
    // Staggered differencing

    if (this->xstart > 1) {
      // More than one guard cell, so set pp and mm values
      // This allows higher-order methods to be used
      BOUT_OMP(parallel)
      {
        stencil s;
        BOUT_FOR_INNER(i, this->getRegion3D(region_str)) {
          s.mm = var[i.xmm()];
          s.m = var[i.xm()];
          s.c = var[i];
          s.p = var[i.xp()];
          s.pp = var[i.xpp()];

          if ((inloc == CELL_CENTRE) && (outloc == CELL_XLOW)) {
            // Producing a stencil centred around a lower X value
            s.pp = s.p;
            s.p = s.c;
          } else if (inloc == CELL_XLOW) {
            // Stencil centred around a cell centre
            s.mm = s.m;
            s.m = s.c;
          }

          result[i] = func(s);
        }
      }
    } else {
      // Only one guard cell, so no pp or mm values
      BOUT_OMP(parallel)
      {
        stencil s;
        BOUT_FOR_INNER(i, this->getRegion3D(region_str)) {
          s.m = var[i.xm()];
          s.c = var[i];
          s.p = var[i.xp()];

          if (outloc == CELL_XLOW) {
            // Producing a stencil centred around a lower X value
            s.pp = s.p;
            s.p = s.c;
          } else if (inloc == CELL_XLOW) {
            // Stencil centred around a cell centre
            s.mm = s.m;
            s.m = s.c;
          }

          result[i] = func(s);
        }
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
        BOUT_FOR_INNER(i, this->getRegion3D(region_str)) {
          s.mm = var[i.xmm()];
          s.m = var[i.xm()];
          s.c = var[i];
          s.p = var[i.xp()];
          s.pp = var[i.xpp()];

          result[i] = func(s);
        }
      }
    } else {
      // Only one guard cell, so no pp or mm values
      BOUT_OMP(parallel)
      {
        stencil s;
        BOUT_FOR_INNER(i, this->getRegion3D(region_str)) {
          s.m = var[i.xm()];
          s.c = var[i];
          s.p = var[i.xp()];

          result[i] = func(s);
        }
      }
    }
  }

  result.setLocation(outloc);

#if CHECK > 0
  // Mark boundaries as invalid
  result.bndry_xin = result.bndry_xout = result.bndry_yup = result.bndry_ydown = false;
#endif

  return result;
}

// Y derivative

const Field2D Mesh::applyYdiff(const Field2D &var, Mesh::deriv_func func, CELL_LOC outloc,
                               REGION region) {
  ASSERT1(this == var.getMesh());
  // Check that the input variable has data
  ASSERT1(var.isAllocated());

  CELL_LOC inloc = var.getLocation();
  if (outloc == CELL_DEFAULT)
    outloc = inloc;
  // Allowed staggers:
  ASSERT1(outloc == inloc);

  if (var.getNy() == 1) {
    auto tmp = Field2D(0., this);
    tmp.setLocation(var.getLocation());
    return tmp;
  }

  /// Convert REGION enum to a Region string identifier
  const auto region_str = REGION_STRING(region);
  
  Field2D result(this);
  result.allocate(); // Make sure data allocated

  if (this->ystart > 1) {
    // More than one guard cell, so set pp and mm values
    // This allows higher-order methods to be used
    BOUT_OMP(parallel)
    {
      stencil s;
      BOUT_FOR_INNER(i, this->getRegion2D(region_str)) {
        s.mm = var[i.ymm()];
        s.m = var[i.ym()];
        s.c = var[i];
        s.p = var[i.yp()];
        s.pp = var[i.ypp()];

        result[i] = func(s);
      }
    }
  } else {
    BOUT_OMP(parallel)
    {
      stencil s;
      BOUT_FOR_INNER(i, this->getRegion2D(region_str)) {
        s.m = var[i.ym()];
        s.c = var[i];
        s.p = var[i.yp()];

        result[i] = func(s);
      }
    }
  }

  result.setLocation(outloc);

#if CHECK > 0
  // Mark boundaries as invalid
  result.bndry_yup = result.bndry_ydown = false;
#endif

  return result;
}

const Field3D Mesh::applyYdiff(const Field3D &var, Mesh::deriv_func func, CELL_LOC outloc,
                               REGION region) {
  ASSERT1(this == var.getMesh());
  // Check that the input variable has data
  ASSERT1(var.isAllocated());
  // Cell location of the input field
  CELL_LOC inloc = var.getLocation();
  if (outloc == CELL_DEFAULT)
    outloc = inloc;
  // Allowed staggers:
  ASSERT1(outloc == inloc || (outloc == CELL_CENTRE && inloc == CELL_YLOW) ||
          (outloc == CELL_YLOW && inloc == CELL_CENTRE));

  if (var.getNy() == 1) {
    auto tmp = Field3D(0., this);
    tmp.setLocation(var.getLocation());
    return tmp;
  }

  /// Convert REGION enum to a Region string identifier
  const auto region_str = REGION_STRING(region);
  
  Field3D result(this);
  result.allocate(); // Make sure data allocated

  if (var.hasYupYdown() && ((&var.yup() != &var) || (&var.ydown() != &var))) {
    // Field "var" has distinct yup and ydown fields which
    // will be used to calculate a derivative along
    // the magnetic field

    if (this->StaggerGrids && (outloc != inloc)) {
      // Staggered differencing

      BOUT_OMP(parallel)
      {
        stencil s;
        BOUT_FOR_INNER(i, this->getRegion3D(region_str)) {

          // Set stencils
          s.m = var.ydown()[i.ym()];
          s.c = var[i];
          s.p = var.yup()[i.yp()];

          if (outloc == CELL_YLOW) {
            // Producing a stencil centred around a lower Y value
            s.pp = s.p;
            s.p = s.c;
          } else if (inloc == CELL_YLOW) {
            // Stencil centred around a cell centre
            s.mm = s.m;
            s.m = s.c;
          }

          result[i] = func(s);
        }
      }
    } else {
      // Non-staggered
      BOUT_OMP(parallel)
      {
        stencil s;
        BOUT_FOR_INNER(i, this->getRegion3D(region_str)) {
          // Set stencils
          s.m = var.ydown()[i.ym()];
          s.c = var[i];
          s.p = var.yup()[i.yp()];

          result[i] = func(s);
        }
      }
    }
  } else {
    // var has no yup/ydown fields, so we need to shift into field-aligned coordinates

    Field3D var_fa = this->toFieldAligned(var);

    if (this->StaggerGrids && (outloc != inloc)) {
      // Staggered differencing

      if (this->ystart > 1) {
        // More than one guard cell, so set pp and mm values
        // This allows higher-order methods to be used
        BOUT_OMP(parallel) {
          stencil s;
          BOUT_FOR_INNER(i, this->getRegion3D(region_str)) {
            // Set stencils
            s.mm = var_fa[i.ymm()];
            s.m = var_fa[i.ym()];
            s.c = var_fa[i];
            s.p = var_fa[i.yp()];
            s.pp = var_fa[i.ypp()];

            if (outloc == CELL_YLOW) {
              // Producing a stencil centred around a lower Y value
              s.pp = s.p;
              s.p = s.c;
            } else if (inloc == CELL_YLOW) {
              // Stencil centred around a cell centre
              s.mm = s.m;
              s.m = s.c;
            }

            result[i] = func(s);
          }
        }
      } else {
        // Only one guard cell, so no pp or mm values
        BOUT_OMP(parallel) {
          stencil s;
          BOUT_FOR_INNER(i, this->getRegion3D(region_str)) {
            // Set stencils
            s.m = var_fa[i.ym()];
            s.c = var_fa[i];
            s.p = var_fa[i.yp()];

            if (outloc == CELL_YLOW) {
              // Producing a stencil centred around a lower Y value
              s.pp = s.p;
              s.p = s.c;
            } else if (inloc == CELL_YLOW) {
              // Stencil centred around a cell centre
              s.mm = s.m;
              s.m = s.c;
            }

            result[i] = func(s);
          }
        }
      }
    } else {
      // Non-staggered differencing

      if (this->ystart > 1) {
        // More than one guard cell, so set pp and mm values
        // This allows higher-order methods to be used
        BOUT_OMP(parallel) {
          stencil s;
          BOUT_FOR_INNER(i, this->getRegion3D(region_str)) {
            // Set stencils
            s.mm = var_fa[i.ymm()];
            s.m = var_fa[i.ym()];
            s.c = var_fa[i];
            s.p = var_fa[i.yp()];
            s.pp = var_fa[i.ypp()];

            result[i] = func(s);
          }
        }
      } else {
        // Only one guard cell, so no pp or mm values
        BOUT_OMP(parallel) {
          stencil s;
          BOUT_FOR_INNER(i, this->getRegion3D(region_str)) {
            // Set stencils
            s.m = var_fa[i.ym()];
            s.c = var_fa[i];
            s.p = var_fa[i.yp()];

            result[i] = func(s);
          }
        }
      }
    }

    // Shift result back

    result = this->fromFieldAligned(result);
  }

  result.setLocation(outloc);

#if CHECK > 0
  // Mark boundaries as invalid
  result.bndry_xin = result.bndry_xout = result.bndry_yup = result.bndry_ydown = false;
#endif

  return result;
}

// Z derivative

const Field3D Mesh::applyZdiff(const Field3D &var, Mesh::deriv_func func, CELL_LOC outloc,
                               REGION region) {
  ASSERT1(this == var.getMesh());
  // Check that the input variable has data
  ASSERT1(var.isAllocated());
  CELL_LOC inloc = var.getLocation();
  // Allowed staggers:
  if (outloc == CELL_DEFAULT)
    outloc = inloc;
  ASSERT1(outloc == inloc);

  if (var.getNz() == 1) {
    auto tmp = Field3D(0., this);
    tmp.setLocation(var.getLocation());
    return tmp;
  }

  /// Convert REGION enum to a Region string identifier
  const auto region_str = REGION_STRING(region);

  Field3D result(this);
  result.allocate(); // Make sure data allocated

  // Check that the input variable has data
  ASSERT1(var.isAllocated());

  BOUT_OMP(parallel)
  {
    stencil s;
    BOUT_FOR_INNER(i, this->getRegion3D(region_str)) {
      s.mm = var[i.zmm()];
      s.m = var[i.zm()];
      s.c = var[i];
      s.p = var[i.zp()];
      s.pp = var[i.zpp()];

      result[i] = func(s);
    }
  }

  result.setLocation(outloc);

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
  // Allowed staggers:
  if (outloc == CELL_DEFAULT)
    outloc = inloc;
  ASSERT1(outloc == inloc || (outloc == CELL_CENTRE && inloc == CELL_XLOW) ||
          (outloc == CELL_XLOW && inloc == CELL_CENTRE));

  Field3D result(this);

  if (this->StaggerGrids && (outloc != inloc)) {
    // Shifting in X. Centre -> Xlow, or Xlow -> Centre

    func = sfDDX;                // Set default
    table = FirstStagDerivTable; // Set table for others
  }

  if (method != DIFF_DEFAULT) {
    // Lookup function
    func = lookupFunc(table, method);
    if (func == nullptr)
      throw BoutException("Cannot use FFT for X derivatives");
  }

  result = applyXdiff(f, func, outloc, region);

  return result;
}

const Field2D Mesh::indexDDX(const Field2D &f, CELL_LOC outloc,
                             DIFF_METHOD method, REGION region) {
  ASSERT1(outloc == CELL_DEFAULT || outloc == f.getLocation());
  ASSERT1(method == DIFF_DEFAULT);
  return applyXdiff(f, fDDX, f.getLocation(), region);
}

////////////// Y DERIVATIVE /////////////////

const Field3D Mesh::indexDDY(const Field3D &f, CELL_LOC outloc, DIFF_METHOD method,
                             REGION region) {
  Mesh::deriv_func func = fDDY; // Set to default function
  DiffLookup *table = FirstDerivTable;

  CELL_LOC inloc = f.getLocation(); // Input location
  // Allowed staggers:
  if (outloc == CELL_DEFAULT)
    outloc = inloc;
  ASSERT1(outloc == inloc || (outloc == CELL_CENTRE && inloc == CELL_YLOW) ||
          (outloc == CELL_YLOW && inloc == CELL_CENTRE));

  Field3D result(this);

  if (this->StaggerGrids && (outloc != inloc)) {
    // Shifting in Y. Centre -> Ylow, or Ylow -> Centre
    func = sfDDY;                // Set default
    table = FirstStagDerivTable; // Set table for others
  }

  if (method != DIFF_DEFAULT) {
    // Lookup function
    func = lookupFunc(table, method);
    if (func == nullptr)
      throw BoutException("Cannot use FFT for Y derivatives");
  }

  result = applyYdiff(f, func, outloc, region);

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
  Mesh::deriv_func func = fDDZ; // Set to default function
  DiffLookup *table = FirstDerivTable;

  CELL_LOC inloc = f.getLocation(); // Input location
  if (outloc == CELL_DEFAULT)
    outloc = inloc;
  // Allowed staggers:
  ASSERT1(outloc == inloc || (outloc == CELL_CENTRE && inloc == CELL_ZLOW) ||
          (outloc == CELL_ZLOW && inloc == CELL_CENTRE));

  Field3D result(this);

  if (this->StaggerGrids && (outloc != inloc)) {
    // Shifting in Z. Centre -> Zlow, or Zlow -> Centre
    func = sfDDZ;                // Set default
    table = FirstStagDerivTable; // Set table for others
  }

  if (method != DIFF_DEFAULT) {
    // Lookup function
    func = lookupFunc(table, method);
  }

  if (func == nullptr) {
    // Use FFT

    BoutReal shift = 0.; // Shifting result in Z?
    if (this->StaggerGrids && (outloc != inloc)) {
      if (outloc == CELL_ZLOW) {
        // Shifting down - multiply by exp(-0.5*i*k*dz)
        shift = -1.;
        throw BoutException("Not tested - probably broken");
      } else {
        // Shifting up
        shift = 1.;
        throw BoutException("Not tested - probably broken");
      }
    }

    result.allocate(); // Make sure data allocated

    // Calculate how many Z wavenumbers will be removed
    const int ncz = this->LocalNz;
    int kfilter =
        static_cast<int>(fft_derivs_filter * ncz / 2); // truncates, rounding down
    if (kfilter < 0)
      kfilter = 0;
    if (kfilter > (ncz / 2))
      kfilter = ncz / 2;
    const int kmax = ncz / 2 - kfilter; // Up to and including this wavenumber index

    const auto region_str = REGION_STRING(region);

    // Only allow a whitelist of regions for now
    ASSERT2(region_str == "RGN_ALL" || region_str == "RGN_NOBNDRY" ||
            region_str == "RGN_NOX" || region_str == "RGN_NOY");

    BOUT_OMP(parallel)
    {
      Array<dcomplex> cv(ncz / 2 + 1);
      const BoutReal kwaveFac = TWOPI / ncz;

      // Note we lookup a 2D region here even though we're operating on a Field3D
      // as we only want to loop over {x, y} and then handle z differently. The
      // Region<Ind2D> blocks are constructed for elements contiguous assuming nz=1,
      // as that isn't the case for Field3D (in general) this shouldn't be expected
      // to vectorise (not that it would anyway) but it should still OpenMP parallelise
      // ok.
      // With this in mind we could perhaps avoid the use of the BOUT_FOR_INNER macro
      // here,
      // but should be ok for now.
      BOUT_FOR_INNER(i, mesh->getRegion2D(region_str)) {
        auto i3D = mesh->ind2Dto3D(i, 0);
        rfft(&f[i3D], ncz, cv.begin()); // Forward FFT

        for (int jz = 0; jz <= kmax; jz++) {
          const BoutReal kwave = jz * kwaveFac; // wave number is 1/[rad]

          cv[jz] *= dcomplex(0, kwave);
          if (shift)
            cv[jz] *= exp(Im * (shift * kwave));
        }
        for (int jz = kmax + 1; jz <= ncz / 2; jz++) {
          cv[jz] = 0.0;
        }

        irfft(cv.begin(), ncz, &result[i3D]); // Reverse FFT
      }
    }

#if CHECK > 0
    // Mark boundaries as invalid
    result.bndry_xin = false;
    result.bndry_xout = false;
    result.bndry_yup = false;
    result.bndry_ydown = false;
#endif

    result.setLocation(outloc);

  } else {
    // All other (non-FFT) functions
    result = applyZdiff(f, func, outloc, region);
  }

  return result;
}

const Field2D Mesh::indexDDZ(const Field2D &f, CELL_LOC UNUSED(outloc),
                             DIFF_METHOD UNUSED(method), REGION UNUSED(region)) {
  ASSERT1(this == f.getMesh());
  auto tmp = Field2D(0., this);
  tmp.setLocation(f.getLocation());
  return tmp;
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
  if (outloc == CELL_DEFAULT)
    outloc = inloc;
  // Allowed staggers:
  ASSERT1(outloc == inloc || (outloc == CELL_CENTRE && inloc == CELL_XLOW) ||
          (outloc == CELL_XLOW && inloc == CELL_CENTRE));

  ASSERT1(this == f.getMesh());

  Field3D result(this);

  if (StaggerGrids && (outloc != inloc)) {
    // Shifting in X. Centre -> Xlow, or Xlow -> Centre
    func = sfD2DX2;               // Set default
    table = SecondStagDerivTable; // Set table for others
  }

  if (method != DIFF_DEFAULT) {
    // Lookup function
    func = lookupFunc(table, method);
    if (func == nullptr)
      throw BoutException("Cannot use FFT for X derivatives");
  }

  result = applyXdiff(f, func, outloc, region);

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

  ASSERT1(this == f.getMesh());

  CELL_LOC inloc = f.getLocation(); // Input location
  // Allowed staggers:
  if (outloc == CELL_DEFAULT)
    outloc = inloc;
  ASSERT1(outloc == inloc || (outloc == CELL_CENTRE && inloc == CELL_YLOW) ||
          (outloc == CELL_YLOW && inloc == CELL_CENTRE));

  Field3D result(this);

  if (StaggerGrids && (outloc != inloc)) {
    // Shifting in Y. Centre -> Ylow, or Ylow -> Centre
    func = sfD2DY2;               // Set default
    table = SecondStagDerivTable; // Set table for others
  }

  if (method != DIFF_DEFAULT) {
    // Lookup function
    func = lookupFunc(table, method);
    if (func == nullptr)
      throw BoutException("Cannot use FFT for Y derivatives");
  }

  result = applyYdiff(f, func, outloc, region);

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

  ASSERT1(this == f.getMesh());

  CELL_LOC inloc = f.getLocation(); // Input location
  // Allowed staggers:
  if (outloc == CELL_DEFAULT)
    outloc = inloc;
  ASSERT1(outloc == inloc || (outloc == CELL_CENTRE && inloc == CELL_ZLOW) ||
          (outloc == CELL_ZLOW && inloc == CELL_CENTRE));

  Field3D result(this);

  if (StaggerGrids && (outloc != inloc)) {
    // Shifting in Z. Centre -> Zlow, or Zlow -> Centre
    func = sfD2DZ2;               // Set default
    table = SecondStagDerivTable; // Set table for others
  }

  if (method != DIFF_DEFAULT) {
    // Lookup function
    func = lookupFunc(table, method);
  }

  if (func == nullptr) {
    // Use FFT

    BoutReal shift = 0.; // Shifting result in Z?
    if (StaggerGrids && (outloc != inloc)) {
      if (outloc == CELL_ZLOW) {
        // Shifting down - multiply by exp(-0.5*i*k*dz)
        throw BoutException("Not tested - probably broken");
      } else {
        // Shifting up
        throw BoutException("Not tested - probably broken");
      }
    }

    result.allocate(); // Make sure data allocated

    // No filtering in 2nd derivative method
    const int ncz = this->LocalNz;
    const int kmax = ncz / 2; // Up to and including this wavenumber index

    const auto region_str = REGION_STRING(region);

    // Only allow a whitelist of regions for now
    ASSERT2(region_str == "RGN_ALL" || region_str == "RGN_NOBNDRY" ||
            region_str == "RGN_NOX" || region_str == "RGN_NOY");

    BOUT_OMP(parallel) {
      Array<dcomplex> cv(ncz / 2 + 1);
      const BoutReal kwaveFac = TWOPI / ncz;

      // Note we lookup a 2D region here even though we're operating on a Field3D
      // as we only want to loop over {x, y} and then handle z differently. The
      // Region<Ind2D> blocks are constructed for elements contiguous assuming nz=1,
      // as that isn't the case for Field3D (in general) this shouldn't be expected
      // to vectorise (not that it would anyway) but it should still OpenMP parallelise
      // ok.
      // With this in mind we could perhaps avoid the use of the BOUT_FOR_INNER macro
      // here,
      // but should be ok for now.
      BOUT_FOR_INNER(i, mesh->getRegion2D(region_str)) {
        auto i3D = mesh->ind2Dto3D(i, 0);

        rfft(&f[i3D], ncz, cv.begin()); // Forward FFT

        for (int jz = 0; jz <= kmax; jz++) {
          const BoutReal kwave = jz * kwaveFac; // wave number is 1/[rad]

          cv[jz] *= -kwave * kwave;
          if (shift)
            cv[jz] *= exp(0.5 * Im * (shift * kwave));
        }
        for (int jz = kmax + 1; jz <= ncz / 2; jz++) {
          cv[jz] = 0.0;
        }

        irfft(cv.begin(), ncz, &result[i3D]); // Reverse FFT
      }
    }

#if CHECK > 0
    // Mark boundaries as invalid
    result.bndry_xin = false;
    result.bndry_xout = false;
    result.bndry_yup = false;
    result.bndry_ydown = false;
#endif

    result.setLocation(outloc);

  } else {
    // All other (non-FFT) functions
    result = applyZdiff(f, func, outloc, region);
  }

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
                               DIFF_METHOD UNUSED(method), REGION UNUSED(region)) {
  ASSERT1(outloc == CELL_DEFAULT || outloc == f.getLocation());
  auto tmp = Field2D(0., this);
  tmp.setLocation(f.getLocation());
  return tmp;
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

  CELL_LOC inloc = f.getLocation();
  CELL_LOC vloc = v.getLocation();
  if (outloc == CELL_DEFAULT)
    outloc = inloc;
  // Allowed staggers:
  ASSERT1(outloc == inloc && inloc == vloc);

  Mesh::upwind_func func = fVDDX;

  if (method != DIFF_DEFAULT) {
    // Lookup function
    func = lookupFunc(UpwindTable, method);
  }

  ASSERT1(this->xstart > 0); // Need at least one guard cell
  ASSERT1(this == f.getMesh());
  ASSERT1(this == v.getMesh());

  /// Convert REGION enum to a Region string identifier
  const auto region_str = REGION_STRING(region);

  Field2D result(this);
  result.allocate(); // Make sure data allocated

  if (this->xstart > 1) {
    // Two or more guard cells
    BOUT_OMP(parallel) {
      stencil s;
      BOUT_FOR_INNER(i, this->getRegion2D(region_str)) {
        s.mm = f[i.xmm()];
        s.m = f[i.xm()];
        s.c = f[i];
        s.p = f[i.xp()];
        s.pp = f[i.xpp()];

        result[i] = func(v[i], s);
      }
    }

  } else {
    // Only one guard cell
    BOUT_OMP(parallel) {
      stencil s;
      BOUT_FOR_INNER(i, this->getRegion2D(region_str)) {
        s.m = f[i.xm()];
        s.c = f[i];
        s.p = f[i.xp()];

        result[i] = func(v[i], s);
      }
    }
  }

#if CHECK > 0
  // Mark boundaries as invalid
  result.bndry_xin = result.bndry_xout = false;
#endif

  result.setLocation(outloc);

  return result;
}

/// General version for 3D objects.
/// 2D objects passed as input will result in copying
const Field3D Mesh::indexVDDX(const Field3D &v, const Field3D &f, CELL_LOC outloc,
                              DIFF_METHOD method, REGION region) {
  TRACE("Mesh::indexVDDX(Field3D, Field3D)");

  ASSERT1(this->xstart > 0); // Need at least one guard cell
  ASSERT1(this == v.getMesh());
  ASSERT1(this == f.getMesh());

  CELL_LOC vloc = v.getLocation();
  CELL_LOC inloc = f.getLocation(); // Input location
  if (outloc == CELL_DEFAULT)
    outloc = inloc;
  // Allowed staggers:
  ASSERT1(outloc == inloc &&
          ((vloc == inloc) || (vloc == CELL_CENTRE && inloc == CELL_XLOW) ||
           (vloc == CELL_XLOW && inloc == CELL_CENTRE)));

  Field3D result(this);
  result.allocate(); // Make sure data allocated

  /// Convert REGION enum to a Region string identifier
  const auto region_str = REGION_STRING(region);

  if (StaggerGrids && (vloc != inloc)) {
    // Staggered grids enabled, and velocity at different location to value

    Mesh::flux_func func = sfVDDX;
    DiffLookup *table = UpwindTable;

    // V staggered w.r.t. variable
    func = sfVDDX;
    table = UpwindStagTable;

    if (method != DIFF_DEFAULT) {
      // Lookup function
      func = lookupFunc(table, method);
    }

    // Note: The velocity stencil contains only (mm, m, p, pp)
    // v.p is v at +1/2, v.m is at -1/2 relative to the field f

    if (this->xstart > 1) {
      // Two or more guard cells

      if (vloc == CELL_XLOW) {
        BOUT_OMP(parallel) {
          stencil fs, vs;
          BOUT_FOR_INNER(i, this->getRegion3D(region_str)) {
            fs.mm = f[i.xmm()];
            fs.m = f[i.xm()];
            fs.c = f[i];
            fs.p = f[i.xp()];
            fs.pp = f[i.xpp()];

            vs.mm = v[i.xm()];
            vs.m = v[i];
            vs.p = v[i.xp()];
            vs.pp = v[i.xpp()];

            result[i] = func(vs, fs);
          }
        }

      } else {
        BOUT_OMP(parallel) {
          stencil fs, vs;
          BOUT_FOR_INNER(i, this->getRegion3D(region_str)) {
            fs.mm = f[i.xmm()];
            fs.m = f[i.xm()];
            fs.c = f[i];
            fs.p = f[i.xp()];
            fs.pp = f[i.xpp()];

            vs.mm = v[i.xmm()];
            vs.m = v[i.xm()];
            vs.p = v[i];
            vs.pp = v[i.xp()];

            result[i] = func(vs, fs);
          }
        }
      }
    } else {
      // One guard cell

      if (vloc == CELL_XLOW) {
        BOUT_OMP(parallel) {
          stencil fs, vs;
          BOUT_FOR_INNER(i, this->getRegion3D(region_str)) {
            fs.m = f[i.xm()];
            fs.c = f[i];
            fs.p = f[i.xp()];

            vs.mm = v[i.xm()];
            vs.m = v[i];
            vs.p = v[i.xp()];

            result[i] = func(vs, fs);
          }
        }

      } else {
        BOUT_OMP(parallel) {
          stencil fs, vs;
          BOUT_FOR_INNER(i, this->getRegion3D(region_str)) {
            fs.m = f[i.xm()];
            fs.c = f[i];
            fs.p = f[i.xp()];

            vs.m = v[i.xm()];
            vs.p = v[i];
            vs.pp = v[i.xp()];

            result[i] = func(vs, fs);
          }
        }
      }
    }

  } else {
    // Not staggered
    Mesh::upwind_func func = fVDDX;
    DiffLookup *table = UpwindTable;

    if (method != DIFF_DEFAULT) {
      // Lookup function
      func = lookupFunc(table, method);
    }

    if (this->xstart > 1) {
      // Two or more guard cells
      BOUT_OMP(parallel) {
        stencil fs;
        BOUT_FOR_INNER(i, this->getRegion3D(region_str)) {
          fs.mm = f[i.xmm()];
          fs.m = f[i.xm()];
          fs.c = f[i];
          fs.p = f[i.xp()];
          fs.pp = f[i.xpp()];

          result[i] = func(v[i], fs);
        }
      }
    } else {
      // Only one guard cell
      BOUT_OMP(parallel) {
        stencil fs;
        BOUT_FOR_INNER(i, this->getRegion3D(region_str)) {
          fs.m = f[i.xm()];
          fs.c = f[i];
          fs.p = f[i.xp()];

          result[i] = func(v[i], fs);
        }
      }
    }
  }

  result.setLocation(outloc);

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
  TRACE("Mesh::indexVDDY");

  ASSERT1(this == v.getMesh());
  ASSERT1(this == f.getMesh());

  CELL_LOC vloc = v.getLocation();
  CELL_LOC inloc = f.getLocation(); // Input location
  if (outloc == CELL_DEFAULT)
    outloc = inloc;
  // Allowed staggers:
  ASSERT1(outloc == inloc &&
          ((vloc == inloc) || (vloc == CELL_CENTRE && inloc == CELL_YLOW) ||
           (vloc == CELL_YLOW && inloc == CELL_CENTRE)));

  Field2D result(this);
  result.allocate(); // Make sure data allocated

  if (this->LocalNy == 1){
    result=0;
    result.setLocation(outloc);
    return result;
  }

  ASSERT1(this->ystart > 0); // Must have at least one guard cell

  /// Convert REGION enum to a Region string identifier
  const auto region_str = REGION_STRING(region);

  if (StaggerGrids && (vloc != inloc)) {
    // Staggered grids enabled, and velocity at different location to value

    Mesh::flux_func func = sfVDDY;
    DiffLookup *table = UpwindTable;

    // V staggered w.r.t. variable
    func = sfVDDY;
    table = UpwindStagTable;

    if (method != DIFF_DEFAULT) {
      // Lookup function
      func = lookupFunc(table, method);
    }

    // Note: vs.c not used for staggered differencing
    // vs.m is at i-1/2, vs.p is as i+1/2
    if (vloc == CELL_YLOW) {
      if (this->ystart > 1) {
        // Two or more guard cells
        BOUT_OMP(parallel) {
          stencil fs, vs;
          BOUT_FOR_INNER(i, this->getRegion2D(region_str)) {
            fs.mm = f[i.ymm()];
            fs.m = f[i.ym()];
            fs.c = f[i];
            fs.p = f[i.yp()];
            fs.pp = f[i.ypp()];

            vs.mm = v[i.ym()];
            vs.m = v[i];
            vs.p = v[i.yp()];
            vs.pp = v[i.ypp()];

            result[i] = func(vs, fs);
          }
        }
      } else {
        // Only one guard cell
        BOUT_OMP(parallel) {
          stencil fs, vs;
          BOUT_FOR_INNER(i, this->getRegion2D(region_str)) {
            fs.m = f[i.ym()];
            fs.c = f[i];
            fs.p = f[i.yp()];

            vs.mm = v[i.ym()];
            vs.m = v[i];
            vs.p = v[i.yp()];

            result[i] = func(vs, fs);
          }
        }
      }
    } else {
      if (this->ystart > 1) {
        // Two or more guard cells
        BOUT_OMP(parallel) {
          stencil fs, vs;
          BOUT_FOR_INNER(i, this->getRegion2D(region_str)) {
            fs.mm = f[i.ymm()];
            fs.m = f[i.ym()];
            fs.c = f[i];
            fs.p = f[i.yp()];
            fs.pp = f[i.ypp()];

            vs.mm = v[i.ymm()];
            vs.m = v[i.ym()];
            vs.p = v[i];
            vs.pp = v[i.yp()];

            result[i] = func(vs, fs);
          }
        }
      } else {
        // Only one guard cell
        BOUT_OMP(parallel) {
          stencil fs, vs;
          BOUT_FOR_INNER(i, this->getRegion2D(region_str)) {
            fs.m = f[i.ym()];
            fs.c = f[i];
            fs.p = f[i.yp()];

            vs.m = v[i.ym()];
            vs.p = v[i];
            vs.pp = v[i.yp()];

            result[i] = func(vs, fs);
          }
        }
      }
    }

  } else {
    // Not staggered

    Mesh::upwind_func func = fVDDY;
    DiffLookup *table = UpwindTable;

    if (method != DIFF_DEFAULT) {
      // Lookup function
      func = lookupFunc(table, method);
    }

    if (this->ystart > 1) {
      // Two or more guard cells
      BOUT_OMP(parallel) {
        stencil fs;
        BOUT_FOR_INNER(i, this->getRegion2D(region_str)) {
          fs.mm = f[i.ymm()];
          fs.m = f[i.ym()];
          fs.c = f[i];
          fs.p = f[i.yp()];
          fs.pp = f[i.ypp()];

          result[i] = func(v[i], fs);
        }
      }
    } else {
      // Only one guard cell
      BOUT_OMP(parallel) {
        stencil fs;
        BOUT_FOR_INNER(i, this->getRegion2D(region_str)) {
          fs.m = f[i.ym()];
          fs.c = f[i];
          fs.p = f[i.yp()];

          result[i] = func(v[i], fs);
        }
      }
    }
  }

  result.setLocation(outloc);

#if CHECK > 0
  // Mark boundaries as invalid
  result.bndry_xin = result.bndry_xout = result.bndry_yup = result.bndry_ydown = false;
#endif

  return result;
}

// general case
const Field3D Mesh::indexVDDY(const Field3D &v, const Field3D &f, CELL_LOC outloc,
                              DIFF_METHOD method, REGION region) {
  TRACE("Mesh::indexVDDY(Field3D, Field3D)");

  ASSERT1(this == v.getMesh());
  ASSERT1(this == f.getMesh());

  CELL_LOC vloc = v.getLocation();
  CELL_LOC inloc = f.getLocation(); // Input location
  if (outloc == CELL_DEFAULT)
    outloc = inloc;
  // Allowed staggers:
  ASSERT1(outloc == inloc &&
          ((vloc == inloc) || (vloc == CELL_CENTRE && inloc == CELL_YLOW) ||
           (vloc == CELL_YLOW && inloc == CELL_CENTRE)));

  Field3D result(this);
  result.allocate(); // Make sure data allocated

  if (this->LocalNy == 1){
    result=0;
    result.setLocation(outloc);
    return result;
  }

  ASSERT1(this->ystart > 0); // Need at least one guard cell

  /// Convert REGION enum to a Region string identifier
  const auto region_str = REGION_STRING(region);

  if (StaggerGrids && (vloc != inloc)) {
    // Staggered grids enabled, and velocity at different location to value

    Mesh::flux_func func = sfVDDY;
    DiffLookup *table = UpwindTable;

    // V staggered w.r.t. variable
    func = sfVDDY;
    table = UpwindStagTable;

    if (method != DIFF_DEFAULT) {
      // Lookup function
      func = lookupFunc(table, method);
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
      BOUT_OMP(parallel) {
        stencil vval, fval;
        BOUT_FOR_INNER(i, this->getRegion3D(region_str)) {
          vval.m = v.ydown()[i.ym()];
          vval.c = v[i];
          vval.p = v.yup()[i.yp()];

          fval.m = f.ydown()[i.ym()];
          fval.c = f[i];
          fval.p = f.yup()[i.yp()];

          // Non-centred stencil
          if (inloc == CELL_YLOW) {
            // Producing a stencil centred around a lower Y value
            vval.pp = vval.p;
            vval.p = vval.c;
          } else {
            // Stencil centred around a cell centre
            vval.mm = vval.m;
            vval.m = vval.c;
          }
          result[i] = func(vval, fval);
        }
      }
    } else {
      // Both must shift to field aligned
      Field3D v_fa = this->toFieldAligned(v);
      Field3D f_fa = this->toFieldAligned(f);
      BOUT_OMP(parallel) {
        stencil vval, fval;
        BOUT_FOR_INNER(i, this->getRegion3D(region_str)) {
          vval.mm = v_fa[i.ymm()];
          vval.m = v_fa[i.ym()];
          vval.c = v_fa[i];
          vval.p = v_fa[i.yp()];
          vval.pp = v_fa[i.ypp()];

          fval.mm = f_fa[i.ymm()];
          fval.m = f_fa[i.ym()];
          fval.c = f[i];
          fval.p = f_fa[i.yp()];
          fval.pp = f_fa[i.ypp()];

          // Non-centred stencil
          if (inloc == CELL_YLOW) {
            // Producing a stencil centred around a lower Y value
            vval.pp = vval.p;
            vval.p = vval.c;
          } else {
            // Stencil centred around a cell centre
            vval.mm = vval.m;
            vval.m = vval.c;
          }
          result[i] = func(vval, fval);
        }
      }
    }
  } else {
    // Non-staggered case

    Mesh::upwind_func func = fVDDY;
    DiffLookup *table = UpwindTable;

    if (method != DIFF_DEFAULT) {
      // Lookup function
      func = lookupFunc(table, method);
    }

    if (f.hasYupYdown() && ((&f.yup() != &f) || (&f.ydown() != &f))) {
      // f has yup and ydown fields which are distinct
      const Field3D f_yup = f.yup();
      const Field3D f_ydown = f.ydown();
      BOUT_OMP(parallel) {
        stencil fs;
        BOUT_FOR_INNER(i, this->getRegion3D(region_str)) {
          fs.m = f_ydown[i.ym()];
          fs.c = f[i];
          fs.p = f_yup[i.yp()];

          result[i] = func(v[i], fs);
        }
      }
    } else {
      // Not using yup/ydown fields, so first transform to field-aligned coordinates
      Field3D f_fa = this->toFieldAligned(f);
      Field3D v_fa = this->toFieldAligned(v);

      if (this->ystart > 1) {
        BOUT_OMP(parallel) {
          stencil fs;
          BOUT_FOR_INNER(i, this->getRegion3D(region_str)) {
            fs.mm = f_fa[i.ymm()];
            fs.m = f_fa[i.ym()];
            fs.c = f_fa[i];
            fs.p = f_fa[i.yp()];
            fs.pp = f_fa[i.ypp()];

            result[i] = func(v_fa[i], fs);
          }
        }
      } else {
        BOUT_OMP(parallel) {
          stencil fs;
          BOUT_FOR_INNER(i, this->getRegion3D(region_str)) {
            fs.m = f_fa[i.ym()];
            fs.c = f_fa[i];
            fs.p = f_fa[i.yp()];

            result[i] = func(v_fa[i], fs);
          }
        }
      }
      // Shift result back
      result = this->fromFieldAligned(result);
    }
  }

  result.setLocation(outloc);

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
  TRACE("Mesh::indexVDDZ");

  ASSERT1(this == v.getMesh());
  ASSERT1(this == f.getMesh());

  CELL_LOC vloc = v.getLocation();
  CELL_LOC inloc = f.getLocation(); // Input location
  if (outloc == CELL_DEFAULT)
    outloc = inloc;
  // Allowed staggers:
  ASSERT1(outloc == inloc &&
          ((vloc == inloc) || (vloc == CELL_CENTRE && inloc == CELL_ZLOW) ||
           (vloc == CELL_ZLOW && inloc == CELL_CENTRE)));

  Field3D result(this);
  result.allocate(); // Make sure data allocated

  /// Convert REGION enum to a Region string identifier
  const auto region_str = REGION_STRING(region);

  if (StaggerGrids && (vloc != inloc)) {
    // Staggered grids enabled, and velocity at different location to value

    Mesh::flux_func func = sfVDDZ;
    DiffLookup *table = UpwindTable;

    // V staggered w.r.t. variable
    func = sfVDDZ;
    table = UpwindStagTable;

    if (method != DIFF_DEFAULT) {
      // Lookup function
      func = lookupFunc(table, method);
    }

    BOUT_OMP(parallel) {
      stencil vval, fval;
      BOUT_FOR_INNER(i, this->getRegion3D(region_str)) {
        fval.mm = f[i.zmm()];
        fval.m = f[i.zm()];
        fval.c = f[i];
        fval.p = f[i.zp()];
        fval.pp = f[i.zpp()];

        vval.mm = v[i.zmm()];
        vval.m = v[i.zm()];
        vval.c = v[i];
        vval.p = v[i.zp()];
        vval.pp = v[i.zpp()];

        if (inloc == CELL_ZLOW) {
          // Producing a stencil centred around a lower Z value
          vval.pp = vval.p;
          vval.p = vval.c;

        } else {
          // Stencil centred around a cell centre
          vval.mm = vval.m;
          vval.m = vval.c;
        }
        result[i] = func(vval, fval);
      }
    }
  } else {
    Mesh::upwind_func func = fVDDZ;
    DiffLookup *table = UpwindTable;

    if (method != DIFF_DEFAULT) {
      // Lookup function
      func = lookupFunc(table, method);
    }

    BOUT_OMP(parallel) {
      stencil fval;
      BOUT_FOR_INNER(i, this->getRegion3D(region_str)) {
        fval.mm = f[i.zmm()];
        fval.m = f[i.zm()];
        fval.c = f[i];
        fval.p = f[i.zp()];
        fval.pp = f[i.zpp()];

        result[i] = func(v[i], fval);
      }
    }
  }

  result.setLocation(outloc);

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

  ASSERT1(this->xstart > 0); // Need at least one guard cell

  if (outloc == CELL_DEFAULT)
    outloc = f.getLocation();

  if ((method == DIFF_SPLIT) || ((method == DIFF_DEFAULT) && (fFDDX == nullptr))) {
    // Split into an upwind and a central differencing part
    // d/dx(v*f) = v*d/dx(f) + f*d/dx(v)
    return indexVDDX(v, f, outloc, DIFF_DEFAULT) + interp_to(f, outloc) * indexDDX(v, outloc);
  }

  ASSERT1(outloc == f.getLocation() && v.getLocation() == f.getLocation());

  Mesh::flux_func func = fFDDX;
  if (method != DIFF_DEFAULT) {
    // Lookup function
    func = lookupFunc(FluxTable, method);
  }

  Field2D result(this);
  result.allocate(); // Make sure data allocated

  ASSERT1(this == v.getMesh());
  ASSERT1(this == f.getMesh());

  /// Convert REGION enum to a Region string identifier
  const auto region_str = REGION_STRING(region);

  if (this->xstart > 1) {
    // Two or more guard cells
    BOUT_OMP(parallel) {
      stencil fs, vs;
      BOUT_FOR_INNER(i, this->getRegion2D(region_str)) {
        fs.mm = f[i.xmm()];
        fs.m = f[i.xm()];
        fs.c = f[i];
        fs.p = f[i.xp()];
        fs.pp = f[i.xpp()];

        vs.mm = v[i.xmm()];
        vs.m = v[i.xm()];
        vs.c = v[i];
        vs.p = v[i.xp()];
        vs.pp = v[i.xpp()];

        result[i] = func(vs, fs);
      }
    }
  } else {
    // Only one guard cell
    BOUT_OMP(parallel) {
      stencil fs, vs;
      BOUT_FOR_INNER(i, this->getRegion2D(region_str)) {
        fs.m = f[i.xm()];
        fs.c = f[i];
        fs.p = f[i.xp()];

        vs.m = v[i.xm()];
        vs.c = v[i];
        vs.p = v[i.xp()];

        result[i] = func(vs, fs);
      }
    }
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

  ASSERT1(this->xstart > 0); // Need at least one guard cell

  CELL_LOC vloc = v.getLocation();
  CELL_LOC inloc = f.getLocation(); // Input location
  if (outloc == CELL_DEFAULT)
    outloc = inloc;

  if ((method == DIFF_SPLIT) || ((method == DIFF_DEFAULT) && (fFDDX == nullptr))) {
    // Split into an upwind and a central differencing part
    // d/dx(v*f) = v*d/dx(f) + f*d/dx(v)
    return indexVDDX(v, f, outloc, DIFF_DEFAULT) + indexDDX(v, outloc, DIFF_DEFAULT) * interp_to(f, outloc);
  }

  ASSERT1(this == f.getMesh());
  ASSERT1(this == v.getMesh());

  // Allowed staggers:
  ASSERT1(outloc == inloc &&
          ((vloc == inloc) || (vloc == CELL_CENTRE && inloc == CELL_XLOW) ||
           (vloc == CELL_XLOW && inloc == CELL_CENTRE)));

  Mesh::flux_func func = fFDDX;
  DiffLookup *table = FluxTable;

  if (StaggerGrids && (vloc != inloc)) {
    // V staggered w.r.t. variable
    func = sfFDDX;
    table = FluxStagTable;
  }

  if (method != DIFF_DEFAULT) {
    // Lookup function
    func = lookupFunc(table, method);
  }

  Field3D result(this);
  result.allocate(); // Make sure data allocated

  /// Convert REGION enum to a Region string identifier
  const auto region_str = REGION_STRING(region);

  if (this->xstart > 1) {
    // Two or more guard cells
    if (StaggerGrids && vloc != inloc) {
      if (inloc == CELL_XLOW) {
        // Producing a stencil centred around a lower X value
        BOUT_OMP(parallel) {
          stencil fs, vs;
          BOUT_FOR_INNER(i, this->getRegion3D(region_str)) {
            // Location of f always the same as the output
            fs.mm = f[i.xmm()];
            fs.m = f[i.xm()];
            fs.c = f[i];
            fs.p = f[i.xp()];
            fs.pp = f[i.xpp()];

            // Note: Location in diffloc
            vs.mm = v[i.xmm()];
            vs.m = v[i.xm()];
            vs.p = v[i];
            vs.pp = v[i.xp()];

            result[i] = func(vs, fs);
          }
        }
      } else {
        // Stencil centred around a cell centre
        BOUT_OMP(parallel) {
          stencil fs, vs;
          BOUT_FOR_INNER(i, this->getRegion3D(region_str)) {
            // Location of f always the same as the output
            fs.mm = f[i.xmm()];
            fs.m = f[i.xm()];
            fs.c = f[i];
            fs.p = f[i.xp()];
            fs.pp = f[i.xpp()];

            vs.mm = v[i.xm()];
            vs.m = v[i];
            vs.p = v[i.xp()];
            vs.pp = v[i.xpp()];

            result[i] = func(vs, fs);
          }
        }
      }
    } else {
      // Non-staggered, two or more guard cells
      BOUT_OMP(parallel) {
        stencil fs, vs;
        BOUT_FOR_INNER(i, this->getRegion3D(region_str)) {
          // Location of f always the same as the output
          fs.mm = f[i.xmm()];
          fs.m = f[i.xm()];
          fs.c = f[i];
          fs.p = f[i.xp()];
          fs.pp = f[i.xpp()];

          // Note: Location in diffloc
          vs.mm = v[i.xmm()];
          vs.m = v[i.xm()];
          vs.c = v[i];
          vs.p = v[i.xp()];
          vs.pp = v[i.xpp()];

          result[i] = func(vs, fs);
        }
      }
    }
  } else {
    // One guard cell
    if (StaggerGrids && vloc != inloc) {
      if (inloc == CELL_XLOW) {
        // Producing a stencil centred around a lower X value
        BOUT_OMP(parallel) {
          stencil fs, vs;
          BOUT_FOR_INNER(i, this->getRegion3D(region_str)) {
            // Location of f always the same as the output
            fs.m = f[i.xm()];
            fs.c = f[i];
            fs.p = f[i.xp()];

            // Note: Location in diffloc
            vs.m = v[i.xm()];
            vs.p = v[i];
            vs.pp = v[i.xp()];

            result[i] = func(vs, fs);
          }
        }
      } else {
        // Stencil centred around a cell centre
        BOUT_OMP(parallel) {
          stencil fs, vs;
          BOUT_FOR_INNER(i, this->getRegion3D(region_str)) {
            // Location of f always the same as the output
            fs.m = f[i.xm()];
            fs.c = f[i];
            fs.p = f[i.xp()];

            vs.mm = v[i.xm()];
            vs.m = v[i];
            vs.p = v[i.xp()];

            result[i] = func(vs, fs);
          }
        }
      }
    } else {
      // Non-staggered, one guard cell
      BOUT_OMP(parallel) {
        stencil fs, vs;
        BOUT_FOR_INNER(i, this->getRegion3D(region_str)) {
          // Location of f always the same as the output
          fs.m = f[i.xm()];
          fs.c = f[i];
          fs.p = f[i.xp()];

          vs.m = v[i.xm()];
          vs.c = v[i];
          vs.p = v[i.xp()];

          result[i] = func(vs, fs);
        }
      }
    }
  }

  result.setLocation(outloc);

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

  ASSERT1(this->ystart > 0); // Need at least one guard cell
  ASSERT1(this == v.getMesh());
  ASSERT1(this == f.getMesh());

  if (outloc == CELL_DEFAULT)
    outloc = f.getLocation();

  if ((method == DIFF_SPLIT) || ((method == DIFF_DEFAULT) && (fFDDY == nullptr))) {
    // Split into an upwind and a central differencing part
    // d/dx(v*f) = v*d/dx(f) + f*d/dx(v)
    return indexVDDY(v, f, outloc, DIFF_DEFAULT) + interp_to(f, outloc) * indexDDY(v, outloc);
  }

  ASSERT1(outloc == f.getLocation() && v.getLocation() == f.getLocation());

  Mesh::flux_func func = fFDDY;
  if (method != DIFF_DEFAULT) {
    // Lookup function
    func = lookupFunc(FluxTable, method);
  }

  Field2D result(this);
  result.allocate(); // Make sure data allocated
  result.setLocation(f.getLocation());

  /// Convert REGION enum to a Region string identifier
  const auto region_str = REGION_STRING(region);

  if (this->ystart > 1) {
    // Two or more guard cells
    BOUT_OMP(parallel) {
      stencil fs, vs;
      BOUT_FOR_INNER(i, this->getRegion2D(region_str)) {
        fs.mm = f[i.ymm()];
        fs.m = f[i.ym()];
        fs.c = f[i];
        fs.p = f[i.yp()];
        fs.pp = f[i.ypp()];

        vs.mm = v[i.ymm()];
        vs.m = v[i.ym()];
        vs.c = v[i];
        vs.p = v[i.yp()];
        vs.pp = v[i.ypp()];

        result[i] = func(vs, fs);
      }
    }

  } else {
    // Only one guard cell
    BOUT_OMP(parallel) {
      stencil fs, vs;
      BOUT_FOR_INNER(i, this->getRegion2D(region_str)) {
        fs.m = f[i.ym()];
        fs.c = f[i];
        fs.p = f[i.yp()];

        vs.m = v[i.ym()];
        vs.c = v[i];
        vs.p = v[i.yp()];

        result[i] = func(vs, fs);
      }
    }
  }

  result.setLocation(outloc);

#if CHECK > 0
  // Mark boundaries as invalid
  result.bndry_xin = result.bndry_xout = false;
#endif

  return result;
}

const Field3D Mesh::indexFDDY(const Field3D &v, const Field3D &f, CELL_LOC outloc,
                              DIFF_METHOD method, REGION region) {
  TRACE("Mesh::indexFDDY");

  CELL_LOC vloc = v.getLocation();
  CELL_LOC inloc = f.getLocation(); // Input location
  if (outloc == CELL_DEFAULT)
    outloc = inloc;

  if ((method == DIFF_SPLIT) || ((method == DIFF_DEFAULT) && (fFDDY == nullptr))) {
    // Split into an upwind and a central differencing part
    // d/dx(v*f) = v*d/dx(f) + f*d/dx(v)
    return indexVDDY(v, f, outloc, DIFF_DEFAULT) + indexDDY(v, outloc, DIFF_DEFAULT) * interp_to(f, outloc);
  }
  Mesh::flux_func func = fFDDY;
  DiffLookup *table = FluxTable;

  // Allowed staggers:
  ASSERT1(outloc == inloc &&
          ((vloc == inloc) || (vloc == CELL_CENTRE && inloc == CELL_YLOW) ||
           (vloc == CELL_YLOW && inloc == CELL_CENTRE)));

  if (StaggerGrids && (vloc != inloc)) {
    // V staggered w.r.t. variable
    func = sfFDDY;
    table = FluxStagTable;
  }

  if (method != DIFF_DEFAULT) {
    // Lookup function
    func = lookupFunc(table, method);
  }

  if (func == nullptr) {
    // To catch when no function
    return indexVDDY(v, f, outloc, DIFF_DEFAULT) + indexDDY(v, outloc, DIFF_DEFAULT) * interp_to(f, outloc);
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

  /// Convert REGION enum to a Region string identifier
  const auto region_str = REGION_STRING(region);

  if (vUseUpDown && fUseUpDown) {
    // Both v and f have up/down fields
    BOUT_OMP(parallel) {
      stencil fval, vval;
      BOUT_FOR_INNER(i, this->getRegion3D(region_str)) {
        fval.m = f.ydown()[i.ym()];
        fval.c = f[i];
        fval.p = f.yup()[i.yp()];

        vval.m = v.ydown()[i.ym()];
        vval.c = v[i];
        vval.p = v.yup()[i.yp()];

        if (StaggerGrids && (inloc != vloc)) {
          // Non-centred stencil
          if (inloc == CELL_YLOW) {
            // Producing a stencil centred around a lower Y value
            vval.pp = vval.p;
            vval.p = vval.c;
          } else {
            // Stencil centred around a cell centre
            vval.mm = vval.m;
            vval.m = vval.c;
          }
        }
        result[i] = func(vval, fval);
      }
    }
  } else {
    // Both must shift to field aligned
    Field3D v_fa = this->toFieldAligned(v);
    Field3D f_fa = this->toFieldAligned(f);
    BOUT_OMP(parallel) {
      stencil fval, vval;
      BOUT_FOR_INNER(i, this->getRegion3D(region_str)) {
        fval.mm = f_fa[i.ymm()];
        fval.m = f_fa[i.ym()];
        fval.c = f_fa[i];
        fval.p = f_fa[i.yp()];
        fval.pp = f_fa[i.ypp()];

        vval.mm = v_fa[i.ymm()];
        vval.m = v_fa[i.ym()];
        vval.c = v_fa[i];
        vval.p = v_fa[i.yp()];
        vval.pp = v_fa[i.ypp()];

        if (StaggerGrids && (inloc != vloc)) {
          // Non-centred stencil
          if (inloc == CELL_YLOW) {
            // Producing a stencil centred around a lower Y value
            vval.pp = vval.p;
            vval.p = vval.c;
          } else {
            // Stencil centred around a cell centre
            vval.mm = vval.m;
            vval.m = vval.c;
          }
        }
        result[i] = func(vval, fval);
      }
    }
  }

  result.setLocation(outloc);

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

  CELL_LOC vloc = v.getLocation();
  CELL_LOC inloc = f.getLocation(); // Input location
  if (outloc == CELL_DEFAULT)
    outloc = inloc;

  if ((method == DIFF_SPLIT) || ((method == DIFF_DEFAULT) && (fFDDZ == nullptr))) {
    // Split into an upwind and a central differencing part
    // d/dx(v*f) = v*d/dx(f) + f*d/dx(v)
    return indexVDDZ(v, f, outloc, DIFF_DEFAULT) +
           indexDDZ(v, outloc, DIFF_DEFAULT, true) * interp_to(f, outloc);
  }

  Mesh::flux_func func = fFDDZ;
  DiffLookup *table = FluxTable;

  // Allowed staggers:
  ASSERT1(outloc == inloc &&
          ((vloc == inloc) || (vloc == CELL_CENTRE && inloc == CELL_ZLOW) ||
           (vloc == CELL_ZLOW && inloc == CELL_CENTRE)));

  if (StaggerGrids && (vloc != inloc)) {
    // V staggered w.r.t. variable
    func = sfFDDZ;
    table = FluxStagTable;
  }

  if (method != DIFF_DEFAULT) {
    // Lookup function
    func = lookupFunc(table, method);
  }

  ASSERT1(this == v.getMesh());
  ASSERT1(this == f.getMesh());

  Field3D result(this);
  result.allocate(); // Make sure data allocated

  /// Convert REGION enum to a Region string identifier
  const auto region_str = REGION_STRING(region);

  BOUT_OMP(parallel) {
    stencil vval, fval;
    BOUT_FOR_INNER(i, this->getRegion3D(region_str)) {
      fval.mm = f[i.zmm()];
      fval.m = f[i.zm()];
      fval.c = f[i];
      fval.p = f[i.zp()];
      fval.pp = f[i.zpp()];

      vval.mm = v[i.zmm()];
      vval.m = v[i.zm()];
      vval.c = v[i];
      vval.p = v[i.zp()];
      vval.pp = v[i.zpp()];

      if (StaggerGrids && (inloc != vloc)) {
        // Non-centred stencil

        if (inloc == CELL_ZLOW) {
          // Producing a stencil centred around a lower Z value
          vval.pp = vval.p;
          vval.p = vval.c;
        } else {
          // Stencil centred around a cell centre
          vval.mm = vval.m;
          vval.m = vval.c;
        }
      }
      result[i] = func(vval, fval);
    }
  }
  result.setLocation(outloc);

#if CHECK > 0
  // Mark boundaries as invalid
  result.bndry_xin = result.bndry_xout = result.bndry_yup = result.bndry_ydown = false;
#endif

  return result;
}
