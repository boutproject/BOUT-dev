/**************************************************************************
 * Basic derivative methods
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

#include "globals.h"
#include "derivs.h"
#include "stencils.h"
#include "utils.h"
#include "fft.h"
#include "interpolation.h"

#include <math.h>
#include <string.h>
#include <stdlib.h>

typedef BoutReal (*deriv_func)(stencil &); // f
typedef BoutReal (*upwind_func)(stencil &, stencil &); // v, f

/*******************************************************************************
 * Basic derivative methods.
 * All expect to have an input grid cell at the same location as the output
 * Hence convert cell centred values -> centred values, or left -> left
 *******************************************************************************/

const BoutReal WENO_SMALL = 1.0e-8; // Small number for WENO schemes

////////////////////// FIRST DERIVATIVES /////////////////////

/// central, 2nd order
BoutReal DDX_C2(stencil &f)
{
  return 0.5*(f.p - f.m);
}

/// central, 4th order
BoutReal DDX_C4(stencil &f)
{
  return (8.*f.p - 8.*f.m + f.mm - f.pp)/12.;
}

/// Central WENO method, 2nd order (reverts to 1st order near shocks)
BoutReal DDX_CWENO2(stencil &f)
{
  BoutReal isl, isr, isc; // Smoothness indicators
  BoutReal al, ar, ac, sa; // Un-normalised weights
  BoutReal dl, dr, dc; // Derivatives using different stencils
  
  dc = 0.5*(f.p - f.m);
  dl = f.c - f.m;
  dr = f.p - f.c;
  
  isl = SQ(dl);
  isr = SQ(dr);
  isc = (13./3.)*SQ(f.p - 2.*f.c + f.m) + 0.25*SQ(f.p-f.m);
  
  al = 0.25/SQ(WENO_SMALL + isl);
  ar = 0.25/SQ(WENO_SMALL + isr);
  ac = 0.5/SQ(WENO_SMALL + isc);
  sa = al + ar + ac;

  return (al*dl + ar*dr + ac*dc)/sa;
}

///////////////////// SECOND DERIVATIVES ////////////////////

/// Second derivative: Central, 2nd order
BoutReal D2DX2_C2(stencil &f)
{
  return f.p + f.m - 2.*f.c;
}

/// Second derivative: Central, 4th order
BoutReal D2DX2_C4(stencil &f)
{
  return (-f.pp + 16.*f.p - 30.*f.c + 16.*f.m - f.mm)/12.;
}

//////////////////////// UPWIND METHODS ///////////////////////

/// Upwinding: Central, 2nd order
BoutReal VDDX_C2(stencil &v, stencil &f)
{
  return v.c*0.5*(f.p - f.m);
}

/// Upwinding: Central, 4th order
BoutReal VDDX_C4(stencil &v, stencil &f)
{
  return v.c*(8.*f.p - 8.*f.m + f.mm - f.pp)/12.;
}

/// upwind, 1st order
BoutReal VDDX_U1(stencil &v, stencil &f)
{
  return v.c>=0.0 ? v.c*(f.c - f.m): v.c*(f.p - f.c);

  /*
  // Velocity at lower end
  BoutReal vs = 0.5*(v.m + v.c);
  BoutReal result = vs >= 0.0 ? vs * f.m : vs * f.c;
  // and at upper 
  vs = 0.5*(v.c + v.p);
  result -= vs >= 0.0 ? vs * f.c : vs * f.p;

  return result;
  */
}

/// upwind, 4th order
BoutReal VDDX_U4(stencil &v, stencil &f)
{
  return v.c >= 0.0 ? v.c*(4.*f.p - 12.*f.m + 2.*f.mm + 6.*f.c)/12.
    : v.c*(-4.*f.m + 12.*f.p - 2.*f.pp - 6.*f.c)/12.;
}

/// Van Leer limiter. Used in TVD code
BoutReal VANLEER(BoutReal r) 
{
  return r + fabs(r)/(1.0 + fabs(r));
}

/// TVD upwinding (2nd order)
/// WARNING WARNING : THIS TVD IMPLEMENTATION DOESN'T WORK PROPERLY
BoutReal VDDX_TVD(BoutReal vc, BoutReal vm, BoutReal vp, BoutReal fc, BoutReal fm, BoutReal fp, BoutReal fmm, BoutReal fpp)
{
  BoutReal denom, res, ri, ri1, ri_1, fluxLeft, fluxRight;

  if (vc>=0.0){

    /*-smoothness indicator*/
    denom=fp-fc;
    if (fabs(denom) < 1e-20) denom=1e-20;
    ri=(fc-fm)/denom;

    denom=fc-fm;
    if (fabs(denom) < 1e-20) denom=1e-20;
    ri_1=(fm-fmm)/denom;

    /*;-nonlinear TVD flux*/
    fluxRight=fc + 0.5*VANLEER(ri)*(fp-fc);
    fluxLeft=fm + 0.5*VANLEER(ri_1)*(fc-fm);
    res = vc*(fluxRight-fluxLeft); /*divide by dx outside*/
    
  }
  else{

    /*;-smoothness indicator*/
    denom=fm-fc;
    if (fabs(denom) < 1e-20) denom=1e-20;
    ri=(fc-fp)/denom;
    
    denom=(fc-fp);
    if (fabs(denom) < 1e-20) denom=1e-20;
    ri1=(fp-fpp)/denom;

    /*;-nonlinear TVD flux*/
    fluxRight=(fp - 0.5*VANLEER(ri1)*(fp-fc));
    fluxLeft=(fc - 0.5*VANLEER(ri)*(fc-fm));
    res = vc*(fluxRight-fluxLeft); /*divide by dx outside*/    
  }
  return res;
}

/// 3rd-order WENO scheme
BoutReal VDDX_WENO3(stencil &v, stencil &f)
{
  BoutReal deriv, w, r;

  if(v.c > 0.0) {
    // Left-biased stencil
    
    r = (WENO_SMALL + SQ(f.c - 2.0*f.m + f.mm)) / (WENO_SMALL + SQ(f.p - 2.0*f.c + f.m));
    w = 1.0 / (1.0 + 2.0*r*r);

    deriv = 0.5*(f.p - f.m) - 0.5*w*(-f.mm + 3.*f.m - 3.*f.c + f.p);
    
  }else {
    // Right-biased
    
    r = (WENO_SMALL + SQ(f.pp - 2.0*f.p + f.c)) / (WENO_SMALL + SQ(f.p - 2.0*f.c + f.m));
    w = 1.0 / (1.0 + 2.0*r*r);
    
    deriv = 0.5*(f.p - f.m) - 0.5*w*( -f.m + 3.*f.c - 3.*f.p + f.pp );
  }

  return v.c*deriv;
}

/// 3rd-order CWENO. Uses the upwinding code and split flux
BoutReal DDX_CWENO3(stencil &f)
{
  BoutReal a, ma = fabs(f.c);
  // Split flux
  a = fabs(f.m); if(a > ma) ma = a;
  a = fabs(f.p); if(a > ma) ma = a;
  a = fabs(f.mm); if(a > ma) ma = a;
  a = fabs(f.pp); if(a > ma) ma = a;
  
  stencil sp, vp, sm, vm;
  
  vp.c = 0.5; vm.c = -0.5;
  
  sp = f + ma;
  sm = ma - f;
  
  return VDDX_WENO3(vp, sp) + VDDX_WENO3(vm, sm);
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
BoutReal DDX_C2_stag(stencil &f)
{
  return f.p - f.m;
}

BoutReal DDX_C4_stag(stencil &f)
{
  //output.write("HERE");
  return ( 27.*(f.p - f.m) - (f.pp - f.mm) ) / 24.;
}

/////////////////////// SECOND DERIVATIVES //////////////////////
// Map Centre -> Low or Low -> Centre

BoutReal D2DX2_C4_stag(stencil &f)
{
  return ( f.pp + f.mm - f.p - f.m ) / 2.;
}

/////////////////////////// UPWINDING ///////////////////////////
// Map (Low, Centre) -> Centre  or (Centre, Low) -> Low
// Hence v contains only (mm, m, p, pp) fields whilst f has 'c' too
//
// v.p is v at +1/2, v.m is at -1/2

BoutReal VDDX_U1_stag(stencil &v, stencil &f)
{
  // Lower cell boundary
  BoutReal result = (v.m >= 0) ? v.m * f.m : v.m * f.c;
  
  // Upper cell boundary
  result -= (v.p >= 0) ? v.p * f.c : v.p * f.p;

  return result;
}

/*******************************************************************************
 * Lookup tables of functions. Map between names, codes and functions
 *******************************************************************************/

/// Translate between DIFF_METHOD codes, and functions
struct DiffLookup {
  DIFF_METHOD method;
  deriv_func func;     // Single-argument differencing function
  upwind_func up_func; // Upwinding function
};

/// Translate between short names, long names and DIFF_METHOD codes
struct DiffNameLookup {
  DIFF_METHOD method;
  const char* label; // Short name
  const char* name;  // Long name
};

/// Differential function name/code lookup
static DiffNameLookup DiffNameTable[] = { {DIFF_U1, "U1", "First order upwinding"},
					  {DIFF_C2, "C2", "Second order central"},
					  {DIFF_W2, "W2", "Second order WENO"},
					  {DIFF_W3, "W3", "Third order WENO"},
					  {DIFF_C4, "C4", "Fourth order central"},
					  {DIFF_U4, "U4", "Fourth order upwinding"},
					  {DIFF_FFT, "FFT", "FFT"},
					  {DIFF_DEFAULT}}; // Use to terminate the list

/// First derivative lookup table
static DiffLookup FirstDerivTable[] = { {DIFF_C2, DDX_C2,     NULL},
					{DIFF_W2, DDX_CWENO2, NULL},
					{DIFF_W3, DDX_CWENO3, NULL},
					{DIFF_C4, DDX_C4,     NULL},
					{DIFF_FFT, NULL,      NULL},
					{DIFF_DEFAULT}};

/// Second derivative lookup table
static DiffLookup SecondDerivTable[] = { {DIFF_C2, D2DX2_C2, NULL},
					{DIFF_C4, D2DX2_C4, NULL},
					{DIFF_FFT, NULL,    NULL},
					{DIFF_DEFAULT}};

/// Upwinding functions lookup table
static DiffLookup UpwindTable[] = { {DIFF_U1, NULL, VDDX_U1},
				    {DIFF_C2, NULL, VDDX_C2},
				    {DIFF_U4, NULL, VDDX_U4},
				    {DIFF_W3, NULL, VDDX_WENO3},
				    {DIFF_C4, NULL, VDDX_C4},
				    {DIFF_DEFAULT}};

/// First staggered derivative lookup
static DiffLookup FirstStagDerivTable[] = { {DIFF_C2, DDX_C2_stag, NULL}, 
					    {DIFF_C4, DDX_C4_stag, NULL},
					    {DIFF_DEFAULT}};

/// Second staggered derivative lookup
static DiffLookup SecondStagDerivTable[] = { {DIFF_C4, D2DX2_C4_stag, NULL},
					     {DIFF_DEFAULT}};

/// Upwinding staggered lookup
static DiffLookup UpwindStagTable[] = { {DIFF_U1, NULL, VDDX_U1_stag},
					{DIFF_DEFAULT} };

/*******************************************************************************
 * Routines to use the above tables to map between function codes, names
 * and pointers
 *******************************************************************************/

deriv_func lookupFunc(DiffLookup* table, DIFF_METHOD method)
{
  int i = 0;
  do {
    if(table[i].method == method)
      return table[i].func;
    i++;
  }while(table[i].method != DIFF_DEFAULT);
  // Not found in list. Return the first 
  
  return table[0].func;
}

upwind_func lookupUpwindFunc(DiffLookup* table, DIFF_METHOD method)
{
  int i = 0;
  do {
    if(table[i].method == method)
      return table[i].up_func;
    i++;
  }while(table[i].method != DIFF_DEFAULT);
  // Not found in list. Return the first 
  
  return table[0].up_func;
}

/// Test if a given DIFF_METHOD exists in a table
bool isImplemented(DiffLookup* table, DIFF_METHOD method)
{
  int i = 0;
  do {
    if(table[i].method == method)
      return true;
    i++;
  }while(table[i].method != DIFF_DEFAULT);
  
  return false;
}

/// This function is used during initialisation only (i.e. doesn't need to be particularly fast)
/// Returns DIFF_METHOD, rather than function so can be applied to central and upwind tables
DIFF_METHOD lookupFunc(DiffLookup *table, const string &label)
{
  DIFF_METHOD matchtype; // code which matches just the first letter ('C', 'U' or 'W')

  if(label.empty())
    return table[0].method;

  matchtype = DIFF_DEFAULT;
  int typeind;
  
  // Loop through the name lookup table
  int i = 0;
  do {
    if((toupper(DiffNameTable[i].label[0]) == toupper(label[0])) && isImplemented(table, DiffNameTable[i].method)) {
      matchtype = DiffNameTable[i].method;
      typeind = i;
      
      if(strcasecmp(label.c_str(), DiffNameTable[i].label) == 0) {// Whole match
	return matchtype;
      }
    }
    i++;
  }while(DiffNameTable[i].method != DIFF_DEFAULT);
  
  // No exact match, so return matchtype.

  if(matchtype == DIFF_DEFAULT) {
    // No type match either. Return the first value in the table
    matchtype = table[0].method;
    output << " No match for '" << label << "' -> ";
  }else
    output << " Type match for '" << label << "' ->";

  return matchtype;
}

void printFuncName(DIFF_METHOD method)
{
  // Find this entry

  int i = 0;
  do {
    if(DiffNameTable[i].method == method) {
      output.write(" %s (%s)\n", DiffNameTable[i].name, DiffNameTable[i].label);
      return;
    }
    i++;
  }while(DiffNameTable[i].method != DIFF_DEFAULT);

  // None
  output.write(" == INVALID DIFFERENTIAL METHOD ==\n");
}

/*******************************************************************************
 * Default functions
 *
 *
 *******************************************************************************/

// Central -> Central (or Left -> Left) functions
deriv_func fDDX, fDDY, fDDZ;     // Differencing methods for each dimension
deriv_func fD2DX2, fD2DY2, fD2DZ2;  // second differential operators
upwind_func fVDDX, fVDDY, fVDDZ; // Upwind functions in the three directions

// Central -> Left (or Left -> Central) functions
deriv_func sfDDX, sfDDY, sfDDZ;
deriv_func sfD2DX2, sfD2DY2, sfD2DZ2;
upwind_func sfVDDX, sfVDDY, sfVDDZ;

/*******************************************************************************
 * Initialisation
 *******************************************************************************/

/// Set the derivative method, given a table and option name
void derivs_set(DiffLookup *table, const char* name, deriv_func &f)
{
/*  char *label = options.getString(name);*/
  string label;
  options.get<string>(name, label, "");

  DIFF_METHOD method = lookupFunc(table, label); // Find the function
  printFuncName(method); // Print differential function name
  f = lookupFunc(table, method); // Find the function pointer
}

void derivs_set(DiffLookup *table, const char* name, upwind_func &f)
{
/*  char *label = options.getString(name);*/
  string label;
  options.get<string>(name, label, "");

  DIFF_METHOD method = lookupFunc(table, label); // Find the function
  printFuncName(method); // Print differential function name
  f = lookupUpwindFunc(table, method);
}

/// Initialise the derivative methods. Must be called before any derivatives are used
int derivs_init()
{
#ifdef CHECK
  msg_stack.push("Initialising derivatives");
#endif

  /// NOTE: StaggerGrids is also in Mesh, but derivs_init needs to come before Mesh
  bool StaggerGrids;
  options.setSection("");
  OPTION(StaggerGrids,   false);

  output.write("Setting X differencing methods\n");
  options.setSection("ddx");
  output.write("\tFirst       : ");
  derivs_set(FirstDerivTable, "first",  fDDX);
  if(StaggerGrids) {
    output.write("\tStag. First : ");
    derivs_set(FirstStagDerivTable, "first",  sfDDX);
  }
  output.write("\tSecond      : ");
  derivs_set(SecondDerivTable, "second", fD2DX2);
  if(StaggerGrids) {
    output.write("\tStag. Second: ");
    derivs_set(SecondDerivTable, "second", sfD2DX2);
  }
  output.write("\tUpwind      : ");
  derivs_set(UpwindTable,     "upwind", fVDDX);
  if(StaggerGrids) {
    output.write("\tStag. Upwind: ");
    derivs_set(UpwindTable,     "upwind", sfVDDX);
  }
  
  if((fDDX == NULL) || (fD2DX2 == NULL)) {
    output.write("\t***Error: FFT cannot be used in X\n");
    return 1;
  }

  output.write("Setting Y differencing methods\n");
  options.setSection("ddy");
  output.write("\tFirst       : ");
  derivs_set(FirstDerivTable, "first",  fDDY);
  if(StaggerGrids) {
    output.write("\tStag. First : ");
    derivs_set(FirstStagDerivTable, "first",  sfDDY);
  }
  output.write("\tSecond      : ");
  derivs_set(SecondDerivTable, "second", fD2DY2);
  if(StaggerGrids) {
    output.write("\tStag. Second: ");
    derivs_set(SecondDerivTable, "second", sfD2DY2);
  }
  output.write("\tUpwind      : ");
  derivs_set(UpwindTable,     "upwind", fVDDY);
  if(StaggerGrids) {
    output.write("\tStag. Upwind: ");
    derivs_set(UpwindTable,     "upwind", sfVDDY);
  }

  if((fDDY == NULL) || (fD2DY2 == NULL)) {
    output.write("\tError: FFT cannot be used in Y\n");
    return 1;
  }

  output.write("Setting Z differencing methods\n");
  options.setSection("ddz");
  output.write("\tFirst       : ");
  derivs_set(FirstDerivTable, "first",  fDDZ);
  if(StaggerGrids) {
    output.write("\tStag. First : ");
    derivs_set(FirstStagDerivTable, "first",  sfDDZ);
  }
  output.write("\tSecond      : ");
  derivs_set(SecondDerivTable, "second", fD2DZ2);
  if(StaggerGrids) {
    output.write("\tStag. Second: ");
    derivs_set(SecondDerivTable, "second", sfD2DZ2);
  }
  output.write("\tUpwind      : ");
  derivs_set(UpwindTable,     "upwind", fVDDZ);
  if(StaggerGrids) {
    output.write("\tStag. Upwind: ");
    derivs_set(UpwindTable,     "upwind", sfVDDZ);
  }

#ifdef CHECK
  msg_stack.pop();
#endif

  return 0;
}

/*******************************************************************************
 * Apply differential operators. These are fairly brain-dead functions
 * which apply a derivative function to a field (sort of like map). Decisions
 * of what to apply are made in the DDX,DDY and DDZ functions lower down.
 *
 * loc  is the cell location of the result
 *******************************************************************************/

// X derivative

const Field2D applyXdiff(const Field2D &var, deriv_func func, const Field2D &dd, CELL_LOC loc = CELL_DEFAULT)
{
  Field2D result;
  result.allocate(); // Make sure data allocated

  bindex bx;
  stencil s;

  BoutReal **r = result.getData();

  start_index(&bx, RGN_NOX);
  do {
    var.setXStencil(s, bx, loc);
    r[bx.jx][bx.jy] = func(s) / dd[bx.jx][bx.jy];
  }while(next_index2(&bx));

#ifdef CHECK
  // Mark boundaries as invalid
  result.bndry_xin = result.bndry_xout = result.bndry_yup = result.bndry_ydown = false;
#endif

  return result;
}

const Field3D applyXdiff(const Field3D &var, deriv_func func, const Field2D &dd, CELL_LOC loc = CELL_DEFAULT)
{
  Field3D result;
  result.allocate(); // Make sure data allocated
  
  stencil s;

  Field3D vs = var;
  if(mesh->ShiftXderivs && (mesh->ShiftOrder == 0)) {
    // Shift in Z using FFT
    vs = var.shiftZ(true); // Shift into BoutReal space
  }
  
  bindex bx;
  BoutReal ***r = result.getData();
  
  start_index(&bx, RGN_NOX);
  do {
    vs.setXStencil(s, bx, loc);
    r[bx.jx][bx.jy][bx.jz] = func(s) / dd[bx.jx][bx.jy];
  }while(next_index3(&bx));
  
  if(mesh->ShiftXderivs && (mesh->ShiftOrder == 0))
    result = result.shiftZ(false); // Shift back

#ifdef CHECK
  // Mark boundaries as invalid
  result.bndry_xin = result.bndry_xout = result.bndry_yup = result.bndry_ydown = false;
#endif

  return result;
}

// Y derivative

const Field2D applyYdiff(const Field2D &var, deriv_func func, const Field2D &dd, CELL_LOC loc = CELL_DEFAULT)
{
  Field2D result;
  result.allocate(); // Make sure data allocated
  BoutReal **r = result.getData();
  
  bindex bx;
  stencil s;

  start_index(&bx, RGN_NOY);
  do{
    var.setYStencil(s, bx, loc);
    r[bx.jx][bx.jy] = func(s) / dd[bx.jx][bx.jy];
  }while(next_index2(&bx));
  
#ifdef CHECK
  // Mark boundaries as invalid
  result.bndry_yup = result.bndry_ydown = false;
#endif

  return result;
}

const Field3D applyYdiff(const Field3D &var, deriv_func func, const Field2D &dd, CELL_LOC loc = CELL_DEFAULT)
{
  Field3D result;
  result.allocate(); // Make sure data allocated
  BoutReal ***r = result.getData();
  
  stencil s;
  bindex bx;
  start_index(&bx, RGN_NOY);
  do {
    var.setYStencil(s, bx, loc);
    
    r[bx.jx][bx.jy][bx.jz] = func(s) / dd[bx.jx][bx.jy];
  
#ifdef CHECK
      if(!finite(r[bx.jx][bx.jy][bx.jz])) {
	msg_stack.push("At [%d][%d][%d]: %e, %e, %e, %e, %e",
		       bx.jx, bx.jy, bx.jz, 
		       s.mm, s.m, s.c, s.p, s.pp);
	bout_error("Non-finite value\n");
      }
#endif
  }while(next_index3(&bx));

#ifdef CHECK
  // Mark boundaries as invalid
  result.bndry_xin = result.bndry_xout = result.bndry_yup = result.bndry_ydown = false;
#endif

  return result;
}

// Z derivative

const Field3D applyZdiff(const Field3D &var, deriv_func func, BoutReal dd, CELL_LOC loc = CELL_DEFAULT)
{
  Field3D result;
  result.allocate(); // Make sure data allocated
  BoutReal ***r = result.getData();
  
  bindex bx;
  stencil s;

  start_index(&bx, RGN_NOZ);
  do {
    var.setZStencil(s, bx, loc);
    r[bx.jx][bx.jy][bx.jz] = func(s) / dd;
  }while(next_index3(&bx));

  return result;
}

/*******************************************************************************
 * First central derivatives
 *******************************************************************************/

////////////// X DERIVATIVE /////////////////

const Field3D DDX(const Field3D &f, CELL_LOC outloc, DIFF_METHOD method)
{
  deriv_func func = fDDX; // Set to default function
  DiffLookup *table = FirstDerivTable;
  
  CELL_LOC inloc = f.getLocation(); // Input location
  CELL_LOC diffloc = inloc; // Location of differential result

  Field3D result;

  if(mesh->StaggerGrids && (outloc == CELL_DEFAULT)) {
    // Take care of CELL_DEFAULT case
    outloc = diffloc; // No shift (i.e. same as no stagger case)
  }

  if(mesh->StaggerGrids && (outloc != inloc)) {
    // Shifting to a new location
    
    if(((inloc == CELL_CENTRE) && (outloc == CELL_XLOW)) ||
       ((inloc == CELL_XLOW) && (outloc == CELL_CENTRE))) {
      // Shifting in X. Centre -> Xlow, or Xlow -> Centre
      
      func = sfDDX; // Set default
      table = FirstStagDerivTable; // Set table for others
      diffloc = (inloc == CELL_CENTRE) ? CELL_XLOW : CELL_CENTRE;
      
    }else {
      // A more complicated shift. Get a result at cell centre, then shift.
      if(inloc == CELL_XLOW) {
	// Shifting
	
	func = sfDDX; // Set default
	table = FirstStagDerivTable; // Set table for others
	diffloc = CELL_CENTRE;

      }else if(inloc != CELL_CENTRE) {
	// Interpolate then (centre -> centre) then interpolate
	return DDX(interp_to(f, CELL_CENTRE), outloc, method);
      }
    }
  }
  
  if(method != DIFF_DEFAULT) {
    // Lookup function
    func = lookupFunc(table, method);
    if(func == NULL)
      bout_error("Cannot use FFT for X derivatives");
  }
  
  result = applyXdiff(f, func, mesh->dx, diffloc);
  result.setLocation(diffloc); // Set the result location

  result = interp_to(result, outloc); // Interpolate if necessary
  
  if(mesh->ShiftXderivs && mesh->IncIntShear) {
    // Using BOUT-06 style shifting
    result += mesh->IntShiftTorsion * DDZ(f, outloc);
  }

  return result;
}

const Field3D DDX(const Field3D &f, DIFF_METHOD method, CELL_LOC outloc)
{
  return DDX(f, outloc, method);
}

const Field3D DDX(const Field3D &f, DIFF_METHOD method)
{
  return DDX(f, CELL_DEFAULT, method);
}

const Field2D DDX(const Field2D &f)
{
  return applyXdiff(f, fDDX, mesh->dx);
}

////////////// Y DERIVATIVE /////////////////

const Field3D DDY(const Field3D &f, CELL_LOC outloc, DIFF_METHOD method)
{
  deriv_func func = fDDY; // Set to default function
  DiffLookup *table = FirstDerivTable;
  
  CELL_LOC inloc = f.getLocation(); // Input location
  CELL_LOC diffloc = inloc; // Location of differential result

  Field3D result;

  if(mesh->StaggerGrids && (outloc == CELL_DEFAULT)) {
    // Take care of CELL_DEFAULT case
    outloc = diffloc; // No shift (i.e. same as no stagger case)
  }

  if(mesh->StaggerGrids && (outloc != inloc)) {
    // Shifting to a new location
    
    //output.write("\nSHIFTING %s -> %s\n", strLocation(inloc), strLocation(outloc));

    if(((inloc == CELL_CENTRE) && (outloc == CELL_YLOW)) ||
       ((inloc == CELL_YLOW) && (outloc == CELL_CENTRE))) {
      // Shifting in Y. Centre -> Ylow, or Ylow -> Centre
      
      //output.write("SHIFT");

      func = sfDDY; // Set default
      table = FirstStagDerivTable; // Set table for others
      diffloc = (inloc == CELL_CENTRE) ? CELL_YLOW : CELL_CENTRE;
      
    }else {
      // A more complicated shift. Get a result at cell centre, then shift.
      if(inloc == CELL_YLOW) {
	// Shifting
	
	func = sfDDY; // Set default
	table = FirstStagDerivTable; // Set table for others
	diffloc = CELL_CENTRE;

      }else if(inloc != CELL_CENTRE) {
	// Interpolate to centre then call DDY again
	return DDY(interp_to(f, CELL_CENTRE), outloc, method);
      }
    }
  }
  
  if(method != DIFF_DEFAULT) {
    // Lookup function
    func = lookupFunc(table, method);
    if(func == NULL)
      bout_error("Cannot use FFT for Y derivatives");
  }
  
  result = applyYdiff(f, func, mesh->dy, diffloc);
  //output.write("SETTING LOC %s\n", strLocation(diffloc));
  result.setLocation(diffloc); // Set the result location

  return interp_to(result, outloc); // Interpolate if necessary
}

const Field3D DDY(const Field3D &f, DIFF_METHOD method, CELL_LOC outloc)
{
  return DDY(f, outloc, method);
}

const Field3D DDY(const Field3D &f, DIFF_METHOD method)
{
  return DDY(f, CELL_DEFAULT, method);
}

const Field2D DDY(const Field2D &f)
{
  return applyYdiff(f, fDDY, mesh->dy);
}

////////////// Z DERIVATIVE /////////////////

const Field3D DDZ(const Field3D &f, CELL_LOC outloc, DIFF_METHOD method, bool inc_xbndry)
{
  deriv_func func = fDDZ; // Set to default function
  DiffLookup *table = FirstDerivTable;
 
  CELL_LOC inloc = f.getLocation(); // Input location
  CELL_LOC diffloc = inloc; // Location of differential result

  if(mesh->StaggerGrids && (outloc == CELL_DEFAULT)) {
    // Take care of CELL_DEFAULT case
    outloc = diffloc; // No shift (i.e. same as no stagger case)
  }

  if(mesh->StaggerGrids && (outloc != inloc)) {
    // Shifting to a new location
    
    if(((inloc == CELL_CENTRE) && (outloc == CELL_ZLOW)) ||
       ((inloc == CELL_ZLOW) && (outloc == CELL_CENTRE))) {
      // Shifting in Z. Centre -> Zlow, or Zlow -> Centre
      
      func = sfDDZ; // Set default
      table = FirstStagDerivTable; // Set table for others
      diffloc = (inloc == CELL_CENTRE) ? CELL_ZLOW : CELL_CENTRE;
      
    }else {
      // A more complicated shift. Get a result at cell centre, then shift.
      if(inloc == CELL_ZLOW) {
	// Shifting
	
	func = sfDDZ; // Set default
	table = FirstStagDerivTable; // Set table for others
	diffloc = CELL_CENTRE;

      }else if(inloc != CELL_CENTRE) {
	// Interpolate then (centre -> centre) then interpolate
	return DDZ(interp_to(f, CELL_CENTRE), outloc, method);
      }
    }
  }
  
  if(method != DIFF_DEFAULT) {
    // Lookup function
    func = lookupFunc(table, method);
  }

  Field3D result;

  if(func == NULL) {
    // Use FFT
    
    BoutReal shift = 0.; // Shifting result in Z?
    if(mesh->StaggerGrids) {
      if((inloc == CELL_CENTRE) && (diffloc == CELL_ZLOW)) {
	// Shifting down - multiply by exp(-0.5*i*k*dz)
	shift = -1.;
      }else if((inloc == CELL_ZLOW) && (diffloc == CELL_CENTRE)) {
	// Shifting up
	shift = 1.;
      }
    }

    result.allocate(); // Make sure data allocated

    static dcomplex *cv = (dcomplex*) NULL;
    int jx, jy, jz;
    BoutReal kwave;
    BoutReal flt;

    int xge = mesh->xstart, xlt = mesh->xend+1;
    if(inc_xbndry) { // Include x boundary region (for mixed XZ derivatives)
      xge = 0;
      xlt = mesh->ngx;
    }
    
    int ncz = mesh->ngz-1;

    if(cv == (dcomplex*) NULL)
      cv = new dcomplex[ncz/2 + 1];

    for(jx=xge;jx<xlt;jx++) {
      for(jy=0;jy<mesh->ngy;jy++) {
	
	rfft(f[jx][jy], ncz, cv); // Forward FFT

	for(jz=0;jz<=ncz/2;jz++) {
	  kwave=jz*2.0*PI/mesh->zlength; // wave number is 1/[rad]

	  if (jz>0.4*ncz) flt=1e-10; else flt=1.0;
	  cv[jz] *= dcomplex(0.0, kwave) * flt;
	  if(mesh->StaggerGrids)
	    cv[jz] *= exp(Im * (shift * kwave * mesh->dz));
	}
	
	irfft(cv, ncz, result[jx][jy]); // Reverse FFT

	result[jx][jy][ncz] = result[jx][jy][0];
	
      }
    }
    
#ifdef CHECK
    // Mark boundaries as invalid
    result.bndry_xin = result.bndry_xout = result.bndry_yup = result.bndry_ydown = false;
#endif
    
  }else {
    // All other (non-FFT) functions 
    result = applyZdiff(f, func, mesh->dz);
  }
  
  result.setLocation(diffloc);

  return interp_to(result, outloc);
}

const Field3D DDZ(const Field3D &f, DIFF_METHOD method, CELL_LOC outloc, bool inc_xbndry)
{
  return DDZ(f, outloc, method, inc_xbndry);
}

const Field3D DDZ(const Field3D &f, DIFF_METHOD method, bool inc_xbndry)
{
  return DDZ(f, CELL_DEFAULT, method, inc_xbndry);
}

const Field3D DDZ(const Field3D &f, bool inc_xbndry)
{
  return DDZ(f, CELL_DEFAULT, DIFF_DEFAULT, inc_xbndry);
}

const Field2D DDZ(const Field2D &f)
{
  Field2D result;
  result = 0.0;
  return result;
}

const Vector3D DDZ(const Vector3D &v, CELL_LOC outloc, DIFF_METHOD method)
{
  Vector3D result;

  result.covariant = v.covariant;

  result.x = DDZ(v.x, outloc, method);
  result.y = DDZ(v.y, outloc, method);
  result.z = DDZ(v.z, outloc, method);

  return result;
}

const Vector3D DDZ(const Vector3D &v, DIFF_METHOD method, CELL_LOC outloc)
{
  return DDZ(v, outloc, method);
}

const Vector2D DDZ(const Vector2D &v)
{
  Vector2D result;

  result.covariant = v.covariant;

  result.x = 0.;
  result.y = 0.;
  result.z = 0.;

  return result;
}

/*******************************************************************************
 * 2nd derivative
 *******************************************************************************/

////////////// X DERIVATIVE /////////////////

const Field3D D2DX2(const Field3D &f, CELL_LOC outloc, DIFF_METHOD method)
{
  deriv_func func = fD2DX2; // Set to default function
  DiffLookup *table = SecondDerivTable;
  
  CELL_LOC inloc = f.getLocation(); // Input location
  CELL_LOC diffloc = inloc; // Location of differential result
  
  Field3D result;
  
  if(mesh->StaggerGrids && (outloc == CELL_DEFAULT)) {
    // Take care of CELL_DEFAULT case
    outloc = diffloc; // No shift (i.e. same as no stagger case)
  }

  if(mesh->StaggerGrids && (outloc != inloc)) {
    // Shifting to a new location
    
    if(((inloc == CELL_CENTRE) && (outloc == CELL_XLOW)) ||
       ((inloc == CELL_XLOW) && (outloc == CELL_CENTRE))) {
      // Shifting in X. Centre -> Xlow, or Xlow -> Centre
      
      func = sfD2DX2; // Set default
      table = SecondStagDerivTable; // Set table for others
      diffloc = (inloc == CELL_CENTRE) ? CELL_XLOW : CELL_CENTRE;
      
    }else {
      // A more complicated shift. Get a result at cell centre, then shift.
      if(inloc == CELL_XLOW) {
	// Shifting
	
	func = sfD2DX2; // Set default
	table = SecondStagDerivTable; // Set table for others
	diffloc = CELL_CENTRE;

      }else if(inloc != CELL_CENTRE) {
	// Interpolate then (centre -> centre) then interpolate
	return D2DX2(interp_to(f, CELL_CENTRE), outloc, method);
      }
    }
  }

  if(method != DIFF_DEFAULT) {
    // Lookup function
    func = lookupFunc(table, method);
    if(func == NULL)
      bout_error("Cannot use FFT for X derivatives");
  }

  result = applyXdiff(f, func, (mesh->dx*mesh->dx));
  result.setLocation(diffloc);

  if(non_uniform) {
    // Correction for non-uniform mesh
    result += applyXdiff(f, fDDX, mesh->d2x); // CHECK THIS!
  }

  result = interp_to(result, outloc);

  if(mesh->ShiftXderivs && mesh->IncIntShear) {
    mesh->IncIntShear = false; // So DDX doesn't try to include I again
    // Add I^2 d^2/dz^2 term
    result += mesh->IntShiftTorsion^2 * D2DZ2(f, outloc);
    // Mixed derivative
    result += 2.*mesh->IntShiftTorsion * D2DXDZ(f);
    // DDZ term
    result += DDX(mesh->IntShiftTorsion) * DDZ(f, outloc);
    mesh->IncIntShear = true;
  }
  
  return result;
}

const Field3D D2DX2(const Field3D &f, DIFF_METHOD method, CELL_LOC outloc)
{
  return D2DX2(f, outloc, method);
}

const Field2D D2DX2(const Field2D &f)
{
  Field2D result;

  result = applyXdiff(f, fD2DX2, (mesh->dx*mesh->dx));

  if(non_uniform) {
    // Correction for non-uniform mesh
    result += applyXdiff(f, fDDX, mesh->d2x);
  }
  
  return(result);
}

////////////// Y DERIVATIVE /////////////////

const Field3D D2DY2(const Field3D &f, CELL_LOC outloc, DIFF_METHOD method)
{
  deriv_func func = fD2DY2; // Set to default function
  DiffLookup *table = SecondDerivTable;
  
  CELL_LOC inloc = f.getLocation(); // Input location
  CELL_LOC diffloc = inloc; // Location of differential result
  
  Field3D result;
  
  if(mesh->StaggerGrids && (outloc == CELL_DEFAULT)) {
    // Take care of CELL_DEFAULT case
    outloc = diffloc; // No shift (i.e. same as no stagger case)
  }

  if(mesh->StaggerGrids && (outloc != inloc)) {
    // Shifting to a new location
    
    if(((inloc == CELL_CENTRE) && (outloc == CELL_YLOW)) ||
       ((inloc == CELL_YLOW) && (outloc == CELL_CENTRE))) {
      // Shifting in Y. Centre -> Ylow, or Ylow -> Centre
      
      func = sfD2DY2; // Set default
      table = SecondStagDerivTable; // Set table for others
      diffloc = (inloc == CELL_CENTRE) ? CELL_YLOW : CELL_CENTRE;
      
    }else {
      // A more complicated shift. Get a result at cell centre, then shift.
      if(inloc == CELL_YLOW) {
	// Shifting
	
	func = sfD2DY2; // Set default
	table = SecondStagDerivTable; // Set table for others
	diffloc = CELL_CENTRE;

      }else if(inloc != CELL_CENTRE) {
	// Interpolate then (centre -> centre) then interpolate
	return D2DY2(interp_to(f, CELL_CENTRE), outloc, method);
      }
    }
  }

  if(method != DIFF_DEFAULT) {
    // Lookup function
    func = lookupFunc(table, method);
    if(func == NULL)
      bout_error("Cannot use FFT for Y derivatives");
  }

  result = applyYdiff(f, func, (mesh->dy*mesh->dy));
  result.setLocation(diffloc);

  if(non_uniform) {
    // Correction for non-uniform mesh
    result += applyYdiff(f, fDDY, mesh->d2y);
  }

  return interp_to(result, outloc);
}

const Field3D D2DY2(const Field3D &f, DIFF_METHOD method, CELL_LOC outloc)
{
  return D2DY2(f, outloc, method);
}

const Field2D D2DY2(const Field2D &f)
{
  Field2D result;

  result = applyYdiff(f, fD2DY2, (mesh->dy*mesh->dy));

  if(non_uniform) {
    // Correction for non-uniform mesh
    result += applyYdiff(f, fDDY, mesh->d2y);
  }

  return(result);
}

////////////// Z DERIVATIVE /////////////////

const Field3D D2DZ2(const Field3D &f, CELL_LOC outloc, DIFF_METHOD method)
{
  deriv_func func = fD2DZ2; // Set to default function
  DiffLookup *table = SecondDerivTable;
  
  CELL_LOC inloc = f.getLocation(); // Input location
  CELL_LOC diffloc = inloc; // Location of differential result
  
  Field3D result;

  if(mesh->StaggerGrids && (outloc == CELL_DEFAULT)) {
    // Take care of CELL_DEFAULT case
    outloc = diffloc; // No shift (i.e. same as no stagger case)
  }

  if(mesh->StaggerGrids && (outloc != inloc)) {
    // Shifting to a new location
    
    if(((inloc == CELL_CENTRE) && (outloc == CELL_ZLOW)) ||
       ((inloc == CELL_ZLOW) && (outloc == CELL_CENTRE))) {
      // Shifting in Z. Centre -> Zlow, or Zlow -> Centre
      
      func = sfD2DZ2; // Set default
      table = SecondStagDerivTable; // Set table for others
      diffloc = (inloc == CELL_CENTRE) ? CELL_ZLOW : CELL_CENTRE;
      
    }else {
      // A more complicated shift. Get a result at cell centre, then shift.
      if(inloc == CELL_ZLOW) {
	// Shifting
	
	func = sfD2DZ2; // Set default
	table = SecondStagDerivTable; // Set table for others
	diffloc = CELL_CENTRE;

      }else if(inloc != CELL_CENTRE) {
	// Interpolate then (centre -> centre) then interpolate
	return D2DZ2(interp_to(f, CELL_CENTRE), outloc, method);
      }
    }
  }

  if(method != DIFF_DEFAULT) {
    // Lookup function
    func = lookupFunc(table, method);
  }

  if(func == NULL) {
    // Use FFT

    BoutReal shift = 0.; // Shifting result in Z?
    if(mesh->StaggerGrids) {
      if((inloc == CELL_CENTRE) && (diffloc == CELL_ZLOW)) {
	// Shifting down - multiply by exp(-0.5*i*k*dz)
	shift = -1.;
      }else if((inloc == CELL_ZLOW) && (diffloc == CELL_CENTRE)) {
	// Shifting up
	shift = 1.;
      }
    }

    BoutReal flt;
    
    result.allocate(); // Make sure data allocated
    static dcomplex *cv = (dcomplex*) NULL;
    int jx, jy, jz;
    BoutReal kwave;

    int ncz = mesh->ngz-1;
    
    if(cv == (dcomplex*) NULL)
      cv = new dcomplex[ncz/2 + 1];

    for(jx=mesh->xstart;jx<=mesh->xend;jx++) {
      for(jy=mesh->ystart;jy<=mesh->yend;jy++) {

	rfft(f[jx][jy], ncz, cv); // Forward FFT
	
	for(jz=0;jz<=ncz/2;jz++) {
	  kwave=jz*2.0*PI/mesh->zlength; // wave number is 1/[rad]
	  
	  if (jz>0.4*ncz) flt=1e-10; else flt=1.0;

	  cv[jz] *= -SQ(kwave) * flt;
	  if(mesh->StaggerGrids)
	    cv[jz] *= exp(Im * (shift * kwave * mesh->dz));
	}

	irfft(cv, ncz, result[jx][jy]); // Reverse FFT
	
	result[jx][jy][ncz] = result[jx][jy][0];
      }
    }

#ifdef CHECK
    // Mark boundaries as invalid
    result.bndry_xin = result.bndry_xout = result.bndry_yup = result.bndry_ydown = false;
#endif

  }else {
    // All other (non-FFT) functions
    
    result = applyZdiff(f, func, SQ(mesh->dz));
  }

  result.setLocation(diffloc);

  return interp_to(result, outloc);
}

const Field3D D2DZ2(const Field3D &f, DIFF_METHOD method, CELL_LOC outloc)
{
  return D2DZ2(f, outloc, method);
}


const Field2D D2DZ2(const Field2D &f)
{
  Field2D result;
  result = 0.0;
  return result;
}

/*******************************************************************************
 * Mixed derivatives
 *******************************************************************************/

/// X-Z mixed derivative
const Field3D D2DXDZ(const Field3D &f)
{
  Field3D result;
  
  // Take derivative in Z, including in X boundaries. Then take derivative in X
  // Maybe should average results of DDX(DDZ) and DDZ(DDX)?
  result = DDX(DDZ(f, true));
  
  return result;
}

/*******************************************************************************
 * Advection schemes
 * 
 * Jan 2009  - Re-written to use Set*Stencil routines
 *******************************************************************************/

////////////// X DERIVATIVE /////////////////

/// Special case where both arguments are 2D. Output location ignored for now
const Field2D VDDX(const Field2D &v, const Field2D &f, CELL_LOC outloc, DIFF_METHOD method)
{
  upwind_func func = fVDDX;

  if(method != DIFF_DEFAULT) {
    // Lookup function
    func = lookupUpwindFunc(UpwindTable, method);
  }

  Field2D result;
  result.allocate(); // Make sure data allocated
  BoutReal **d = result.getData();

  bindex bx;
  stencil vs, fs;
  start_index(&bx);
  do {
    f.setXStencil(fs, bx);
    v.setXStencil(vs, bx);
    
    d[bx.jx][bx.jy] = func(vs, fs) / mesh->dx[bx.jx][bx.jy];
  }while(next_index2(&bx));

#ifdef CHECK
  // Mark boundaries as invalid
  result.bndry_xin = result.bndry_xout = false;
#endif

  return result;
}

const Field2D VDDX(const Field2D &v, const Field2D &f, DIFF_METHOD method)
{
  return VDDX(v, f, CELL_DEFAULT, method);
}

/// General version for 2 or 3-D objects
const Field3D VDDX(const Field &v, const Field &f, CELL_LOC outloc, DIFF_METHOD method)
{
  upwind_func func = fVDDX;
  DiffLookup *table = UpwindTable;

  CELL_LOC vloc = v.getLocation();
  CELL_LOC inloc = f.getLocation(); // Input location
  CELL_LOC diffloc = inloc; // Location of differential result
  
  if(mesh->StaggerGrids && (outloc == CELL_DEFAULT)) {
    // Take care of CELL_DEFAULT case
    outloc = diffloc; // No shift (i.e. same as no stagger case)
  }
  
  if(mesh->StaggerGrids && (vloc != inloc)) {
    // Staggered grids enabled, and velocity at different location to value
    if(vloc == CELL_XLOW) {
      // V staggered w.r.t. variable
      func = sfVDDX;
      table = UpwindStagTable;
      diffloc = CELL_CENTRE;
    }else if((vloc == CELL_CENTRE) && (inloc == CELL_XLOW)) {
      // Shifted
      func = sfVDDX;
      table = UpwindStagTable;
      diffloc = CELL_XLOW;
    }else {
      // More complicated. Deciding what to do here isn't straightforward
      // For now, interpolate velocity to the same location as f.

      // Should be able to do something like:
      //return VDDX(interp_to(v, inloc), f, outloc, method);
      
      // Instead, pretend it's been shifted FIX THIS
      diffloc = vloc;
    }
  }

  if(method != DIFF_DEFAULT) {
    // Lookup function
    func = lookupUpwindFunc(table, method);
  }

  /// Clone inputs (for shifting)
  Field *vp = v.clone();
  Field *fp = f.clone();
  
  if(mesh->ShiftXderivs && (mesh->ShiftOrder == 0)) {
    // Shift in Z using FFT if needed
    vp->shiftToReal(true);
    fp->shiftToReal(true);
  }
  
  Field3D result;
  result.allocate(); // Make sure data allocated
  BoutReal ***d = result.getData();

  bindex bx;
  stencil vval, fval;
  
  start_index(&bx);
  do {
    vp->setXStencil(vval, bx, diffloc);
    fp->setXStencil(fval, bx); // Location is always the same as input
    
    d[bx.jx][bx.jy][bx.jz] = func(vval, fval) / mesh->dx[bx.jx][bx.jy];
  }while(next_index3(&bx));
  
  if(mesh->ShiftXderivs && (mesh->ShiftOrder == 0))
    result = result.shiftZ(false); // Shift back
  
  result.setLocation(inloc);

#ifdef CHECK
  // Mark boundaries as invalid
  result.bndry_xin = result.bndry_xout = result.bndry_yup = result.bndry_ydown = false;
#endif
  
  // Delete clones
  delete vp;
  delete fp;
  
  return interp_to(result, outloc);
}

const Field3D VDDX(const Field &v, const Field &f, DIFF_METHOD method, CELL_LOC outloc)
{
  return VDDX(v, f, outloc, method);
}

////////////// Y DERIVATIVE /////////////////

// special case where both are 2D
const Field2D VDDY(const Field2D &v, const Field2D &f, CELL_LOC outloc, DIFF_METHOD method)
{
  upwind_func func = fVDDY;
  DiffLookup *table = UpwindTable;

  CELL_LOC vloc = v.getLocation();
  CELL_LOC inloc = f.getLocation(); // Input location
  CELL_LOC diffloc = inloc; // Location of differential result
  
  if(mesh->StaggerGrids && (outloc == CELL_DEFAULT)) {
    // Take care of CELL_DEFAULT case
    outloc = diffloc; // No shift (i.e. same as no stagger case)
  }

  if(mesh->StaggerGrids && (vloc != inloc)) {
    // Staggered grids enabled, and velocity at different location to value
    if(vloc == CELL_YLOW) {
      // V staggered w.r.t. variable
      func = sfVDDY;
      table = UpwindStagTable;
      diffloc = CELL_CENTRE;
    }else if((vloc == CELL_CENTRE) && (inloc == CELL_YLOW)) {
      // Shifted
      func = sfVDDY;
      table = UpwindStagTable;
      diffloc = CELL_YLOW;
    }else {
      // More complicated. Deciding what to do here isn't straightforward
      // For now, interpolate velocity to the same location as f.

      // Should be able to do something like:
      //return VDDY(interp_to(v, inloc), f, outloc, method);
      
      // Instead, pretend it's been shifted FIX THIS
      diffloc = vloc;
    }
  }

  if(method != DIFF_DEFAULT) {
    // Lookup function
    func = lookupUpwindFunc(table, method);
  }

  bindex bx;
  stencil fval, vval;
  
  Field2D result;
  result.allocate(); // Make sure data allocated
  BoutReal **d = result.getData();

  start_index(&bx);
  do {
    f.setYStencil(fval, bx);
    v.setYStencil(vval, bx, diffloc);
    d[bx.jx][bx.jy] = func(vval,fval)/mesh->dy[bx.jx][bx.jy];
  }while(next_index2(&bx));

  result.setLocation(inloc);
  
#ifdef CHECK
  // Mark boundaries as invalid
  result.bndry_xin = result.bndry_xout = result.bndry_yup = result.bndry_ydown = false;
#endif

  return interp_to(result, outloc);
}

const Field2D VDDY(const Field2D &v, const Field2D &f, DIFF_METHOD method)
{
  return VDDY(v, f, CELL_DEFAULT, method);
}

// general case
const Field3D VDDY(const Field &v, const Field &f, CELL_LOC outloc, DIFF_METHOD method)
{
  upwind_func func = fVDDY;
  DiffLookup *table = UpwindTable;

  CELL_LOC vloc = v.getLocation();
  CELL_LOC inloc = f.getLocation(); // Input location
  CELL_LOC diffloc = inloc; // Location of differential result
  
  if(mesh->StaggerGrids && (outloc == CELL_DEFAULT)) {
    // Take care of CELL_DEFAULT case
    outloc = diffloc; // No shift (i.e. same as no stagger case)
  }

  if(mesh->StaggerGrids && (vloc != inloc)) {
    // Staggered grids enabled, and velocity at different location to value
    if(vloc == CELL_YLOW) {
      // V staggered w.r.t. variable
      func = sfVDDY;
      table = UpwindStagTable;
      diffloc = CELL_CENTRE;
    }else if((vloc == CELL_CENTRE) && (inloc == CELL_YLOW)) {
      // Shifted
      func = sfVDDY;
      table = UpwindStagTable;
      diffloc = CELL_YLOW;
    }else {
      // More complicated. Deciding what to do here isn't straightforward
      // For now, interpolate velocity to the same location as f.

      // Should be able to do something like:
      //return VDDY(interp_to(v, inloc), f, outloc, method);
      
      // Instead, pretend it's been shifted FIX THIS
      diffloc = vloc;
    }
  }
  
  if(method != DIFF_DEFAULT) {
    // Lookup function
    func = lookupUpwindFunc(table, method);
  }
  bindex bx;
  stencil vval, fval;
  
  Field3D result;
  result.allocate(); // Make sure data allocated
  BoutReal ***d = result.getData();

  start_index(&bx);
  do {
    v.setYStencil(vval, bx, diffloc);
    f.setYStencil(fval, bx);
    
    d[bx.jx][bx.jy][bx.jz] = func(vval, fval)/mesh->dy[bx.jx][bx.jy];
  }while(next_index3(&bx));

  result.setLocation(inloc);

#ifdef CHECK
  // Mark boundaries as invalid
  result.bndry_xin = result.bndry_xout = result.bndry_yup = result.bndry_ydown = false;
#endif

  return interp_to(result, outloc);
}

const Field3D VDDY(const Field &v, const Field &f, DIFF_METHOD method, CELL_LOC outloc)
{
  return VDDY(v, f, outloc, method);
}

////////////// Z DERIVATIVE /////////////////

// special case where both are 2D
const Field2D VDDZ(const Field2D &v, const Field2D &f)
{
  Field2D result;
  result = 0.0;
  return result;
}

// Note that this is zero because no compression is included
const Field2D VDDZ(const Field3D &v, const Field2D &f)
{
  Field2D result;
  result = 0.0;
  return result;
}

// general case
const Field3D VDDZ(const Field &v, const Field &f, CELL_LOC outloc, DIFF_METHOD method)
{
  upwind_func func = fVDDZ;
  DiffLookup *table = UpwindTable;

  CELL_LOC vloc = v.getLocation();
  CELL_LOC inloc = f.getLocation(); // Input location
  CELL_LOC diffloc = inloc; // Location of differential result
  
  if(mesh->StaggerGrids && (outloc == CELL_DEFAULT)) {
    // Take care of CELL_DEFAULT case
    outloc = diffloc; // No shift (i.e. same as no stagger case)
  }

  if(mesh->StaggerGrids && (vloc != inloc)) {
    // Staggered grids enabled, and velocity at different location to value
    if(vloc == CELL_ZLOW) {
      // V staggered w.r.t. variable
      func = sfVDDZ;
      table = UpwindStagTable;
      diffloc = CELL_CENTRE;
    }else if((vloc == CELL_CENTRE) && (inloc == CELL_ZLOW)) {
      // Shifted
      func = sfVDDZ;
      table = UpwindStagTable;
      diffloc = CELL_ZLOW;
    }else {
      // More complicated. Deciding what to do here isn't straightforward
      // For now, interpolate velocity to the same location as f.

      // Should be able to do something like:
      //return VDDY(interp_to(v, inloc), f, outloc, method);
      
      // Instead, pretend it's been shifted FIX THIS
      diffloc = vloc;
    }
  }

  if(method != DIFF_DEFAULT) {
    // Lookup function
    func = lookupUpwindFunc(table, method);
  }

  bindex bx;
  stencil vval, fval;
  
  Field3D result;
  result.allocate(); // Make sure data allocated
  BoutReal ***d = result.getData();
  
  start_index(&bx);
  do {
    v.setZStencil(vval, bx, diffloc);
    f.setZStencil(fval, bx);
    
    d[bx.jx][bx.jy][bx.jz] = func(vval, fval)/mesh->dz;
  }while(next_index3(&bx));

  result.setLocation(inloc);

#ifdef CHECK
  // Mark boundaries as invalid
  result.bndry_xin = result.bndry_xout = result.bndry_yup = result.bndry_ydown = false;
#endif

  return interp_to(result, outloc);
}

const Field3D VDDZ(const Field &v, const Field &f, DIFF_METHOD method, CELL_LOC outloc)
{
  return VDDZ(v, f, outloc, method);
}
