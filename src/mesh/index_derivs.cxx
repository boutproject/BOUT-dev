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

#include <globals.hxx>
#include <derivs.hxx>
#include <stencils.hxx>
#include <utils.hxx>
#include <fft.hxx>
#include <interpolation.hxx>
#include <bout/constants.hxx>
#include <msg_stack.hxx>

#include <cmath>
#include <string.h>
#include <stdlib.h>

#include <output.hxx>

#include <bout/mesh.hxx>

/*******************************************************************************
 * Limiters
 *******************************************************************************/

/// Van Leer limiter. Used in TVD code
BoutReal VANLEER(BoutReal r) {
  return r + fabs(r)/(1.0 + fabs(r));
}

// Superbee limiter
BoutReal SUPERBEE(BoutReal r) {
  return BOUTMAX(0.0, BOUTMIN(2.*r, 1.0), BOUTMIN(r, 2.));
}

/*******************************************************************************
 * Basic derivative methods.
 * All expect to have an input grid cell at the same location as the output
 * Hence convert cell centred values -> centred values, or left -> left
 *******************************************************************************/

const BoutReal WENO_SMALL = 1.0e-8; // Small number for WENO schemes

////////////////////// FIRST DERIVATIVES /////////////////////

/// central, 2nd order
BoutReal DDX_C2(stencil &f) {
  return 0.5*(f.p - f.m);
}

/// forward, 2nd order
Mesh::boundary_derivs_pair DDX_F2(forward_stencil &f) {
  Mesh::boundary_derivs_pair result;
  result.inner = -1.5*f.c+2.*f.p-0.5*f.p2;
  result.outer = -1.5*f.m+2.*f.c-0.5*f.p;
  return result;
}

/// backward 2nd order
Mesh::boundary_derivs_pair DDX_B2(backward_stencil &f) {
  Mesh::boundary_derivs_pair result;
  result.inner = 0.5*f.m2-2.*f.m+1.5*f.c;
  result.outer = 0.5*f.m-2.*f.c+1.5*f.p;
  return result;
}

/// central, 4th order
BoutReal DDX_C4(stencil &f) {
  return (8.*f.p - 8.*f.m + f.mm - f.pp)/12.;
}

/// forward, 4th order
Mesh::boundary_derivs_pair DDX_F4(forward_stencil &f) {
  Mesh::boundary_derivs_pair result;
  result.inner = -1./4.*f.m-5./6.*f.c+3./2.*f.p-1./2.*f.p2+1./12.*f.p3; // uncentred (forward-biased) derivative
  result.outer = -25./12.*f.m+4.*f.c-3.*f.p+4./3.*f.p2-1./4.*f.p3; // forward derivative
  return result;
}

/// backward, 4th order
Mesh::boundary_derivs_pair DDX_B4(backward_stencil &f) {
  Mesh::boundary_derivs_pair result;
  result.inner = 1./4.*f.p+5./6.*f.c-3./2.*f.m+1./2.*f.m2-1./12.*f.m3; // uncentred (forward-biased) derivative
  result.outer = 25./12.*f.p-4.*f.c+3.*f.m-4./3.*f.m2+1./4.*f.m3; // forward derivative
  return result;
}

/// Central WENO method, 2nd order (reverts to 1st order near shocks)
BoutReal DDX_CWENO2(stencil &f) {
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

// Smoothing 2nd order derivative
BoutReal DDX_S2(stencil &f) {
  
  // 4th-order differencing
  BoutReal result = (8.*f.p - 8.*f.m + f.mm - f.pp)/12.;

  result += SIGN(f.c)*(f.pp - 4.*f.p + 6.*f.c - 4.*f.m + f.mm)/12.;
  
  return result;
}

///////////////////// SECOND DERIVATIVES ////////////////////

/// Second derivative: Central, 2nd order
BoutReal D2DX2_C2(stencil &f) {
  return f.p + f.m - 2.*f.c;
}

/// Second derivative: Forward, 2nd order
Mesh::boundary_derivs_pair D2DX2_F2(forward_stencil &f) {
  Mesh::boundary_derivs_pair result;
  result.inner = 2.*f.c-5.*f.p+4.*f.p2-3.*f.p3;
  result.outer = 2.*f.m-5.*f.c+4.*f.p-3.*f.p2;
  return result;
}

/// Second derivative: Backward, 2nd order
Mesh::boundary_derivs_pair D2DX2_B2(backward_stencil &f) {
  Mesh::boundary_derivs_pair result;
  result.inner = -2.*f.c+5.*f.m-4.*f.m2+3.*f.m3;
  result.outer = -2.*f.p+5.*f.c-4.*f.m+3.*f.m2;
  return result;
}

/// Second derivative: Central, 4th order
BoutReal D2DX2_C4(stencil &f) {
  return (-f.pp + 16.*f.p - 30.*f.c + 16.*f.m - f.mm)/12.;
}

/// Second derivatives: Forward, 4th order
Mesh::boundary_derivs_pair D2DX2_F4(forward_stencil &f) {
  Mesh::boundary_derivs_pair result;
  result.inner = 5./6.*f.m-5./4.*f.c-1./3.*f.p+7./6.*f.p2-1./2.*f.p3+1./12.*f.p4; // uncentred (forward-biased) derivative
  result.outer = 15./4.*f.m-77./6.*f.c+107./6.*f.p-13.*f.p2+61./12.*f.p3-5./6.*f.p4; // forward derivative
  return result;
}

/// Second derivatives: Backward, 4th order
Mesh::boundary_derivs_pair D2DX2_B4(backward_stencil &f) {
  Mesh::boundary_derivs_pair result;
  result.inner = 5./6.*f.p-5./4.*f.c-1./3.*f.m+7./6.*f.m2-1./2.*f.m3+1./12.*f.m4; // uncentred (backward-biased) derivative
  result.outer = 15./4.*f.p-77./6.*f.c+107./6.*f.m-13.*f.m2+61./12.*f.m3-5./6.*f.m4; // backward derivative
  return result;
}

//////////////////////// UPWIND METHODS ///////////////////////

/// Upwinding: Central, 2nd order
BoutReal VDDX_C2(BoutReal vc, stencil &f) {
  return vc*0.5*(f.p - f.m);
}

/// Upwinding: Central, 4th order
BoutReal VDDX_C4(BoutReal vc, stencil &f) {
  return vc*(8.*f.p - 8.*f.m + f.mm - f.pp)/12.;
}

/// upwind, 1st order
BoutReal VDDX_U1(BoutReal vc, stencil &f) {
  return vc>=0.0 ? vc*(f.c - f.m): vc*(f.p - f.c);
}

/// upwind, 2nd order
BoutReal VDDX_U2(BoutReal vc, stencil &f) {
  return vc>=0.0 ? vc*(1.5*f.c - 2.0*f.m + 0.5*f.mm): vc*(-0.5*f.pp + 2.0*f.p - 1.5*f.c);
}

/// upwind, 4th order
BoutReal VDDX_U4(BoutReal vc, stencil &f) {
  return vc >= 0.0 ? vc*(4.*f.p - 12.*f.m + 2.*f.mm + 6.*f.c)/12.
    : vc*(-4.*f.m + 12.*f.p - 2.*f.pp - 6.*f.c)/12.;
}

/// 3rd-order WENO scheme
BoutReal VDDX_WENO3(BoutReal vc, stencil &f) {
  BoutReal deriv, w, r;

  if(vc > 0.0) {
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

  return vc*deriv;
}

/// 3rd-order CWENO. Uses the upwinding code and split flux
BoutReal DDX_CWENO3(stencil &f) {
  BoutReal a, ma = fabs(f.c);
  // Split flux
  a = fabs(f.m); if(a > ma) ma = a;
  a = fabs(f.p); if(a > ma) ma = a;
  a = fabs(f.mm); if(a > ma) ma = a;
  a = fabs(f.pp); if(a > ma) ma = a;
  
  stencil sp, vp, sm, vm;
  
  sp = f + ma;
  sm = ma - f;
  
  return VDDX_WENO3(0.5, sp) + VDDX_WENO3(-0.5, sm);
}

//////////////////////// FLUX METHODS ///////////////////////

BoutReal FDDX_U1(stencil &v, stencil &f) {
  // Velocity at lower end
  BoutReal vs = 0.5*(v.m + v.c);
  BoutReal result = (vs >= 0.0) ? vs * f.m : vs * f.c;
  // and at upper 
  vs = 0.5*(v.c + v.p);
  result -= (vs >= 0.0) ? vs * f.c : vs * f.p;

  return result;
}

BoutReal FDDX_C2(stencil &v, stencil &f) {
  return 0.5*(v.p*f.p - v.m*f.m);
}

BoutReal FDDX_C4(stencil &v, stencil &f) {
  return (8.*v.p*f.p - 8.*v.m*f.m + v.mm*f.mm - v.pp*f.pp)/12.;
}

/// Non-oscillatory, containing No free parameters and Dissipative (NND) scheme
/// http://arxiv.org/abs/1010.4135v1
BoutReal FDDX_NND(stencil &v, stencil &f) {
  // f{+-} i
  BoutReal fp = 0.5*(v.c + fabs(v.c))*f.c;
  BoutReal fm = 0.5*(v.c - fabs(v.c))*f.c;
  
  // f{+-} i+1
  BoutReal fp1 = 0.5*(v.p + fabs(v.p))*f.p;
  BoutReal fm1 = 0.5*(v.p - fabs(v.p))*f.p;
  
  // f{+-} i+2
  BoutReal fm2 = 0.5*(v.pp - fabs(v.pp))*f.pp;

  // f{+-} i-1
  BoutReal fp_1 = 0.5*(v.m + fabs(v.m))*f.m;
  BoutReal fm_1 = 0.5*(v.m - fabs(v.m))*f.m;
  
  // f{+-} i-2
  BoutReal fp_2 = 0.5*(v.mm + fabs(v.mm))*f.mm;

  // f^{LR} {i+1/2}
  BoutReal flp = fp  + 0.5*MINMOD(fp1 - fp, fp - fp_1);
  BoutReal frp = fm1 - 0.5*MINMOD(fm1 - fm, fm2 - fm1);
  
  // f^{LR} {i-1/2}
  BoutReal flm = fp_1  + 0.5*MINMOD(fp - fp_1, fp_1 - fp_2);
  BoutReal frm = fm - 0.5*MINMOD(fm - fm_1, fm1 - fm);
    
  // h{+-}
  BoutReal hp = flp + frp;
  BoutReal hm = flm + frm; 
  
  return hp - hm;
}

//////////////////////// MUSCL scheme ///////////////////////

void DDX_KT_LR(const stencil &f, BoutReal &fLp, BoutReal &fRp, BoutReal &fLm, BoutReal &fRm) {
  // Limiter functions
  BoutReal phi   = SUPERBEE( (f.c - f.m) / (f.p - f.c) );
  BoutReal phi_m = SUPERBEE( (f.m - f.mm) / (f.c - f.m) );
  BoutReal phi_p = SUPERBEE( (f.p - f.c) / (f.pp - f.p) );

  fLp = f.c + 0.5*phi*(f.p - f.c);
  fRp = f.p - 0.5*phi_p*(f.pp - f.p);

  fLm = f.m + 0.5*phi_m*(f.c - f.m);
  fRm = f.c - 0.5*phi*(f.p - f.c);
}

// du/dt = d/dx(f)  with maximum local velocity Vmax
BoutReal DDX_KT(const stencil &f, const stencil &u, const BoutReal Vmax) {
  BoutReal uLp, uRp, uLm, uRm;
  BoutReal fLp, fRp, fLm, fRm;
  
  DDX_KT_LR(u, uLp, uRp, uLm, uRm);
  DDX_KT_LR(f, fLp, fRp, fLm, fRm);
  
  BoutReal Fm = 0.5*( fRm + fLm - Vmax*(uRm - uLm));
  BoutReal Fp = 0.5*( fRp + fLp - Vmax*(uRp - uLp));
  
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
BoutReal DDX_C2_stag(stencil &f) {
  return f.p - f.m;
}

Mesh::boundary_derivs_pair DDX_F2_stag(forward_stencil &f) {
  Mesh::boundary_derivs_pair result;
  result.inner = -2.*f.c+3*f.p-f.p2;
  result.outer = -2.*f.m+3*f.c-f.p;
  return result;
}

Mesh::boundary_derivs_pair DDX_B2_stag(backward_stencil &f) {
  Mesh::boundary_derivs_pair result;
  result.inner = 2.*f.c-3*f.m+f.m2;
  result.outer = 2.*f.p-3*f.c+f.m;
  return result;
}

BoutReal DDX_C4_stag(stencil &f) {
  return ( 27.*(f.p - f.m) - (f.pp - f.mm) ) / 24.;
}

Mesh::boundary_derivs_pair DDX_F4_stag(forward_stencil &f) {
  Mesh::boundary_derivs_pair result;
  result.inner = -11./12.*f.m+17./24.*f.c+3./8.*f.p-5./24.*f.p2+1./24.*f.p3; // uncentred (forward-biased) derivative
  result.outer = -31./8*f.m+229./24.*f.c-75./8.*f.p+37./8.*f.p2-11./12.*f.p3; // forward derivative
  return result;
}

Mesh::boundary_derivs_pair DDX_B4_stag(backward_stencil &f) {
  Mesh::boundary_derivs_pair result;
  result.inner = 11./12.*f.p-17./24.*f.c-3./8.*f.m+5./24.*f.m2-1./24.*f.m3; // uncentred (backward-biased) derivative
  result.outer = 31./8*f.p-229./24.*f.c+75./8.*f.m-37./8.*f.m2+11./12.*f.m3; // backward derivative
  return result;
}

/////////////////////// SECOND DERIVATIVES //////////////////////
// Map Centre -> Low or Low -> Centre

Mesh::boundary_derivs_pair D2DX2_F2_stag(forward_stencil &f) {
  Mesh::boundary_derivs_pair result;
  result.inner = 5./2.*f.c-13./2.*f.p+11./2.*f.p2-3./2.*f.p3;
  result.outer = 5./2.*f.m-13./2.*f.c+11./2.*f.p-3./2.*f.p2;
  return result;
}

Mesh::boundary_derivs_pair D2DX2_B2_stag(backward_stencil &f) {
  Mesh::boundary_derivs_pair result;
  result.inner = 5./2.*f.c-13./2.*f.m+11./2.*f.m2-3./2.*f.m3;
  result.outer = 5./2.*f.p-13./2.*f.c+11./2.*f.m-3./2.*f.m2;
  return result;
}

BoutReal D2DX2_C4_stag(stencil &f) {
  return ( f.pp + f.mm - f.p - f.m ) / 2.;
}

Mesh::boundary_derivs_pair D2DX2_F4_stag(forward_stencil &f) {
  Mesh::boundary_derivs_pair result;
  result.inner = 95./48.*f.m-269./48.*f.c+49./8.*f.p-85./24.*f.p2+59./48.*f.p3-3./16.*f.p4;
  result.outer = 301./48.*f.m-377./16.*f.c+865./24.*f.p-683./24.*f.p2+187./16.*f.p3-95./48.*f.p4;
  return result;
}

Mesh::boundary_derivs_pair D2DX2_B4_stag(backward_stencil &f) {
  Mesh::boundary_derivs_pair result;
  result.inner = 95./48.*f.p-269./48.*f.c+49./8.*f.m-85./24.*f.m2+59./48.*f.m3-3./16.*f.m4;
  result.outer = 301./48.*f.p-377./16.*f.c+865./24.*f.m-683./24.*f.m2+187./16.*f.m3-95./48.*f.m4;
  return result;
}

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
  result -= f.c*(v.p - v.m);
  
  return result;
}

BoutReal VDDX_U2_stag(stencil &v, stencil &f) {
  BoutReal result;
  
  if (v.p>0 && v.m>0) {
    // Extrapolate v to centre from below, use 2nd order backward difference on f
    result = (1.5*v.m - .5*v.mm) * (.5*f.mm - 2.*f.m + 1.5*f.c);
  }
  else if (v.p<0 && v.m<0) {
    // Extrapolate v to centre from above, use 2nd order forward difference on f
    result = (1.5*v.p - .5*v.pp) * (-1.5*f.c + 2.*f.p - .5*f.pp);
  }
  else {
    // Velocity changes sign, hence is almost zero: use centred interpolation/differencing
    result = .25 * (v.p + v.m) * (f.p - f.m);
  }
  
  return result;
}

BoutReal VDDX_C2_stag(stencil &v, stencil &f) {
  // Result is needed at location of f: interpolate v to f's location and take an unstaggered derivative of f
  return 0.5*(v.p+v.m) * 0.5*(f.p - f.m);
}

BoutReal VDDX_C4_stag(stencil &v, stencil &f) {
  // Result is needed at location of f: interpolate v to f's location and take an unstaggered derivative of f
  return (9.*(v.m + v.p) - v.mm - v.pp)/16. * (8.*f.p - 8.*f.m + f.mm - f.pp)/12.;
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

  return result;
}

/*******************************************************************************
 * Lookup tables of functions. Map between names, codes and functions
 *******************************************************************************/

/// Translate between DIFF_METHOD codes, and functions
struct DiffLookup {
  DIFF_METHOD method;
  Mesh::deriv_func func;     // Single-argument differencing function
  Mesh::inner_boundary_deriv_func inner_boundary_func; // Differencing function using forward derivatives
  Mesh::outer_boundary_deriv_func outer_boundary_func; // Differencing function using backward derivatives
  Mesh::upwind_func up_func; // Upwinding function
  Mesh::flux_func fl_func; // Flux function
  Mesh::inner_boundary_upwind_func inner_boundary_up_func; // Upwinding function using forward derivatives
  Mesh::outer_boundary_upwind_func outer_boundary_up_func; // Upwinding function using backward derivatives
};

/// Translate between short names, long names and DIFF_METHOD codes
struct DiffNameLookup {
  DIFF_METHOD method;
  const char* label; // Short name
  const char* name;  // Long name
};

/// Differential function name/code lookup
static DiffNameLookup DiffNameTable[] = { {DIFF_U1, "U1", "First order upwinding"},
					  {DIFF_U2, "U2", "Second order upwinding"},
					  {DIFF_C2, "C2", "Second order central"},
					  {DIFF_W2, "W2", "Second order WENO"},
					  {DIFF_W3, "W3", "Third order WENO"},
					  {DIFF_C4, "C4", "Fourth order central"},
					  {DIFF_U4, "U4", "Fourth order upwinding"},
                      {DIFF_S2, "S2", "Smoothing 2nd order"},
					  {DIFF_FFT, "FFT", "FFT"},
                      {DIFF_NND, "NND", "NND"},
                      {DIFF_SPLIT, "SPLIT", "Split into upwind and central"},
					  {DIFF_DEFAULT, NULL, NULL}}; // Use to terminate the list

/// First derivative lookup table
static DiffLookup FirstDerivTable[] = { {DIFF_C2, DDX_C2,     DDX_F2, DDX_B2, NULL, NULL, NULL, NULL},
					{DIFF_W2, DDX_CWENO2, DDX_F2, DDX_B2, NULL, NULL, NULL, NULL},
					{DIFF_W3, DDX_CWENO3, DDX_F4, DDX_B4, NULL, NULL, NULL, NULL},
					{DIFF_C4, DDX_C4,     DDX_F4, DDX_B4, NULL, NULL, NULL, NULL},
                                        {DIFF_S2, DDX_S2,     NULL,   NULL,   NULL, NULL, NULL, NULL},
					{DIFF_FFT, NULL,      NULL,   NULL,   NULL, NULL, NULL, NULL},
					{DIFF_DEFAULT, NULL,  NULL,   NULL,   NULL, NULL, NULL, NULL}};

/// Second derivative lookup table
static DiffLookup SecondDerivTable[] = { {DIFF_C2, D2DX2_C2, D2DX2_F2, D2DX2_B2, NULL, NULL, NULL, NULL},
					 {DIFF_C4, D2DX2_C4, D2DX2_F4, D2DX2_B4, NULL, NULL, NULL, NULL},
					 {DIFF_FFT, NULL,    NULL,     NULL,     NULL, NULL, NULL, NULL},
					 {DIFF_DEFAULT, NULL,NULL,     NULL,     NULL, NULL, NULL, NULL}};

/// Upwinding functions lookup table
static DiffLookup UpwindTable[] = { {DIFF_U1, NULL, NULL, NULL, VDDX_U1, NULL, NULL, NULL},
					{DIFF_U2, NULL, NULL, NULL, VDDX_U2, NULL, NULL, NULL}, 
				    {DIFF_C2, NULL, NULL, NULL, VDDX_C2, NULL, NULL, NULL},
				    {DIFF_U4, NULL, NULL, NULL, VDDX_U4, NULL, NULL, NULL},
				    {DIFF_W3, NULL, NULL, NULL, VDDX_WENO3, NULL, NULL, NULL},
				    {DIFF_C4, NULL, NULL, NULL, VDDX_C4, NULL, NULL, NULL},
				    {DIFF_DEFAULT, NULL, NULL, NULL, NULL, NULL, NULL, NULL}};

/// Flux functions lookup table
static DiffLookup FluxTable[] = { {DIFF_SPLIT, NULL, NULL, NULL, NULL, NULL, NULL, NULL},
                                  {DIFF_U1, NULL, NULL, NULL, NULL, FDDX_U1, NULL, NULL},
                                  {DIFF_C2, NULL, NULL, NULL, NULL, FDDX_C2, NULL, NULL},
                                  {DIFF_C4, NULL, NULL, NULL, NULL, FDDX_C4, NULL, NULL},
                                  {DIFF_NND, NULL, NULL, NULL, NULL, FDDX_NND, NULL, NULL},
                                  {DIFF_DEFAULT, NULL, NULL, NULL, NULL, NULL, NULL, NULL}};

/// First staggered derivative lookup
static DiffLookup FirstStagDerivTable[] = { {DIFF_C2, DDX_C2_stag, DDX_F2_stag, DDX_B2_stag, NULL, NULL, NULL, NULL}, 
					    {DIFF_C4, DDX_C4_stag, DDX_F4_stag, DDX_B4_stag, NULL, NULL, NULL, NULL},
					    {DIFF_DEFAULT, NULL, NULL, NULL, NULL, NULL, NULL, NULL}};

/// Second staggered derivative lookup
static DiffLookup SecondStagDerivTable[] = { {DIFF_C4, D2DX2_C4_stag, D2DX2_F4_stag, D2DX2_B4_stag, NULL, NULL, NULL, NULL},
					     {DIFF_DEFAULT, NULL, NULL, NULL, NULL, NULL, NULL, NULL}};

/// Upwinding staggered lookup
static DiffLookup UpwindStagTable[] = { {DIFF_U1, NULL, NULL, NULL, NULL, VDDX_U1_stag, NULL, NULL},
					{DIFF_U2, NULL, NULL, NULL, NULL, VDDX_U2_stag, NULL, NULL},
					{DIFF_C2, NULL, NULL, NULL, NULL, VDDX_C2_stag, NULL, NULL},
					{DIFF_C4, NULL, NULL, NULL, NULL, VDDX_C4_stag, NULL, NULL},
					{DIFF_DEFAULT, NULL, NULL, NULL, NULL, NULL, NULL, NULL} };

/// Flux staggered lookup
static DiffLookup FluxStagTable[] = { {DIFF_SPLIT, NULL, NULL, NULL, NULL, NULL, NULL, NULL},
                                      {DIFF_U1, NULL, NULL, NULL, NULL, FDDX_U1_stag, NULL, NULL},
                                      {DIFF_DEFAULT, NULL, NULL, NULL, NULL, NULL, NULL, NULL}};

/*******************************************************************************
 * Routines to use the above tables to map between function codes, names
 * and pointers
 *******************************************************************************/

Mesh::deriv_func lookupFunc(DiffLookup* table, DIFF_METHOD method) {
  int i = 0;
  do {
    if(table[i].method == method)
      return table[i].func;
    i++;
  }while(table[i].method != DIFF_DEFAULT);
  // Not found in list. Return the first 
  
  return table[0].func;
}

Mesh::inner_boundary_deriv_func lookupInnerBoundaryFunc(DiffLookup* table, DIFF_METHOD method) {
  int i = 0;
  do {
    if(table[i].method == method)
      return table[i].inner_boundary_func;
    i++;
  }while(table[i].method != DIFF_DEFAULT);
  // Not found in list. Return the first 
  
  return table[0].inner_boundary_func;
}

Mesh::outer_boundary_deriv_func lookupOuterBoundaryFunc(DiffLookup* table, DIFF_METHOD method) {
  int i = 0;
  do {
    if(table[i].method == method)
      return table[i].outer_boundary_func;
    i++;
  }while(table[i].method != DIFF_DEFAULT);
  // Not found in list. Return the first 
  
  return table[0].outer_boundary_func;
}

Mesh::upwind_func lookupUpwindFunc(DiffLookup* table, DIFF_METHOD method) {
  int i = 0;
  do {
    if(table[i].method == method)
      return table[i].up_func;
    i++;
  }while(table[i].method != DIFF_DEFAULT);
  // Not found in list. Return the first 
  
  return table[0].up_func;
}

Mesh::flux_func lookupFluxFunc(DiffLookup* table, DIFF_METHOD method) {
  int i = 0;
  do {
    if(table[i].method == method)
      return table[i].fl_func;
    i++;
  }while(table[i].method != DIFF_DEFAULT);
  // Not found in list. Return the first 
  
  return table[0].fl_func;
}

Mesh::inner_boundary_upwind_func lookupInnerBoundaryUpwindFunc(DiffLookup* table, DIFF_METHOD method) {
  int i = 0;
  do {
    if(table[i].method == method)
      return table[i].inner_boundary_up_func;
    i++;
  }while(table[i].method != DIFF_DEFAULT);
  // Not found in list. Return the first 
  
  return table[0].inner_boundary_up_func;
}

Mesh::outer_boundary_upwind_func lookupOuterBoundaryUpwindFunc(DiffLookup* table, DIFF_METHOD method) {
  int i = 0;
  do {
    if(table[i].method == method)
      return table[i].outer_boundary_up_func;
    i++;
  }while(table[i].method != DIFF_DEFAULT);
  // Not found in list. Return the first 
  
  return table[0].outer_boundary_up_func;
}

/// Test if a given DIFF_METHOD exists in a table
bool isImplemented(DiffLookup* table, DIFF_METHOD method) {
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
DIFF_METHOD lookupFunc(DiffLookup *table, const string &label) {

  if(label.empty())
    return table[0].method;

  // Loop through the name lookup table
  for (int i = 0; DiffNameTable[i].method != DIFF_DEFAULT ; ++i ){
    if(strcasecmp(label.c_str(), DiffNameTable[i].label) == 0) {// Whole match
      return  DiffNameTable[i].method;
    }
  }
  
  // No exact match, so throw
  std::string avail{};
  for (int i = 0; DiffNameTable[i].method != DIFF_DEFAULT ; ++i ){
    avail += DiffNameTable[i].label;
    avail += "\n";
  }
  throw BoutException("Unknown option %s.\nAvailable options are:\n%s",label.c_str(),avail.c_str());
}

void printFuncName(DIFF_METHOD method) {
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
Mesh::deriv_func fDDX, fDDY, fDDZ;        ///< Differencing methods for each dimension
Mesh::deriv_func fD2DX2, fD2DY2, fD2DZ2;  ///< second differential operators
Mesh::upwind_func fVDDX, fVDDY, fVDDZ;    ///< Upwind functions in the three directions
Mesh::flux_func fFDDX, fFDDY, fFDDZ;    ///< Default flux functions
Mesh::inner_boundary_deriv_func fDDX_in, fDDY_in; ///< Differencing methods in the inner boundaries
Mesh::outer_boundary_deriv_func fDDX_out, fDDY_out; ///< Differencing methods in the outer boundaries
Mesh::inner_boundary_deriv_func fD2DX2_in, fD2DY2_in; ///< Second derivative methods in the inner boundaries
Mesh::outer_boundary_deriv_func fD2DX2_out, fD2DY2_out; ///< Second derivative methods in the outer boundaries
Mesh::inner_boundary_upwind_func fVDDX_in, fVDDY_in;    ///< Upwind functions in the inner boundaries
Mesh::outer_boundary_upwind_func fVDDX_out, fVDDY_out;    ///< Upwind functions in the outer boundaries
Mesh::inner_boundary_upwind_func fFDDX_in, fFDDY_in;    ///< Default flux functions in the inner boundaries
Mesh::outer_boundary_upwind_func fFDDX_out, fFDDY_out;    ///< Default flux functions in the outer boundaries

// Central -> Left (or Left -> Central) functions
Mesh::deriv_func sfDDX, sfDDY, sfDDZ;
Mesh::deriv_func sfD2DX2, sfD2DY2, sfD2DZ2;
Mesh::flux_func sfVDDX, sfVDDY, sfVDDZ;
Mesh::flux_func sfFDDX, sfFDDY, sfFDDZ;
Mesh::inner_boundary_deriv_func sfDDX_in, sfDDY_in;
Mesh::outer_boundary_deriv_func sfDDX_out, sfDDY_out;
Mesh::inner_boundary_deriv_func sfD2DX2_in, sfD2DY2_in;
Mesh::outer_boundary_deriv_func sfD2DX2_out, sfD2DY2_out;
Mesh::inner_boundary_upwind_func sfVDDX_in, sfVDDY_in;
Mesh::outer_boundary_upwind_func sfVDDX_out, sfVDDY_out;
Mesh::inner_boundary_upwind_func sfFDDX_in, sfFDDY_in;
Mesh::outer_boundary_upwind_func sfFDDX_out, sfFDDY_out;

/*******************************************************************************
 * Initialisation
 *******************************************************************************/

/// Set the derivative method, given a table and option name
void derivs_set(Options *options, DiffLookup *table, const char* name, Mesh::deriv_func &f) {
  string label;
  options->get(name, label, "", false);

  DIFF_METHOD method = lookupFunc(table, label); // Find the function
  printFuncName(method); // Print differential function name
  f = lookupFunc(table, method); // Find the function pointer
}

void derivs_set(Options *options, DiffLookup *table, const char* name, Mesh::upwind_func &f) {
  string label;
  options->get(name, label, "", false);

  DIFF_METHOD method = lookupFunc(table, label); // Find the function
  printFuncName(method); // Print differential function name
  f = lookupUpwindFunc(table, method);
}

void derivs_set(Options *options, DiffLookup *table, const char* name, Mesh::flux_func &f) {
  string label;
  options->get(name, label, "", false);

  DIFF_METHOD method = lookupFunc(table, label); // Find the function
  printFuncName(method); // Print differential function name
  f = lookupFluxFunc(table, method);
}

/// Set the derivative methods including for boundaries, given a table and option name
void derivs_set(Options *options, DiffLookup *table, const char* name, Mesh::deriv_func &f, Mesh::inner_boundary_deriv_func &f_in, Mesh::outer_boundary_deriv_func &f_out) {
  string label;
  options->get(name, label, "", false);

  DIFF_METHOD method = lookupFunc(table, label); // Find the function
  printFuncName(method); // Print differential function name
  f = lookupFunc(table, method); // Find the function pointers
  f_in = lookupInnerBoundaryFunc(table, method);
  f_out = lookupOuterBoundaryFunc(table, method);
}

void derivs_set(Options *options, DiffLookup *table, const char* name, Mesh::upwind_func &f, Mesh::inner_boundary_upwind_func &f_in, Mesh::outer_boundary_upwind_func &f_out) {
  string label;
  options->get(name, label, "", false);

  DIFF_METHOD method = lookupFunc(table, label); // Find the function
  printFuncName(method); // Print differential function name
  f = lookupUpwindFunc(table, method);
  f_in = lookupInnerBoundaryUpwindFunc(table, method);
  f_out = lookupOuterBoundaryUpwindFunc(table, method);
}

void derivs_set(Options *options, DiffLookup *table, const char* name, Mesh::flux_func &f, Mesh::inner_boundary_upwind_func &f_in, Mesh::outer_boundary_upwind_func &f_out) {
  string label;
  options->get(name, label, "", false);

  DIFF_METHOD method = lookupFunc(table, label); // Find the function
  printFuncName(method); // Print differential function name
  f = lookupFluxFunc(table, method);
  f_in = lookupInnerBoundaryUpwindFunc(table, method);
  f_out = lookupOuterBoundaryUpwindFunc(table, method);
}

/// Initialise derivatives from options
void derivs_initialise(Options *options, bool StaggerGrids,
                 Mesh::deriv_func &fdd, Mesh::deriv_func &sfdd, 
                 Mesh::deriv_func &fd2d, Mesh::deriv_func &sfd2d, 
                 Mesh::upwind_func &fu, Mesh::flux_func &sfu,
                 Mesh::flux_func &ff, Mesh::flux_func &sff) {
  output.write("\tFirst       : ");
  derivs_set(options, FirstDerivTable, "first",  fdd);
  if(StaggerGrids) {
    output.write("\tStag. First : ");
    derivs_set(options, FirstStagDerivTable, "first",  sfdd);
  }
  output.write("\tSecond      : ");
  derivs_set(options, SecondDerivTable, "second", fd2d);
  if(StaggerGrids) {
    output.write("\tStag. Second: ");
    derivs_set(options, SecondStagDerivTable, "second", sfd2d);
  }
  output.write("\tUpwind      : ");
  derivs_set(options, UpwindTable,     "upwind", fu);
  if(StaggerGrids) {
    output.write("\tStag. Upwind: ");
    derivs_set(options, UpwindStagTable,     "upwind", sfu);
  }
  output.write("\tFlux        : ");
  derivs_set(options, FluxTable,     "flux", ff);
  if(StaggerGrids) {
    output.write("\tStag. Flux  : ");
    derivs_set(options, FluxStagTable,     "flux", sff);
  }
}

void derivs_initialise(Options *options, bool StaggerGrids,
                 Mesh::deriv_func &fdd, Mesh::deriv_func &sfdd, 
                 Mesh::deriv_func &fd2d, Mesh::deriv_func &sfd2d, 
                 Mesh::upwind_func &fu, Mesh::flux_func &sfu,
                 Mesh::flux_func &ff, Mesh::flux_func &sff,
		 Mesh::inner_boundary_deriv_func &inner_fdd, Mesh::inner_boundary_deriv_func &inner_sfdd,
		 Mesh::outer_boundary_deriv_func &outer_fdd, Mesh::outer_boundary_deriv_func &outer_sfdd,
		 Mesh::inner_boundary_deriv_func &inner_fd2d, Mesh::inner_boundary_deriv_func &inner_sfd2d,
		 Mesh::outer_boundary_deriv_func &outer_fd2d, Mesh::outer_boundary_deriv_func &outer_sfd2d, 
                 Mesh::inner_boundary_upwind_func &inner_fu, Mesh::inner_boundary_upwind_func &inner_sfu,
                 Mesh::outer_boundary_upwind_func &outer_fu, Mesh::outer_boundary_upwind_func &outer_sfu,
                 Mesh::inner_boundary_upwind_func &inner_ff, Mesh::inner_boundary_upwind_func &inner_sff,
                 Mesh::outer_boundary_upwind_func &outer_ff, Mesh::outer_boundary_upwind_func &outer_sff
		) {
  output.write("\tFirst       : ");
  derivs_set(options, FirstDerivTable, "first",  fdd, inner_fdd, outer_fdd);
  if(StaggerGrids) {
    output.write("\tStag. First : ");
    derivs_set(options, FirstStagDerivTable, "first",  sfdd, inner_sfdd, outer_sfdd);
  }
  output.write("\tSecond      : ");
  derivs_set(options, SecondDerivTable, "second", fd2d, inner_fd2d, outer_fd2d);
  if(StaggerGrids) {
    output.write("\tStag. Second: ");
    derivs_set(options, SecondStagDerivTable, "second", sfd2d, inner_sfd2d, outer_sfd2d);
  }
  output.write("\tUpwind      : ");
  derivs_set(options, UpwindTable,     "upwind", fu, inner_fu, outer_fu);
  if(StaggerGrids) {
    output.write("\tStag. Upwind: ");
    derivs_set(options, UpwindStagTable,     "upwind", sfu, inner_sfu, outer_sfu);
  }
  output.write("\tFlux        : ");
  derivs_set(options, FluxTable,     "flux", ff, inner_ff, outer_ff);
  if(StaggerGrids) {
    output.write("\tStag. Flux  : ");
    derivs_set(options, FluxStagTable,     "flux", sff, inner_sff, outer_sff);
  }
}

/// Initialise the derivative methods. Must be called before any derivatives are used
void Mesh::derivs_init(Options* options) {
#ifdef CHECK
  msg_stack.push("Initialising derivatives");
#endif

  output.write("Setting X differencing methods\n");
  derivs_initialise(options->getSection("ddx"), 
              StaggerGrids,
              fDDX, sfDDX, 
              fD2DX2, sfD2DX2,
              fVDDX, sfVDDX,
              fFDDX, sfFDDX,
	      fDDX_in, sfDDX_in,
	      fDDX_out, sfDDX_out,
	      fD2DX2_in, sfD2DX2_in,
	      fD2DX2_out, sfD2DX2_out,
              fVDDX_in, sfVDDX_in,
              fVDDX_out, sfVDDX_out,
              fFDDX_in, sfFDDX_in,
              fFDDX_out, sfFDDX_out
 	    );
  
  if((fDDX == NULL) || (fD2DX2 == NULL))
    throw BoutException("FFT cannot be used in X\n");
  
  output.write("Setting Y differencing methods\n");
  derivs_initialise(options->getSection("ddy"), 
              StaggerGrids,
              fDDY, sfDDY, 
              fD2DY2, sfD2DY2,
              fVDDY, sfVDDY,
              fFDDY, sfFDDY,
	      fDDY_in, sfDDY_in,
	      fDDY_out, sfDDY_out,
	      fD2DY2_in, sfD2DY2_in,
	      fD2DY2_out, sfD2DY2_out,
              fVDDX_in, sfVDDX_in,
              fVDDX_out, sfVDDX_out,
              fFDDX_in, sfFDDX_in,
              fFDDX_out, sfFDDX_out);
  
  if((fDDY == NULL) || (fD2DY2 == NULL))
    throw BoutException("FFT cannot be used in Y\n");
  
  output.write("Setting Z differencing methods\n");
  derivs_initialise(options->getSection("ddz"), 
              StaggerGrids,
              fDDZ, sfDDZ, 
              fD2DZ2, sfD2DZ2,
              fVDDZ, sfVDDZ,
              fFDDZ, sfFDDZ);
  
#ifdef CHECK
  msg_stack.pop();
#endif
}

/*******************************************************************************
 * Apply differential operators. These are fairly brain-dead functions
 * which apply a derivative function to a field (sort of like map). Decisions
 * of what to apply are made in the DDX,DDY and DDZ functions lower down.
 *
 * loc  is the cell location of the result
 *******************************************************************************/

// X derivative

const Field2D Mesh::applyXdiff(const Field2D &var, Mesh::deriv_func func, Mesh::inner_boundary_deriv_func func_in, Mesh::outer_boundary_deriv_func func_out, CELL_LOC loc) {
  if (var.getNx() == 1){
    return 0.;
  }
  Field2D result;
  result.allocate(); // Make sure data allocated

  bindex bx;

  start_index(&bx, RGN_NOX);
  stencil s;
  do {
    var.setXStencil(s, bx, loc);
    result(bx.jx,bx.jy) = func(s);
  }while(next_index2(&bx));

#ifdef CHECK
  // Mark boundaries as invalid
  result.bndry_xin = result.bndry_xout = result.bndry_yup = result.bndry_ydown = false;
#endif
  
  if (mesh->freeboundary_ydown) {
    for (RangeIterator it=mesh->iterateBndryLowerY(); !it.isDone(); it++)
      for (bx.jy=mesh->ystart-1; bx.jy>=0; bx.jy--) {
	bx.jx=it.ind;
	calc_index(&bx);
	var.setXStencil(s, bx, loc);
	result(bx.jx,bx.jy) = func(s);
      }
    #ifdef CHECK
      result.bndry_ydown = true;
    #endif
  }
  if (mesh->freeboundary_yup) {
    for (RangeIterator it=mesh->iterateBndryUpperY(); !it.isDone(); it++)
      for (bx.jy=mesh->yend+1; bx.jy<mesh->LocalNy; bx.jy++) {
	bx.jx=it.ind;
	calc_index(&bx);
	var.setXStencil(s, bx, loc);
	result(bx.jx,bx.jy) = func(s);
      }
    #ifdef CHECK
      result.bndry_yup = true;
    #endif
  }
  if (mesh->freeboundary_xin && mesh->firstX() && !mesh->periodicX) {
    forward_stencil fs;
    Mesh::boundary_derivs_pair funcs_pair;
    bx.jx=mesh->xstart-1;
    for (bx.jy=mesh->ystart; bx.jy<=mesh->yend; bx.jy++) {
      calc_index(&bx);
      var.setXStencil(fs, bx, loc);
      funcs_pair = func_in(fs);
      result(bx.jx,bx.jy) = funcs_pair.inner;
      result(bx.jxm,bx.jy) = funcs_pair.outer;
    }
    #ifdef CHECK
      result.bndry_xin = true;
    #endif
  }
  if (mesh->freeboundary_xout && mesh->lastX() && !mesh->periodicX) {
    backward_stencil bs;
    Mesh::boundary_derivs_pair funcs_pair;
    bx.jx=mesh->xend+1;
    for (bx.jy=mesh->ystart; bx.jy<=mesh->yend; bx.jy++) {
      calc_index(&bx);
      var.setXStencil(bs, bx, loc);
      funcs_pair = func_out(bs);
      result(bx.jx,bx.jy) = funcs_pair.inner;
      result(bx.jxp,bx.jy) = funcs_pair.outer;
    }
    #ifdef CHECK
      result.bndry_xout = true;
    #endif
  }

  return result;
}

const Field3D Mesh::applyXdiff(const Field3D &var, Mesh::deriv_func func, Mesh::inner_boundary_deriv_func func_in, Mesh::outer_boundary_deriv_func func_out, CELL_LOC loc) {
  if (var.getNx() == 1) {
    return 0.;
  }
  
  // Check that the input variable has data
  ASSERT1(var.isAllocated());

  Field3D result;
  result.allocate(); // Make sure data allocated
  
  if (mesh->StaggerGrids && 
      (loc != CELL_DEFAULT) && (loc != var.getLocation())) {
    // Staggered differencing

    CELL_LOC location = var.getLocation();
    
    if (mesh->xstart > 1) {
      // More than one guard cell, so set pp and mm values
      // This allows higher-order methods to be used
      for(const auto &i : result.region(RGN_NOX)) {
        stencil s;
        s.c = var[i];
        s.p = var[i.xp()];
        s.m = var[i.xm()];
        s.pp = var[i.offset(2,0,0)];
        s.mm = var[i.offset(-2,0,0)];
        
        if ((location == CELL_CENTRE) && (loc == CELL_XLOW)) {
          // Producing a stencil centred around a lower X value
          s.pp = s.p;
          s.p  = s.c;
        } else if (location == CELL_XLOW) {
          // Stencil centred around a cell centre
          s.mm = s.m;
          s.m  = s.c;
        }

        result[i] = func(s);
      }
    } else {
      // Only one guard cell, so no pp or mm values
      for(const auto &i : result.region(RGN_NOX)) {
        stencil s;
        s.c = var[i];
        s.p = var[i.xp()];
        s.m = var[i.xm()];
        s.pp = nan("");
        s.mm = nan("");
        
        if ((location == CELL_CENTRE) && (loc == CELL_XLOW)) {
          // Producing a stencil centred around a lower X value
          s.pp = s.p;
          s.p  = s.c;
        } else if (location == CELL_XLOW) {
          // Stencil centred around a cell centre
          s.mm = s.m;
          s.m  = s.c;
        }
        
        result[i] = func(s);
      }
    }
    
  } else {
    // Non-staggered differencing
    
    if (mesh->xstart > 1) {
      // More than one guard cell, so set pp and mm values
      // This allows higher-order methods to be used
      for(const auto &i : result.region(RGN_NOX)) {
        stencil s;
        s.c = var[i];
        s.p = var[i.xp()];
        s.m = var[i.xm()];
        s.pp = var[i.offset(2,0,0)];
        s.mm = var[i.offset(-2,0,0)];
        
        result[i] = func(s);
      }
    } else {
      // Only one guard cell, so no pp or mm values
      for(const auto &i : result.region(RGN_NOX)) {
        stencil s;
        s.c = var[i];
        s.p = var[i.xp()];
        s.m = var[i.xm()];
        s.pp = nan("");
        s.mm = nan("");
        
        result[i] = func(s);
      }
    }
  }
  
  bindex bx;
  stencil s;
  
#ifdef CHECK
  // Mark boundaries as invalid
  result.bndry_xin = result.bndry_xout = result.bndry_yup = result.bndry_ydown = false;
#endif
  
  if (mesh->freeboundary_ydown) {
    for (RangeIterator it=mesh->iterateBndryLowerY(); !it.isDone(); it++)
      for (bx.jy=mesh->ystart-1; bx.jy>=0; bx.jy--)
	for (bx.jz=0; bx.jz<mesh->LocalNz; bx.jz++) {
	  bx.jx=it.ind;
	  calc_index(&bx);
	  var.setXStencil(s, bx, loc);
	  result(bx.jx,bx.jy,bx.jz) = func(s);
	}
    #ifdef CHECK
      result.bndry_ydown = true;
    #endif
  }
  if (mesh->freeboundary_yup) {
    for (RangeIterator it=mesh->iterateBndryUpperY(); !it.isDone(); it++)
      for (bx.jy=mesh->yend+1; bx.jy<mesh->LocalNy; bx.jy++)
	for (bx.jz=0; bx.jz<mesh->LocalNz; bx.jz++) {
	  bx.jx=it.ind;
	  calc_index(&bx);
	  var.setXStencil(s, bx, loc);
	  result(bx.jx,bx.jy,bx.jz) = func(s);
	}
    #ifdef CHECK
      result.bndry_yup = true;
    #endif
  }
  if (mesh->freeboundary_xin && mesh->firstX() && !mesh->periodicX) {
    forward_stencil fs;
    Mesh::boundary_derivs_pair funcs_pair;
    bx.jx=mesh->xstart-1;
    for (bx.jy=mesh->ystart; bx.jy<=mesh->yend; bx.jy++)
      for (bx.jz=0; bx.jz<mesh->LocalNz; bx.jz++) {
	calc_index(&bx);
	var.setXStencil(fs, bx, loc);
	funcs_pair = func_in(fs);
	result(bx.jx,bx.jy,bx.jz) = funcs_pair.inner;
	result(bx.jxm,bx.jy,bx.jz) = funcs_pair.outer;
      }
    #ifdef CHECK
      result.bndry_xin = true;
    #endif
  }
  if (mesh->freeboundary_xout && mesh->lastX() && !mesh->periodicX) {
    backward_stencil bs;
    Mesh::boundary_derivs_pair funcs_pair;
    bx.jx=mesh->xend+1;
    for (bx.jy=mesh->ystart; bx.jy<=mesh->yend; bx.jy++)
      for (bx.jz=0; bx.jz<mesh->LocalNz; bx.jz++) {
	calc_index(&bx);
	var.setXStencil(bs, bx, loc);
	funcs_pair = func_out(bs);
	result(bx.jx,bx.jy,bx.jz) = funcs_pair.inner;
	result(bx.jxp,bx.jy,bx.jz) = funcs_pair.outer;
      }
    #ifdef CHECK
      result.bndry_xout = true;
    #endif
  }
  
  return result;
}

// Y derivative

const Field2D Mesh::applyYdiff(const Field2D &var, Mesh::deriv_func func, Mesh::inner_boundary_deriv_func func_in, Mesh::outer_boundary_deriv_func func_out, CELL_LOC loc) {
  if (var.getNy() == 1) {
    return 0.;
  }

  // Check that the input variable has data
  ASSERT1(var.isAllocated());

  Field2D result;
  result.allocate(); // Make sure data allocated
  
  if (mesh->ystart > 1) {
    // More than one guard cell, so set pp and mm values
    // This allows higher-order methods to be used
    
    for(const auto &i : result.region(RGN_NOBNDRY)) {
      // Set stencils
      stencil s;
      s.c = var[i];
      s.p = var[i.yp()];
      s.m = var[i.ym()];
      s.pp = var[i.offset(0,2,0)];
      s.mm = var[i.offset(0,-2,0)];

      result[i] = func(s);
    }
  } else {
    // Only one guard cell, so no pp or mm values
    for(const auto &i : result.region(RGN_NOBNDRY)) {
      // Set stencils
      stencil s;
      s.c = var[i];
      s.p = var[i.yp()];
      s.m = var[i.ym()];
      s.pp = nan("");
      s.mm = nan("");

      result[i] = func(s);
    }
  }
    
  bindex bx;
  stencil s;
#ifdef CHECK
  // Mark boundaries as invalid
  result.bndry_yup = result.bndry_ydown = false;
#endif

  if (mesh->freeboundary_xin && mesh->firstX() && !mesh->periodicX) {
    for (bx.jx=mesh->xstart-1; bx.jx>=0; bx.jx--)
      for (bx.jy=mesh->ystart; bx.jy<=mesh->ystart; bx.jy++) {
	calc_index(&bx);
	var.setYStencil(s, bx, loc);
	result(bx.jx,bx.jy) = func(s);
      }
    #ifdef CHECK
      result.bndry_xin = true;
    #endif
  }
  if (mesh->freeboundary_xout && mesh->lastX() && !mesh->periodicX) {
    for (bx.jx=mesh->xend+1; bx.jx<mesh->LocalNx; bx.jx++)
      for (bx.jy=mesh->ystart; bx.jy<=mesh->ystart; bx.jy++) {
	calc_index(&bx);
	var.setYStencil(s, bx, loc);
	result(bx.jx,bx.jy) = func(s);
      }
    #ifdef CHECK
      result.bndry_xout = true;
    #endif
  }
  if (mesh->freeboundary_ydown) {
    forward_stencil fs;
    Mesh::boundary_derivs_pair funcs_pair;
    bx.jy=mesh->ystart-1;
    for (RangeIterator it=mesh->iterateBndryLowerY(); !it.isDone(); it++) {
      bx.jx = it.ind;
      calc_index(&bx);
      var.setYStencil(fs, bx, loc);
      funcs_pair = func_in(fs);
      result(bx.jx,bx.jy) = funcs_pair.inner;
      result(bx.jx,bx.jym) = funcs_pair.outer;
    }
    #ifdef CHECK
      result.bndry_ydown = true;
    #endif
  }
  if (mesh->freeboundary_yup) {
    backward_stencil bs;
    Mesh::boundary_derivs_pair funcs_pair;
    bx.jy=mesh->yend+1;
    for (RangeIterator it=mesh->iterateBndryUpperY(); !it.isDone(); it++) {
      bx.jx = it.ind;
      calc_index(&bx);
      var.setYStencil(bs, bx, loc);
      funcs_pair = func_out(bs);
      result(bx.jx,bx.jy) = funcs_pair.inner;
      result(bx.jx,bx.jyp) = funcs_pair.outer;
    }
    #ifdef CHECK
      result.bndry_yup = true;
    #endif
  }

  return result;
}

const Field3D Mesh::applyYdiff(const Field3D &var, Mesh::deriv_func func, Mesh::inner_boundary_deriv_func func_in, Mesh::outer_boundary_deriv_func func_out, CELL_LOC loc) {
  if (var.getNy() == 1){
    return 0.;
  }

  // Check that the input variable has data
  ASSERT1(var.isAllocated());

  Field3D result;
  result.allocate(); // Make sure data allocated
  
  if (var.hasYupYdown() && 
      ( (&var.yup() != &var) || (&var.ydown() != &var))) {
    // Field "var" has distinct yup and ydown fields which
    // will be used to calculate a derivative along 
    // the magnetic field
    
    if (mesh->StaggerGrids && (loc != CELL_DEFAULT) && (loc != var.getLocation())) {
      // Staggered differencing

      // Cell location of the input field
      CELL_LOC location = var.getLocation();
      
      for(const auto &i : result.region(RGN_NOBNDRY)) {
        // Set stencils
        stencil s;
        s.c = var[i];
        s.p = var.yup()[i.yp()];
        s.m = var.ydown()[i.ym()];
        s.pp = nan("");
        s.mm = nan("");
        
        if ((location == CELL_CENTRE) && (loc == CELL_YLOW)) {
          // Producing a stencil centred around a lower Y value
          s.pp = s.p;
          s.p  = s.c;
        } else if(location == CELL_YLOW) {
          // Stencil centred around a cell centre
          s.mm = s.m;
          s.m  = s.c;
        }

        result[i] = func(s);
      }
    } else {
      // Non-staggered
      for(const auto &i : result.region(RGN_NOBNDRY)) {
        // Set stencils
        stencil s;
        s.c = var[i];
        s.p = var.yup()[i.yp()];
        s.m = var.ydown()[i.ym()];
        s.pp = nan("");
        s.mm = nan("");
        
        result[i] = func(s);
      }
    }
  }else {
    // var has no yup/ydown fields, so we need to shift into field-aligned coordinates
    
    Field3D var_fa = mesh->toFieldAligned(var);
    
    if (mesh->StaggerGrids && (loc != CELL_DEFAULT) && (loc != var.getLocation())) {
      // Staggered differencing
      
      // Cell location of the input field
      CELL_LOC location = var.getLocation();
      
      if (mesh->ystart > 1) {
        // More than one guard cell, so set pp and mm values
        // This allows higher-order methods to be used
        for(const auto &i : result.region(RGN_NOBNDRY)) {
          // Set stencils
          stencil s;
          s.c = var_fa[i];
          s.p = var_fa[i.yp()];
          s.m = var_fa[i.ym()];
          s.pp = var_fa[i.offset(0,2,0)];
          s.mm = var_fa[i.offset(0,-2,0)];
          
          if ((location == CELL_CENTRE) && (loc == CELL_YLOW)) {
            // Producing a stencil centred around a lower Y value
            s.pp = s.p;
            s.p  = s.c;
          } else if(location == CELL_YLOW) {
            // Stencil centred around a cell centre
            s.mm = s.m;
            s.m  = s.c;
          }
          
          result[i] = func(s);
        }
      } else {
        // Only one guard cell, so no pp or mm values
        for(const auto &i : result.region(RGN_NOBNDRY)) {
          // Set stencils
          stencil s;
          s.c = var_fa[i];
          s.p = var_fa[i.yp()];
          s.m = var_fa[i.ym()];
          s.pp = nan("");
          s.mm = nan("");
          
          if ((location == CELL_CENTRE) && (loc == CELL_YLOW)) {
            // Producing a stencil centred around a lower Y value
            s.pp = s.p;
            s.p  = s.c;
          } else if(location == CELL_YLOW) {
            // Stencil centred around a cell centre
            s.mm = s.m;
            s.m  = s.c;
          }
          
          result[i] = func(s);
        }
      }
      
    } else {
      // Non-staggered differencing
      
      if (mesh->ystart > 1) {
        // More than one guard cell, so set pp and mm values
        // This allows higher-order methods to be used
        for(const auto &i : result.region(RGN_NOBNDRY)) {
          // Set stencils
          stencil s;
          s.c = var_fa[i];
          s.p = var_fa[i.yp()];
          s.m = var_fa[i.ym()];
          s.pp = var_fa[i.offset(0,2,0)];
          s.mm = var_fa[i.offset(0,-2,0)];
          
          result[i] = func(s);
        }
      } else {
        // Only one guard cell, so no pp or mm values
        for(const auto &i : result.region(RGN_NOBNDRY)) {
          // Set stencils
          stencil s;
          s.c = var_fa[i];
          s.p = var_fa[i.yp()];
          s.m = var_fa[i.ym()];
          s.pp = nan("");
          s.mm = nan("");
          
          result[i] = func(s);
        }
      }
    }
    
    // Shift result back
    
    result = mesh->fromFieldAligned(result);
  }
#ifdef CHECK
  // Mark boundaries as invalid
  result.bndry_xin = result.bndry_xout = result.bndry_yup = result.bndry_ydown = false;
#endif
  
  bindex bx;
  if (mesh->freeboundary_xin && mesh->firstX() && !mesh->periodicX) {
    for (bx.jx=mesh->xstart-1; bx.jx>=0; bx.jx--)
      for (bx.jy=mesh->ystart; bx.jy<=mesh->ystart; bx.jy++)
	for (bx.jz=0; bx.jz<mesh->LocalNz; bx.jz++) {
          stencil s;
	  calc_index(&bx);
	  var.setYStencil(s, bx, loc);
	  result(bx.jx,bx.jy,bx.jz) = func(s);
	}
    #ifdef CHECK
      result.bndry_xin = true;
    #endif
  }
  if (mesh->freeboundary_xout && mesh->lastX() && !mesh->periodicX) {
    for (bx.jx=mesh->xend+1; bx.jx<mesh->LocalNx; bx.jx++)
      for (bx.jy=mesh->ystart; bx.jy<=mesh->ystart; bx.jy++)
	for (bx.jz=0; bx.jz<mesh->LocalNz; bx.jz++) {
          stencil s;
	  calc_index(&bx);
	  var.setYStencil(s, bx, loc);
	  result(bx.jx,bx.jy,bx.jz) = func(s);
	}
    #ifdef CHECK
      result.bndry_xout = true;
    #endif
  }
  if (mesh->freeboundary_ydown) {
    forward_stencil fs;
    Mesh::boundary_derivs_pair funcs_pair;
    bx.jy=mesh->ystart-1;
    for (RangeIterator it=mesh->iterateBndryLowerY(); !it.isDone(); it++)
      for (bx.jz=0; bx.jz<mesh->LocalNz; bx.jz++) {
	bx.jx = it.ind;
        stencil s;
	calc_index(&bx);
	var.setYStencil(fs, bx, loc);
	funcs_pair = func_in(fs);
	result(bx.jx,bx.jy,bx.jz) = funcs_pair.inner;
	result(bx.jx,bx.jym,bx.jz) = funcs_pair.outer;
      }
    #ifdef CHECK
      result.bndry_ydown = true;
    #endif
  }
  if (mesh->freeboundary_yup) {
    backward_stencil bs;
    Mesh::boundary_derivs_pair funcs_pair;
    bx.jy=mesh->yend+1;
    for (RangeIterator it=mesh->iterateBndryUpperY(); !it.isDone(); it++)
      for (bx.jz=0; bx.jz<mesh->LocalNz; bx.jz++) {
	bx.jx = it.ind;
        stencil s;
	calc_index(&bx);
	var.setYStencil(bs, bx, loc);
	funcs_pair = func_out(bs);
	result(bx.jx,bx.jy,bx.jz) = funcs_pair.inner;
	result(bx.jx,bx.jyp,bx.jz) = funcs_pair.outer;
      }
    #ifdef CHECK
      result.bndry_yup = true;
    #endif
  }

  return result;
}

// Z derivative

const Field3D Mesh::applyZdiff(const Field3D &var, Mesh::deriv_func func, CELL_LOC loc) {
  Field3D result;
  result.allocate(); // Make sure data allocated
  
  // Check that the input variable has data
  ASSERT1(var.isAllocated());
  
  for(const auto &i : result.region(RGN_ALL)) {
    stencil s;
    s.c = var[i];
    s.p = var[i.zp()];
    s.m = var[i.zm()];
    s.pp = var[i.offset(0,0,2)];
    s.mm = var[i.offset(0,0,-2)];
    
    result[i] = func(s);
  }

  bindex bx;
  stencil s;

  if (mesh->freeboundary_xin && mesh->firstX() && !mesh->periodicX) {
    for (bx.jx=mesh->xstart-1; bx.jx>=0; bx.jx--)
      for (bx.jy=mesh->ystart; bx.jy<=mesh->ystart; bx.jy++)
	for (bx.jz=0; bx.jz<mesh->LocalNz; bx.jz++) {
	  calc_index(&bx);
	  var.setZStencil(s, bx, loc);
	  result(bx.jx,bx.jy,bx.jz) = func(s);
	}
    #ifdef CHECK
      result.bndry_xin = true;
    #endif
  }

  if (mesh->freeboundary_xout && mesh->lastX() && !mesh->periodicX) {
    for (bx.jx=mesh->xend+1; bx.jx<mesh->LocalNx; bx.jx++)
      for (bx.jy=mesh->ystart; bx.jy<=mesh->ystart; bx.jy++)
	for (bx.jz=0; bx.jz<mesh->LocalNz; bx.jz++) {
	  calc_index(&bx);
	  var.setZStencil(s, bx, loc);
	  result(bx.jx,bx.jy,bx.jz) = func(s);
	}
    #ifdef CHECK
      result.bndry_xout = true;
    #endif
  }

  if (mesh->freeboundary_ydown) {
    for (RangeIterator it=mesh->iterateBndryLowerY(); !it.isDone(); it++)
      for (bx.jy=mesh->ystart-1; bx.jy>=0; bx.jy--)
	for (bx.jz=0; bx.jz<mesh->LocalNz; bx.jz++) {
	  bx.jx = it.ind;
	  calc_index(&bx);
	  var.setZStencil(s, bx, loc);
	  result(bx.jx,bx.jy,bx.jz) = func(s);
	}
    #ifdef CHECK
      result.bndry_ydown = true;
    #endif
  }

  if (mesh->freeboundary_yup) {
    for (RangeIterator it=mesh->iterateBndryUpperY(); !it.isDone(); it++)
      for (bx.jy=mesh->yend; bx.jy<mesh->LocalNy; bx.jy++)
	for (bx.jz=0; bx.jz<mesh->LocalNz; bx.jz++) {
	  bx.jx = it.ind;
	  calc_index(&bx);
	  var.setZStencil(s, bx, loc);
	  result(bx.jx,bx.jy,bx.jz) = func(s);
	}
    #ifdef CHECK
      result.bndry_yup = true;
    #endif
  }
  
  return result;
}

/*******************************************************************************
 * First central derivatives
 *******************************************************************************/

////////////// X DERIVATIVE /////////////////

const Field3D Mesh::indexDDX(const Field3D &f, CELL_LOC outloc, DIFF_METHOD method) {
  Mesh::deriv_func func = fDDX; // Set to default function
  Mesh::inner_boundary_deriv_func func_in = fDDX_in;
  Mesh::outer_boundary_deriv_func func_out = fDDX_out;
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
      func_in = sfDDX_in;
      func_out = sfDDX_out;
      table = FirstStagDerivTable; // Set table for others
      diffloc = (inloc == CELL_CENTRE) ? CELL_XLOW : CELL_CENTRE;
      
    }else {
      // A more complicated shift. Get a result at cell centre, then shift.
      if(inloc == CELL_XLOW) {
	// Shifting
	
	func = sfDDX; // Set default
	func_in = sfDDX_in;
	func_out = sfDDX_out;
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
    func_in = lookupInnerBoundaryFunc(table, method);
    func_out = lookupOuterBoundaryFunc(table, method);
    if(func == NULL)
      throw BoutException("Cannot use FFT for X derivatives");
  }
  
  result = applyXdiff(f, func, func_in, func_out, diffloc);
  result.setLocation(diffloc); // Set the result location

  result = interp_to(result, outloc); // Interpolate if necessary

  return result;
}

const Field2D Mesh::indexDDX(const Field2D &f) {
  return applyXdiff(f, fDDX, fDDX_in, fDDX_out);
}

////////////// Y DERIVATIVE /////////////////

const Field3D Mesh::indexDDY(const Field3D &f, CELL_LOC outloc, DIFF_METHOD method) {
  Mesh::deriv_func func = fDDY; // Set to default function
  Mesh::inner_boundary_deriv_func func_in = fDDY_in;
  Mesh::outer_boundary_deriv_func func_out = fDDY_out;
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
      func_in = sfDDY_in;
      func_out = sfDDY_out;
      table = FirstStagDerivTable; // Set table for others
      diffloc = (inloc == CELL_CENTRE) ? CELL_YLOW : CELL_CENTRE;
      
    }else {
      // A more complicated shift. Get a result at cell centre, then shift.
      if(inloc == CELL_YLOW) {
	// Shifting
	
	func = sfDDY; // Set default
	func_in = sfDDY_in;
	func_out = sfDDY_out;
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
    func_in = lookupInnerBoundaryFunc(table, method);
    func_out = lookupOuterBoundaryFunc(table, method);
    if(func == NULL)
      throw BoutException("Cannot use FFT for Y derivatives");
  }
  
  result = applyYdiff(f, func, func_in, func_out, diffloc);
  //output.write("SETTING LOC %s -> %s\n", strLocation(diffloc), strLocation(outloc));
  result.setLocation(diffloc); // Set the result location
  
  return interp_to(result, outloc); // Interpolate if necessary
}

const Field2D Mesh::indexDDY(const Field2D &f) {
  return applyYdiff(f, fDDY, fDDY_in, fDDY_out);
}

////////////// Z DERIVATIVE /////////////////

const Field3D Mesh::indexDDZ(const Field3D &f, CELL_LOC outloc, DIFF_METHOD method, bool inc_xbndry) {
  Mesh::deriv_func func = fDDZ; // Set to default function
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

    int ncz = mesh->LocalNz;
    
#ifndef _OPENMP
    static dcomplex *cv = (dcomplex*) NULL;
#else
    static dcomplex *globalcv;
    static int nthreads = 0;
#endif 

    #pragma omp parallel
    {
#ifndef _OPENMP
      // Serial, so can have a single static array
      if(cv == (dcomplex*) NULL)
        cv = new dcomplex[ncz/2 + 1];
#else
      // Parallel, so allocate a separate array for each thread
      
      int th_id = omp_get_thread_num(); // thread ID
        int n_th = omp_get_num_threads();
      if(th_id == 0) {
        if(nthreads < n_th) {
          // Allocate memory in thread zero
          if(nthreads > 0)
            delete[] globalcv;
          globalcv = new dcomplex[n_th*(ncz/2 + 1)];
          nthreads = n_th;
        }
      }
      // Wait for memory to be allocated
      #pragma omp barrier
      
      dcomplex *cv = globalcv + th_id*(ncz/2 + 1); // Separate array for each thread
#endif
      int xs = mesh->xstart;
      int xe = mesh->xend;
      int ys = mesh->ystart;
      int ye = mesh->yend;
      if(inc_xbndry) { // Include x boundary region (for mixed XZ derivatives)
        xs = 0;
        xe = mesh->LocalNx-1;
      }
      if (mesh->freeboundary_xin && mesh->firstX() && !mesh->periodicX)
        xs = 0;
      if (mesh->freeboundary_xout && mesh->lastX() && !mesh->periodicX)
        xe = mesh->LocalNx-1;
      if (mesh->freeboundary_ydown)
        ys = 0;
      if (mesh->freeboundary_yup)
        ye = mesh->LocalNy-1;
      #pragma omp for
      for(int jx=xs;jx<=xe;jx++) {
        for(int jy=ys;jy<=ye;jy++) {
          rfft(f(jx, jy), ncz, cv); // Forward FFT
          
        for(int jz=0;jz<=ncz/2;jz++) {
            BoutReal kwave=jz*2.0*PI/ncz; // wave number is 1/[rad]
            
          BoutReal flt;
          if (jz>0.4*ncz) flt=1e-10; else flt=1.0;
          cv[jz] *= dcomplex(0.0, kwave) * flt;
          if(mesh->StaggerGrids)
            cv[jz] *= exp(Im * (shift * kwave));
        }
          
        irfft(cv, ncz, result(jx,jy)); // Reverse FFT
      }
    }
    }
    // End of parallel section
    
#ifdef CHECK
    // Mark boundaries as invalid
    if (mesh->freeboundary_xin) result.bndry_xin = true;
    else result.bndry_xin = false;
    if (mesh->freeboundary_xout) result.bndry_xout = true;
    else result.bndry_xout = false;
    if (mesh->freeboundary_yup) result.bndry_yup = true;
    else result.bndry_yup = false;
    if (mesh->freeboundary_ydown) result.bndry_ydown = true;
    else result.bndry_ydown = false;
#endif
    
  }else {
    // All other (non-FFT) functions 
    result = applyZdiff(f, func);
  }
  
  result.setLocation(diffloc);

  return interp_to(result, outloc);
}

const Field2D Mesh::indexDDZ(const Field2D &UNUSED(f)) {
  Field2D result;
  result = 0.0;
  return result;
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
const Field3D Mesh::indexD2DX2(const Field3D &f, CELL_LOC outloc, DIFF_METHOD method) {
  Mesh::deriv_func func = fD2DX2; // Set to default function
  Mesh::inner_boundary_deriv_func func_in = fD2DX2_in;
  Mesh::outer_boundary_deriv_func func_out = fD2DX2_out;
  DiffLookup *table = SecondDerivTable;
  
  CELL_LOC inloc = f.getLocation(); // Input location
  CELL_LOC diffloc = inloc; // Location of differential result
  
  Field3D result;
  
  if(StaggerGrids && (outloc == CELL_DEFAULT)) {
    // Take care of CELL_DEFAULT case
    outloc = diffloc; // No shift (i.e. same as no stagger case)
  }

  if(StaggerGrids && (outloc != inloc)) {
    // Shifting to a new location
    
    if(((inloc == CELL_CENTRE) && (outloc == CELL_XLOW)) ||
       ((inloc == CELL_XLOW) && (outloc == CELL_CENTRE))) {
      // Shifting in X. Centre -> Xlow, or Xlow -> Centre
      
      func = sfD2DX2; // Set default
      func_in = sfD2DX2_in;
      func_out = sfD2DX2_out;
      table = SecondStagDerivTable; // Set table for others
      diffloc = (inloc == CELL_CENTRE) ? CELL_XLOW : CELL_CENTRE;
      
    }else {
      // A more complicated shift. Get a result at cell centre, then shift.
      if(inloc == CELL_XLOW) {
	// Shifting
	
	func = sfD2DX2; // Set default
	func_in = sfD2DX2_in;
	func_out = sfD2DX2_out;
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
    func_in = lookupInnerBoundaryFunc(table, method);
    func_out = lookupOuterBoundaryFunc(table, method);
    if(func == NULL)
      throw BoutException("Cannot use FFT for X derivatives");
  }
  
  result = applyXdiff(f, func, func_in, func_out);
  result.setLocation(diffloc);
  
  result = interp_to(result, outloc);

  /*
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
  */
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
const Field2D Mesh::indexD2DX2(const Field2D &f) {
  Field2D result;
  
  result = applyXdiff(f, fD2DX2, fD2DX2_in, fD2DX2_out);
  
  return result;
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
const Field3D Mesh::indexD2DY2(const Field3D &f, CELL_LOC outloc, DIFF_METHOD method) {
  Mesh::deriv_func func = fD2DY2; // Set to default function
  Mesh::inner_boundary_deriv_func func_in = fD2DY2_in;
  Mesh::outer_boundary_deriv_func func_out = fD2DY2_out;
  DiffLookup *table = SecondDerivTable;
  
  CELL_LOC inloc = f.getLocation(); // Input location
  CELL_LOC diffloc = inloc; // Location of differential result
  
  Field3D result;
  
  if(StaggerGrids && (outloc == CELL_DEFAULT)) {
    // Take care of CELL_DEFAULT case
    outloc = diffloc; // No shift (i.e. same as no stagger case)
  }

  if(StaggerGrids && (outloc != inloc)) {
    // Shifting to a new location
    
    if(((inloc == CELL_CENTRE) && (outloc == CELL_YLOW)) ||
       ((inloc == CELL_YLOW) && (outloc == CELL_CENTRE))) {
      // Shifting in Y. Centre -> Ylow, or Ylow -> Centre
      
      func = sfD2DY2; // Set default
      func_in = sfD2DY2_in;
      func_out = sfD2DY2_out;
      table = SecondStagDerivTable; // Set table for others
      diffloc = (inloc == CELL_CENTRE) ? CELL_YLOW : CELL_CENTRE;
      
    }else {
      // A more complicated shift. Get a result at cell centre, then shift.
      if(inloc == CELL_YLOW) {
	// Shifting
	
	func = sfD2DY2; // Set default
	func_in = sfD2DY2_in;
	func_out = sfD2DY2_out;
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
    func_in = lookupInnerBoundaryFunc(table, method);
    func_out = lookupOuterBoundaryFunc(table, method);
    if(func == NULL)
      throw BoutException("Cannot use FFT for Y derivatives");
  }
  
  result = applyYdiff(f, func, func_in, func_out);
  result.setLocation(diffloc);

  return interp_to(result, outloc);
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
const Field2D Mesh::indexD2DY2(const Field2D &f) {
  return applyYdiff(f, fD2DY2, fD2DY2_in, fD2DY2_out);
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
const Field3D Mesh::indexD2DZ2(const Field3D &f, CELL_LOC outloc, DIFF_METHOD method, bool inc_xbndry) {
  Mesh::deriv_func func = fD2DZ2; // Set to default function
  DiffLookup *table = SecondDerivTable;
  
  CELL_LOC inloc = f.getLocation(); // Input location
  CELL_LOC diffloc = inloc; // Location of differential result
  
  Field3D result;

  if(StaggerGrids && (outloc == CELL_DEFAULT)) {
    // Take care of CELL_DEFAULT case
    outloc = diffloc; // No shift (i.e. same as no stagger case)
  }

  if(StaggerGrids && (outloc != inloc)) {
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
    if(StaggerGrids) {
      if((inloc == CELL_CENTRE) && (diffloc == CELL_ZLOW)) {
	// Shifting down - multiply by exp(-0.5*i*k*dz)
	shift = -1.;
      }else if((inloc == CELL_ZLOW) && (diffloc == CELL_CENTRE)) {
	// Shifting up
	shift = 1.;
      }
    }
    
    result.allocate(); // Make sure data allocated

    int ncz = mesh->LocalNz;
    
    static dcomplex *cv = (dcomplex*) NULL;
    
    // Serial, so can have a single static array
    if(cv == (dcomplex*) NULL)
      cv = new dcomplex[ncz/2 + 1];

    int xs = mesh->xstart;
    int xe = mesh->xend;
    int ys = mesh->ystart;
    int ye = mesh->yend;
    if(inc_xbndry) { // Include x boundary region (for mixed XZ derivatives)
      xs = 0;
      xe = mesh->LocalNx-1;
    }
    if (mesh->freeboundary_xin && mesh->firstX() && !mesh->periodicX)
      xs = 0;
    if (mesh->freeboundary_xout && mesh->lastX() && !mesh->periodicX)
      xe = mesh->LocalNx-1;
    if (mesh->freeboundary_ydown)
      ys = 0;
    if (mesh->freeboundary_yup)
      ye = mesh->LocalNy-1;
      
    for(int jx=xs;jx<=xe;jx++) {
      for(int jy=ys;jy<=ye;jy++) {
          
	rfft(f(jx,jy), ncz, cv); // Forward FFT
	
	for(int jz=0;jz<=ncz/2;jz++) {
	  BoutReal kwave=jz*2.0*PI/ncz; // wave number is 1/[rad]

	  cv[jz] *= -SQ(kwave);
	  if(StaggerGrids)
	    cv[jz] *= exp(0.5*Im * (shift * kwave));
	}

	irfft(cv, ncz, result(jx,jy)); // Reverse FFT
      }
    }

#ifdef CHECK
    // Mark boundaries as invalid
    if (mesh->freeboundary_xin) result.bndry_xin = true;
    else result.bndry_xin = false;
    if (mesh->freeboundary_xout) result.bndry_xout = true;
    else result.bndry_xout = false;
    if (mesh->freeboundary_yup) result.bndry_yup = true;
    else result.bndry_yup = false;
    if (mesh->freeboundary_ydown) result.bndry_ydown = true;
    else result.bndry_ydown = false;
#endif

  }
  else {
    // All other (non-FFT) functions
    result = applyZdiff(f, func);
  }

  result.setLocation(diffloc);

  return interp_to(result, outloc);
}

/*******************************************************************************
 * Fourth derivatives
 *******************************************************************************/

BoutReal D4DX4_C2(stencil &f) {
  return (f.pp - 4.*f.p + 6.*f.c - 4.*f.m + f.mm);
}

Mesh::boundary_derivs_pair D4D4_F2(forward_stencil &f) {
  Mesh::boundary_derivs_pair result;
  result.inner = 2.*f.m-9.*f.c+16.*f.p-14.*f.p2+6.*f.p3-f.p4;
  result.outer = 3.*f.m-14.*f.c+26.*f.p-24.*f.p2+11.*f.p3-2.*f.p4;
  return result;
}

Mesh::boundary_derivs_pair D4D4_B2(backward_stencil &f) {
  Mesh::boundary_derivs_pair result;
  result.inner = 2.*f.p-9.*f.c+16.*f.m-14.*f.m2+6.*f.m3-f.m4;
  result.outer = 3.*f.p-14.*f.c+26.*f.m-24.*f.m2+11.*f.m3-2.*f.m4;
  return result;
}

const Field3D Mesh::indexD4DX4(const Field3D &f) {
  return applyXdiff(f, D4DX4_C2, D4D4_F2, D4D4_B2);
}

const Field2D Mesh::indexD4DX4(const Field2D &f) {
  return applyXdiff(f, D4DX4_C2, D4D4_F2, D4D4_B2);
}

const Field3D Mesh::indexD4DY4(const Field3D &f) {
  return applyYdiff(f, D4DX4_C2, D4D4_F2, D4D4_B2);
}

const Field2D Mesh::indexD4DY4(const Field2D &f) {
  return applyYdiff(f, D4DX4_C2, D4D4_F2, D4D4_B2);
}

const Field3D Mesh::indexD4DZ4(const Field3D &f) {
  return applyZdiff(f, D4DX4_C2);
}

/*******************************************************************************
 * Mixed derivatives
 *******************************************************************************/


/*******************************************************************************
 * Advection schemes
 * 
 * Jan 2009  - Re-written to use Set*Stencil routines
 *******************************************************************************/

////////////// X DERIVATIVE /////////////////

/// Special case where both arguments are 2D. Output location ignored for now
const Field2D Mesh::indexVDDX(const Field2D &v, const Field2D &f, CELL_LOC UNUSED(outloc), DIFF_METHOD method) {
  Mesh::upwind_func func = fVDDX;

  if(method != DIFF_DEFAULT) {
    // Lookup function
    func = lookupUpwindFunc(UpwindTable, method);
  }

  Field2D result;
  result.allocate(); // Make sure data allocated

  bindex bx;
  stencil vs, fs;
  start_index(&bx);
  do {
    f.setXStencil(fs, bx);
    
    result(bx.jx,bx.jy) = func(v[{bx.jx,bx.jy,0}], fs);
  }while(next_index2(&bx));

#ifdef CHECK
  // Mark boundaries as invalid
  result.bndry_xin = result.bndry_xout = false;
#endif

  return result;
}

/// General version for 2 or 3-D objects
const Field3D Mesh::indexVDDX(const Field &v, const Field &f, CELL_LOC outloc, DIFF_METHOD method) {
  TRACE("Mesh::indexVDDX(Field, Field)");
  
  Field3D result;
  result.allocate(); // Make sure data allocated
  
  CELL_LOC vloc = v.getLocation();
  CELL_LOC inloc = f.getLocation(); // Input location
  CELL_LOC diffloc = inloc; // Location of differential result
  
  if(StaggerGrids && (outloc == CELL_DEFAULT)) {
    // Take care of CELL_DEFAULT case
    outloc = diffloc; // No shift (i.e. same as no stagger case)
  }
  
  if(StaggerGrids && (vloc != inloc)) {
    // Staggered grids enabled, and velocity at different location to value
    
    Mesh::flux_func func = sfVDDX;
    DiffLookup *table = UpwindTable;
    
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
      
      throw BoutException("Unhandled shift in Mesh::indexVDDX");
    }
    if(method != DIFF_DEFAULT) {
      // Lookup function
      func = lookupFluxFunc(table, method);
    }
    
    bindex bx;
    start_index(&bx);
    stencil vval, fval;
    do {
      v.setXStencil(vval, bx, diffloc);
      f.setXStencil(fval, bx); // Location is always the same as input
      
      result(bx.jx, bx.jy, bx.jz) = func(vval, fval);
    }while(next_index3(&bx));
    
  }else {
    // Not staggered
    Mesh::upwind_func func = fVDDX;
    DiffLookup *table = UpwindTable;
    
    if(method != DIFF_DEFAULT) {
      // Lookup function
      func = lookupUpwindFunc(table, method);
    }
    
    bindex bx;
    start_index(&bx);
    stencil vval, fval;
    do {
      f.setXStencil(fval, bx); // Location is always the same as input
      result(bx.jx, bx.jy, bx.jz) = func(v[{bx.jx, bx.jy, bx.jz}], fval);
    }while(next_index3(&bx));
  }
  
  result.setLocation(inloc);

#ifdef CHECK
  // Mark boundaries as invalid
  result.bndry_xin = result.bndry_xout = result.bndry_yup = result.bndry_ydown = false;
#endif
  
  return interp_to(result, outloc);
}

////////////// Y DERIVATIVE /////////////////

// special case where both are 2D
const Field2D Mesh::indexVDDY(const Field2D &v, const Field2D &f, CELL_LOC outloc, DIFF_METHOD method) {
  TRACE("Mesh::indexVDDY");
  
  Field2D result;
  result.allocate(); // Make sure data allocated

  CELL_LOC vloc = v.getLocation();
  CELL_LOC inloc = f.getLocation(); // Input location
  CELL_LOC diffloc = inloc; // Location of differential result
  
  if(StaggerGrids && (outloc == CELL_DEFAULT)) {
    // Take care of CELL_DEFAULT case
    outloc = diffloc; // No shift (i.e. same as no stagger case)
  }

  if(StaggerGrids && (vloc != inloc)) {
    // Staggered grids enabled, and velocity at different location to value

    Mesh::flux_func func = sfVDDY;
    DiffLookup *table = UpwindTable;
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
      
      // Instead, pretend it's been shifted
      throw BoutException("Unhandled shift in indexVDDY(Field2D, Field2D)");
    }
    
    if(method != DIFF_DEFAULT) {
      // Lookup function
      func = lookupFluxFunc(table, method);
    }
    
    bindex bx;
    stencil fval, vval;
    
    start_index(&bx);
    do {
      f.setYStencil(fval, bx);
      v.setYStencil(vval, bx, diffloc);
      result(bx.jx, bx.jy) = func(vval,fval);
    }while(next_index2(&bx));
  }else {
    Mesh::upwind_func func = fVDDY;
    DiffLookup *table = UpwindTable;
    
    if(method != DIFF_DEFAULT) {
      // Lookup function
      func = lookupUpwindFunc(table, method);
    }
    bindex bx;
    stencil fval, vval;
    start_index(&bx);
    do {
      f.setYStencil(fval, bx);
      result(bx.jx, bx.jy) = func(v[{bx.jx, bx.jy, 0}],fval);
    }while(next_index2(&bx));
    
  }
  result.setLocation(inloc);
    
#ifdef CHECK
  // Mark boundaries as invalid
  result.bndry_xin = result.bndry_xout = result.bndry_yup = result.bndry_ydown = false;
#endif

  return interp_to(result, outloc);
}

// general case
const Field3D Mesh::indexVDDY(const Field &v, const Field &f, CELL_LOC outloc, DIFF_METHOD method) {
  TRACE("Mesh::indexVDDY(Field, Field)");
  
  Field3D result;
  result.allocate(); // Make sure data allocated

  CELL_LOC vloc = v.getLocation();
  CELL_LOC inloc = f.getLocation(); // Input location
  CELL_LOC diffloc = inloc; // Location of differential result
  
  if(StaggerGrids && (outloc == CELL_DEFAULT)) {
    // Take care of CELL_DEFAULT case
    outloc = diffloc; // No shift (i.e. same as no stagger case)
  }

  if(StaggerGrids && (vloc != inloc)) {
    // Staggered grids enabled, and velocity at different location to value
    
    Mesh::flux_func func = sfVDDY;
    DiffLookup *table = UpwindTable;
  
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
      
      throw BoutException("Unhandled shift in VDDY(Field, Field)");
    }

    if(method != DIFF_DEFAULT) {
      // Lookup function
      func = lookupFluxFunc(table, method);
    }
    
    bindex bx;
    start_index(&bx);
    stencil vval, fval;
    do {
      v.setYStencil(vval, bx, diffloc);
      f.setYStencil(fval, bx);
      
      result(bx.jx, bx.jy, bx.jz) = func(vval, fval);
    }while(next_index3(&bx));
    
  }else {
    Mesh::upwind_func func = fVDDY;
    DiffLookup *table = UpwindTable;
    
    if(method != DIFF_DEFAULT) {
      // Lookup function
      func = lookupUpwindFunc(table, method);
    }
  
    bindex bx;
    start_index(&bx);
    stencil vval, fval;
    do {
      f.setYStencil(fval, bx);
      
      result(bx.jx, bx.jy, bx.jz) = func(v[{bx.jx, bx.jy, bx.jz}], fval);
    }while(next_index3(&bx));
  }
  
  result.setLocation(inloc);
    
#ifdef CHECK
  // Mark boundaries as invalid
  result.bndry_xin = result.bndry_xout = result.bndry_yup = result.bndry_ydown = false;
#endif

  return interp_to(result, outloc);
}

////////////// Z DERIVATIVE /////////////////

// general case
const Field3D Mesh::indexVDDZ(const Field &v, const Field &f, CELL_LOC outloc, DIFF_METHOD method) {
  TRACE("Mesh::indexVDDZ");
  
  Field3D result;
  result.allocate(); // Make sure data allocated

  CELL_LOC vloc = v.getLocation();
  CELL_LOC inloc = f.getLocation(); // Input location
  CELL_LOC diffloc = inloc; // Location of differential result
  
  if(StaggerGrids && (outloc == CELL_DEFAULT)) {
    // Take care of CELL_DEFAULT case
    outloc = diffloc; // No shift (i.e. same as no stagger case)
  }

  if(StaggerGrids && (vloc != inloc)) {
    // Staggered grids enabled, and velocity at different location to value

    Mesh::flux_func func = sfVDDZ;
    DiffLookup *table = UpwindTable;
    
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
      
      throw BoutException("Unhandled shift in indexVDDZ");
    }

    if(method != DIFF_DEFAULT) {
      // Lookup function
      func = lookupFluxFunc(table, method);
    }

    bindex bx;
    start_index(&bx);
    stencil vval, fval;
    do {
      v.setZStencil(vval, bx, diffloc);
      f.setZStencil(fval, bx);
      
      result(bx.jx, bx.jy, bx.jz) = func(vval, fval);
    }while(next_index3(&bx));
    
  }else {
    Mesh::upwind_func func = fVDDZ;
    DiffLookup *table = UpwindTable;

    if(method != DIFF_DEFAULT) {
      // Lookup function
      func = lookupUpwindFunc(table, method);
    }
    
    bindex bx;
    start_index(&bx);
    stencil vval, fval;
    do {
      f.setZStencil(fval, bx);
      result(bx.jx, bx.jy, bx.jz) = func(v[{bx.jx, bx.jy, bx.jz}], fval);
    }while(next_index3(&bx));
  }
  
  result.setLocation(inloc);
  
#ifdef CHECK
  // Mark boundaries as invalid
  result.bndry_xin = result.bndry_xout = result.bndry_yup = result.bndry_ydown = false;
#endif

  return interp_to(result, outloc);
}

/*******************************************************************************
 * Flux conserving schemes
 *******************************************************************************/

const Field2D Mesh::indexFDDX(const Field2D &v, const Field2D &f, CELL_LOC outloc, DIFF_METHOD method) {
  TRACE("Mesh::::indexFDDX(Field2D, Field2D)");
  
  if( (method == DIFF_SPLIT) || ((method == DIFF_DEFAULT) && (fFDDX == NULL)) ) {
    // Split into an upwind and a central differencing part
    // d/dx(v*f) = v*d/dx(f) + f*d/dx(v)
    return indexVDDX(v, f, outloc, DIFF_DEFAULT) + f * indexDDX(v);
  }
  
  Mesh::flux_func func = fFDDX;
  if(method != DIFF_DEFAULT) {
    // Lookup function
    func = lookupFluxFunc(FluxTable, method);
  }
  Field2D result;
  result.allocate(); // Make sure data allocated

  bindex bx;
  stencil vs, fs;
  start_index(&bx);
  do {
    f.setXStencil(fs, bx);
    v.setXStencil(vs, bx);
    
    result(bx.jx, bx.jy) = func(vs, fs);
  }while(next_index2(&bx));

#ifdef CHECK
  // Mark boundaries as invalid
  result.bndry_xin = result.bndry_xout = false;
#endif
  
  return result;
}

const Field3D Mesh::indexFDDX(const Field3D &v, const Field3D &f, CELL_LOC outloc, DIFF_METHOD method) {
  TRACE("Mesh::indexFDDX");
  
  if( (method == DIFF_SPLIT) || ((method == DIFF_DEFAULT) && (fFDDX == NULL)) ) {
    // Split into an upwind and a central differencing part
    // d/dx(v*f) = v*d/dx(f) + f*d/dx(v)
    return indexVDDX(v, f, outloc, DIFF_DEFAULT) + indexDDX(v, outloc, DIFF_DEFAULT) * f;
  }
  
  Mesh::flux_func func = fFDDX;
  DiffLookup *table = FluxTable;

  CELL_LOC vloc = v.getLocation();
  CELL_LOC inloc = f.getLocation(); // Input location
  CELL_LOC diffloc = inloc; // Location of differential result
  
  if(StaggerGrids && (outloc == CELL_DEFAULT)) {
    // Take care of CELL_DEFAULT case
    outloc = diffloc; // No shift (i.e. same as no stagger case)
  }
  
  if(StaggerGrids && (vloc != inloc)) {
    // Staggered grids enabled, and velocity at different location to value
    if(vloc == CELL_XLOW) {
      // V staggered w.r.t. variable
      func = sfFDDX;
      table = FluxStagTable;
      diffloc = CELL_CENTRE;
    }else if((vloc == CELL_CENTRE) && (inloc == CELL_XLOW)) {
      // Shifted
      func = sfFDDX;
      table = FluxStagTable;
      diffloc = CELL_XLOW;
    }else {
      // More complicated. Deciding what to do here isn't straightforward
      // For now, interpolate velocity to the same location as f.

      // Should be able to do something like:
      //return VDDX(interp_to(v, inloc), f, outloc, method);
      
      throw BoutException("Unhandled shift in indexFDDX");
    }
  }

  if(method != DIFF_DEFAULT) {
    // Lookup function
    func = lookupFluxFunc(table, method);
  }
  
  Field3D result;
  result.allocate(); // Make sure data allocated

  bindex bx;
  stencil vval, fval;
  
  start_index(&bx);
  do {
    v.setXStencil(vval, bx, diffloc);
    f.setXStencil(fval, bx); // Location is always the same as input
    
    result(bx.jx,bx.jy,bx.jz) = func(vval, fval);
  }while(next_index3(&bx));
  
  result.setLocation(inloc);

#ifdef CHECK
  // Mark boundaries as invalid
  result.bndry_xin = result.bndry_xout = result.bndry_yup = result.bndry_ydown = false;
#endif
  
  return interp_to(result, outloc);
}

/////////////////////////////////////////////////////////////////////////

const Field2D Mesh::indexFDDY(const Field2D &v, const Field2D &f, CELL_LOC outloc, DIFF_METHOD method) {
  TRACE("Mesh::indexFDDY(Field2D, Field2D)");
  
  if( (method == DIFF_SPLIT) || ((method == DIFF_DEFAULT) && (fFDDY == NULL)) ) {
    // Split into an upwind and a central differencing part
    // d/dx(v*f) = v*d/dx(f) + f*d/dx(v)
    return indexVDDY(v, f, outloc, DIFF_DEFAULT) + f * indexDDY(v);
  }
  
  Mesh::flux_func func = fFDDY;
  if(method != DIFF_DEFAULT) {
    // Lookup function
    func = lookupFluxFunc(FluxTable, method);
  }
  Field2D result;
  result.allocate(); // Make sure data allocated

  bindex bx;
  stencil vs, fs;
  start_index(&bx);
  do {
    f.setYStencil(fs, bx);
    v.setYStencil(vs, bx);
    
    result(bx.jx, bx.jy) = func(vs, fs);
  }while(next_index2(&bx));

#ifdef CHECK
  // Mark boundaries as invalid
  result.bndry_xin = result.bndry_xout = false;
#endif
  
  return result;
}

const Field3D Mesh::indexFDDY(const Field3D &v, const Field3D &f, CELL_LOC outloc, DIFF_METHOD method) {
  TRACE("Mesh::indexFDDY");
  
  if( (method == DIFF_SPLIT) || ((method == DIFF_DEFAULT) && (fFDDY == NULL)) ) {
    // Split into an upwind and a central differencing part
    // d/dx(v*f) = v*d/dx(f) + f*d/dx(v)
    return indexVDDY(v, f, outloc, DIFF_DEFAULT) + indexDDY(v, outloc, DIFF_DEFAULT) * f;
  }
  Mesh::flux_func func = fFDDY;
  DiffLookup *table = FluxTable;

  CELL_LOC vloc = v.getLocation();
  CELL_LOC inloc = f.getLocation(); // Input location
  CELL_LOC diffloc = inloc; // Location of differential result
  
  if(StaggerGrids && (outloc == CELL_DEFAULT)) {
    // Take care of CELL_DEFAULT case
    outloc = diffloc; // No shift (i.e. same as no stagger case)
  }
  
  if(StaggerGrids && (vloc != inloc)) {
    // Staggered grids enabled, and velocity at different location to value    
    if(vloc == CELL_YLOW) {
      // V staggered w.r.t. variable
      func = sfFDDY;
      table = FluxStagTable;
      diffloc = CELL_CENTRE;
    }else if((vloc == CELL_CENTRE) && (inloc == CELL_YLOW)) {
      // Shifted
      func = sfFDDY;
      table = FluxStagTable;
      diffloc = CELL_YLOW;
    }else {
      // More complicated. Deciding what to do here isn't straightforward
      // For now, interpolate velocity to the same location as f.
      
      // Should be able to do something like:
      //return VDDX(interp_to(v, inloc), f, outloc, method);
      
      throw BoutException("Unhandled shift in indexFDDY");
    }
  }

  if(method != DIFF_DEFAULT) {
    // Lookup function
    func = lookupFluxFunc(table, method);
  }
  
  if(func == NULL) {
    // To catch when no function
    return indexVDDY(v, f, outloc, DIFF_DEFAULT) + indexDDY(v, outloc, DIFF_DEFAULT) * f;
  }

  Field3D result;
  result.allocate(); // Make sure data allocated

  bindex bx;
  stencil vval, fval;
  
  start_index(&bx);
  do {
    v.setYStencil(vval, bx, diffloc);
    f.setYStencil(fval, bx); // Location is always the same as input
    
    result(bx.jx, bx.jy, bx.jz) = func(vval, fval);

  }while(next_index3(&bx));

  result.setLocation(inloc);

#ifdef CHECK
  // Mark boundaries as invalid
  result.bndry_xin = result.bndry_xout = result.bndry_yup = result.bndry_ydown = false;
#endif

  return interp_to(result, outloc);
}

/////////////////////////////////////////////////////////////////////////

const Field3D Mesh::indexFDDZ(const Field3D &v, const Field3D &f, CELL_LOC outloc, DIFF_METHOD method) {
  TRACE("Mesh::indexFDDZ(Field3D, Field3D)");
  if( (method == DIFF_SPLIT) || ((method == DIFF_DEFAULT) && (fFDDZ == NULL)) ) {
    // Split into an upwind and a central differencing part
    // d/dx(v*f) = v*d/dx(f) + f*d/dx(v)
    return indexVDDZ(v, f, outloc, DIFF_DEFAULT) + indexDDZ(v, outloc, DIFF_DEFAULT,true) * f;
  }
  
  Mesh::flux_func func = fFDDZ;
  DiffLookup *table = FluxTable;

  CELL_LOC vloc = v.getLocation();
  CELL_LOC inloc = f.getLocation(); // Input location
  CELL_LOC diffloc = inloc; // Location of differential result
  
  if(StaggerGrids && (outloc == CELL_DEFAULT)) {
    // Take care of CELL_DEFAULT case
    outloc = diffloc; // No shift (i.e. same as no stagger case)
  }
  
  if(StaggerGrids && (vloc != inloc)) {
    // Staggered grids enabled, and velocity at different location to value
    if(vloc == CELL_ZLOW) {
      // V staggered w.r.t. variable
      func = sfFDDZ;
      table = FluxStagTable;
      diffloc = CELL_CENTRE;
    }else if((vloc == CELL_CENTRE) && (inloc == CELL_ZLOW)) {
      // Shifted
      func = sfFDDZ;
      table = FluxStagTable;
      diffloc = CELL_ZLOW;
    }else {
      // More complicated. Deciding what to do here isn't straightforward
      // For now, interpolate velocity to the same location as f.

      // Should be able to do something like:
      //return VDDX(interp_to(v, inloc), f, outloc, method);
      
      throw BoutException("Unhandled shift in indexFDDZ");
    }
  }

  if(method != DIFF_DEFAULT) {
    // Lookup function
    func = lookupFluxFunc(table, method);
  }
  
  Field3D result;
  result.allocate(); // Make sure data allocated

  bindex bx;
  stencil vval, fval;
  
  start_index(&bx);
  do {
    v.setZStencil(vval, bx, diffloc);
    f.setZStencil(fval, bx); // Location is always the same as input
    
    result(bx.jx, bx.jy, bx.jz) = func(vval, fval);
  }while(next_index3(&bx));
  
  result.setLocation(inloc);

#ifdef CHECK
  // Mark boundaries as invalid
  result.bndry_xin = result.bndry_xout = result.bndry_yup = result.bndry_ydown = false;
#endif
  
  return interp_to(result, outloc);
}
