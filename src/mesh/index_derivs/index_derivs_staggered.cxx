#include <bout/index_derivs.hxx>

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

////////////////////////////////////////////////////////////////////////////////
/// Upwind methods
////////////////////////////////////////////////////////////////////////////////
REGISTER_UPWIND_DERIVATIVE_STAGGERED(VDDX_U1_stag, "U1", 1, DERIV::Upwind) {
  // Lower cell boundary
  BoutReal result = (v.m >= 0) ? v.m * f.m : v.m * f.c;

  // Upper cell boundary
  result -= (v.p >= 0) ? v.p * f.c : v.p * f.p;
  result *= -1;

  // result is now d/dx(v*f), but want v*d/dx(f) so subtract f*d/dx(v)
  result -= f.c * (v.p - v.m);
  return result;
}

REGISTER_UPWIND_DERIVATIVE_STAGGERED(VDDX_U2_stag, "U2", 2, DERIV::Upwind) {
  // Calculate d(v*f)/dx = (v*f)[i+1/2] - (v*f)[i-1/2]

  // Upper cell boundary
  BoutReal result =
      (v.p >= 0.) ? v.p * (1.5 * f.c - 0.5 * f.m) : v.p * (1.5 * f.p - 0.5 * f.pp);

  // Lower cell boundary
  result -= (v.m >= 0.) ? v.m * (1.5 * f.m - 0.5 * f.mm) : v.m * (1.5 * f.c - 0.5 * f.p);

  // result is now d/dx(v*f), but want v*d/dx(f) so subtract f*d/dx(v)
  result -= f.c * (v.p - v.m);

  return result;
}

REGISTER_UPWIND_DERIVATIVE_STAGGERED(VDDX_C2_stag, "C2", 1, DERIV::Upwind) {
  // Result is needed at location of f: interpolate v to f's location and take an
  // unstaggered derivative of f
  return 0.5 * (v.p + v.m) * 0.5 * (f.p - f.m);
}

REGISTER_UPWIND_DERIVATIVE_STAGGERED(VDDX_C4_stag, "C4", 2, DERIV::Upwind) {
  // Result is needed at location of f: interpolate v to f's location and take an
  // unstaggered derivative of f
  return (9. * (v.m + v.p) - v.mm - v.pp) / 16. * (8. * f.p - 8. * f.m + f.mm - f.pp)
         / 12.;
}

////////////////////////////////////////////////////////////////////////////////
/// Flux methods
////////////////////////////////////////////////////////////////////////////////
REGISTER_FLUX_DERIVATIVE_STAGGERED(FDDX_U1_stag, "U1", 1, DERIV::Flux) {
  // Lower cell boundary
  BoutReal result = (v.m >= 0) ? v.m * f.m : v.m * f.c;

  // Upper cell boundary
  result -= (v.p >= 0) ? v.p * f.c : v.p * f.p;

  return -result;
}

REGISTER_FLUX_DERIVATIVE_STAGGERED(FDDX_U2_stag, "U2", 2, DERIV::Flux) {
  // Calculate d(v*f)/dx = (v*f)[i+1/2] - (v*f)[i-1/2]

  // Upper cell boundary
  BoutReal result =
      (v.p >= 0.) ? v.p * (1.5 * f.c - 0.5 * f.m) : v.p * (1.5 * f.p - 0.5 * f.pp);

  // Lower cell boundary
  result -= (v.m >= 0.) ? v.m * (1.5 * f.m - 0.5 * f.mm) : v.m * (1.5 * f.c - 0.5 * f.p);

  return result;
}

