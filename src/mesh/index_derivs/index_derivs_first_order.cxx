#include <bout/index_derivs.hxx>

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

