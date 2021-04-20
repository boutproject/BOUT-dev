#include <bout/index_derivs.hxx>

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

