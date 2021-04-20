#include <bout/index_derivs.hxx>

//////////////////////////////
//--- Fourth order derivatives
//////////////////////////////
REGISTER_STANDARD_DERIVATIVE(D4DX4_C2, "C2", 2, DERIV::StandardFourth) {
  return (f.pp - 4. * f.p + 6. * f.c - 4. * f.m + f.mm);
}

