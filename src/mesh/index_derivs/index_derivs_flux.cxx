#include <bout/index_derivs.hxx>

////////////////////////////////////////////////////////////////////////////////
/// Flux non-staggered methods
///
/// Basic derivative methods.
/// All expect to have an input grid cell at the same location as the output
/// Hence convert cell centred values -> centred values, or left -> left
///
////////////////////////////////////////////////////////////////////////////////

REGISTER_FLUX_DERIVATIVE(FDDX_U1, "U1", 1, DERIV::Flux) { // No vec

  // Velocity at lower end
  BoutReal vs = 0.5 * (v.m + v.c);
  BoutReal result = (vs >= 0.0) ? vs * f.m : vs * f.c;
  // and at upper
  vs = 0.5 * (v.c + v.p);
  // Existing form doesn't vectorise due to branching
  result -= (vs >= 0.0) ? vs * f.c : vs * f.p;
  return -result;

  // Alternative form would but may involve more operations
  const auto vSplit = vUpDown(vs);
  return result - std::get<0>(vSplit) * f.c + std::get<1>(vSplit) * f.p;
}

REGISTER_FLUX_DERIVATIVE(FDDX_U2, "U2", 2, DERIV::Flux) { // No vec

  // Velocity at lower end
  BoutReal vs = 0.5 * (v.m + v.c);
  BoutReal result = (vs >= 0.0) ? vs * (1.5*f.m - 0.5*f.mm) : vs * (1.5*f.c - 0.5*f.p);
  // and at upper
  vs = 0.5 * (v.c + v.p);
  // Existing form doesn't vectorise due to branching
  result -= (vs >= 0.0) ? vs * (1.5*f.c - 0.5*f.m) : vs * (1.5*f.p - 0.5*f.pp);
  return -result;
}

REGISTER_FLUX_DERIVATIVE(FDDX_C2, "C2", 2, DERIV::Flux) {
  return 0.5 * (v.p * f.p - v.m * f.m);
}

REGISTER_FLUX_DERIVATIVE(FDDX_C4, "C4", 2, DERIV::Flux) {
  return (8. * v.p * f.p - 8. * v.m * f.m + v.mm * f.mm - v.pp * f.pp) / 12.;
}
