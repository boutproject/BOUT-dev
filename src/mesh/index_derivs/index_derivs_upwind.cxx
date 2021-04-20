#include <bout/index_derivs.hxx>

////////////////////////////////////////////////////////////////////////////////
/// Upwind non-staggered methods
///
/// Basic derivative methods.
/// All expect to have an input grid cell at the same location as the output
/// Hence convert cell centred values -> centred values, or left -> left
///
////////////////////////////////////////////////////////////////////////////////

/// Upwinding: Central, 2nd order
REGISTER_UPWIND_DERIVATIVE(VDDX_C2, "C2", 1, DERIV::Upwind) {
  return vc * 0.5 * (f.p - f.m);
}

/// Upwinding: Central, 4th order
REGISTER_UPWIND_DERIVATIVE(VDDX_C4, "C4", 2, DERIV::Upwind) {
  return vc * (8. * f.p - 8. * f.m + f.mm - f.pp) / 12.;
}

/// upwind, 1st order
REGISTER_UPWIND_DERIVATIVE(VDDX_U1, "U1", 1, DERIV::Upwind) { // No vec
  // Existing form doesn't vectorise due to branching
  return vc >= 0.0 ? vc * (f.c - f.m) : vc * (f.p - f.c);
  // Alternative form would but may involve more operations
  const auto vSplit = vUpDown(vc);
  return (std::get<0>(vSplit) * (f.p - f.c) + std::get<1>(vSplit) * (f.c - f.m));
}

/// upwind, 2nd order
REGISTER_UPWIND_DERIVATIVE(VDDX_U2, "U2", 2, DERIV::Upwind) { // No vec
  // Existing form doesn't vectorise due to branching
  return vc >= 0.0 ? vc * (1.5 * f.c - 2.0 * f.m + 0.5 * f.mm)
                   : vc * (-0.5 * f.pp + 2.0 * f.p - 1.5 * f.c);
  // Alternative form would but may involve more operations
  const auto vSplit = vUpDown(vc);
  return (std::get<0>(vSplit) * (1.5 * f.c - 2.0 * f.m + 0.5 * f.mm)
          + std::get<1>(vSplit) * (-0.5 * f.pp + 2.0 * f.p - 1.5 * f.c));
}

/// upwind, 3rd order
REGISTER_UPWIND_DERIVATIVE(VDDX_U3, "U3", 2, DERIV::Upwind) { // No vec
  // Existing form doesn't vectorise due to branching
  return vc >= 0.0 ? vc * (4. * f.p - 12. * f.m + 2. * f.mm + 6. * f.c) / 12.
                   : vc * (-4. * f.m + 12. * f.p - 2. * f.pp - 6. * f.c) / 12.;
  // Alternative form would but may involve more operations
  const auto vSplit = vUpDown(vc);
  return (std::get<0>(vSplit) * (4. * f.p - 12. * f.m + 2. * f.mm + 6. * f.c)
          + std::get<1>(vSplit) * (-4. * f.m + 12. * f.p - 2. * f.pp - 6. * f.c))
         / 12.;
}

/// 3rd-order WENO scheme
REGISTER_UPWIND_DERIVATIVE(VDDX_WENO3, "W3", 2, DERIV::Upwind) { // No vec
  BoutReal deriv, w, r;
  // Existing form doesn't vectorise due to branching

  if (vc > 0.0) {
    // Left-biased stencil

    r = (WENO_SMALL + SQ(f.c - 2.0 * f.m + f.mm))
        / (WENO_SMALL + SQ(f.p - 2.0 * f.c + f.m));

    deriv = (-f.mm + 3. * f.m - 3. * f.c + f.p);

  } else {
    // Right-biased

    r = (WENO_SMALL + SQ(f.pp - 2.0 * f.p + f.c))
        / (WENO_SMALL + SQ(f.p - 2.0 * f.c + f.m));

    deriv = (-f.m + 3. * f.c - 3. * f.p + f.pp);
  }

  w = 1.0 / (1.0 + 2.0 * r * r);
  deriv = 0.5 * ((f.p - f.m) - w * deriv);

  return vc * deriv;
}

///-----------------------------------------------------------------
/// 3rd-order CWENO. Uses the upwinding code and split flux
REGISTER_STANDARD_DERIVATIVE(DDX_CWENO3, "W3", 2, DERIV::Standard) {
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

  const VDDX_WENO3 upwindOp{};
  return upwindOp(0.5, sp) + upwindOp(-0.5, sm);
}

