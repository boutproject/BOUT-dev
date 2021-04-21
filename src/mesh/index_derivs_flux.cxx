// These derivatives split out from index_derivs.cxx in an effort to
// stop certain compilers from choking when parsing that file

#include <bout/index_derivs.hxx>

#include "bout/mesh.hxx"
#include "fft.hxx"
#include "interpolation.hxx"

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

/////////////////////////////////////////////////////////////////////////////////////
/// Here's an example of defining and registering a custom method that doesn't fit
/// into the standard stencil based approach.
// /////////////////////////////////////////////////////////////////////////////////

#if BOUT_HAS_FFTW
class FFTDerivativeType {
public:
  template <DIRECTION direction, STAGGER stagger, int nGuards, typename T>
  void standard(const T& var, T& result, const std::string& region) const {
    AUTO_TRACE();
    ASSERT2(meta.derivType == DERIV::Standard)
    ASSERT2(var.getMesh()->getNguard(direction) >= nGuards);
    ASSERT2(direction == DIRECTION::Z); // Only in Z for now
    ASSERT2(stagger == STAGGER::None);  // Staggering not currently supported
    ASSERT2(bout::utils::is_Field3D<T>::value); // Should never need to call this with Field2D

    auto* theMesh = var.getMesh();

    // Calculate how many Z wavenumbers will be removed
    const int ncz = theMesh->getNpoints(direction);

    int kfilter = static_cast<int>(theMesh->fft_derivs_filter * ncz
                                   / 2); // truncates, rounding down
    if (kfilter < 0)
      kfilter = 0;
    if (kfilter > (ncz / 2))
      kfilter = ncz / 2;
    const int kmax = ncz / 2 - kfilter; // Up to and including this wavenumber index

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
      BOUT_FOR_INNER(i, theMesh->getRegion2D(region)) {
        auto i3D = theMesh->ind2Dto3D(i, 0);
        rfft(&var[i3D], ncz, cv.begin()); // Forward FFT

        for (int jz = 0; jz <= kmax; jz++) {
          const BoutReal kwave = jz * kwaveFac; // wave number is 1/[rad]
          cv[jz] *= dcomplex(0, kwave);
        }
        for (int jz = kmax + 1; jz <= ncz / 2; jz++) {
          cv[jz] = 0.0;
        }

        irfft(cv.begin(), ncz, &result[i3D]); // Reverse FFT
      }
    }
  }

  template <DIRECTION direction, STAGGER stagger, int nGuards, typename T>
  void upwindOrFlux(const T& UNUSED(vel), const T& UNUSED(var), T& UNUSED(result),
                    const std::string& UNUSED(region)) const {
    AUTO_TRACE();
    throw BoutException("The FFT METHOD isn't available in upwind/Flux");
  }
  static constexpr metaData meta{"FFT", 0, DERIV::Standard};
};
constexpr metaData FFTDerivativeType::meta;

class FFT2ndDerivativeType {
public:
  template <DIRECTION direction, STAGGER stagger, int nGuards, typename T>
  void standard(const T& var, T& result, const std::string& region) const {
    AUTO_TRACE();
    ASSERT2(meta.derivType == DERIV::StandardSecond);
    ASSERT2(var.getMesh()->getNguard(direction) >= nGuards);
    ASSERT2(direction == DIRECTION::Z); // Only in Z for now
    ASSERT2(stagger == STAGGER::None);  // Staggering not currently supported
    ASSERT2(bout::utils::is_Field3D<T>::value); // Should never need to call this with Field2D

    auto* theMesh = var.getMesh();

    // Calculate how many Z wavenumbers will be removed
    const int ncz = theMesh->getNpoints(direction);
    const int kmax = ncz / 2;

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
      BOUT_FOR_INNER(i, theMesh->getRegion2D(region)) {
        auto i3D = theMesh->ind2Dto3D(i, 0);
        rfft(&var[i3D], ncz, cv.begin()); // Forward FFT

        for (int jz = 0; jz <= kmax; jz++) {
          const BoutReal kwave = jz * kwaveFac; // wave number is 1/[rad]
          cv[jz] *= -kwave * kwave;
        }
        for (int jz = kmax + 1; jz <= ncz / 2; jz++) {
          cv[jz] = 0.0;
        }

        irfft(cv.begin(), ncz, &result[i3D]); // Reverse FFT
      }
    }
  }

  template <DIRECTION direction, STAGGER stagger, int nGuards, typename T>
  void upwindOrFlux(const T& UNUSED(vel), const T& UNUSED(var), T& UNUSED(result),
                    const std::string& UNUSED(region)) const {
    AUTO_TRACE();
    throw BoutException("The FFT METHOD isn't available in upwind/Flux");
  }
  static constexpr metaData meta{"FFT", 0, DERIV::StandardSecond};
};
constexpr metaData FFT2ndDerivativeType::meta;

produceCombinations<Set<WRAP_ENUM(DIRECTION, Z)>, Set<WRAP_ENUM(STAGGER, None)>,
                    Set<TypeContainer<Field3D>>,
                    Set<FFTDerivativeType, FFT2ndDerivativeType>>
    registerFFTDerivative(registerMethod{});
#endif

class SplitFluxDerivativeType {
public:
  template <DIRECTION direction, STAGGER stagger, int nGuards, typename T>
  void standard(const T&, T&, const std::string) const {
    AUTO_TRACE();
    throw BoutException("The SPLIT method isn't available for standard");
  }

  template <DIRECTION direction, STAGGER stagger, int nGuards, typename T>
  void upwindOrFlux(const T& vel, const T& var, T& result, const std::string region) const {
    AUTO_TRACE();
    // Split into an upwind and a central differencing part
    // d/dx(v*f) = v*d/dx(f) + f*d/dx(v)
    result = bout::derivatives::index::flowDerivative<T, direction, DERIV::Upwind>(
        vel, var, result.getLocation(), "DEFAULT", region);
    result += bout::derivatives::index::standardDerivative<T, direction, DERIV::Standard>(
                  vel, result.getLocation(), "DEFAULT", region)
              * interp_to(var, result.getLocation());
  }
  static constexpr metaData meta{"SPLIT", 2, DERIV::Flux};
};
constexpr metaData SplitFluxDerivativeType::meta;

produceCombinations<Set<WRAP_ENUM(DIRECTION, X), WRAP_ENUM(DIRECTION, Y),
                        WRAP_ENUM(DIRECTION, YOrthogonal), WRAP_ENUM(DIRECTION, Z)>,
                    Set<WRAP_ENUM(STAGGER, None), WRAP_ENUM(STAGGER, C2L),
                        WRAP_ENUM(STAGGER, L2C)>,
                    Set<TypeContainer<Field3D>, TypeContainer<Field2D>>,
                    Set<SplitFluxDerivativeType>>
    registerSplitDerivative(registerMethod{});
