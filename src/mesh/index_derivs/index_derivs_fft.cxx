#include <bout/index_derivs.hxx>

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
