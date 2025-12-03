#ifndef BOUT_XZMONOTONICHERMITESPLINE_HXX
#define BOUT_XZMONOTONICHERMITESPLINE_HXX

#include "hermite_spline_xz.hxx"

#include <bout/bout_types.hxx>

/// Monotonic Hermite spline interpolator
///
/// Similar to XZHermiteSpline, so uses most of the same code.
/// Forces the interpolated result to be in the range of the
/// neighbouring cell values. This prevents unphysical overshoots,
/// but also degrades accuracy near maxima and minima.
/// Perhaps should only impose near boundaries, since that is where
/// problems most obviously occur.
///
/// You can control how tight the clipping to the range of the neighbouring cell
/// values through ``rtol`` and ``atol``:
///
///     diff = (max_of_neighours - min_of_neighours) * rtol + atol
///
/// and the interpolated value is instead clipped to the range
/// ``[min_of_neighours - diff, max_of_neighours + diff]``
class XZMonotonicHermiteSpline : public XZHermiteSpline {
  /// Absolute tolerance for clipping
  BoutReal atol = 0.0;
  /// Relative tolerance for clipping
  BoutReal rtol = 1.0;

public:
  XZMonotonicHermiteSpline(Mesh* mesh = nullptr, Options* options = nullptr)
      : XZHermiteSpline(0, mesh),
        atol{(*options)["atol"]
                 .doc("Absolute tolerance for clipping overshoot")
                 .withDefault(0.0)},
        rtol{(*options)["rtol"]
                 .doc("Relative tolerance for clipping overshoot")
                 .withDefault(1.0)} {
    if (localmesh->getNXPE() > 1) {
      throw BoutException("Do not support MPI splitting in X");
    }
  }
  XZMonotonicHermiteSpline(int y_offset = 0, Mesh* mesh = nullptr)
      : XZHermiteSpline(y_offset, mesh) {
    if (localmesh->getNXPE() > 1) {
      throw BoutException("Do not support MPI splitting in X");
    }
  }
  XZMonotonicHermiteSpline(const BoutMask& mask, int y_offset = 0, Mesh* mesh = nullptr)
      : XZHermiteSpline(mask, y_offset, mesh) {
    if (localmesh->getNXPE() > 1) {
      throw BoutException("Do not support MPI splitting in X");
    }
  }

  using XZHermiteSpline::interpolate;
  /// Interpolate using precalculated weights.
  /// This function is called by the other interpolate functions
  /// in the base class XZHermiteSpline.
  Field3D interpolate(const Field3D& f,
                      const std::string& region = "RGN_NOBNDRY") const override;
};

namespace {
const RegisterXZInterpolation<XZMonotonicHermiteSpline>
    registerinterpmonotonichermitespline{"monotonichermitespline"};
}

#endif // BOUT_XZMONOTONICHERMITESPLINE_HXX
