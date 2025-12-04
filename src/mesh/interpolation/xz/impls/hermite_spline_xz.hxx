#ifndef BOUT_XZHERMITESPLINE_HXX
#define BOUT_XZHERMITESPLINE_HXX

#include "bout/interpolation_xz.hxx"

#include "bout/assert.hxx"
#include "bout/bout_types.hxx"
#include "bout/mask.hxx"
#include "bout/paralleltransform.hxx"
#include "bout/region.hxx"
#include "bout/utils.hxx"
#include <bout/field3d.hxx>

#include <cmath>
#include <string>
#include <vector>

/// Hermite spline interpolation in XZ
///
/// Does not support MPI splitting in X
class XZHermiteSpline : public XZInterpolation {
  Tensor<Ind3D> i_corner; // index of bottom-left grid point
  Tensor<int> k_corner;   // z-index of bottom-left grid point

  // Basis functions for cubic Hermite spline interpolation
  //    see http://en.wikipedia.org/wiki/Cubic_Hermite_spline
  // The h00 and h01 basis functions are applied to the function itself
  // and the h10 and h11 basis functions are applied to its derivative
  // along the interpolation direction.

  Field3D h00_x;
  Field3D h01_x;
  Field3D h10_x;
  Field3D h11_x;
  Field3D h00_z;
  Field3D h01_z;
  Field3D h10_z;
  Field3D h11_z;

protected:
  struct InterpolatePointResult {
    Ind3D ic;
    Ind3D iczp;
    Ind3D icxp;
    Ind3D icxpzp;
    BoutReal result{};
  };

  // Interpolate a single point at index `i`
  //
  // Protected so this can be reused between `XZHermiteSpline` and `XZMonotonicHermiteSpline`
  auto interpolate_point(const Field3D& f, const Field3D& fx, const Field3D& fz,
                         const Field3D& fxz, const Ind3D& i) const
      -> InterpolatePointResult {
    const auto ic = i_corner[i];
    const auto iczp = ic.zp();
    const auto icxp = ic.xp();
    const auto icxpzp = iczp.xp();

    // Interpolate f in X at Z
    const BoutReal f_z = (f[ic] * h00_x[i]) + (f[icxp] * h01_x[i]) + (fx[ic] * h10_x[i])
                         + (fx[icxp] * h11_x[i]);

    // Interpolate f in X at Z+1
    const BoutReal f_zp1 = (f[iczp] * h00_x[i]) + (f[icxpzp] * h01_x[i])
                           + (fx[iczp] * h10_x[i]) + (fx[icxpzp] * h11_x[i]);

    // Interpolate fz in X at Z
    const BoutReal fz_z = (fz[ic] * h00_x[i]) + (fz[icxp] * h01_x[i])
                          + (fxz[ic] * h10_x[i]) + (fxz[icxp] * h11_x[i]);

    // Interpolate fz in X at Z+1
    const BoutReal fz_zp1 = (fz[iczp] * h00_x[i]) + (fz[icxpzp] * h01_x[i])
                            + (fxz[iczp] * h10_x[i]) + (fxz[icxpzp] * h11_x[i]);

    // Interpolate in Z
    const BoutReal result =
        (+f_z * h00_z[i]) + (f_zp1 * h01_z[i]) + (fz_z * h10_z[i]) + (fz_zp1 * h11_z[i]);

    ASSERT2(std::isfinite(result) || i.x() < localmesh->xstart
            || i.x() > localmesh->xend);

    return {ic, iczp, icxp, icxpzp, result};
  }

public:
  XZHermiteSpline(Mesh* mesh = nullptr, [[maybe_unused]] Options* options = nullptr)
      : XZHermiteSpline(0, mesh) {}
  XZHermiteSpline(int y_offset = 0, Mesh* mesh = nullptr);
  XZHermiteSpline(const BoutMask& mask, int y_offset = 0, Mesh* mesh = nullptr)
      : XZHermiteSpline(y_offset, mesh) {
    setRegion(regionFromMask(mask, localmesh));
  }
  ~XZHermiteSpline() override = default;

  void calcWeights(const Field3D& delta_x, const Field3D& delta_z,
                   const std::string& region = "RGN_NOBNDRY") override;
  void calcWeights(const Field3D& delta_x, const Field3D& delta_z, const BoutMask& mask,
                   const std::string& region = "RGN_NOBNDRY") override;

  // Use precalculated weights
  Field3D interpolate(const Field3D& f,
                      const std::string& region = "RGN_NOBNDRY") const override;
  // Calculate weights and interpolate
  Field3D interpolate(const Field3D& f, const Field3D& delta_x, const Field3D& delta_z,
                      const std::string& region = "RGN_NOBNDRY") override;
  Field3D interpolate(const Field3D& f, const Field3D& delta_x, const Field3D& delta_z,
                      const BoutMask& mask,
                      const std::string& region = "RGN_NOBNDRY") override;
  std::vector<ParallelTransform::PositionsAndWeights>
  getWeightsForYApproximation(int i, int j, int k, int yoffset) override;
};

namespace {
const RegisterXZInterpolation<XZHermiteSpline> registerinterphermitespline{
    "hermitespline"};
}

#endif // BOUT_XZHERMITESPLINE_HXX
