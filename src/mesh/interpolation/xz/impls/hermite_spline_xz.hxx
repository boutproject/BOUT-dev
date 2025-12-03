#ifndef BOUT_XZHERMITESPLINE_HXX
#define BOUT_XZHERMITESPLINE_HXX

#include "bout/interpolation_xz.hxx"

#include "bout/mask.hxx"
#include "bout/paralleltransform.hxx"
#include "bout/region.hxx"
#include "bout/utils.hxx"
#include <bout/field3d.hxx>

#include <string>
#include <vector>

/// Hermite spline interpolation in XZ
///
/// Does not support MPI splitting in X
class XZHermiteSpline : public XZInterpolation {
protected:
  /// This is protected rather than private so that it can be
  /// extended and used by HermiteSplineMonotonic

  Tensor<Ind3D> i_corner; // index of bottom-left grid point
  Tensor<int> k_corner;                           // z-index of bottom-left grid point

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
