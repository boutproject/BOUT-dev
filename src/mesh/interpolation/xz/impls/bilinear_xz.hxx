#ifndef BOUT_XZBILINEAR_HXX
#define BOUT_XZBILINEAR_HXX

#include "bout/field3d.hxx"
#include "bout/interpolation_xz.hxx"
#include "bout/mask.hxx"
#include "bout/utils.hxx"

#include <string>

/// XZBilinear interpolation class
///
/// Does not support MPI splitting in X.
class XZBilinear : public XZInterpolation {
  Tensor<int> i_corner; // x-index of bottom-left grid point
  Tensor<int> k_corner; // z-index of bottom-left grid point

  Field3D w0, w1, w2, w3;

public:
  XZBilinear(Mesh* mesh = nullptr, [[maybe_unused]] Options* options = nullptr)
      : XZBilinear(0, mesh) {}
  XZBilinear(int y_offset = 0, Mesh* mesh = nullptr);
  XZBilinear(const BoutMask& mask, int y_offset = 0, Mesh* mesh = nullptr)
      : XZBilinear(y_offset, mesh) {
    setRegion(regionFromMask(mask, localmesh));
  }

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
};

namespace {
const RegisterXZInterpolation<XZBilinear> registerinterpbilinear{"bilinear"};
} // namespace

#endif // BOUT_XZBILINEAR_HXX
