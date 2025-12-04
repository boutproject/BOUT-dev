#ifndef BOUT_XZLAGRANGE_HXX
#define BOUT_XZLAGRANGE_HXX

#include "bout/field3d.hxx"
#include "bout/interpolation_xz.hxx"
#include "bout/mask.hxx"
#include "bout/utils.hxx"
#include <string>

/// XZLagrange4pt interpolation class
///
/// Does not support MPI splitting in X
class XZLagrange4pt : public XZInterpolation {
  Tensor<int> i_corner; // x-index of bottom-left grid point
  Tensor<int> k_corner; // z-index of bottom-left grid point

  Field3D t_x, t_z;

public:
  XZLagrange4pt(Mesh* mesh = nullptr, [[maybe_unused]] Options* options = nullptr)
      : XZLagrange4pt(0, mesh) {}
  XZLagrange4pt(int y_offset = 0, Mesh* mesh = nullptr);
  XZLagrange4pt(const BoutMask& mask, int y_offset = 0, Mesh* mesh = nullptr)
      : XZLagrange4pt(y_offset, mesh) {
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
const RegisterXZInterpolation<XZLagrange4pt> registerinterplagrange4pt{"lagrange4pt"};
}

#endif // BOUT_XZLAGRANGE_HXX
