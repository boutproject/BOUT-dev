#ifndef BOUT_XZHERMITESPLINE_HXX
#define BOUT_XZHERMITESPLINE_HXX

#include "bout/interpolation_xz.hxx"

#include <bout/build_defines.hxx>
#include <bout/bout_types.hxx>
#include <bout/field3d.hxx>

#include <vector>

#if BOUT_HAS_PETSC
#include <bout/petsclib.hxx>
#endif

class XZHermiteSpline : public XZInterpolation {
protected:
  /// This is protected rather than private so that it can be
  /// extended and used by HermiteSplineMonotonic

  Tensor<SpecificInd<IND_TYPE::IND_3D>> i_corner; // index of bottom-left grid point
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

  std::vector<Field3D> newWeights;

#if BOUT_HAS_PETSC
  PetscLib* petsclib;
  bool isInit{false};
  Mat petscWeights;
  Vec rhs, result;
#endif

public:
  XZHermiteSpline(Mesh* mesh = nullptr, [[maybe_unused]] Options* options = nullptr)
      : XZHermiteSpline(0, mesh) {}
  XZHermiteSpline(int y_offset = 0, Mesh* mesh = nullptr);
  XZHermiteSpline(const BoutMask& mask, int y_offset = 0, Mesh* mesh = nullptr)
      : XZHermiteSpline(y_offset, mesh) {
    setRegion(regionFromMask(mask, localmesh));
  }
  ~XZHermiteSpline() {
#if BOUT_HAS_PETSC
    if (isInit) {
      MatDestroy(&petscWeights);
      VecDestroy(&rhs);
      VecDestroy(&result);
      isInit = false;
      delete petsclib;
    }
#endif
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
  std::vector<ParallelTransform::PositionsAndWeights>
  getWeightsForYApproximation(int i, int j, int k, int yoffset) override;
};

namespace {
const RegisterXZInterpolation<XZHermiteSpline> registerinterphermitespline{"hermitespline"};
}

#endif // BOUT_XZHERMITESPLINE_HXX
