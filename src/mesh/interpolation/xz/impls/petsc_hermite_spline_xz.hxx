#ifndef BOUT_XZPETSCHERMITESPLINE_HXX
#define BOUT_XZPETSCHERMITESPLINE_HXX

#include <bout/interpolation_xz.hxx>
#include <bout/build_defines.hxx>

#if not BOUT_HAS_PETSC
namespace {
  const XZInterpolationFactory::RegisterUnavailableInFactory registerunavailablepetschermitespline("petschermitespline", "BOUT++ was not configured with PETSc");
}
#else

#include <bout/bout_types.hxx>
#include <bout/field3d.hxx>
#include <bout/petsclib.hxx>
#include <bout/region.hxx>

#include <vector>

class XZPetscHermiteSpline : public XZInterpolation {
  Tensor<Ind3D> i_corner; //< index of bottom-left grid point
  Tensor<int> k_corner;   //< z-index of bottom-left grid point

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

  PetscLib petsclib;
  bool isInit{false};
  Mat petscWeights;
  Vec rhs, result;

public:
  XZPetscHermiteSpline(Mesh* mesh = nullptr, [[maybe_unused]] Options* options = nullptr)
      : XZPetscHermiteSpline(0, mesh) {}
  XZPetscHermiteSpline(int y_offset = 0, Mesh* mesh = nullptr);
  XZPetscHermiteSpline(const BoutMask& mask, int y_offset = 0, Mesh* mesh = nullptr)
      : XZPetscHermiteSpline(y_offset, mesh) {
    setRegion(regionFromMask(mask, localmesh));
  }
  ~XZPetscHermiteSpline() override {
    if (isInit) {
      MatDestroy(&petscWeights);
      VecDestroy(&rhs);
      VecDestroy(&result);
      isInit = false;
    }
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
const RegisterXZInterpolation<XZPetscHermiteSpline> registerinterppetschermitespline{"petschermitespline"};
}

#endif

#endif // BOUT_XZPETSCHERMITESPLINE_HXX
