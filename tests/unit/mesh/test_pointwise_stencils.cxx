#include "gtest/gtest.h"

#include "bout/constants.hxx"
#include "bout/difops.hxx"
#include "bout/fv_ops.hxx"
#include "bout/globals.hxx"
#include "bout/smoothing.hxx"

#include "fake_mesh_fixture.hxx"

#include <cmath>

namespace {
template <int NZ, int MYG>
struct StencilGrid {
  static constexpr int nz = NZ;
  static constexpr int myg = MYG;
};

template <class Grid>
class PointwiseStencilsTest
    : public FakeMeshFixture_tmpl<9, 9, Grid::nz, false, Grid::myg> {};

using StencilGrids =
    ::testing::Types<StencilGrid<2, 1>, StencilGrid<7, 1>, StencilGrid<7, 2>>;
TYPED_TEST_SUITE(PointwiseStencilsTest, StencilGrids);

TYPED_TEST(PointwiseStencilsTest, SmoothingQuadraticIncludingGuards) {
  Field3D f{0.0};
  for (const auto& i : f.getRegion("RGN_ALL")) {
    f[i] = i.x() * i.x() + i.y() * i.y() + 4.0 * i.z();
  }
  f.splitParallelSlices();
  f.yup() = f + 4.0;
  f.ydown() = f + 8.0;

  const auto sx = smooth_x(f);
  const auto sy = smooth_y(f);
  const auto sxy = smoothXY(f);
  for (const auto& i : f.getRegion("RGN_ALL")) {
    // A 1-2-1 filter adds half the quadratic coefficient in its direction.
    const bool x_edge = i.x() == 0 || i.x() == 8;
    const bool y_edge = i.y() == 0 || i.y() == 8;
    EXPECT_DOUBLE_EQ(sx[i], f[i] + (x_edge ? 0.0 : 0.5));
    // The two distinct parallel slices add a further (4 + 8)/4.
    EXPECT_DOUBLE_EQ(sy[i], f[i] + (y_edge ? 0.0 : 3.5));
    if (i.x() >= 2 && i.x() <= 6 && i.y() >= 2 && i.y() <= 6) {
      EXPECT_DOUBLE_EQ(sxy[i], f[i] + 0.75);
    }
  }
}

TYPED_TEST(PointwiseStencilsTest, AveragingBroadcastsIncludeGuards) {
  Field2D f2{0.0};
  for (const auto& i : f2.getRegion("RGN_ALL")) {
    f2[i] = i.x() + i.y();
  }
  Field3D f{f2};
  for (const auto& i : f.getRegion("RGN_ALL")) {
    f[i] += 4.0 * i.z();
  }
  const auto ax2 = averageX(f2);
  const auto ay2 = averageY(f2);
  const auto ax = averageX(f);
  const auto ay = averageY(f);
  // Both interior lines are symmetric about index 4, for either guard width.
  for (const auto& i : f2.getRegion("RGN_ALL")) {
    EXPECT_DOUBLE_EQ(ax2[i], 4.0 + i.y());
    EXPECT_DOUBLE_EQ(ay2[i], i.x() + 4.0);
  }
  for (const auto& i : f.getRegion("RGN_ALL")) {
    EXPECT_DOUBLE_EQ(ax[i], 4.0 + i.y() + 4.0 * i.z());
    EXPECT_DOUBLE_EQ(ay[i], i.x() + 4.0 + 4.0 * i.z());
  }
}

} // namespace
