#include "gtest/gtest.h"

#include "bout/boutexception.hxx"
//#include "bout/boundary_region_iter.hxx"
#include "bout/yboundary_regions.hxx"

#include "fake_mesh_fixture.hxx"
#include "test_extras.hxx"

using YBTest = FakeMeshFixture;

using bout::globals::mesh;

TEST_F(YBTest, dirichlet_o2) {
  Field3D test = 1.0;
  YBoundary sheath(YBndryType::all, nullptr, *mesh);
  sheath.iter([&](auto& point) { point.dirichlet_o2(test, 2); });
  EXPECT_TRUE(IsFieldEqual(test, 3.0, "RGN_YGUARDS"));
}

// Not yet implemented for FA
// TEST_F(YBTest, dirichlet_o3) {
//   Field3D test = 1.0;
//   YBoundary sheath(YBndryType::all, nullptr, *mesh);
//   sheath.iter([&](auto& point) { point.dirichlet_o3(test, 2); });
//   EXPECT_TRUE(IsFieldEqual(test, 3.0, "RGN_YGUARDS"));
// }

TEST_F(YBTest, interpolate_boundary_o2) {
  Field3D test = 1.0;
  YBoundary sheath(YBndryType::all, nullptr, *mesh);
  sheath.iter(
      [&](auto& point) { EXPECT_DOUBLE_EQ(point.interpolate_boundary_o2(test), 1.0); });
}
