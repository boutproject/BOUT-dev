#include "gtest/gtest.h"

#include "bout/boutexception.hxx"
//#include "bout/boundary_region_iter.hxx"
#include "bout/yboundary_regions.hxx"

#include "fake_mesh_fixture.hxx"

using YBTest = FakeMeshFixture;

using bout::globals::mesh;

TEST_F(YBTest, dirichlet_o2) {

  Field3D test = 0.0;

  YBoundary sheath(YBndryType::all, nullptr, *mesh);

  sheath.iter([&](auto& point) { point.dirichlet_o2(test, 1); });
  for (int x = mesh->xstart; x <= mesh->xend; ++x) {
    for (int z = mesh->zstart; z <= mesh->zend; ++z) {
      EXPECT_DOUBLE_EQ(test(x, 0, z), 2.0);
      EXPECT_DOUBLE_EQ(test(x, mesh->yend + 1, z), 2.0);
    }
  }
}
