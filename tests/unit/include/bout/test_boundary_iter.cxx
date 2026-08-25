#include "gtest/gtest.h"

#include "bout/boundary_region.hxx"
#include "bout/boundary_region_iter.hxx"

#include "fake_mesh_fixture.hxx"

using BoundaryIteratorTest = FakeMeshFixture;

using namespace bout::boundary;

TEST_F(BoundaryIteratorTest, Contains) {
  BoundaryRegionFCI bndry{"test", BndryLoc::yup, +1, bout::globals::mesh};

  bndry.add_point(1, 0, 3, 0.0, 0.0, 0.0, 0.0, 0, 0);
  bndry.add_point(1, 1, 3, 0.0, 0.0, 0.0, 0.0, 0, 0);
  bndry.add_point(1, 2, 3, 0.0, 0.0, 0.0, 0.0, 0, 0);

  ASSERT_TRUE(bndry.contains(1, 2, 3));
  ASSERT_FALSE(bndry.contains(1, 3, 3));
  ASSERT_FALSE(bndry.contains(3, 2, 1));
}

TEST_F(BoundaryIteratorTest, SetValid) {
  BoundaryRegionFCI bndry{"test", BndryLoc::yup, +1, bout::globals::mesh};

  bndry.add_point(1, 2, 3, 0.0, 0.0, 0.0, 0.0, 0, 0);

  for (auto& point : bndry) {
    point.setValid(6);
  }

  const auto first = *bndry.begin();
  ASSERT_EQ(first.valid(), 6);
}
