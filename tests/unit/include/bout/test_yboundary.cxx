#include "gtest/gtest.h"

#include "bout/boundary_common.hxx"
#include "bout/boundary_region.hxx"
#include "bout/boundary_region_iter.hxx"
#include "bout/bout_types.hxx"
#include "bout/field3d.hxx"
#include "bout/region.hxx"
#include "bout/utils.hxx"
#include "bout/yboundary_regions.hxx"

#include "fake_mesh_fixture.hxx"
#include "test_extras.hxx"

using YBTest = FakeMeshFixture_tmpl<4, 5, 7>;

using bout::boundary::BoundaryIterator;
using bout::boundary::YBoundary;
using bout::globals::mesh;

TEST_F(YBTest, dirichlet_o2_rgn) {
  dynamic_cast<FakeMesh*>(mesh)->createBoundaries();
  Field3D test = 1.0;
  YBoundary sheath(YBndryType::all, nullptr, *mesh);
  sheath.iter([&](const BoundaryIterator auto& point) { dirichlet_o2(point, test, 2); });
  EXPECT_TRUE(IsFieldEqual(test, 3.0, "RGN_YGUARDS"));
}

TEST_F(YBTest, neumann_o1_rgn) {
  dynamic_cast<FakeMesh*>(mesh)->createBoundaries();
  Field3D test = 1.0;
  YBoundary sheath(YBndryType::all, nullptr, *mesh);
  sheath.iter([&](const BoundaryIterator auto& point) { neumann_o1(point, test, 1); });
  EXPECT_TRUE(IsFieldEqual(test, 2.0, "RGN_YGUARDS"));
}

TEST_F(YBTest, neumann_o2_rgn) {
  dynamic_cast<FakeMesh*>(mesh)->createBoundaries();
  Field3D test = 1.0;
  YBoundary sheath(YBndryType::all, nullptr, *mesh);
  sheath.iter([&](const BoundaryIterator auto& point) { neumann_o2(point, test, 1); });
  EXPECT_TRUE(IsFieldEqual(test, 3.0, "RGN_YGUARDS"));
}

TEST_F(YBTest, neumann_o3_rgn) {
  dynamic_cast<FakeMesh*>(mesh)->createBoundaries();
  Field3D test = 1.0;
  YBoundary sheath(YBndryType::all, nullptr, *mesh);
  sheath.iter([&](const BoundaryIterator auto& point) { neumann_o3(point, test, 1); });
  EXPECT_TRUE(IsFieldEqual(test, 0.0, "RGN_YGUARDS"));
}

TEST_F(YBTest, bndry_size) {
  dynamic_cast<FakeMesh*>(mesh)->createBoundaries();
  YBoundary sheath(YBndryType::all, nullptr, *mesh);
  int sum = 0;
  sheath.iter([&]([[maybe_unused]] const BoundaryIterator auto& point) { sum++; });
  EXPECT_EQ(sum, 28);
}

TEST_F(YBTest, interpolate_boundary_o2) {
  Field3D test = 1.0;
  YBoundary sheath(YBndryType::all, nullptr, *mesh);
  sheath.iter([&](const BoundaryIterator auto& point) {
    EXPECT_DOUBLE_EQ(interpolate_boundary_o2(point, test), 1.0);
  });
}

TEST_F(YBTest, interpolate_boundary_o2_const) {
  const Field3D test = 1.0;
  YBoundary sheath(YBndryType::all, nullptr, *mesh);
  sheath.iter([&](const BoundaryIterator auto& point) {
    EXPECT_DOUBLE_EQ(interpolate_boundary_o2(point, test), 1.0);
  });
}

TEST_F(YBTest, extrapolate_boundary_o2) {
  Field3D test = 1.0;
  YBoundary sheath(YBndryType::all, nullptr, *mesh);
  sheath.iter([&](const BoundaryIterator auto& point) {
    EXPECT_DOUBLE_EQ(extrapolate_boundary_o2(point, test), 1.0);
  });
}

TEST_F(YBTest, dirichlet_o1) {
  Field3D test = 1.0;
  YBoundary sheath(YBndryType::all, nullptr, *mesh);
  sheath.iter([&](const BoundaryIterator auto& point) {
    dirichlet_o1(point, test, 2.0);
    EXPECT_DOUBLE_EQ(interpolate_boundary_o2(point, test), 1.5);
  });
}

TEST_F(YBTest, dirichlet_o2) {
  Field3D test = 1.0;
  YBoundary sheath(YBndryType::all, nullptr, *mesh);
  sheath.iter([&](const BoundaryIterator auto& point) {
    dirichlet_o2(point, test, 2.0);
    EXPECT_DOUBLE_EQ(interpolate_boundary_o2(point, test), 2.0);
  });
}

TEST_F(YBTest, dirichlet_o3) {
  Field3D test = 1.0;
  YBoundary sheath(YBndryType::all, nullptr, *mesh);
  sheath.iter([&](const BoundaryIterator auto& point) {
    dirichlet_o3(point, test, 2.0);
    EXPECT_DOUBLE_EQ(interpolate_boundary_o2(point, test), 7. / 3.);
  });
}

TEST_F(YBTest, extrapolate_boundary_free) {
  Field3D test = makeField<Field3D>([&](auto& i) { return SQ(i.y() - 2) + 1; }, mesh);
  YBoundary sheath(YBndryType::all, nullptr, *mesh);
  using namespace bout::boundary;
  sheath.iter([&](const BoundaryIterator auto& point) {
    EXPECT_DOUBLE_EQ(
        extrapolate_boundary_free(point, test, BoundaryFreeExtrapolation::limited), 3);
    EXPECT_DOUBLE_EQ(
        extrapolate_boundary_free(point, test, BoundaryFreeExtrapolation::linear), 3.5);
    EXPECT_DOUBLE_EQ(
        extrapolate_boundary_free(point, test, BoundaryFreeExtrapolation::exponential),
        5);
  });
}

TEST_F(YBTest, set_free) {
  Field3D test = makeField<Field3D>([&](auto& i) { return SQ(i.y() - 2) + 1; }, mesh);
  YBoundary sheath(YBndryType::all, nullptr, *mesh);
  sheath.iter([&](const BoundaryIterator auto& point) {
    set_free(point, test, bout::boundary::BoundaryFreeExtrapolation::limited);
    EXPECT_DOUBLE_EQ(interpolate_boundary_o2(point, test), 3);
  });
}

TEST_F(YBTest, interpolate_boundary_o2_square) {
  Field3D test = makeField<Field3D>([&](auto& i) { return SQ(i.y() - 2); }, mesh);
  YBoundary sheath(YBndryType::all, nullptr, *mesh);
  sheath.iter([&](const BoundaryIterator auto& point) {
    EXPECT_DOUBLE_EQ(interpolate_boundary_o2(point, test), 2.5);
  });
}

TEST_F(YBTest, limit_at_least) {
  Field3D test = 1.0;
  YBoundary sheath(YBndryType::all, nullptr, *mesh);
  sheath.iter([&](const BoundaryIterator auto& point) {
    limit_at_least(point, test, 2.0);
    EXPECT_DOUBLE_EQ(interpolate_boundary_o2(point, test), 1.5);
  });
}

TEST_F(YBTest, extrapolate_grad_o2_square) {
  Field3D test = makeField<Field3D>([&](auto& i) { return SQ(i.y() - 2); }, mesh);
  YBoundary sheath(YBndryType::all, nullptr, *mesh);
  sheath.iter([&](const BoundaryIterator auto& point) {
    EXPECT_DOUBLE_EQ(extrapolate_grad_o2(point, test), 1);
  });
}

TEST_F(YBTest, extrapolate_next_o1_square) {
  Field3D test = makeField<Field3D>([&](auto& i) { return SQ(i.y() - 2); }, mesh);
  YBoundary sheath(YBndryType::all, nullptr, *mesh);
  sheath.iter([&](const BoundaryIterator auto& point) {
    EXPECT_DOUBLE_EQ(extrapolate_next_o1(point, test), 1);
  });
}

TEST_F(YBTest, extrapolate_next_o2_square) {
  Field3D test = makeField<Field3D>([&](auto& i) { return SQ(i.y() - 2); }, mesh);
  YBoundary sheath(YBndryType::all, nullptr, *mesh);
  sheath.iter([&](const BoundaryIterator auto& point) {
    EXPECT_DOUBLE_EQ(extrapolate_next_o2(point, test), 2);
  });
}

TEST_F(YBTest, next_square) {
  Field3D test = makeField<Field3D>([&](auto& i) { return SQ(i.y() - 2); }, mesh);
  YBoundary sheath(YBndryType::all, nullptr, *mesh);
  sheath.iter(
      [&](const BoundaryIterator auto& point) { EXPECT_DOUBLE_EQ(point.next(test), 4); });
}

TEST_F(YBTest, current_square) {
  Field3D test = makeField<Field3D>([&](auto& i) { return SQ(i.y() - 2); }, mesh);
  YBoundary sheath(YBndryType::all, nullptr, *mesh);
  sheath.iter([&](const BoundaryIterator auto& point) {
    EXPECT_DOUBLE_EQ(point.current(test), 1);
  });
}

TEST_F(YBTest, prev_square) {
  Field3D test = makeField<Field3D>([&](auto& i) { return SQ(i.y() - 2); }, mesh);
  YBoundary sheath(YBndryType::all, nullptr, *mesh);
  sheath.iter(
      [&](const BoundaryIterator auto& point) { EXPECT_DOUBLE_EQ(point.prev(test), 0); });
}

TEST_F(YBTest, next_const) {
  const Field3D test = makeField<Field3D>([&](auto& i) { return SQ(i.y() - 2); }, mesh);
  YBoundary sheath(YBndryType::all, nullptr, *mesh);
  sheath.iter(
      [&](const BoundaryIterator auto& point) { EXPECT_DOUBLE_EQ(point.next(test), 4); });
}

TEST_F(YBTest, current_const) {
  const Field3D test = makeField<Field3D>([&](auto& i) { return SQ(i.y() - 2); }, mesh);
  YBoundary sheath(YBndryType::all, nullptr, *mesh);
  sheath.iter([&](const BoundaryIterator auto& point) {
    EXPECT_DOUBLE_EQ(point.current(test), 1);
  });
}

TEST_F(YBTest, prev_const) {
  const Field3D test = makeField<Field3D>([&](auto& i) { return SQ(i.y() - 2); }, mesh);
  YBoundary sheath(YBndryType::all, nullptr, *mesh);
  sheath.iter(
      [&](const BoundaryIterator auto& point) { EXPECT_DOUBLE_EQ(point.prev(test), 0); });
}

TEST_F(YBTest, getAt_square) {
  Field3D test = makeField<Field3D>([&](auto& i) { return SQ(i.y() - 2); }, mesh);
  YBoundary sheath(YBndryType::all, nullptr, *mesh);
  sheath.iter([&](const BoundaryIterator auto& point) {
    EXPECT_DOUBLE_EQ(point.at(test, 0), 4);
    EXPECT_DOUBLE_EQ(point.at(test, 1), 1);
    EXPECT_DOUBLE_EQ(point.at(test, 2), 0);
  });
}

TEST_F(YBTest, at_func) {
  Field3D test = makeField<Field3D>([&](auto& i) { return i.y() - 2; }, mesh);
  YBoundary sheath(YBndryType::all, nullptr, *mesh);
  auto square = [&]([[maybe_unused]] int yo, Ind3D ind) -> BoutReal {
    return SQ(test[ind]);
  };
  sheath.iter([&](const BoundaryIterator auto& point) {
    EXPECT_DOUBLE_EQ(point.at(square, 0), 4);
    EXPECT_DOUBLE_EQ(point.at(square, 1), 1);
    EXPECT_DOUBLE_EQ(point.at(square, 2), 0);
  });
}
