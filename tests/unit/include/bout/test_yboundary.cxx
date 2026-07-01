#include "gtest/gtest.h"

#include "bout/boutexception.hxx"
#include "bout/output_bout_types.hxx"
#include "bout/yboundary_regions.hxx"

#include "fake_mesh_fixture.hxx"
#include "test_extras.hxx"

using YBTest = FakeMeshFixture_tmpl<4, 5, 7>;

using bout::globals::mesh;

TEST_F(YBTest, dirichlet_o2_rgn) {
  dynamic_cast<FakeMesh*>(mesh)->createBoundaries();
  Field3D test = 1.0;
  YBoundary sheath(YBndryType::all, nullptr, *mesh);
  sheath.iter([&](auto& point) { point.dirichlet_o2(test, 2); });
  EXPECT_TRUE(IsFieldEqual(test, 3.0, "RGN_YGUARDS"));
}

TEST_F(YBTest, bndry_size) {
  dynamic_cast<FakeMesh*>(mesh)->createBoundaries();
  YBoundary sheath(YBndryType::all, nullptr, *mesh);
  int sum = 0;
  sheath.iter([&](auto& point) { sum++; });
  EXPECT_EQ(sum, 28);
}

TEST_F(YBTest, interpolate_boundary_o2) {
  Field3D test = 1.0;
  YBoundary sheath(YBndryType::all, nullptr, *mesh);
  sheath.iter(
      [&](auto& point) { EXPECT_DOUBLE_EQ(point.interpolate_boundary_o2(test), 1.0); });
}

TEST_F(YBTest, extrapolate_boundary_o2) {
  Field3D test = 1.0;
  YBoundary sheath(YBndryType::all, nullptr, *mesh);
  sheath.iter(
      [&](auto& point) { EXPECT_DOUBLE_EQ(point.extrapolate_boundary_o2(test), 1.0); });
}

TEST_F(YBTest, dirichlet_o1) {
  Field3D test = 1.0;
  YBoundary sheath(YBndryType::all, nullptr, *mesh);
  sheath.iter([&](auto& point) {
    point.dirichlet_o1(test, 2.0);
    EXPECT_DOUBLE_EQ(point.interpolate_boundary_o2(test), 1.5);
  });
}

TEST_F(YBTest, dirichlet_o2) {
  Field3D test = 1.0;
  YBoundary sheath(YBndryType::all, nullptr, *mesh);
  sheath.iter([&](auto& point) {
    point.dirichlet_o2(test, 2.0);
    EXPECT_DOUBLE_EQ(point.interpolate_boundary_o2(test), 2.0);
  });
}

TEST_F(YBTest, dirichlet_o3) {
  Field3D test = 1.0;
  YBoundary sheath(YBndryType::all, nullptr, *mesh);
  sheath.iter([&](auto& point) {
    point.dirichlet_o3(test, 2.0);
    EXPECT_DOUBLE_EQ(point.interpolate_boundary_o2(test), 7. / 3.);
  });
}

TEST_F(YBTest, extrapolate_boundary_free) {
  Field3D test = makeField<Field3D>([&](auto& i) { return SQ(i.y() - 2) + 1; }, mesh);
  YBoundary sheath(YBndryType::all, nullptr, *mesh);
  sheath.iter([&](auto& point) {
    EXPECT_DOUBLE_EQ(point.extrapolate_boundary_free(
                         test, bout::boundary::BoundaryFreeExtrapolation::limited),
                     3);
    EXPECT_DOUBLE_EQ(point.extrapolate_boundary_free(
                         test, bout::boundary::BoundaryFreeExtrapolation::linear),
                     3.5);
    EXPECT_DOUBLE_EQ(point.extrapolate_boundary_free(
                         test, bout::boundary::BoundaryFreeExtrapolation::exponential),
                     5);
  });
}

TEST_F(YBTest, set_free) {
  Field3D test = makeField<Field3D>([&](auto& i) { return SQ(i.y() - 2) + 1; }, mesh);
  YBoundary sheath(YBndryType::all, nullptr, *mesh);
  sheath.iter([&](auto& point) {
    point.set_free(test, bout::boundary::BoundaryFreeExtrapolation::limited);
    EXPECT_DOUBLE_EQ(point.interpolate_boundary_o2(test), 3);
  });
}

TEST_F(YBTest, interpolate_boundary_o2_square) {
  Field3D test = makeField<Field3D>([&](auto& i) { return SQ(i.y() - 2); }, mesh);
  YBoundary sheath(YBndryType::all, nullptr, *mesh);
  sheath.iter(
      [&](auto& point) { EXPECT_DOUBLE_EQ(point.interpolate_boundary_o2(test), 2.5); });
}

TEST_F(YBTest, extrapolate_boundary_o2_square) {
  Field3D test = makeField<Field3D>([&](auto& i) { return SQ(i.y() - 2); }, mesh);
  YBoundary sheath(YBndryType::all, nullptr, *mesh);
  sheath.iter(
      [&](auto& point) { EXPECT_DOUBLE_EQ(point.extrapolate_boundary_o2(test), 1.5); });
}

TEST_F(YBTest, extrapolate_grad_o2_square) {
  Field3D test = makeField<Field3D>([&](auto& i) { return SQ(i.y() - 2); }, mesh);
  YBoundary sheath(YBndryType::all, nullptr, *mesh);
  sheath.iter([&](auto& point) { EXPECT_DOUBLE_EQ(point.extrapolate_grad_o2(test), 1); });
}

TEST_F(YBTest, extrapolate_next_o1_square) {
  Field3D test = makeField<Field3D>([&](auto& i) { return SQ(i.y() - 2); }, mesh);
  YBoundary sheath(YBndryType::all, nullptr, *mesh);
  sheath.iter([&](auto& point) { EXPECT_DOUBLE_EQ(point.extrapolate_next_o1(test), 1); });
}

TEST_F(YBTest, extrapolate_next_o2_square) {
  Field3D test = makeField<Field3D>([&](auto& i) { return SQ(i.y() - 2); }, mesh);
  YBoundary sheath(YBndryType::all, nullptr, *mesh);
  sheath.iter([&](auto& point) { EXPECT_DOUBLE_EQ(point.extrapolate_next_o2(test), 2); });
}

TEST_F(YBTest, next_square) {
  Field3D test = makeField<Field3D>([&](auto& i) { return SQ(i.y() - 2); }, mesh);
  YBoundary sheath(YBndryType::all, nullptr, *mesh);
  sheath.iter([&](auto& point) { EXPECT_DOUBLE_EQ(point.next(test), 4); });
}

TEST_F(YBTest, current_square) {
  Field3D test = makeField<Field3D>([&](auto& i) { return SQ(i.y() - 2); }, mesh);
  YBoundary sheath(YBndryType::all, nullptr, *mesh);
  sheath.iter([&](auto& point) { EXPECT_DOUBLE_EQ(point.current(test), 1); });
}

TEST_F(YBTest, prev_square) {
  Field3D test = makeField<Field3D>([&](auto& i) { return SQ(i.y() - 2); }, mesh);
  YBoundary sheath(YBndryType::all, nullptr, *mesh);
  sheath.iter([&](auto& point) { EXPECT_DOUBLE_EQ(point.prev(test), 0); });
}

TEST_F(YBTest, getAt_square) {
  Field3D test = makeField<Field3D>([&](auto& i) { return SQ(i.y() - 2); }, mesh);
  YBoundary sheath(YBndryType::all, nullptr, *mesh);
  sheath.iter([&](auto& point) {
    EXPECT_DOUBLE_EQ(point.getAt(test, 0), 4);
    EXPECT_DOUBLE_EQ(point.getAt(test, 1), 1);
    EXPECT_DOUBLE_EQ(point.getAt(test, 2), 0);
  });
}

TEST_F(YBTest, getAt_func) {
  Field3D test = makeField<Field3D>([&](auto& i) { return i.y() - 2; }, mesh);
  YBoundary sheath(YBndryType::all, nullptr, *mesh);
  auto square = [&](int yo, Ind3D ind) -> double { return SQ(test[ind]); };
  sheath.iter([&](auto& point) {
    EXPECT_DOUBLE_EQ(point.getAt(square, 0), 4);
    EXPECT_DOUBLE_EQ(point.getAt(square, 1), 1);
    EXPECT_DOUBLE_EQ(point.getAt(square, 2), 0);
  });
}
