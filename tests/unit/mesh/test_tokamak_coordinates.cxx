#include "gtest/gtest.h"
#include <cmath>

#include "bout/build_defines.hxx"
#include "bout/constants.hxx"
#include "bout/mesh.hxx"
#include "bout/output.hxx"
#include "bout/tokamak_coordinates.hxx"

#include "fake_mesh.hxx"
#include "fake_mesh_fixture.hxx"
#include "test_extras.hxx"

using bout::globals::mesh;

using TokamakCoordinatesTestFCI = FakeMeshFixtureFCI;

using TokamakCoordinatesTest = FakeMeshFixture;

// Non-FCI tokamak coordinates should return the full component-wise
// normalisation contract, including dx/J factors and no aggregate `g`.
TEST_F(TokamakCoordinatesTest, Normalisation) {
  const auto norm = bout::TokamakOrFCIMetricNormaliser(mesh, 2.0, 3.0);

  EXPECT_DOUBLE_EQ(*norm.g11, 1. / 36.);
  EXPECT_DOUBLE_EQ(*norm.g22, 9.);
  EXPECT_DOUBLE_EQ(*norm.g33, 9.);
  EXPECT_DOUBLE_EQ(*norm.g12, 1. / 2.);
  EXPECT_DOUBLE_EQ(*norm.g13, 1. / 2.);
  EXPECT_DOUBLE_EQ(*norm.g23, 9.);
  EXPECT_DOUBLE_EQ(*norm.dx, 18.);
  EXPECT_DOUBLE_EQ(*norm.J, 1.5);
  EXPECT_DOUBLE_EQ(*norm.Bxy, 2.);
  EXPECT_FALSE(norm.g.has_value());
}

// FCI coordinates should use the aggregate `g`/`J` contract and leave the
// component-wise metric normalisation factors unset.
TEST_F(TokamakCoordinatesTestFCI, Normalisation) {
  const auto norm = bout::TokamakOrFCIMetricNormaliser(mesh, 2.0, 3.0);

  EXPECT_DOUBLE_EQ(*norm.g, 9.);
  EXPECT_DOUBLE_EQ(*norm.J, 27.);
  EXPECT_DOUBLE_EQ(*norm.Bxy, 2.);
  EXPECT_FALSE(norm.g11.has_value());
  EXPECT_FALSE(norm.g22.has_value());
  EXPECT_FALSE(norm.g33.has_value());
  EXPECT_FALSE(norm.g12.has_value());
  EXPECT_FALSE(norm.g13.has_value());
  EXPECT_FALSE(norm.g23.has_value());
  EXPECT_FALSE(norm.dx.has_value());
}
