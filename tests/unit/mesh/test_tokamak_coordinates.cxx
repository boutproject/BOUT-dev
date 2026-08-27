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

class TokamakCoordinatesTest : public FakeMeshFixture {
public:
  using FieldMetric = Coordinates::FieldMetric;
  WithQuietOutput info{output_info};
  WithQuietOutput warn{output_warn};
  WithQuietOutput progress{output_progress};
};

constexpr BoutReal default_dz{TWOPI / TokamakCoordinatesTest::nz};

TEST_F(TokamakCoordinatesTest, Normalisation) {
  const auto norm = bout::TokamakOrFCIMetricNormaliser(mesh, 2.0, 3.0);

  EXPECT_DOUBLE_EQ(*norm.g11, 1. / 36.);
  EXPECT_DOUBLE_EQ(*norm.g22, 9.);
  EXPECT_DOUBLE_EQ(*norm.g33, 9.);
  EXPECT_DOUBLE_EQ(*norm.g12, 1. / 2.);
  EXPECT_DOUBLE_EQ(*norm.g13, 1. / 2.);
  EXPECT_DOUBLE_EQ(*norm.g23, 9.);
  EXPECT_DOUBLE_EQ(*norm.Bxy, 2.);
}

TEST_F(TokamakCoordinatesTestFCI, Normalisation) {
  const auto norm = bout::TokamakOrFCIMetricNormaliser(mesh, 2.0, 3.0);

  EXPECT_DOUBLE_EQ(*norm.g, 9.);
  EXPECT_DOUBLE_EQ(*norm.J, 27.);
  EXPECT_DOUBLE_EQ(*norm.Bxy, 2.);
}
