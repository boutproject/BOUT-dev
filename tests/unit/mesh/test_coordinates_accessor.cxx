#include "gtest/gtest.h"

#include "bout/build_defines.hxx"
#include "bout/coordinates_accessor.hxx"

#include "fake_mesh_fixture.hxx"

using bout::globals::mesh;

using CoordinatesAccessorTest = FakeMeshFixture;
using FieldMetric = Coordinates::FieldMetric;

TEST_F(CoordinatesAccessorTest, CreateAccessor) {
  CoordinatesAccessor acc(mesh->getCoordinates());

  // Basic sanity checks
  EXPECT_TRUE(acc.mesh_nz == mesh->LocalNz);
  EXPECT_TRUE(acc.data != nullptr);
  EXPECT_FLOAT_EQ(mesh->getCoordinates()->dx(0, 0, 0), acc.dx(0));
}

TEST_F(CoordinatesAccessorTest, CreateTwoAccessors) {
  CoordinatesAccessor acc(mesh->getCoordinates());

  // Create another accessor from the same Coordinates
  CoordinatesAccessor acc2(mesh->getCoordinates());

  // Should share the same data pointer
  EXPECT_EQ(acc.data, acc2.data);
}

TEST_F(CoordinatesAccessorTest, ClearTwice) {
  CoordinatesAccessor::clear(); // Should clear any left over

  EXPECT_EQ(CoordinatesAccessor::clear(), 0);
}

TEST_F(CoordinatesAccessorTest, ClearOne) {
  CoordinatesAccessor::clear(); // Should clear any left over

  CoordinatesAccessor acc(mesh->getCoordinates());

  EXPECT_EQ(CoordinatesAccessor::clear(), 1);
}

TEST_F(CoordinatesAccessorTest, ClearOneReused) {
  CoordinatesAccessor::clear(); // Should clear any left over

  CoordinatesAccessor acc(mesh->getCoordinates());
  CoordinatesAccessor acc2(mesh->getCoordinates());

  EXPECT_EQ(CoordinatesAccessor::clear(), 1);
}

TEST_F(CoordinatesAccessorTest, ClearBoth) {
  CoordinatesAccessor::clear(); // Should clear any left over

  // Make a new Coordinates object to access
  Coordinates coords{mesh,
                     FieldMetric{1.0},  // dx
                     FieldMetric{1.0},  // dy
                     FieldMetric{1.0},  // dz
                     FieldMetric{1.0},  // J
                     FieldMetric{1.0},  // Bxy
                     FieldMetric{1.0},  // g11
                     FieldMetric{1.0},  // g22
                     FieldMetric{1.0},  // g33
                     FieldMetric{0.0},  // g12
                     FieldMetric{0.0},  // g13
                     FieldMetric{0.0},  // g23
                     FieldMetric{1.0},  // g_11
                     FieldMetric{1.0},  // g_22
                     FieldMetric{1.0},  // g_23
                     FieldMetric{0.0},  // g_12
                     FieldMetric{0.0},  // g_13
                     FieldMetric{0.0},  // g_23
                     FieldMetric{0.0},  // ShiftTorsion
                     FieldMetric{0.0}}; // IntShiftTorsion
  // Need to set geometry information
  coords.G1 = coords.G2 = coords.G3 = 0.2;
  coords.non_uniform = true;
  coords.d1_dx = coords.d1_dy = coords.d1_dz = 0.1;
#if BOUT_USE_METRIC_3D
  coords.Bxy.splitParallelSlices();
  coords.Bxy.yup() = coords.Bxy.ydown() = coords.Bxy;
#endif

  CoordinatesAccessor acc(mesh->getCoordinates());
  CoordinatesAccessor acc2(&coords); // Different from previous Coordinates

  // Clear both
  EXPECT_EQ(CoordinatesAccessor::clear(), 2);
}

TEST_F(CoordinatesAccessorTest, ClearOneTwo) {
  CoordinatesAccessor::clear(); // Should clear any left over

  // Make a new Coordinates object to access
  Coordinates coords{mesh,
                     FieldMetric{1.0},  // dx
                     FieldMetric{1.0},  // dy
                     FieldMetric{1.0},  // dz
                     FieldMetric{1.0},  // J
                     FieldMetric{1.0},  // Bxy
                     FieldMetric{1.0},  // g11
                     FieldMetric{1.0},  // g22
                     FieldMetric{1.0},  // g33
                     FieldMetric{0.0},  // g12
                     FieldMetric{0.0},  // g13
                     FieldMetric{0.0},  // g23
                     FieldMetric{1.0},  // g_11
                     FieldMetric{1.0},  // g_22
                     FieldMetric{1.0},  // g_23
                     FieldMetric{0.0},  // g_12
                     FieldMetric{0.0},  // g_13
                     FieldMetric{0.0},  // g_23
                     FieldMetric{0.0},  // ShiftTorsion
                     FieldMetric{0.0}}; // IntShiftTorsion
  // Need to set geometry information
  coords.G1 = coords.G2 = coords.G3 = 0.2;
  coords.non_uniform = true;
  coords.d1_dx = coords.d1_dy = coords.d1_dz = 0.1;
#if BOUT_USE_METRIC_3D
  coords.Bxy.splitParallelSlices();
  coords.Bxy.yup() = coords.Bxy.ydown() = coords.Bxy;
#endif

  CoordinatesAccessor acc(mesh->getCoordinates());
  CoordinatesAccessor acc2(&coords); // Different from previous Coordinates

  // clear first one
  EXPECT_EQ(CoordinatesAccessor::clear(mesh->getCoordinates()), 1);
  // then the second one
  EXPECT_EQ(CoordinatesAccessor::clear(), 1);
}

TEST_F(CoordinatesAccessorTest, ClearTwoOneNone) {
  CoordinatesAccessor::clear(); // Should clear any left over

  // Make a new Coordinates object to access
  Coordinates coords{mesh,
                     FieldMetric{1.0},  // dx
                     FieldMetric{1.0},  // dy
                     FieldMetric{1.0},  // dz
                     FieldMetric{1.0},  // J
                     FieldMetric{1.0},  // Bxy
                     FieldMetric{1.0},  // g11
                     FieldMetric{1.0},  // g22
                     FieldMetric{1.0},  // g33
                     FieldMetric{0.0},  // g12
                     FieldMetric{0.0},  // g13
                     FieldMetric{0.0},  // g23
                     FieldMetric{1.0},  // g_11
                     FieldMetric{1.0},  // g_22
                     FieldMetric{1.0},  // g_23
                     FieldMetric{0.0},  // g_12
                     FieldMetric{0.0},  // g_13
                     FieldMetric{0.0},  // g_23
                     FieldMetric{0.0},  // ShiftTorsion
                     FieldMetric{0.0}}; // IntShiftTorsion
  // Need to set geometry information
  coords.G1 = coords.G2 = coords.G3 = 0.2;
  coords.non_uniform = true;
  coords.d1_dx = coords.d1_dy = coords.d1_dz = 0.1;
#if BOUT_USE_METRIC_3D
  coords.Bxy.splitParallelSlices();
  coords.Bxy.yup() = coords.Bxy.ydown() = coords.Bxy;
#endif

  CoordinatesAccessor acc(mesh->getCoordinates());
  CoordinatesAccessor acc2(&coords); // Different from previous Coordinates

  // clear second one
  EXPECT_EQ(CoordinatesAccessor::clear(&coords), 1);
  // then the second one
  EXPECT_EQ(CoordinatesAccessor::clear(mesh->getCoordinates()), 1);
  // Should be none left
  EXPECT_EQ(CoordinatesAccessor::clear(), 0);
}
