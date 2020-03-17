#include "gtest/gtest.h"

#include "bout/coordinates.hxx"
#include "bout/mesh.hxx"
#include "output.hxx"

#include "test_extras.hxx"

/// Global mesh
namespace bout {
namespace globals {
extern Mesh* mesh;
}
} // namespace bout

using bout::globals::mesh;

using CoordinatesTest = FakeMeshFixture;

TEST_F(CoordinatesTest, ZLength) {
  Coordinates coords{
      mesh,         Field2D{1.0}, Field2D{1.0}, BoutReal{1.0}, Field2D{1.0}, Field2D{1.0},
      Field2D{1.0}, Field2D{1.0}, Field2D{1.0}, Field2D{0.0},  Field2D{0.0}, Field2D{0.0},
      Field2D{1.0}, Field2D{1.0}, Field2D{1.0}, Field2D{0.0},  Field2D{0.0}, Field2D{0.0},
      Field2D{0.0}, Field2D{0.0}, false};

  EXPECT_DOUBLE_EQ(coords.zlength(), 7.0);
}

TEST_F(CoordinatesTest, Jacobian) {
  Coordinates coords{
      mesh,         Field2D{1.0}, Field2D{1.0}, BoutReal{1.0}, Field2D{1.0}, Field2D{1.0},
      Field2D{1.0}, Field2D{1.0}, Field2D{1.0}, Field2D{0.0},  Field2D{0.0}, Field2D{0.0},
      Field2D{1.0}, Field2D{1.0}, Field2D{1.0}, Field2D{0.0},  Field2D{0.0}, Field2D{0.0},
      Field2D{0.0}, Field2D{0.0}, false};

  EXPECT_NO_THROW(coords.jacobian());

  EXPECT_TRUE(IsFieldEqual(coords.J, 1.0));
  EXPECT_TRUE(IsFieldEqual(coords.Bxy, 1.0));
}

/// To do generalise these tests
#ifndef COORDINATES_USE_3D
TEST_F(CoordinatesTest, CalcContravariant) {
  Coordinates coords{
      mesh,         Field2D{1.0}, Field2D{1.0}, BoutReal{1.0}, Field2D{0.0}, Field2D{0.0},
      Field2D{1.0}, Field2D{1.0}, Field2D{1.0}, Field2D{0.0},  Field2D{0.0}, Field2D{0.0},
      Field2D{0.0}, Field2D{0.0}, Field2D{0.0}, Field2D{0.0},  Field2D{0.0}, Field2D{0.0},
      Field2D{0.0}, Field2D{0.0}, false};

  output_info.disable();
  coords.calcCovariant();
  output_info.enable();

  EXPECT_TRUE(IsFieldEqual(coords.g_11, 1.0));
  EXPECT_TRUE(IsFieldEqual(coords.g_22, 1.0));
  EXPECT_TRUE(IsFieldEqual(coords.g_33, 1.0));
  EXPECT_TRUE(IsFieldEqual(coords.g_12, 0.0));
  EXPECT_TRUE(IsFieldEqual(coords.g_13, 0.0));
  EXPECT_TRUE(IsFieldEqual(coords.g_23, 0.0));
}

TEST_F(CoordinatesTest, CalcCovariant) {
  Coordinates coords{
      mesh,         Field2D{1.0}, Field2D{1.0}, BoutReal{1.0}, Field2D{0.0}, Field2D{0.0},
      Field2D{0.0}, Field2D{0.0}, Field2D{0.0}, Field2D{0.0},  Field2D{0.0}, Field2D{0.0},
      Field2D{1.0}, Field2D{1.0}, Field2D{1.0}, Field2D{0.0},  Field2D{0.0}, Field2D{0.0},
      Field2D{0.0}, Field2D{0.0}, false};

  output_info.disable();
  coords.calcContravariant();
  output_info.enable();

  EXPECT_TRUE(IsFieldEqual(coords.g11, 1.0));
  EXPECT_TRUE(IsFieldEqual(coords.g22, 1.0));
  EXPECT_TRUE(IsFieldEqual(coords.g33, 1.0));
  EXPECT_TRUE(IsFieldEqual(coords.g12, 0.0));
  EXPECT_TRUE(IsFieldEqual(coords.g13, 0.0));
  EXPECT_TRUE(IsFieldEqual(coords.g23, 0.0));
}
#endif
