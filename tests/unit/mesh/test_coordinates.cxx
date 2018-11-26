#include "gtest/gtest.h"

#include "bout/coordinates.hxx"
#include "bout/mesh.hxx"
#include "output.hxx"

#include "test_extras.hxx"

/// Global mesh
extern Mesh *mesh;

using CoordinatesTest = FakeMeshFixture;

TEST_F(CoordinatesTest, ZLength) {
  FakeCoordinates coords{mesh,         Field2D{1.0}, Field2D{1.0}, BoutReal{1.0},
                         Field2D{1.0}, Field2D{1.0}, Field2D{1.0}, Field2D{1.0},
                         Field2D{1.0}, Field2D{0.0}, Field2D{0.0}, Field2D{0.0},
                         Field2D{1.0}, Field2D{1.0}, Field2D{1.0}, Field2D{0.0},
                         Field2D{0.0}, Field2D{0.0}};

  EXPECT_DOUBLE_EQ(coords.zlength(), 7.0);
}

TEST_F(CoordinatesTest, Jacobian) {
  FakeCoordinates coords{mesh,         Field2D{1.0}, Field2D{1.0}, BoutReal{1.0},
                         Field2D{1.0}, Field2D{1.0}, Field2D{1.0}, Field2D{1.0},
                         Field2D{1.0}, Field2D{0.0}, Field2D{0.0}, Field2D{0.0},
                         Field2D{1.0}, Field2D{1.0}, Field2D{1.0}, Field2D{0.0},
                         Field2D{0.0}, Field2D{0.0}};

  EXPECT_NO_THROW(coords.jacobian());

  EXPECT_TRUE(IsField3DEqualBoutReal(coords.J, 1.0));
  EXPECT_TRUE(IsField3DEqualBoutReal(coords.Bxy, 1.0));
}

TEST_F(CoordinatesTest, CalcContravariant) {
  FakeCoordinates coords{mesh,         Field2D{1.0}, Field2D{1.0}, BoutReal{1.0},
                         Field2D{0.0}, Field2D{0.0}, Field2D{1.0}, Field2D{1.0},
                         Field2D{1.0}, Field2D{0.0}, Field2D{0.0}, Field2D{0.0},
                         Field2D{0.0}, Field2D{0.0}, Field2D{0.0}, Field2D{0.0},
                         Field2D{0.0}, Field2D{0.0}};

  output_info.disable();
  coords.calcCovariant();
  output_info.enable();

  EXPECT_TRUE(IsField2DEqualBoutReal(coords.g_11, 1.0));
  EXPECT_TRUE(IsField2DEqualBoutReal(coords.g_22, 1.0));
  EXPECT_TRUE(IsField2DEqualBoutReal(coords.g_33, 1.0));
  EXPECT_TRUE(IsField2DEqualBoutReal(coords.g_12, 0.0));
  EXPECT_TRUE(IsField2DEqualBoutReal(coords.g_13, 0.0));
  EXPECT_TRUE(IsField2DEqualBoutReal(coords.g_23, 0.0));
}

TEST_F(CoordinatesTest, CalcCovariant) {
  FakeCoordinates coords{mesh,         Field2D{1.0}, Field2D{1.0}, BoutReal{1.0},
                         Field2D{0.0}, Field2D{0.0}, Field2D{0.0}, Field2D{0.0},
                         Field2D{0.0}, Field2D{0.0}, Field2D{0.0}, Field2D{0.0},
                         Field2D{1.0}, Field2D{1.0}, Field2D{1.0}, Field2D{0.0},
                         Field2D{0.0}, Field2D{0.0}};

  output_info.disable();
  coords.calcContravariant();
  output_info.enable();

  EXPECT_TRUE(IsField2DEqualBoutReal(coords.g11, 1.0));
  EXPECT_TRUE(IsField2DEqualBoutReal(coords.g22, 1.0));
  EXPECT_TRUE(IsField2DEqualBoutReal(coords.g33, 1.0));
  EXPECT_TRUE(IsField2DEqualBoutReal(coords.g12, 0.0));
  EXPECT_TRUE(IsField2DEqualBoutReal(coords.g13, 0.0));
  EXPECT_TRUE(IsField2DEqualBoutReal(coords.g23, 0.0));
}
