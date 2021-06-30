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

const BoutReal default_dz{2 * 3.141592653589793 / 7};

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

TEST_F(CoordinatesTest, DefaultConstructor) {
  output_info.disable();
  output_warn.disable();
  Coordinates coords(mesh);
  output_warn.enable();
  output_info.enable();

  EXPECT_TRUE(IsFieldEqual(coords.dx, 1.0));
  EXPECT_TRUE(IsFieldEqual(coords.dy, 1.0));
  EXPECT_DOUBLE_EQ(coords.dz, default_dz);

  EXPECT_TRUE(IsFieldEqual(coords.g11, 1.0));
  EXPECT_TRUE(IsFieldEqual(coords.g22, 1.0));
  EXPECT_TRUE(IsFieldEqual(coords.g33, 1.0));
  EXPECT_TRUE(IsFieldEqual(coords.g12, 0.0));
  EXPECT_TRUE(IsFieldEqual(coords.g13, 0.0));
  EXPECT_TRUE(IsFieldEqual(coords.g23, 0.0));

  EXPECT_TRUE(IsFieldEqual(coords.J, 1.0));
  EXPECT_TRUE(IsFieldEqual(coords.Bxy, 1.0));
}

TEST_F(CoordinatesTest, ConstructWithMeshSpacing) {

  static_cast<FakeMesh*>(bout::globals::mesh)
      ->setGridDataSource(new FakeGridDataSource({{"dx", 2.0}, {"dy", 3.2}, {"dz", 42}}));

  output_info.disable();
  output_warn.disable();
  Coordinates coords(mesh);
  output_warn.enable();
  output_info.enable();

  EXPECT_TRUE(IsFieldEqual(coords.dx, 2.0));
  EXPECT_TRUE(IsFieldEqual(coords.dy, 3.2));
  EXPECT_DOUBLE_EQ(coords.dz, 42.);

  EXPECT_TRUE(IsFieldEqual(coords.g11, 1.0));
  EXPECT_TRUE(IsFieldEqual(coords.g22, 1.0));
  EXPECT_TRUE(IsFieldEqual(coords.g33, 1.0));
  EXPECT_TRUE(IsFieldEqual(coords.g12, 0.0));
  EXPECT_TRUE(IsFieldEqual(coords.g13, 0.0));
  EXPECT_TRUE(IsFieldEqual(coords.g23, 0.0));

  EXPECT_TRUE(IsFieldEqual(coords.J, 1.0));
  EXPECT_TRUE(IsFieldEqual(coords.Bxy, 1.0));
}

TEST_F(CoordinatesTest, SmallMeshSpacing) {
  static_cast<FakeMesh*>(bout::globals::mesh)
      ->setGridDataSource(new FakeGridDataSource({{"dx", 1e-9}}));

  output_info.disable();
  output_warn.disable();
  EXPECT_THROW(Coordinates coords(mesh), BoutException);
  output_warn.enable();
  output_info.enable();
}

TEST_F(CoordinatesTest, ConstructWithDiagonalContravariantMetric) {

  static_cast<FakeMesh*>(bout::globals::mesh)
      ->setGridDataSource(new FakeGridDataSource({{"g11", 2.0}, {"g22", 3.2}, {"g33", 42}}));

  output_info.disable();
  output_warn.disable();
  Coordinates coords(mesh);
  output_warn.enable();
  output_info.enable();

  // Didn't specify grid spacing, so default to 1
  EXPECT_TRUE(IsFieldEqual(coords.dx, 1.0));
  EXPECT_TRUE(IsFieldEqual(coords.dy, 1.0));
  EXPECT_DOUBLE_EQ(coords.dz, default_dz);

  // Diagonal contravariant metric
  EXPECT_TRUE(IsFieldEqual(coords.g11, 2.0));
  EXPECT_TRUE(IsFieldEqual(coords.g22, 3.2));
  EXPECT_TRUE(IsFieldEqual(coords.g33, 42));
  EXPECT_TRUE(IsFieldEqual(coords.g12, 0.0));
  EXPECT_TRUE(IsFieldEqual(coords.g13, 0.0));
  EXPECT_TRUE(IsFieldEqual(coords.g23, 0.0));

  // Covariant metric should be inverse
  // Note: Not calculated in corners
  EXPECT_TRUE(IsFieldEqual(coords.g_11, 1./2.0, "RGN_NOCORNERS"));
  EXPECT_TRUE(IsFieldEqual(coords.g_22, 1./3.2, "RGN_NOCORNERS"));
  EXPECT_TRUE(IsFieldEqual(coords.g_33, 1./42, "RGN_NOCORNERS"));

  EXPECT_TRUE(IsFieldEqual(coords.J, 1. / sqrt(2.0 * 3.2 * 42), "RGN_NOCORNERS"));
  EXPECT_TRUE(IsFieldEqual(coords.Bxy, sqrt(2.0 * 42), "RGN_NOCORNERS", 1e-10));
}

TEST_F(CoordinatesTest, NegativeJacobian) {
  static_cast<FakeMesh*>(bout::globals::mesh)
      ->setGridDataSource(new FakeGridDataSource({{"J", -1.0}}));

  output_info.disable();
  output_warn.disable();
  EXPECT_THROW(Coordinates coords(mesh), BoutException);
  output_warn.enable();
  output_info.enable();
}

TEST_F(CoordinatesTest, NegativeB) {
  static_cast<FakeMesh*>(bout::globals::mesh)
      ->setGridDataSource(new FakeGridDataSource({{"Bxy", -1.0}}));

  output_info.disable();
  output_warn.disable();
  EXPECT_THROW(Coordinates coords(mesh), BoutException);
  output_warn.enable();
  output_info.enable();
}
