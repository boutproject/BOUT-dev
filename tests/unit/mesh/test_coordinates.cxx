#include "gtest/gtest.h"

#include "bout/build_config.hxx"
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

class CoordinatesTest : public FakeMeshFixture {
public:
  using FieldMetric = Coordinates::FieldMetric;
  CoordinatesTest() : FakeMeshFixture() {}
};

TEST_F(CoordinatesTest, ZLength) {
  Coordinates coords{mesh,
                     FieldMetric{1.0}, // dx
                     FieldMetric{1.0}, // dy
                     FieldMetric{1.0}, // dz
                     FieldMetric{1.0}, // J
                     FieldMetric{1.0}, // Bxy
                     FieldMetric{1.0}, // g11
                     FieldMetric{1.0}, // g22
                     FieldMetric{1.0}, // g33
                     FieldMetric{0.0}, // g12
                     FieldMetric{0.0}, // g13
                     FieldMetric{0.0}, // g23
                     FieldMetric{1.0}, // g_11
                     FieldMetric{1.0}, // g_22
                     FieldMetric{1.0}, // g_23
                     FieldMetric{0.0}, // g_12
                     FieldMetric{0.0}, // g_13
                     FieldMetric{0.0}, // g_23
                     FieldMetric{0.0}, // ShiftTorsion
                     FieldMetric{0.0}, // IntShiftTorsion
                     false};           // calculate_geometry

  EXPECT_TRUE(IsFieldEqual(coords.zlength(), 7.0));
}

#if BOUT_USE_METRIC_3D
TEST_F(CoordinatesTest, ZLength3D) {
  auto dz = makeField<FieldMetric>([](const Ind3D& i) -> BoutReal {
    return static_cast<BoutReal>(i.x() + i.y()) / nz;
  });
  auto expected =
      makeField<Field2D>([](const Ind2D& i) -> BoutReal { return i.x() + i.y(); });

  Coordinates coords{mesh,
                     FieldMetric{1.0}, // dx
                     FieldMetric{1.0}, // dy
                     dz,               // dz
                     FieldMetric{1.0}, // J
                     FieldMetric{1.0}, // Bxy
                     FieldMetric{1.0}, // g11
                     FieldMetric{1.0}, // g22
                     FieldMetric{1.0}, // g33
                     FieldMetric{0.0}, // g12
                     FieldMetric{0.0}, // g13
                     FieldMetric{0.0}, // g23
                     FieldMetric{1.0}, // g_11
                     FieldMetric{1.0}, // g_22
                     FieldMetric{1.0}, // g_23
                     FieldMetric{0.0}, // g_12
                     FieldMetric{0.0}, // g_13
                     FieldMetric{0.0}, // g_23
                     FieldMetric{0.0}, // ShiftTorsion
                     FieldMetric{0.0}, // IntShiftTorsion
                     false};           // calculate_geometry

  EXPECT_TRUE(IsFieldEqual(coords.zlength(), expected));
}
#endif

TEST_F(CoordinatesTest, Jacobian) {
  Coordinates coords{mesh,
                     FieldMetric{1.0}, // dx
                     FieldMetric{1.0}, // dy
                     FieldMetric{1.0}, // dz
                     FieldMetric{1.0}, // J
                     FieldMetric{1.0}, // Bxy
                     FieldMetric{1.0}, // g11
                     FieldMetric{1.0}, // g22
                     FieldMetric{1.0}, // g33
                     FieldMetric{0.0}, // g12
                     FieldMetric{0.0}, // g13
                     FieldMetric{0.0}, // g23
                     FieldMetric{1.0}, // g_11
                     FieldMetric{1.0}, // g_22
                     FieldMetric{1.0}, // g_23
                     FieldMetric{0.0}, // g_12
                     FieldMetric{0.0}, // g_13
                     FieldMetric{0.0}, // g_23
                     FieldMetric{0.0}, // ShiftTorsion
                     FieldMetric{0.0}, // IntShiftTorsion
                     false};           // calculate_geometry

  EXPECT_NO_THROW(coords.jacobian());

  EXPECT_TRUE(IsFieldEqual(coords.J, 1.0));
  EXPECT_TRUE(IsFieldEqual(coords.Bxy, 1.0));
}

/// To do generalise these tests
#if not(BOUT_USE_METRIC_3D)
TEST_F(CoordinatesTest, CalcContravariant) {
  Coordinates coords{mesh,
                     FieldMetric{1.0}, // dx
                     FieldMetric{1.0}, // dy
                     FieldMetric{1.0}, // dz
                     FieldMetric{0.0}, // J
                     FieldMetric{0.0}, // Bxy
                     FieldMetric{1.0}, // g11
                     FieldMetric{1.0}, // g22
                     FieldMetric{1.0}, // g33
                     FieldMetric{0.0}, // g12
                     FieldMetric{0.0}, // g13
                     FieldMetric{0.0}, // g23
                     FieldMetric{0.0}, // g_11
                     FieldMetric{0.0}, // g_22
                     FieldMetric{0.0}, // g_23
                     FieldMetric{0.0}, // g_12
                     FieldMetric{0.0}, // g_13
                     FieldMetric{0.0}, // g_23
                     FieldMetric{0.0}, // ShiftTorsion
                     FieldMetric{0.0}, // IntShiftTorsion
                     false};           // calculate_geometry

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
  Coordinates coords{mesh,
                     FieldMetric{1.0}, // dx
                     FieldMetric{1.0}, // dy
                     FieldMetric{1.0}, // dz
                     FieldMetric{0.0}, // J
                     FieldMetric{0.0}, // Bxy
                     FieldMetric{0.0}, // g11
                     FieldMetric{0.0}, // g22
                     FieldMetric{0.0}, // g33
                     FieldMetric{0.0}, // g12
                     FieldMetric{0.0}, // g13
                     FieldMetric{0.0}, // g23
                     FieldMetric{1.0}, // g_11
                     FieldMetric{1.0}, // g_22
                     FieldMetric{1.0}, // g_23
                     FieldMetric{0.0}, // g_12
                     FieldMetric{0.0}, // g_13
                     FieldMetric{0.0}, // g_23
                     FieldMetric{0.0}, // ShiftTorsion
                     FieldMetric{0.0}, // IntShiftTorsion
                     false};           // calculate_geometry

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
