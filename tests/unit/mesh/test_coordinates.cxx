#include "gtest/gtest.h"

#include "bout/build_config.hxx"
#include "bout/constants.hxx"
#include "bout/coordinates.hxx"
#include "bout/mesh.hxx"
#include "bout/output.hxx"

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

constexpr BoutReal default_dz{TWOPI / CoordinatesTest::nz};

TEST_F(CoordinatesTest, ZLength) {
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
  // No call to Coordinates::geometry() needed here

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
                     FieldMetric{1.0},  // dx
                     FieldMetric{1.0},  // dy
                     dz,                // dz
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
  // No call to Coordinates::geometry() needed here

  EXPECT_TRUE(IsFieldEqual(coords.zlength(), expected));
}
#endif

TEST_F(CoordinatesTest, Jacobian) {
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
  // No call to Coordinates::geometry() needed here

  EXPECT_NO_THROW(coords.jacobian());

  EXPECT_TRUE(IsFieldEqual(coords.J(), 1.0));
  EXPECT_TRUE(IsFieldEqual(coords.Bxy(), 1.0));
}

/// To do generalise these tests
// #if not(BOUT_USE_METRIC_3D)
TEST_F(CoordinatesTest, CalcContravariant) {
  Coordinates coords{mesh,
                     FieldMetric{1.0},  // dx
                     FieldMetric{1.0},  // dy
                     FieldMetric{1.0},  // dz
                     FieldMetric{0.0},  // J
                     FieldMetric{0.0},  // Bxy
                     FieldMetric{1.0},  // g11
                     FieldMetric{1.0},  // g22
                     FieldMetric{1.0},  // g33
                     FieldMetric{0.0},  // g12
                     FieldMetric{0.0},  // g13
                     FieldMetric{0.0},  // g23
                     FieldMetric{0.0},  // g_11
                     FieldMetric{0.0},  // g_22
                     FieldMetric{0.0},  // g_23
                     FieldMetric{0.0},  // g_12
                     FieldMetric{0.0},  // g_13
                     FieldMetric{0.0},  // g_23
                     FieldMetric{0.0},  // ShiftTorsion
                     FieldMetric{0.0}}; // IntShiftTorsion
  // No call to Coordinates::geometry() needed here

  output_info.disable();
  coords.calcCovariant();
  output_info.enable();

  EXPECT_TRUE(IsFieldEqual(coords.g_11(), 1.0));
  EXPECT_TRUE(IsFieldEqual(coords.g_22(), 1.0));
  EXPECT_TRUE(IsFieldEqual(coords.g_33(), 1.0));
  EXPECT_TRUE(IsFieldEqual(coords.g_12(), 0.0));
  EXPECT_TRUE(IsFieldEqual(coords.g_13(), 0.0));
  EXPECT_TRUE(IsFieldEqual(coords.g_23(), 0.0));
}

TEST_F(CoordinatesTest, CalcCovariant) {
  Coordinates coords{mesh,
                     FieldMetric{1.0},  // dx
                     FieldMetric{1.0},  // dy
                     FieldMetric{1.0},  // dz
                     FieldMetric{0.0},  // J
                     FieldMetric{0.0},  // Bxy
                     FieldMetric{0.0},  // g11
                     FieldMetric{0.0},  // g22
                     FieldMetric{0.0},  // g33
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
  // No call to Coordinates::geometry() needed here

  output_info.disable();
  coords.calcContravariant();
  output_info.enable();

  EXPECT_TRUE(IsFieldEqual(coords.g11(), 1.0));
  EXPECT_TRUE(IsFieldEqual(coords.g22(), 1.0));
  EXPECT_TRUE(IsFieldEqual(coords.g33(), 1.0));
  EXPECT_TRUE(IsFieldEqual(coords.g12(), 0.0));
  EXPECT_TRUE(IsFieldEqual(coords.g13(), 0.0));
  EXPECT_TRUE(IsFieldEqual(coords.g23(), 0.0));
}
// #endif

TEST_F(CoordinatesTest, DefaultConstructor) {
  output_info.disable();
  output_warn.disable();
  Coordinates coords(mesh);
  output_warn.enable();
  output_info.enable();

  EXPECT_TRUE(IsFieldEqual(coords.dx, 1.0));
  EXPECT_TRUE(IsFieldEqual(coords.dy, 1.0));
  EXPECT_TRUE(IsFieldEqual(coords.dz, default_dz));

  EXPECT_TRUE(IsFieldEqual(coords.g11(), 1.0));
  EXPECT_TRUE(IsFieldEqual(coords.g22(), 1.0));
  EXPECT_TRUE(IsFieldEqual(coords.g33(), 1.0));
  EXPECT_TRUE(IsFieldEqual(coords.g12(), 0.0));
  EXPECT_TRUE(IsFieldEqual(coords.g13(), 0.0));
  EXPECT_TRUE(IsFieldEqual(coords.g23(), 0.0));

  EXPECT_TRUE(IsFieldEqual(coords.J(), 1.0));
  EXPECT_TRUE(IsFieldEqual(coords.Bxy(), 1.0));
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
  EXPECT_TRUE(IsFieldEqual(coords.dz, 42.));

  EXPECT_TRUE(IsFieldEqual(coords.g11(), 1.0));
  EXPECT_TRUE(IsFieldEqual(coords.g22(), 1.0));
  EXPECT_TRUE(IsFieldEqual(coords.g33(), 1.0));
  EXPECT_TRUE(IsFieldEqual(coords.g12(), 0.0));
  EXPECT_TRUE(IsFieldEqual(coords.g13(), 0.0));
  EXPECT_TRUE(IsFieldEqual(coords.g23(), 0.0));

  EXPECT_TRUE(IsFieldEqual(coords.J(), 1.0));
  EXPECT_TRUE(IsFieldEqual(coords.Bxy(), 1.0));
}

TEST_F(CoordinatesTest, SmallMeshSpacing) {
  static_cast<FakeMesh*>(bout::globals::mesh)
      ->setGridDataSource(new FakeGridDataSource({{"dx", 1e-9}}));

  output_info.disable();
  output_warn.disable();
  Coordinates coords(mesh);
  EXPECT_THROW(coords.geometry(), BoutException);
  output_warn.enable();
  output_info.enable();
}

TEST_F(CoordinatesTest, ConstructWithDiagonalContravariantMetric) {

  static_cast<FakeMesh*>(bout::globals::mesh)
      ->setGridDataSource(
          new FakeGridDataSource({{"g11", 2.0}, {"g22", 3.2}, {"g33", 42}}));

  output_info.disable();
  output_warn.disable();
  Coordinates coords(mesh);
  output_warn.enable();
  output_info.enable();

  // Didn't specify grid spacing, so default to 1
  EXPECT_TRUE(IsFieldEqual(coords.dx, 1.0));
  EXPECT_TRUE(IsFieldEqual(coords.dy, 1.0));
  EXPECT_TRUE(IsFieldEqual(coords.dz, default_dz));

  // Diagonal contravariant metric
  EXPECT_TRUE(IsFieldEqual(coords.g11(), 2.0));
  EXPECT_TRUE(IsFieldEqual(coords.g22(), 3.2));
  EXPECT_TRUE(IsFieldEqual(coords.g33(), 42));
  EXPECT_TRUE(IsFieldEqual(coords.g12(), 0.0));
  EXPECT_TRUE(IsFieldEqual(coords.g13(), 0.0));
  EXPECT_TRUE(IsFieldEqual(coords.g23(), 0.0));

  // Covariant metric should be inverse
  // Note: Not calculated in corners
  EXPECT_TRUE(IsFieldEqual(coords.g_11(), 1. / 2.0, "RGN_NOCORNERS"));
  EXPECT_TRUE(IsFieldEqual(coords.g_22(), 1. / 3.2, "RGN_NOCORNERS"));
  EXPECT_TRUE(IsFieldEqual(coords.g_33(), 1. / 42, "RGN_NOCORNERS"));

  EXPECT_TRUE(IsFieldEqual(coords.J(), 1. / sqrt(2.0 * 3.2 * 42), "RGN_NOCORNERS"));
  EXPECT_TRUE(IsFieldEqual(coords.Bxy(), sqrt(2.0 * 42), "RGN_NOCORNERS", 1e-10));
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

TEST_F(CoordinatesTest, GetContravariantMetricTensor) {
  Coordinates coords{mesh,
                     FieldMetric{1.0},  // dx
                     FieldMetric{1.0},  // dy
                     FieldMetric{1.0},  // dz
                     FieldMetric{0.0},  // J
                     FieldMetric{0.0},  // Bxy
                     FieldMetric{1.2},  // g11
                     FieldMetric{2.3},  // g22
                     FieldMetric{3.4},  // g33
                     FieldMetric{4.5},  // g12
                     FieldMetric{5.6},  // g13
                     FieldMetric{6.7},  // g23
                     FieldMetric{1.0},  // g_11
                     FieldMetric{1.0},  // g_22
                     FieldMetric{1.0},  // g_23
                     FieldMetric{0.0},  // g_12
                     FieldMetric{0.0},  // g_13
                     FieldMetric{0.0},  // g_23
                     FieldMetric{0.0},  // ShiftTorsion
                     FieldMetric{0.0}}; // IntShiftTorsion

  EXPECT_TRUE(IsFieldEqual(coords.g11(), 1.2));
  EXPECT_TRUE(IsFieldEqual(coords.g22(), 2.3));
  EXPECT_TRUE(IsFieldEqual(coords.g33(), 3.4));
  EXPECT_TRUE(IsFieldEqual(coords.g12(), 4.5));
  EXPECT_TRUE(IsFieldEqual(coords.g13(), 5.6));
  EXPECT_TRUE(IsFieldEqual(coords.g23(), 6.7));
}

TEST_F(CoordinatesTest, SetContravariantMetricTensor) {
  // Set initial values for the metric tensor in the Coordinates constructor
  Coordinates coords{mesh,
                     FieldMetric{1.0},  // dx
                     FieldMetric{1.0},  // dy
                     FieldMetric{1.0},  // dz
                     FieldMetric{0.0},  // J
                     FieldMetric{0.0},  // Bxy
                     FieldMetric{0.0},  // g11
                     FieldMetric{0.0},  // g22
                     FieldMetric{0.0},  // g33
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

  //  Modify with setter
  auto updated_metric_tensor = MetricTensor(1.7, 2.3, 3.1, 0.9, 5.7, 1.9);
  coords.setContravariantMetricTensor(updated_metric_tensor);

  //  Get values with getter and check they have been modified as expected
  EXPECT_TRUE(IsFieldEqual(coords.g11(), 1.7));
  EXPECT_TRUE(IsFieldEqual(coords.g22(), 2.3));
  EXPECT_TRUE(IsFieldEqual(coords.g33(), 3.1));
  EXPECT_TRUE(IsFieldEqual(coords.g12(), 0.9));
  EXPECT_TRUE(IsFieldEqual(coords.g13(), 5.7));
  EXPECT_TRUE(IsFieldEqual(coords.g23(), 1.9));
}

TEST_F(CoordinatesTest, GetCovariantMetricTensor) {
  Coordinates coords{mesh,
                     FieldMetric{1.0},  // dx
                     FieldMetric{1.0},  // dy
                     FieldMetric{1.0},  // dz
                     FieldMetric{0.0},  // J
                     FieldMetric{0.0},  // Bxy
                     FieldMetric{1.2},  // g11
                     FieldMetric{2.3},  // g22
                     FieldMetric{3.4},  // g33
                     FieldMetric{4.5},  // g12
                     FieldMetric{5.6},  // g13
                     FieldMetric{6.7},  // g23
                     FieldMetric{9.7},  // g_11
                     FieldMetric{7.5},  // g_22
                     FieldMetric{4.7},  // g_23
                     FieldMetric{3.9},  // g_12
                     FieldMetric{1.7},  // g_13
                     FieldMetric{5.3},  // g_23
                     FieldMetric{0.0},  // ShiftTorsion
                     FieldMetric{0.0}}; // IntShiftTorsion

  EXPECT_TRUE(IsFieldEqual(coords.g_11(), 9.7));
  EXPECT_TRUE(IsFieldEqual(coords.g_22(), 7.5));
  EXPECT_TRUE(IsFieldEqual(coords.g_33(), 4.7));
  EXPECT_TRUE(IsFieldEqual(coords.g_12(), 3.9));
  EXPECT_TRUE(IsFieldEqual(coords.g_13(), 1.7));
  EXPECT_TRUE(IsFieldEqual(coords.g_23(), 5.3));
}

TEST_F(CoordinatesTest, SetCovariantMetricTensor) {
  {
    // Set initial values for the metric tensor in the Coordinates constructor
    Coordinates coords{mesh,
                       FieldMetric{1.0},  // dx
                       FieldMetric{1.0},  // dy
                       FieldMetric{1.0},  // dz
                       FieldMetric{0.0},  // J
                       FieldMetric{0.0},  // Bxy
                       FieldMetric{0.0},  // g11
                       FieldMetric{0.0},  // g22
                       FieldMetric{0.0},  // g33
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

    //  Modify with setter
    auto updated_metric_tensor = MetricTensor(1.7, 2.3, 3.1, 0.9, 5.7, 1.9);
    coords.setCovariantMetricTensor(updated_metric_tensor);

    //  Get values with getter and check they have been modified as expected
    EXPECT_TRUE(IsFieldEqual(coords.g_11(), 1.7));
    EXPECT_TRUE(IsFieldEqual(coords.g_22(), 2.3));
    EXPECT_TRUE(IsFieldEqual(coords.g_33(), 3.1));
    EXPECT_TRUE(IsFieldEqual(coords.g_12(), 0.9));
    EXPECT_TRUE(IsFieldEqual(coords.g_13(), 5.7));
    EXPECT_TRUE(IsFieldEqual(coords.g_23(), 1.9));
  }
}