#include "gtest/gtest.h"

#include "bout/build_defines.hxx"
#include "bout/constants.hxx"
#include "bout/coordinates.hxx"
#include "bout/mesh.hxx"
#include "bout/output.hxx"

#include "fake_mesh_fixture.hxx"
#include "test_extras.hxx"

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

  EXPECT_NO_THROW(coords.recalculateJacobian());

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

  coords.setContravariantMetricTensor(
      ContravariantMetricTensor(1.0, 1.0, 1.0, 0.0, 0.0, 0.0));

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

  coords.setCovariantMetricTensor(CovariantMetricTensor(1.0, 1.0, 1.0, 0.0, 0.0, 0.0));

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

  EXPECT_TRUE(IsFieldEqual(coords.dx(), 1.0));
  EXPECT_TRUE(IsFieldEqual(coords.dy(), 1.0));
  EXPECT_TRUE(IsFieldEqual(coords.dz(), default_dz));

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

  EXPECT_TRUE(IsFieldEqual(coords.dx(), 2.0));
  EXPECT_TRUE(IsFieldEqual(coords.dy(), 3.2));
  EXPECT_TRUE(IsFieldEqual(coords.dz(), 42.));

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

  WithQuietOutput quiet_info{output_info};
  WithQuietOutput quiet_warn{output_warn};
  EXPECT_THROW(Coordinates{mesh}, BoutException);
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
  EXPECT_TRUE(IsFieldEqual(coords.dx(), 1.0));
  EXPECT_TRUE(IsFieldEqual(coords.dy(), 1.0));
  EXPECT_TRUE(IsFieldEqual(coords.dz(), default_dz));

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
  auto updated_metric_tensor = ContravariantMetricTensor(1.0, 2.0, 0.4, 1.0, 0.0, 0.2);
  coords.setContravariantMetricTensor(updated_metric_tensor);

  //  Get values with getter and check they have been modified as expected
  EXPECT_TRUE(IsFieldEqual(coords.g11(), 1.0));
  EXPECT_TRUE(IsFieldEqual(coords.g22(), 2.0));
  EXPECT_TRUE(IsFieldEqual(coords.g33(), 0.4));
  EXPECT_TRUE(IsFieldEqual(coords.g12(), 1.0));
  EXPECT_TRUE(IsFieldEqual(coords.g13(), 0.0));
  EXPECT_TRUE(IsFieldEqual(coords.g23(), 0.2));
}

TEST_F(CoordinatesTest, CheckCovariantCalculatedFromContravariant) {

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

  //  Modify contravariant components
  constexpr double g11 = 1.0;
  constexpr double g22 = 1.0;
  constexpr double g33 = 1.0;
  const double g12 = sqrt(3.0) / 4.0; // (std::sqrt only constexpr since C++26)
  constexpr double g13 = 0.5;
  const double g23 = sqrt(3.0) / 4.0; // (std::sqrt only constexpr since C++26)
  auto updated_metric_tensor = ContravariantMetricTensor(g11, g22, g33, g12, g13, g23);
  coords.setContravariantMetricTensor(updated_metric_tensor);

  //  Check that the covariant components have been calculated corrected
  constexpr double expected_g_11 = 13.0 / 9.0;
  constexpr double expected_g_22 = 4.0 / 3.0;
  constexpr double expected_g_33 = 13.0 / 9.0;
  const double expected_g_12 = -2.0 * sqrt(3.0) / 9.0;
  constexpr double expected_g_13 = -5.0 / 9.0;
  const double expected_g_23 = -2.0 * sqrt(3.0) / 9.0;

  EXPECT_TRUE(IsFieldEqual(coords.g_11(), expected_g_11));
  EXPECT_TRUE(IsFieldEqual(coords.g_22(), expected_g_22));
  EXPECT_TRUE(IsFieldEqual(coords.g_33(), expected_g_33));
  EXPECT_TRUE(IsFieldEqual(coords.g_12(), expected_g_12));
  EXPECT_TRUE(IsFieldEqual(coords.g_13(), expected_g_13));
  EXPECT_TRUE(IsFieldEqual(coords.g_23(), expected_g_23));
}

TEST_F(CoordinatesTest, CheckContravariantCalculatedFromCovariant) {

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

  //  Modify covariant components
  constexpr double g_11 = 1.0;
  constexpr double g_22 = 1.0;
  constexpr double g_33 = 1.0;
  const double g_12 = sqrt(3.0) / 4.0; // (std::sqrt only constexpr since C++26)
  constexpr double g_13 = 0.5;
  const double g_23 = sqrt(3.0) / 4.0; // (std::sqrt only constexpr since C++26)
  auto updated_metric_tensor = CovariantMetricTensor(g_11, g_22, g_33, g_12, g_13, g_23);
  coords.setCovariantMetricTensor(updated_metric_tensor);

  //  Check that the contravariant components have been calculated corrected
  constexpr double expected_g11 = 13.0 / 9.0;
  constexpr double expected_g22 = 4.0 / 3.0;
  constexpr double expected_g33 = 13.0 / 9.0;
  const double expected_g12 = -2.0 * sqrt(3.0) / 9.0;
  constexpr double expected_g13 = -5.0 / 9.0;
  const double expected_g23 = -2.0 * sqrt(3.0) / 9.0;

  EXPECT_TRUE(IsFieldEqual(coords.g11(), expected_g11));
  EXPECT_TRUE(IsFieldEqual(coords.g22(), expected_g22));
  EXPECT_TRUE(IsFieldEqual(coords.g33(), expected_g33));
  EXPECT_TRUE(IsFieldEqual(coords.g12(), expected_g12));
  EXPECT_TRUE(IsFieldEqual(coords.g13(), expected_g13));
  EXPECT_TRUE(IsFieldEqual(coords.g23(), expected_g23));
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
    auto updated_metric_tensor = CovariantMetricTensor(1.0, 2.0, 0.4, 1.0, 0.0, 0.2);
    coords.setCovariantMetricTensor(updated_metric_tensor);

    //  Get values with getter and check they have been modified as expected
    EXPECT_TRUE(IsFieldEqual(coords.g_11(), 1.0));
    EXPECT_TRUE(IsFieldEqual(coords.g_22(), 2.0));
    EXPECT_TRUE(IsFieldEqual(coords.g_33(), 0.4));
    EXPECT_TRUE(IsFieldEqual(coords.g_12(), 1.0));
    EXPECT_TRUE(IsFieldEqual(coords.g_13(), 0.0));
    EXPECT_TRUE(IsFieldEqual(coords.g_23(), 0.2));
  }
}

TEST_F(CoordinatesTest, IndexedAccessors) {

  int x = mesh->xstart;
  int y = mesh->ystart;
#if BOUT_USE_METRIC_3D
  int z = mesh->LocalNz;
#endif

  output_info.disable();
  output_warn.disable();
  Coordinates coords(mesh);
  output_warn.enable();
  output_info.enable();

  const auto& dx = coords.dx();
  const auto& dy = coords.dy();
#if BOUT_USE_METRIC_3D
  const auto& dz = coords.dz();
#endif

#if not(BOUT_USE_METRIC_3D)
  const BoutReal expected_dx = dx(x, y);
  const BoutReal expected_dy = dy(x, y);
#else
  const BoutReal expected_dx = dx(x, y, z);
  const BoutReal expected_dy = dy(x, y, z);
  const BoutReal expected_dz = dz(x, y, z);
#endif

#if not(BOUT_USE_METRIC_3D)
  const FieldMetric& actual_dx = coords.dx(x, y);
  const FieldMetric& actual_dy = coords.dy(x, y);
#else
  const Field3D& actual_dx = coords.dx(x, y, z);
  const Field3D& actual_dy = coords.dy(x, y, z);
  const Field3D& actual_dz = coords.dz(x, y, z);
#endif

  EXPECT_EQ(actual_dx, expected_dx);
  EXPECT_EQ(actual_dy, expected_dy);
#if BOUT_USE_METRIC_3D
  EXPECT_EQ(actual_dz, expected_dz);
#endif
}
