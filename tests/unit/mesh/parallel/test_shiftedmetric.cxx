#include "bout/build_defines.hxx"

#include "gtest/gtest.h"

#include "test_extras.hxx"
#include "bout/fft.hxx"

#if BOUT_HAS_FFTW
#include "fake_mesh.hxx"

// The unit tests use the global mesh
using namespace bout::globals;

class ShiftedMetricTest : public ::testing::Test {
public:
  ShiftedMetricTest() {
    WithQuietOutput quiet_info{output_info};
    WithQuietOutput quiet_warn{output_warn};

    delete mesh;
    mesh = new FakeMesh(nx, ny, nz);
    static_cast<FakeMesh*>(mesh)->setCoordinates(nullptr);

    // Use two y-guards to test multiple parallel slices
    mesh->ystart = 2;
    mesh->yend = mesh->LocalNy - 3;

    mesh->createDefaultRegions();

    zShift = Field2D{mesh};

    fillField(zShift, {{1., 2., 3., 4., 5., 6., 7.},
                       {2., 4., 6., 8., 10., 12., 14.},
                       {3., 6., 9., 12., 15., 18., 21.}});

    static_cast<FakeMesh*>(mesh)->setCoordinates(std::make_shared<Coordinates>(
        mesh, Field2D{1.0}, Field2D{1.0}, BoutReal{1.0}, Field2D{1.0}, Field2D{0.0},
        Field2D{1.0}, Field2D{1.0}, Field2D{1.0}, Field2D{0.0}, Field2D{0.0},
        Field2D{0.0}, Field2D{1.0}, Field2D{1.0}, Field2D{1.0}, Field2D{0.0},
        Field2D{0.0}, Field2D{0.0}, Field2D{0.0}, Field2D{0.0}));
    // No call to Coordinates::geometry() needed here

    auto coords = mesh->getCoordinates();
    coords->setParallelTransform(bout::utils::make_unique<ShiftedMetric>(
        *mesh, CELL_CENTRE, zShift, coords->zlength()(0, 0)));

    Field3D input_temp{mesh};

    // input values have been slightly shuffled to ensure that input is not
    // constant in y, as this can hide bugs. Shuffling also means the rows are
    // different by something that is not a shift in the z-direction to ensure
    // that this also cannot hide bugs.
    fillField(input_temp, {{{1., 2., 3., 4., 5.},
                            {2., 1., 3., 4., 5.},
                            {1., 3., 2., 4., 5.},
                            {1., 2., 4., 3., 5.},
                            {1., 2., 3., 5., 4.},
                            {1., 2., 3., 4., 5.},
                            {2., 1., 3., 4., 5.}},

                           {{2., 1., 3., 4., 5.},
                            {1., 3., 2., 4., 5.},
                            {1., 2., 4., 3., 5.},
                            {1., 2., 3., 5., 4.},
                            {1., 2., 3., 4., 5.},
                            {2., 1., 3., 4., 5.},
                            {1., 3., 2., 4., 5.}},

                           {{1., 3., 2., 4., 5.},
                            {1., 2., 4., 3., 5.},
                            {1., 2., 3., 5., 4.},
                            {1., 2., 3., 4., 5.},
                            {2., 1., 3., 4., 5.},
                            {1., 3., 2., 4., 5.},
                            {1., 2., 4., 3., 5.}}});

    input = std::move(input_temp);
  }

  virtual ~ShiftedMetricTest() {
    delete mesh;
    mesh = nullptr;
  }

  static constexpr int nx = 3;
  static constexpr int ny = 7;
  static constexpr int nz = 5;

  Field2D zShift;
  Field3D input;
};

TEST_F(ShiftedMetricTest, ToFieldAligned) {
  Field3D expected{mesh};
  expected.setDirectionY(YDirectionType::Aligned);

  fillField(expected, {{{2., 3., 4., 5., 1.},
                        {3., 4., 5., 2., 1.},
                        {4., 5., 1., 3., 2.},
                        {5., 1., 2., 4., 3.},
                        {1., 2., 3., 5., 4.},
                        {2., 3., 4., 5., 1.},
                        {3., 4., 5., 2., 1.}},

                       {{3., 4., 5., 2., 1.},
                        {5., 1., 3., 2., 4.},
                        {2., 4., 3., 5., 1.},
                        {5., 4., 1., 2., 3.},
                        {1., 2., 3., 4., 5.},
                        {3., 4., 5., 2., 1.},
                        {5., 1., 3., 2., 4.}},

                       {{4., 5., 1., 3., 2.},
                        {2., 4., 3., 5., 1.},
                        {4., 1., 2., 3., 5.},
                        {3., 4., 5., 1., 2.},
                        {2., 1., 3., 4., 5.},
                        {4., 5., 1., 3., 2.},
                        {2., 4., 3., 5., 1.}}});

  Field3D result = toFieldAligned(input);

  EXPECT_TRUE(IsFieldEqual(result, expected, "RGN_ALL", FFTTolerance));
  EXPECT_TRUE(IsFieldEqual(fromFieldAligned(result), input));
  EXPECT_TRUE(areFieldsCompatible(result, expected));
  EXPECT_FALSE(areFieldsCompatible(result, input));
}

TEST_F(ShiftedMetricTest, FromFieldAligned) {
  // reset input.yDirectionType so that fromFieldAligned is not a null
  // operation
  input.setDirectionY(YDirectionType::Aligned);

  Field3D expected{mesh, CELL_CENTRE};
  expected.setDirectionY(YDirectionType::Standard);

  fillField(expected, {{{5., 1., 2., 3., 4.},
                        {4., 5., 2., 1., 3.},
                        {2., 4., 5., 1., 3.},
                        {2., 4., 3., 5., 1.},
                        {1., 2., 3., 5., 4.},
                        {5., 1., 2., 3., 4.},
                        {4., 5., 2., 1., 3.}},

                       {{4., 5., 2., 1., 3.},
                        {3., 2., 4., 5., 1.},
                        {5., 1., 2., 4., 3.},
                        {3., 5., 4., 1., 2.},
                        {1., 2., 3., 4., 5.},
                        {4., 5., 2., 1., 3.},
                        {3., 2., 4., 5., 1.}},

                       {{2., 4., 5., 1., 3.},
                        {5., 1., 2., 4., 3.},
                        {2., 3., 5., 4., 1.},
                        {4., 5., 1., 2., 3.},
                        {2., 1., 3., 4., 5.},
                        {2., 4., 5., 1., 3.},
                        {5., 1., 2., 4., 3.}}});

  Field3D result = fromFieldAligned(input);

  // Loosen tolerance a bit due to FFTs
  EXPECT_TRUE(IsFieldEqual(result, expected, "RGN_ALL", FFTTolerance));
  EXPECT_TRUE(IsFieldEqual(toFieldAligned(result), input, "RGN_ALL", FFTTolerance));
  EXPECT_TRUE(areFieldsCompatible(result, expected));
  EXPECT_FALSE(areFieldsCompatible(result, input));
}

TEST_F(ShiftedMetricTest, FromToFieldAligned) {
  EXPECT_TRUE(IsFieldEqual(fromFieldAligned(toFieldAligned(input)), input, "RGN_ALL",
                           FFTTolerance));
}

TEST_F(ShiftedMetricTest, ToFromFieldAligned) {
  input.setDirectionY(YDirectionType::Aligned);

  EXPECT_TRUE(IsFieldEqual(toFieldAligned(fromFieldAligned(input)), input, "RGN_ALL",
                           FFTTolerance));
}

TEST_F(ShiftedMetricTest, ToFieldAlignedFieldPerp) {
  Field3D expected{mesh};
  expected.setDirectionY(YDirectionType::Aligned);

  fillField(expected, {{{2., 3., 4., 5., 1.},
                        {3., 4., 5., 2., 1.},
                        {4., 5., 1., 3., 2.},
                        {5., 1., 2., 4., 3.},
                        {1., 2., 3., 5., 4.},
                        {2., 3., 4., 5., 1.},
                        {3., 4., 5., 2., 1.}},

                       {{3., 4., 5., 2., 1.},
                        {5., 1., 3., 2., 4.},
                        {2., 4., 3., 5., 1.},
                        {5., 4., 1., 2., 3.},
                        {1., 2., 3., 4., 5.},
                        {3., 4., 5., 2., 1.},
                        {5., 1., 3., 2., 4.}},

                       {{4., 5., 1., 3., 2.},
                        {2., 4., 3., 5., 1.},
                        {4., 1., 2., 3., 5.},
                        {3., 4., 5., 1., 2.},
                        {2., 1., 3., 4., 5.},
                        {4., 5., 1., 3., 2.},
                        {2., 4., 3., 5., 1.}}});

  FieldPerp result = toFieldAligned(sliceXZ(input, 3), "RGN_NOX");

  // Note that the region argument does not do anything for FieldPerp, as
  // FieldPerp does not have a getRegion2D() method. Values are never set in
  // the x-guard or x-boundary cells
  EXPECT_TRUE(IsFieldEqual(result, sliceXZ(expected, 3), "RGN_NOBNDRY", FFTTolerance));
  EXPECT_TRUE(IsFieldEqual(fromFieldAligned(result, "RGN_NOX"), sliceXZ(input, 3),
                           "RGN_NOBNDRY", FFTTolerance));
  EXPECT_TRUE(areFieldsCompatible(result, sliceXZ(expected, 3)));
  EXPECT_FALSE(areFieldsCompatible(result, sliceXZ(input, 3)));
}

TEST_F(ShiftedMetricTest, FromFieldAlignedFieldPerp) {
  // reset input.yDirectionType so that fromFieldAligned is not a null
  // operation
  input.setDirectionY(YDirectionType::Aligned);

  Field3D expected{mesh, CELL_CENTRE};
  expected.setDirectionY(YDirectionType::Standard);

  fillField(expected, {{{5., 1., 2., 3., 4.},
                        {4., 5., 2., 1., 3.},
                        {2., 4., 5., 1., 3.},
                        {2., 4., 3., 5., 1.},
                        {1., 2., 3., 5., 4.},
                        {5., 1., 2., 3., 4.},
                        {4., 5., 2., 1., 3.}},

                       {{4., 5., 2., 1., 3.},
                        {3., 2., 4., 5., 1.},
                        {5., 1., 2., 4., 3.},
                        {3., 5., 4., 1., 2.},
                        {1., 2., 3., 4., 5.},
                        {4., 5., 2., 1., 3.},
                        {3., 2., 4., 5., 1.}},

                       {{2., 4., 5., 1., 3.},
                        {5., 1., 2., 4., 3.},
                        {2., 3., 5., 4., 1.},
                        {4., 5., 1., 2., 3.},
                        {2., 1., 3., 4., 5.},
                        {2., 4., 5., 1., 3.},
                        {5., 1., 2., 4., 3.}}});

  FieldPerp result = fromFieldAligned(sliceXZ(input, 4), "RGN_NOX");

  // Note that the region argument does not do anything for FieldPerp, as
  // FieldPerp does not have a getRegion2D() method. Values are never set in
  // the x-guard or x-boundary cells
  EXPECT_TRUE(IsFieldEqual(result, sliceXZ(expected, 4), "RGN_NOBNDRY", FFTTolerance));
  EXPECT_TRUE(IsFieldEqual(toFieldAligned(result, "RGN_NOX"), sliceXZ(input, 4),
                           "RGN_NOBNDRY", FFTTolerance));
  EXPECT_TRUE(areFieldsCompatible(result, sliceXZ(expected, 4)));
  EXPECT_FALSE(areFieldsCompatible(result, sliceXZ(input, 4)));
}

TEST_F(ShiftedMetricTest, FromToFieldAlignedFieldPerp) {
  // Note that the region argument does not do anything for FieldPerp, as
  // FieldPerp does not have a getRegion2D() method. Values are never set in
  // the x-guard or x-boundary cells
  EXPECT_TRUE(IsFieldEqual(
      fromFieldAligned(toFieldAligned(sliceXZ(input, 2), "RGN_NOX"), "RGN_NOX"),
      sliceXZ(input, 2), "RGN_NOBNDRY", FFTTolerance));
}

TEST_F(ShiftedMetricTest, ToFromFieldAlignedFieldPerp) {
  // Note that the region argument does not do anything for FieldPerp, as
  // FieldPerp does not have a getRegion2D() method. Values are never set in
  // the x-guard or x-boundary cells
  input.setDirectionY(YDirectionType::Aligned);

  EXPECT_TRUE(IsFieldEqual(
      toFieldAligned(fromFieldAligned(sliceXZ(input, 6), "RGN_NOX"), "RGN_NOX"),
      sliceXZ(input, 6), "RGN_NOBNDRY", FFTTolerance));
}

TEST_F(ShiftedMetricTest, CalcParallelSlices) {
  // Actual interesting bit here!
  input.getCoordinates()->getParallelTransform().calcParallelSlices(input);
  // Expected output values

  Field3D expected_up_1{mesh};

  // Note: here zeroes are for values we don't expect to read
  fillField(expected_up_1, {{{0., 0., 0., 0., 0.},
                             {0., 0., 0., 0., 0.},
                             {0., 0., 0., 0., 0.},
                             {2., 4., 3., 5., 1.},
                             {2., 3., 5., 4., 1.},
                             {2., 3., 4., 5., 1.},
                             {0., 0., 0., 0., 0.}},

                            {{0., 0., 0., 0., 0.},
                             {0., 0., 0., 0., 0.},
                             {0., 0., 0., 0., 0.},
                             {3., 5., 4., 1., 2.},
                             {3., 4., 5., 1., 2.},
                             {3., 4., 5., 2., 1.},
                             {0., 0., 0., 0., 0.}},

                            {{0., 0., 0., 0., 0.},
                             {0., 0., 0., 0., 0.},
                             {0., 0., 0., 0., 0.},
                             {4., 5., 1., 2., 3.},
                             {4., 5., 2., 1., 3.},
                             {4., 5., 1., 3., 2.},
                             {0., 0., 0., 0., 0.}}});

  Field3D expected_up_2{mesh};

  fillField(expected_up_2, {{{0., 0., 0., 0., 0.},
                             {0., 0., 0., 0., 0.},
                             {0., 0., 0., 0., 0.},
                             {0., 0., 0., 0., 0.},
                             {3., 5., 4., 1., 2.},
                             {3., 4., 5., 1., 2.},
                             {3., 4., 5., 2., 1.}},

                            {{0., 0., 0., 0., 0.},
                             {0., 0., 0., 0., 0.},
                             {0., 0., 0., 0., 0.},
                             {0., 0., 0., 0., 0.},
                             {5., 1., 2., 3., 4.},
                             {5., 2., 1., 3., 4.},
                             {5., 1., 3., 2., 4.}},

                            {{0., 0., 0., 0., 0.},
                             {0., 0., 0., 0., 0.},
                             {0., 0., 0., 0., 0.},
                             {0., 0., 0., 0., 0.},
                             {1., 3., 4., 5., 2.},
                             {3., 2., 4., 5., 1.},
                             {2., 4., 3., 5., 1.}}});

  Field3D expected_down_1{mesh};

  fillField(expected_down_1, {{{0., 0., 0., 0., 0.},
                               {5., 2., 1., 3., 4.},
                               {5., 1., 3., 2., 4.},
                               {5., 1., 2., 4., 3.},
                               {0., 0., 0., 0., 0.},
                               {0., 0., 0., 0., 0.},
                               {0., 0., 0., 0., 0.}},

                              {{0., 0., 0., 0., 0.},
                               {4., 5., 1., 3., 2.},
                               {3., 5., 1., 2., 4.},
                               {5., 4., 1., 2., 3.},
                               {0., 0., 0., 0., 0.},
                               {0., 0., 0., 0., 0.},
                               {0., 0., 0., 0., 0.}},

                              {{0., 0., 0., 0., 0.},
                               {4., 3., 5., 1., 2.},
                               {3., 5., 4., 1., 2.},
                               {3., 4., 5., 1., 2.},
                               {0., 0., 0., 0., 0.},
                               {0., 0., 0., 0., 0.},
                               {0., 0., 0., 0., 0.}}});

  Field3D expected_down2{mesh};

  fillField(expected_down2, {{{4., 5., 1., 2., 3.},
                              {4., 5., 2., 1., 3.},
                              {4., 5., 1., 3., 2.},
                              {0., 0., 0., 0., 0.},
                              {0., 0., 0., 0., 0.},
                              {0., 0., 0., 0., 0.},
                              {0., 0., 0., 0., 0.}},

                             {{1., 3., 4., 5., 2.},
                              {3., 2., 4., 5., 1.},
                              {2., 4., 3., 5., 1.},
                              {0., 0., 0., 0., 0.},
                              {0., 0., 0., 0., 0.},
                              {0., 0., 0., 0., 0.},
                              {0., 0., 0., 0., 0.}},

                             {{5., 1., 3., 2., 4.},
                              {5., 1., 2., 4., 3.},
                              {4., 1., 2., 3., 5.},
                              {0., 0., 0., 0., 0.},
                              {0., 0., 0., 0., 0.},
                              {0., 0., 0., 0., 0.},
                              {0., 0., 0., 0., 0.}}});

  EXPECT_TRUE(IsFieldEqual(input.ynext(1), expected_up_1, "RGN_YPAR_+1", FFTTolerance));
  EXPECT_TRUE(IsFieldEqual(input.ynext(2), expected_up_2, "RGN_YPAR_+2", FFTTolerance));
  EXPECT_TRUE(
      IsFieldEqual(input.ynext(-1), expected_down_1, "RGN_YPAR_-1", FFTTolerance));
  EXPECT_TRUE(IsFieldEqual(input.ynext(-2), expected_down2, "RGN_YPAR_-2", FFTTolerance));
}
#endif
