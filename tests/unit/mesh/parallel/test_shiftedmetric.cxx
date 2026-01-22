#include "bout/build_defines.hxx"

#include "gtest/gtest.h"
#include <type_traits>

#include "../../src/mesh/parallel/shiftedmetricinterp.hxx"
#include "test_extras.hxx"
#include "bout/fft.hxx"
#include "bout/field2d.hxx"
#include "bout/globals.hxx"
#include "bout/paralleltransform.hxx"
#include "bout/traits.hxx"

#if BOUT_HAS_FFTW
#include "fake_mesh.hxx"

// The unit tests use the global mesh
using namespace bout::globals;

template <class Transform>
class ShiftedMetricTest : public ::testing::Test {
public:
  ShiftedMetricTest() {
    WithQuietOutput quiet_info{output_info};
    WithQuietOutput quiet_warn{output_warn};

    test_mesh.setCoordinates(nullptr);
    mesh = &test_mesh;

    // Use two y-guards to test multiple parallel slices
    mesh->ystart = 2;
    mesh->yend = mesh->LocalNy - 3;

    mesh->createDefaultRegions();

    zShift = Field2D{mesh};

    fillField(zShift, {{1., 2., 3., 4., 5., 6., 7.},
                       {2., 4., 6., 8., 10., 12., 14.},
                       {3., 6., 9., 12., 15., 18., 21.}});

    test_mesh.setCoordinates(std::make_shared<Coordinates>(
        mesh, Field2D{1.0}, Field2D{1.0}, BoutReal{1.0}, Field2D{1.0}, Field2D{0.0},
        Field2D{1.0}, Field2D{1.0}, Field2D{1.0}, Field2D{0.0}, Field2D{0.0},
        Field2D{0.0}, Field2D{1.0}, Field2D{1.0}, Field2D{1.0}, Field2D{0.0},
        Field2D{0.0}, Field2D{0.0}, Field2D{0.0}, Field2D{0.0}));
    // No call to Coordinates::geometry() needed here

    auto* coords = mesh->getCoordinates();
    coords->setParallelTransform(bout::utils::make_unique<Transform>(
        *mesh, CELL_CENTRE, zShift, coords->zlength()(0, 0)));

    Field3D input_temp{mesh};

    // input values have been slightly shuffled to ensure that input is not
    // constant in y, as this can hide bugs. Shuffling also means the rows are
    // different by something that is not a shift in the z-direction to ensure
    // that this also cannot hide bugs.
    fillField(input_temp, {{{-77., 1., 2., 3., 4., 5., -99.},
                            {-77., 2., 1., 3., 4., 5., -99.},
                            {-77., 1., 3., 2., 4., 5., -99.},
                            {-77., 1., 2., 4., 3., 5., -99.},
                            {-77., 1., 2., 3., 5., 4., -99.},
                            {-77., 1., 2., 3., 4., 5., -99.},
                            {-77., 2., 1., 3., 4., 5., -99.}},

                           {{-77., 2., 1., 3., 4., 5., -99.},
                            {-77., 1., 3., 2., 4., 5., -99.},
                            {-77., 1., 2., 4., 3., 5., -99.},
                            {-77., 1., 2., 3., 5., 4., -99.},
                            {-77., 1., 2., 3., 4., 5., -99.},
                            {-77., 2., 1., 3., 4., 5., -99.},
                            {-77., 1., 3., 2., 4., 5., -99.}},

                           {{-77., 1., 3., 2., 4., 5., -99.},
                            {-77., 1., 2., 4., 3., 5., -99.},
                            {-77., 1., 2., 3., 5., 4., -99.},
                            {-77., 1., 2., 3., 4., 5., -99.},
                            {-77., 2., 1., 3., 4., 5., -99.},
                            {-77., 1., 3., 2., 4., 5., -99.},
                            {-77., 1., 2., 4., 3., 5., -99.}}});

    original = copy(input_temp);
    // Make sure Z guards are copied
    mesh->communicate_no_slices(input_temp);
    input = std::move(input_temp);
  }

  ShiftedMetricTest(const ShiftedMetricTest&) = delete;
  ShiftedMetricTest(ShiftedMetricTest&&) = delete;
  ShiftedMetricTest& operator=(const ShiftedMetricTest&) = delete;
  ShiftedMetricTest& operator=(ShiftedMetricTest&&) = delete;
  ~ShiftedMetricTest() override { bout::globals::mesh = nullptr; }

  static constexpr int nx = 3;
  static constexpr int ny = 7;
  static constexpr int nz = 7;

  FakeMesh test_mesh{nx, ny, nz};

  Field2D zShift;
  Field3D input;
  Field3D original;
};

using ShiftedMetricTypes = ::testing::Types<ShiftedMetric, ShiftedMetricInterp>;
TYPED_TEST_SUITE(ShiftedMetricTest, ShiftedMetricTypes);

TYPED_TEST(ShiftedMetricTest, ToFieldAligned) {
  Field3D expected{mesh};
  expected.setDirectionY(YDirectionType::Aligned);

  fillField(expected, {{{-99., 2., 3., 4., 5., 1., -99.},
                        {-99., 3., 4., 5., 2., 1., -99.},
                        {-99., 4., 5., 1., 3., 2., -99.},
                        {-99., 5., 1., 2., 4., 3., -99.},
                        {-99., 1., 2., 3., 5., 4., -99.},
                        {-99., 2., 3., 4., 5., 1., -99.},
                        {-99., 3., 4., 5., 2., 1., -99.}},

                       {{-99., 3., 4., 5., 2., 1., -99.},
                        {-99., 5., 1., 3., 2., 4., -99.},
                        {-99., 2., 4., 3., 5., 1., -99.},
                        {-99., 5., 4., 1., 2., 3., -99.},
                        {-99., 1., 2., 3., 4., 5., -99.},
                        {-99., 3., 4., 5., 2., 1., -99.},
                        {-99., 5., 1., 3., 2., 4., -99.}},

                       {{-99., 4., 5., 1., 3., 2., -99.},
                        {-99., 2., 4., 3., 5., 1., -99.},
                        {-99., 4., 1., 2., 3., 5., -99.},
                        {-99., 3., 4., 5., 1., 2., -99.},
                        {-99., 2., 1., 3., 4., 5., -99.},
                        {-99., 4., 5., 1., 3., 2., -99.},
                        {-99., 2., 4., 3., 5., 1., -99.}}});

  Field3D result = toFieldAligned(this->input);

  EXPECT_TRUE(IsFieldEqual(result, expected, "RGN_NOZ", FFTTolerance));
  EXPECT_TRUE(IsFieldEqual(fromFieldAligned(result), this->input, "RGN_NOZ"));
  EXPECT_TRUE(areFieldsCompatible(result, expected));
  EXPECT_FALSE(areFieldsCompatible(result, this->input));
}

TYPED_TEST(ShiftedMetricTest, FromFieldAligned) {
  // reset input.yDirectionType so that fromFieldAligned is not a null
  // operation
  this->input.setDirectionY(YDirectionType::Aligned);

  Field3D expected{mesh, CELL_CENTRE};
  expected.setDirectionY(YDirectionType::Standard);

  fillField(expected, {{{-99., 5., 1., 2., 3., 4., -99.},
                        {-99., 4., 5., 2., 1., 3., -99.},
                        {-99., 2., 4., 5., 1., 3., -99.},
                        {-99., 2., 4., 3., 5., 1., -99.},
                        {-99., 1., 2., 3., 5., 4., -99.},
                        {-99., 5., 1., 2., 3., 4., -99.},
                        {-99., 4., 5., 2., 1., 3., -99.}},

                       {{-99., 4., 5., 2., 1., 3., -99.},
                        {-99., 3., 2., 4., 5., 1., -99.},
                        {-99., 5., 1., 2., 4., 3., -99.},
                        {-99., 3., 5., 4., 1., 2., -99.},
                        {-99., 1., 2., 3., 4., 5., -99.},
                        {-99., 4., 5., 2., 1., 3., -99.},
                        {-99., 3., 2., 4., 5., 1., -99.}},

                       {{-99., 2., 4., 5., 1., 3., -99.},
                        {-99., 5., 1., 2., 4., 3., -99.},
                        {-99., 2., 3., 5., 4., 1., -99.},
                        {-99., 4., 5., 1., 2., 3., -99.},
                        {-99., 2., 1., 3., 4., 5., -99.},
                        {-99., 2., 4., 5., 1., 3., -99.},
                        {-99., 5., 1., 2., 4., 3., -99.}}});

  Field3D result = fromFieldAligned(this->input);

  // Loosen tolerance a bit due to FFTs
  EXPECT_TRUE(IsFieldEqual(result, expected, "RGN_NOZ", FFTTolerance));
  EXPECT_TRUE(IsFieldEqual(toFieldAligned(result), this->input, "RGN_NOZ", FFTTolerance));
  EXPECT_TRUE(areFieldsCompatible(result, expected));
  EXPECT_FALSE(areFieldsCompatible(result, this->input));
}

TYPED_TEST(ShiftedMetricTest, FromToFieldAligned) {
  EXPECT_TRUE(IsFieldEqual(fromFieldAligned(toFieldAligned(this->input)), this->input,
                           "RGN_NOZ", FFTTolerance));
}

TYPED_TEST(ShiftedMetricTest, ToFromFieldAligned) {
  this->input.setDirectionY(YDirectionType::Aligned);

  EXPECT_TRUE(IsFieldEqual(toFieldAligned(fromFieldAligned(this->input)), this->input,
                           "RGN_NOZ", FFTTolerance));
}

TYPED_TEST(ShiftedMetricTest, ToFieldAlignedFieldPerp) {
  if constexpr (std::is_same_v<TypeParam, ShiftedMetricInterp>) {
    GTEST_SKIP_("Not implemented yet");
  }

  Field3D expected{mesh};
  expected.setDirectionY(YDirectionType::Aligned);

  fillField(expected, {{{-99., 2., 3., 4., 5., 1., -99.},
                        {-99., 3., 4., 5., 2., 1., -99.},
                        {-99., 4., 5., 1., 3., 2., -99.},
                        {-99., 5., 1., 2., 4., 3., -99.},
                        {-99., 1., 2., 3., 5., 4., -99.},
                        {-99., 2., 3., 4., 5., 1., -99.},
                        {-99., 3., 4., 5., 2., 1., -99.}},

                       {{-99., 3., 4., 5., 2., 1., -99.},
                        {-99., 5., 1., 3., 2., 4., -99.},
                        {-99., 2., 4., 3., 5., 1., -99.},
                        {-99., 5., 4., 1., 2., 3., -99.},
                        {-99., 1., 2., 3., 4., 5., -99.},
                        {-99., 3., 4., 5., 2., 1., -99.},
                        {-99., 5., 1., 3., 2., 4., -99.}},

                       {{-99., 4., 5., 1., 3., 2., -99.},
                        {-99., 2., 4., 3., 5., 1., -99.},
                        {-99., 4., 1., 2., 3., 5., -99.},
                        {-99., 3., 4., 5., 1., 2., -99.},
                        {-99., 2., 1., 3., 4., 5., -99.},
                        {-99., 4., 5., 1., 3., 2., -99.},
                        {-99., 2., 4., 3., 5., 1., -99.}}});

  FieldPerp result = toFieldAligned(sliceXZ(this->input, 3), "RGN_NOX");

  // Note that the region argument does not do anything for FieldPerp, as
  // FieldPerp does not have a getRegion2D() method. Values are never set in
  // the x-guard or x-boundary cells
  EXPECT_TRUE(IsFieldEqual(result, sliceXZ(expected, 3), "RGN_NOBNDRY", FFTTolerance));
  EXPECT_TRUE(IsFieldEqual(fromFieldAligned(result, "RGN_NOX"), sliceXZ(this->input, 3),
                           "RGN_NOBNDRY", FFTTolerance));
  EXPECT_TRUE(areFieldsCompatible(result, sliceXZ(expected, 3)));
  EXPECT_FALSE(areFieldsCompatible(result, sliceXZ(this->input, 3)));
}

TYPED_TEST(ShiftedMetricTest, FromFieldAlignedFieldPerp) {
  if constexpr (std::is_same_v<TypeParam, ShiftedMetricInterp>) {
    GTEST_SKIP_("Not implemented yet");
  }

  // reset this->input.yDirectionType so that fromFieldAligned is not a null
  // operation
  this->input.setDirectionY(YDirectionType::Aligned);

  Field3D expected{mesh, CELL_CENTRE};
  expected.setDirectionY(YDirectionType::Standard);

  fillField(expected, {{{-99., 5., 1., 2., 3., 4., -99.},
                        {-99., 4., 5., 2., 1., 3., -99.},
                        {-99., 2., 4., 5., 1., 3., -99.},
                        {-99., 2., 4., 3., 5., 1., -99.},
                        {-99., 1., 2., 3., 5., 4., -99.},
                        {-99., 5., 1., 2., 3., 4., -99.},
                        {-99., 4., 5., 2., 1., 3., -99.}},

                       {{-99., 4., 5., 2., 1., 3., -99.},
                        {-99., 3., 2., 4., 5., 1., -99.},
                        {-99., 5., 1., 2., 4., 3., -99.},
                        {-99., 3., 5., 4., 1., 2., -99.},
                        {-99., 1., 2., 3., 4., 5., -99.},
                        {-99., 4., 5., 2., 1., 3., -99.},
                        {-99., 3., 2., 4., 5., 1., -99.}},

                       {{-99., 2., 4., 5., 1., 3., -99.},
                        {-99., 5., 1., 2., 4., 3., -99.},
                        {-99., 2., 3., 5., 4., 1., -99.},
                        {-99., 4., 5., 1., 2., 3., -99.},
                        {-99., 2., 1., 3., 4., 5., -99.},
                        {-99., 2., 4., 5., 1., 3., -99.},
                        {-99., 5., 1., 2., 4., 3., -99.}}});

  FieldPerp result = fromFieldAligned(sliceXZ(this->input, 4), "RGN_NOX");

  // Note that the region argument does not do anything for FieldPerp, as
  // FieldPerp does not have a getRegion2D() method. Values are never set in
  // the x-guard or x-boundary cells
  EXPECT_TRUE(IsFieldEqual(result, sliceXZ(expected, 4), "RGN_NOBNDRY", FFTTolerance));
  EXPECT_TRUE(IsFieldEqual(toFieldAligned(result, "RGN_NOX"), sliceXZ(this->input, 4),
                           "RGN_NOBNDRY", FFTTolerance));
  EXPECT_TRUE(areFieldsCompatible(result, sliceXZ(expected, 4)));
  EXPECT_FALSE(areFieldsCompatible(result, sliceXZ(this->input, 4)));
}

TYPED_TEST(ShiftedMetricTest, FromToFieldAlignedFieldPerp) {
  if constexpr (std::is_same_v<TypeParam, ShiftedMetricInterp>) {
    GTEST_SKIP_("Not implemented yet");
  }

  // Note that the region argument does not do anything for FieldPerp, as
  // FieldPerp does not have a getRegion2D() method. Values are never set in
  // the x-guard or x-boundary cells
  EXPECT_TRUE(IsFieldEqual(
      fromFieldAligned(toFieldAligned(sliceXZ(this->input, 2), "RGN_NOX"), "RGN_NOX"),
      sliceXZ(this->input, 2), "RGN_NOBNDRY", FFTTolerance));
}

TYPED_TEST(ShiftedMetricTest, ToFromFieldAlignedFieldPerp) {
  if constexpr (std::is_same_v<TypeParam, ShiftedMetricInterp>) {
    GTEST_SKIP_("Not implemented yet");
  }

  // Note that the region argument does not do anything for FieldPerp, as
  // FieldPerp does not have a getRegion2D() method. Values are never set in
  // the x-guard or x-boundary cells
  this->input.setDirectionY(YDirectionType::Aligned);

  EXPECT_TRUE(IsFieldEqual(
      toFieldAligned(fromFieldAligned(sliceXZ(this->input, 6), "RGN_NOX"), "RGN_NOX"),
      sliceXZ(this->input, 6), "RGN_NOBNDRY", FFTTolerance));
}

TYPED_TEST(ShiftedMetricTest, CalcParallelSlices) {
  WithQuietOutput quiet_info{output_info};

  // We don't shift in the guard cells, and the parallel slices are
  // stored offset in y, therefore we need to make new regions that we
  // can compare the expected and actual outputs over

  constexpr bool is_shifted_interp = std::is_same_v<TypeParam, ShiftedMetricInterp>;

  // ShiftedMetricInterp doesn't seem to store interpolated values in _any_ of the
  // guards, so we can only check the interior X point
  const int xstart = is_shifted_interp ? mesh->xstart : 0;
  const int xend = is_shifted_interp ? mesh->xend : mesh->LocalNx - 1;

  // It also means we can't check in the _y_ guards either.
  // TODO(peter): Is this a bug?
  const int yup_1_end = is_shifted_interp ? mesh->yend : mesh->yend + 1;
  const int yup_2_end = is_shifted_interp ? mesh->yend : mesh->yend + 2;

  const int ydown_1_start = is_shifted_interp ? mesh->ystart : mesh->ystart - 1;
  const int ydown_2_start = is_shifted_interp ? mesh->ystart : mesh->ystart - 2;

  mesh->addRegion3D("RGN_YUP",
                    Region<Ind3D>(xstart, xend, mesh->ystart + 1, yup_1_end, mesh->zstart,
                                  mesh->zend, mesh->LocalNy, mesh->LocalNz));
  mesh->addRegion3D("RGN_YUP2",
                    Region<Ind3D>(xstart, xend, mesh->ystart + 2, yup_2_end, mesh->zstart,
                                  mesh->zend, mesh->LocalNy, mesh->LocalNz));

  mesh->addRegion3D("RGN_YDOWN", Region<Ind3D>(xstart, xend, ydown_1_start,
                                               mesh->yend - 1, mesh->zstart, mesh->zend,
                                               mesh->LocalNy, mesh->LocalNz));
  mesh->addRegion3D("RGN_YDOWN2", Region<Ind3D>(xstart, xend, ydown_2_start,
                                                mesh->yend - 2, mesh->zstart, mesh->zend,
                                                mesh->LocalNy, mesh->LocalNz));
  output_info.enable();

  // Actual interesting bit here!
  this->input.getCoordinates()->getParallelTransform().calcParallelSlices(this->input);
  // Expected output values

  Field3D expected_up_1{mesh};

  // Note: here zeroes are for values we don't expect to read
  fillField(expected_up_1, {{{0., 0., 0., 0., 0., 0., 0.},
                             {0., 0., 0., 0., 0., 0., 0.},
                             {0., 0., 0., 0., 0., 0., 0.},
                             {0., 2., 4., 3., 5., 1., 0.},
                             {0., 2., 3., 5., 4., 1., 0.},
                             {0., 2., 3., 4., 5., 1., 0.},
                             {0., 0., 0., 0., 0., 0., 0.}},

                            {{0., 0., 0., 0., 0., 0., 0.},
                             {0., 0., 0., 0., 0., 0., 0.},
                             {0., 0., 0., 0., 0., 0., 0.},
                             {0., 3., 5., 4., 1., 2., 0.},
                             {0., 3., 4., 5., 1., 2., 0.},
                             {0., 3., 4., 5., 2., 1., 0.},
                             {0., 0., 0., 0., 0., 0., 0.}},

                            {{0., 0., 0., 0., 0., 0., 0.},
                             {0., 0., 0., 0., 0., 0., 0.},
                             {0., 0., 0., 0., 0., 0., 0.},
                             {0., 4., 5., 1., 2., 3., 0.},
                             {0., 4., 5., 2., 1., 3., 0.},
                             {0., 4., 5., 1., 3., 2., 0.},
                             {0., 0., 0., 0., 0., 0., 0.}}});

  Field3D expected_up_2{mesh};

  fillField(expected_up_2, {{{0., 0., 0., 0., 0., 0., 0.},
                             {0., 0., 0., 0., 0., 0., 0.},
                             {0., 0., 0., 0., 0., 0., 0.},
                             {0., 0., 0., 0., 0., 0., 0.},
                             {0., 3., 5., 4., 1., 2., 0.},
                             {0., 3., 4., 5., 1., 2., 0.},
                             {0., 3., 4., 5., 2., 1., 0.}},

                            {{0., 0., 0., 0., 0., 0., 0.},
                             {0., 0., 0., 0., 0., 0., 0.},
                             {0., 0., 0., 0., 0., 0., 0.},
                             {0., 0., 0., 0., 0., 0., 0.},
                             {0., 5., 1., 2., 3., 4., 0.},
                             {0., 5., 2., 1., 3., 4., 0.},
                             {0., 5., 1., 3., 2., 4., 0.}},

                            {{0., 0., 0., 0., 0., 0., 0.},
                             {0., 0., 0., 0., 0., 0., 0.},
                             {0., 0., 0., 0., 0., 0., 0.},
                             {0., 0., 0., 0., 0., 0., 0.},
                             {0., 1., 3., 4., 5., 2., 0.},
                             {0., 3., 2., 4., 5., 1., 0.},
                             {0., 2., 4., 3., 5., 1., 0.}}});

  Field3D expected_down_1{mesh};

  fillField(expected_down_1, {{{0., 0., 0., 0., 0., 0., 0.},
                               {0., 5., 2., 1., 3., 4., 0.},
                               {0., 5., 1., 3., 2., 4., 0.},
                               {0., 5., 1., 2., 4., 3., 0.},
                               {0., 0., 0., 0., 0., 0., 0.},
                               {0., 0., 0., 0., 0., 0., 0.},
                               {0., 0., 0., 0., 0., 0., 0.}},

                              {{0., 0., 0., 0., 0., 0., 0.},
                               {0., 4., 5., 1., 3., 2., 0.},
                               {0., 3., 5., 1., 2., 4., 0.},
                               {0., 5., 4., 1., 2., 3., 0.},
                               {0., 0., 0., 0., 0., 0., 0.},
                               {0., 0., 0., 0., 0., 0., 0.},
                               {0., 0., 0., 0., 0., 0., 0.}},

                              {{0., 0., 0., 0., 0., 0., 0.},
                               {0., 4., 3., 5., 1., 2., 0.},
                               {0., 3., 5., 4., 1., 2., 0.},
                               {0., 3., 4., 5., 1., 2., 0.},
                               {0., 0., 0., 0., 0., 0., 0.},
                               {0., 0., 0., 0., 0., 0., 0.},
                               {0., 0., 0., 0., 0., 0., 0.}}});

  Field3D expected_down_2{mesh};

  fillField(expected_down_2, {{{0., 4., 5., 1., 2., 3., 0.},
                               {0., 4., 5., 2., 1., 3., 0.},
                               {0., 4., 5., 1., 3., 2., 0.},
                               {0., 0., 0., 0., 0., 0., 0.},
                               {0., 0., 0., 0., 0., 0., 0.},
                               {0., 0., 0., 0., 0., 0., 0.},
                               {0., 0., 0., 0., 0., 0., 0.}},

                              {{0., 1., 3., 4., 5., 2., 0.},
                               {0., 3., 2., 4., 5., 1., 0.},
                               {0., 2., 4., 3., 5., 1., 0.},
                               {0., 0., 0., 0., 0., 0., 0.},
                               {0., 0., 0., 0., 0., 0., 0.},
                               {0., 0., 0., 0., 0., 0., 0.},
                               {0., 0., 0., 0., 0., 0., 0.}},

                              {{0., 5., 1., 3., 2., 4., 0.},
                               {0., 5., 1., 2., 4., 3., 0.},
                               {0., 4., 1., 2., 3., 5., 0.},
                               {0., 0., 0., 0., 0., 0., 0.},
                               {0., 0., 0., 0., 0., 0., 0.},
                               {0., 0., 0., 0., 0., 0., 0.},
                               {0., 0., 0., 0., 0., 0., 0.}}});

  EXPECT_TRUE(IsFieldEqual(this->input.ynext(1), expected_up_1, "RGN_YUP", FFTTolerance));
  EXPECT_TRUE(
      IsFieldEqual(this->input.ynext(2), expected_up_2, "RGN_YUP2", FFTTolerance));
  EXPECT_TRUE(
      IsFieldEqual(this->input.ynext(-1), expected_down_1, "RGN_YDOWN", FFTTolerance));
  EXPECT_TRUE(
      IsFieldEqual(this->input.ynext(-2), expected_down_2, "RGN_YDOWN2", FFTTolerance));
}
#endif
