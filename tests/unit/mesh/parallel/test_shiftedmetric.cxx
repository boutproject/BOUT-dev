#include "gtest/gtest.h"

#include "fft.hxx"
#include "test_extras.hxx"
#include "bout/paralleltransform.hxx"

// The unit tests use the global mesh
using namespace bout::globals;

namespace bout {
namespace globals {
extern Mesh* mesh;
} // namespace globals
} // namespace bout

class ShiftedMetricTest : public ::testing::Test {
public:
  ShiftedMetricTest() {
    // Delete any existing mesh
    if (mesh != nullptr) {
      delete mesh;
      mesh = nullptr;
    }
    mesh = new FakeMesh(nx, ny, nz);

    // Use two y-guards to test multiple parallel slices
    mesh->ystart = 2;
    mesh->yend = mesh->LocalNy - 3;

    output_info.disable();
    mesh->createDefaultRegions();
    output_info.enable();

    zShift = Field2D{mesh};

    fillField(zShift, {{1., 2., 3., 4., 5., 6., 7.},
                       {2., 4., 6., 8., 10., 12., 14.},
                       {3., 6., 9., 12., 15., 18., 21.}});

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

    dynamic_cast<FakeMesh*>(mesh)->setCoordinates(std::make_shared<Coordinates>(
        mesh, Options::root(), Field2D{1.0}, Field2D{1.0}, BoutReal{1.0}, Field2D{1.0},
        Field2D{0.0}, Field2D{1.0}, Field2D{1.0}, Field2D{1.0}, Field2D{0.0},
        Field2D{0.0}, Field2D{0.0}, Field2D{1.0}, Field2D{1.0}, Field2D{1.0},
        Field2D{0.0}, Field2D{0.0}, Field2D{0.0}, Field2D{0.0}, Field2D{0.0}, false));
  }

  ~ShiftedMetricTest() {
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
  ShiftedMetric shifted{*mesh, zShift};

  Field3D expected{mesh};

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

  EXPECT_TRUE(IsFieldEqual(shifted.toFieldAligned(input, RGN_ALL), expected, "RGN_ALL",
                           FFTTolerance));
}

TEST_F(ShiftedMetricTest, FromFieldAligned) {
  ShiftedMetric shifted{*mesh, zShift};

  Field3D expected{mesh};

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

  // Loosen tolerance a bit due to FFTs
  EXPECT_TRUE(IsFieldEqual(shifted.fromFieldAligned(input, RGN_ALL), expected, "RGN_ALL",
                           FFTTolerance));
}

TEST_F(ShiftedMetricTest, CalcYUpDown) {
  // We don't shift in the guard cells, and the parallel slices are
  // stored offset in y, therefore we need to make new regions that we
  // can compare the expected and actual outputs over
  output_info.disable();
  mesh->addRegion3D("RGN_YUP",
                    Region<Ind3D>(0, mesh->LocalNx - 1, mesh->ystart + 1, mesh->yend + 1,
                                  0, mesh->LocalNz - 1, mesh->LocalNy, mesh->LocalNz));
  mesh->addRegion3D("RGN_YUP2",
                    Region<Ind3D>(0, mesh->LocalNx - 1, mesh->ystart + 2, mesh->yend + 2,
                                  0, mesh->LocalNz - 1, mesh->LocalNy, mesh->LocalNz));

  mesh->addRegion3D("RGN_YDOWN",
                    Region<Ind3D>(0, mesh->LocalNx - 1, mesh->ystart - 1, mesh->yend - 1,
                                  0, mesh->LocalNz - 1, mesh->LocalNy, mesh->LocalNz));
  mesh->addRegion3D("RGN_YDOWN2",
                    Region<Ind3D>(0, mesh->LocalNx - 1, mesh->ystart - 2, mesh->yend - 2,
                                  0, mesh->LocalNz - 1, mesh->LocalNy, mesh->LocalNz));
  output_info.enable();

  // Actual interesting bit here!
  ShiftedMetric shifted{*mesh, zShift};
  shifted.calcYUpDown(input);

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

  EXPECT_TRUE(IsFieldEqual(input.ynext(1), expected_up_1, "RGN_YUP", FFTTolerance));
  EXPECT_TRUE(IsFieldEqual(input.ynext(2), expected_up_2, "RGN_YUP2", FFTTolerance));
  EXPECT_TRUE(IsFieldEqual(input.ynext(-1), expected_down_1, "RGN_YDOWN", FFTTolerance));
  EXPECT_TRUE(IsFieldEqual(input.ynext(-2), expected_down2, "RGN_YDOWN2", FFTTolerance));
}
