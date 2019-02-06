#include "gtest/gtest.h"

#include "bout/paralleltransform.hxx"
#include "test_extras.hxx"

// The unit tests use the global mesh
using namespace bout::globals;

namespace bout{
namespace globals{
extern Mesh *mesh;
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
    output_info.disable();
    mesh->createDefaultRegions();
    output_info.enable();

    zShift = Field2D{mesh};

    fillField(zShift, {{1., 2., 3., 4., 5.}, {1., 2., 3., 4., 5.}, {1., 2., 3., 4., 5.}});

    dynamic_cast<FakeMesh*>(mesh)->setCoordinates(std::make_shared<Coordinates>(
        mesh, Field2D{1.0}, Field2D{1.0}, BoutReal{1.0}, Field2D{1.0}, Field2D{0.0},
        Field2D{1.0}, Field2D{1.0}, Field2D{1.0}, Field2D{0.0}, Field2D{0.0},
        Field2D{0.0}, Field2D{1.0}, Field2D{1.0}, Field2D{1.0}, Field2D{0.0},
        Field2D{0.0}, Field2D{0.0}, Field2D{0.0}, Field2D{0.0}, false));
  }

  ~ShiftedMetricTest() {
    delete mesh;
    mesh = nullptr;
  }

  static constexpr int nx = 3;
  static constexpr int ny = 5;
  static constexpr int nz = 7;

  Field2D zShift;
};

TEST_F(ShiftedMetricTest, ToFieldAligned) {
  ShiftedMetric shifted{*mesh, zShift};

  Field3D input{mesh};

  fillField(input, {{{1., 2., 3., 4., 5., 6., 7.},
                     {1., 2., 3., 4., 5., 6., 7.},
                     {1., 2., 3., 4., 5., 6., 7.},
                     {1., 2., 3., 4., 5., 6., 7.},
                     {1., 2., 3., 4., 5., 6., 7.}},

                    {{1., 2., 3., 4., 5., 6., 7.},
                     {1., 2., 3., 4., 5., 6., 7.},
                     {1., 2., 3., 4., 5., 6., 7.},
                     {1., 2., 3., 4., 5., 6., 7.},
                     {1., 2., 3., 4., 5., 6., 7.}},

                    {{1., 2., 3., 4., 5., 6., 7.},
                     {1., 2., 3., 4., 5., 6., 7.},
                     {1., 2., 3., 4., 5., 6., 7.},
                     {1., 2., 3., 4., 5., 6., 7.},
                     {1., 2., 3., 4., 5., 6., 7.}}});

  Field3D expected{mesh};

  fillField(expected, {{{2., 3., 4., 5., 6., 7., 1.},
                        {3., 4., 5., 6., 7., 1., 2.},
                        {4., 5., 6., 7., 1., 2., 3.},
                        {5., 6., 7., 1., 2., 3., 4.},
                        {6., 7., 1., 2., 3., 4., 5.}},

                       {{2., 3., 4., 5., 6., 7., 1.},
                        {3., 4., 5., 6., 7., 1., 2.},
                        {4., 5., 6., 7., 1., 2., 3.},
                        {5., 6., 7., 1., 2., 3., 4.},
                        {6., 7., 1., 2., 3., 4., 5.}},

                       {{2., 3., 4., 5., 6., 7., 1.},
                        {3., 4., 5., 6., 7., 1., 2.},
                        {4., 5., 6., 7., 1., 2., 3.},
                        {5., 6., 7., 1., 2., 3., 4.},
                        {6., 7., 1., 2., 3., 4., 5.}}});

  EXPECT_TRUE(IsFieldEqual(shifted.toFieldAligned(input), expected));
}

TEST_F(ShiftedMetricTest, FromFieldAligned) {
  ShiftedMetric shifted{*mesh, zShift};

  Field3D input{mesh};

  fillField(input, {{{1., 2., 3., 4., 5., 6., 7.},
                     {1., 2., 3., 4., 5., 6., 7.},
                     {1., 2., 3., 4., 5., 6., 7.},
                     {1., 2., 3., 4., 5., 6., 7.},
                     {1., 2., 3., 4., 5., 6., 7.}},

                    {{1., 2., 3., 4., 5., 6., 7.},
                     {1., 2., 3., 4., 5., 6., 7.},
                     {1., 2., 3., 4., 5., 6., 7.},
                     {1., 2., 3., 4., 5., 6., 7.},
                     {1., 2., 3., 4., 5., 6., 7.}},

                    {{1., 2., 3., 4., 5., 6., 7.},
                     {1., 2., 3., 4., 5., 6., 7.},
                     {1., 2., 3., 4., 5., 6., 7.},
                     {1., 2., 3., 4., 5., 6., 7.},
                     {1., 2., 3., 4., 5., 6., 7.}}});

  Field3D expected{mesh};

  fillField(expected, {{{7., 1., 2., 3., 4., 5., 6.},
                        {6., 7., 1., 2., 3., 4., 5.},
                        {5., 6., 7., 1., 2., 3., 4.},
                        {4., 5., 6., 7., 1., 2., 3.},
                        {3., 4., 5., 6., 7., 1., 2.}},

                       {{7., 1., 2., 3., 4., 5., 6.},
                        {6., 7., 1., 2., 3., 4., 5.},
                        {5., 6., 7., 1., 2., 3., 4.},
                        {4., 5., 6., 7., 1., 2., 3.},
                        {3., 4., 5., 6., 7., 1., 2.}},

                       {{7., 1., 2., 3., 4., 5., 6.},
                        {6., 7., 1., 2., 3., 4., 5.},
                        {5., 6., 7., 1., 2., 3., 4.},
                        {4., 5., 6., 7., 1., 2., 3.},
                        {3., 4., 5., 6., 7., 1., 2.}}});

  // Loosen tolerance a bit due to FFTs
  EXPECT_TRUE(IsFieldEqual(shifted.fromFieldAligned(input), expected, "RGN_ALL", 1.e-12));
}

TEST_F(ShiftedMetricTest, CalcYUpDown) {
  output_info.disable();
  auto region_yup = mesh->getRegion("RGN_NOY");
  region_yup.periodicShift(ShiftedMetricTest::nz,
                           ShiftedMetricTest::ny * ShiftedMetricTest::nz);
  mesh->addRegion("RGN_YUP", region_yup);

  auto region_ydown = mesh->getRegion("RGN_NOY");
  region_ydown.periodicShift(-ShiftedMetricTest::nz,
                             ShiftedMetricTest::ny * ShiftedMetricTest::nz);
  mesh->addRegion("RGN_YDOWN", region_ydown);
  output_info.enable();

  ShiftedMetric shifted{*mesh, zShift};

  Field3D input{mesh};

  fillField(input, {{{1., 2., 3., 4., 5., 6., 7.},
                     {1., 2., 3., 4., 5., 6., 7.},
                     {1., 2., 3., 4., 5., 6., 7.},
                     {1., 2., 3., 4., 5., 6., 7.},
                     {1., 2., 3., 4., 5., 6., 7.}},

                    {{1., 2., 3., 4., 5., 6., 7.},
                     {1., 2., 3., 4., 5., 6., 7.},
                     {1., 2., 3., 4., 5., 6., 7.},
                     {1., 2., 3., 4., 5., 6., 7.},
                     {1., 2., 3., 4., 5., 6., 7.}},

                    {{1., 2., 3., 4., 5., 6., 7.},
                     {1., 2., 3., 4., 5., 6., 7.},
                     {1., 2., 3., 4., 5., 6., 7.},
                     {1., 2., 3., 4., 5., 6., 7.},
                     {1., 2., 3., 4., 5., 6., 7.}}});

  shifted.calcYUpDown(input);

  Field3D expected_up{mesh};

  fillField(expected_up, {{{0., 0., 0., 0., 0., 0., 0.},
                           {2., 3., 4., 5., 6., 7., 1.},
                           {2., 3., 4., 5., 6., 7., 1.},
                           {2., 3., 4., 5., 6., 7., 1.},
                           {2., 3., 4., 5., 6., 7., 1.}},

                          {{0., 0., 0., 0., 0., 0., 0.},
                           {2., 3., 4., 5., 6., 7., 1.},
                           {2., 3., 4., 5., 6., 7., 1.},
                           {2., 3., 4., 5., 6., 7., 1.},
                           {2., 3., 4., 5., 6., 7., 1.}},

                          {{0., 0., 0., 0., 0., 0., 0.},
                           {2., 3., 4., 5., 6., 7., 1.},
                           {2., 3., 4., 5., 6., 7., 1.},
                           {2., 3., 4., 5., 6., 7., 1.},
                           {2., 3., 4., 5., 6., 7., 1.}}});

  Field3D expected_down{mesh};

  fillField(expected_down, {{{7., 1., 2., 3., 4., 5., 6.},
                             {7., 1., 2., 3., 4., 5., 6.},
                             {7., 1., 2., 3., 4., 5., 6.},
                             {7., 1., 2., 3., 4., 5., 6.},
                             {0., 0., 0., 0., 0., 0., 0.}},

                            {{7., 1., 2., 3., 4., 5., 6.},
                             {7., 1., 2., 3., 4., 5., 6.},
                             {7., 1., 2., 3., 4., 5., 6.},
                             {7., 1., 2., 3., 4., 5., 6.},
                             {0., 0., 0., 0., 0., 0., 0.}},

                            {{7., 1., 2., 3., 4., 5., 6.},
                             {7., 1., 2., 3., 4., 5., 6.},
                             {7., 1., 2., 3., 4., 5., 6.},
                             {7., 1., 2., 3., 4., 5., 6.},
                             {0., 0., 0., 0., 0., 0., 0.}}});

  EXPECT_TRUE(IsFieldEqual(input.yup(), expected_up, "RGN_YUP"));
  EXPECT_TRUE(IsFieldEqual(input.ydown(), expected_down, "RGN_YDOWN"));
}
