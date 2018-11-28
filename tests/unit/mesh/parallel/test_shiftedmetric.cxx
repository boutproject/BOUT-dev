#include "gtest/gtest.h"

#include "bout/paralleltransform.hxx"
#include "test_extras.hxx"

extern Mesh* mesh;

using ShiftedMetricTest = FakeMeshFixture;

TEST_F(ShiftedMetricTest, ToFieldAligned) {
  Field2D zShift{mesh};

  fillField(zShift, {{1., 1., 1., 1., 1.}, {1., 1., 1., 1., 1.}, {1., 1., 1., 1., 1.}});

  dynamic_cast<FakeMesh*>(mesh)->setCoordinates(std::make_shared<FakeCoordinates>(
      mesh, Field2D{1.0}, Field2D{1.0}, BoutReal{1.0}, Field2D{1.0}, Field2D{0.0},
      Field2D{1.0}, Field2D{1.0}, Field2D{1.0}, Field2D{0.0}, Field2D{0.0}, Field2D{0.0},
      Field2D{1.0}, Field2D{1.0}, Field2D{1.0}, Field2D{0.0}, Field2D{0.0},
      Field2D{0.0}));

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

  Field3D output{mesh};

  fillField(output, {{{2., 3., 4., 5., 6., 7., 1.},
                      {2., 3., 4., 5., 6., 7., 1.},
                      {2., 3., 4., 5., 6., 7., 1.},
                      {2., 3., 4., 5., 6., 7., 1.},
                      {2., 3., 4., 5., 6., 7., 1.}},

                     {{2., 3., 4., 5., 6., 7., 1.},
                      {2., 3., 4., 5., 6., 7., 1.},
                      {2., 3., 4., 5., 6., 7., 1.},
                      {2., 3., 4., 5., 6., 7., 1.},
                      {2., 3., 4., 5., 6., 7., 1.}},

                     {{2., 3., 4., 5., 6., 7., 1.},
                      {2., 3., 4., 5., 6., 7., 1.},
                      {2., 3., 4., 5., 6., 7., 1.},
                      {2., 3., 4., 5., 6., 7., 1.},
                      {2., 3., 4., 5., 6., 7., 1.}}});

  EXPECT_TRUE(IsField3DEqualField3D(shifted.toFieldAligned(input), output));
}
