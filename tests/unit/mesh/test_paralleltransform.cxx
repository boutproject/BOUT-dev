#include "gtest/gtest.h"

#include "test_extras.hxx"
#include "bout/paralleltransform.hxx"

extern Mesh* mesh;

using ParallelTransformTest = FakeMeshFixture;

TEST_F(ParallelTransformTest, IdentityCalcYUpDown) {

  ParallelTransformIdentity transform{};

  Field3D field{1.0};

  transform.calcYUpDown(field);

  EXPECT_TRUE(IsField3DEqualBoutReal(field.yup(), 1.0));
  EXPECT_TRUE(IsField3DEqualBoutReal(field.ydown(), 1.0));
}

TEST_F(ParallelTransformTest, IdentityCalcYUpDownTwoSlices) {

  ParallelTransformIdentity transform{};

  mesh->ystart = 2;

  Field3D field{1.0};

  transform.calcYUpDown(field);

  EXPECT_TRUE(IsField3DEqualBoutReal(field.yup(0), 1.0));
  EXPECT_TRUE(IsField3DEqualBoutReal(field.yup(1), 1.0));

  EXPECT_TRUE(IsField3DEqualBoutReal(field.ydown(0), 1.0));
  EXPECT_TRUE(IsField3DEqualBoutReal(field.ydown(1), 1.0));
}
