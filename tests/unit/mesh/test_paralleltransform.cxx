#include "gtest/gtest.h"

#include "test_extras.hxx"
#include "bout/paralleltransform.hxx"

namespace bout {
namespace globals {
extern Mesh* mesh;
}
} // namespace bout

using ParallelTransformTest = FakeMeshFixture;

TEST_F(ParallelTransformTest, IdentityCalcParallelSlices) {

  ParallelTransformIdentity transform{*bout::globals::mesh};

  Field3D field{1.0};

  transform.calcParallelSlices(field);

  EXPECT_TRUE(IsFieldEqual(field.yup(), 1.0));
  EXPECT_TRUE(IsFieldEqual(field.ydown(), 1.0));
}

TEST_F(ParallelTransformTest, IdentityCalcTwoParallelSlices) {

  ParallelTransformIdentity transform{*bout::globals::mesh};

  bout::globals::mesh->ystart = 2;

  Field3D field{1.0};

  transform.calcParallelSlices(field);

  EXPECT_TRUE(IsFieldEqual(field.yup(0), 1.0));
  EXPECT_TRUE(IsFieldEqual(field.yup(1), 1.0));

  EXPECT_TRUE(IsFieldEqual(field.ydown(0), 1.0));
  EXPECT_TRUE(IsFieldEqual(field.ydown(1), 1.0));
}
