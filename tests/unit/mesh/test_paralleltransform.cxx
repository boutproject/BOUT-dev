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

TEST_F(ParallelTransformTest, IdentityToFieldAligned) {

  ParallelTransformIdentity transform{*bout::globals::mesh};

  Field3D field{1.0};

  Field3D result = transform.toFieldAligned(field, RGN_ALL);

  EXPECT_TRUE(IsFieldEqual(result, 1.0));
  EXPECT_TRUE(result.getDirectionY() == YDirectionType::Aligned);
}

TEST_F(ParallelTransformTest, IdentityFromFieldAligned) {

  ParallelTransformIdentity transform{*bout::globals::mesh};

  Field3D field{1.0};
  field.setDirectionY(YDirectionType::Aligned);

  Field3D result = transform.fromFieldAligned(field, RGN_ALL);

  EXPECT_TRUE(IsFieldEqual(result, 1.0));
  EXPECT_TRUE(result.getDirectionY() == YDirectionType::Standard);
}

TEST_F(ParallelTransformTest, IdentityToFieldAlignedFieldPerp) {

  ParallelTransformIdentity transform{*bout::globals::mesh};

  FieldPerp field{1.0};
  field.setIndex(2);

  FieldPerp result = transform.toFieldAligned(field, RGN_ALL);

  EXPECT_TRUE(IsFieldEqual(result, 1.0));
  EXPECT_TRUE(result.getDirectionY() == YDirectionType::Aligned);
}

TEST_F(ParallelTransformTest, IdentityFromFieldAlignedFieldPerp) {

  ParallelTransformIdentity transform{*bout::globals::mesh};

  FieldPerp field{1.0};
  field.setIndex(2);
  field.setDirectionY(YDirectionType::Aligned);

  FieldPerp result = transform.fromFieldAligned(field, RGN_ALL);

  EXPECT_TRUE(IsFieldEqual(result, 1.0));
  EXPECT_TRUE(result.getDirectionY() == YDirectionType::Standard);
}
