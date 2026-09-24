#include "gtest/gtest.h"

#include "test_extras.hxx"
#include "bout/field2d.hxx"
#include "bout/field3d.hxx"

#include "fake_mesh_fixture.hxx"

#include <type_traits>

using IfElseTest = FakeMeshFixture;

TEST_F(IfElseTest, Field2DChoosesSelectedBranch) {
  const Field2D lhs{makeField<Field2D>(
      [](const Ind2D& i) { return static_cast<BoutReal>(i.x() + i.y()); })};
  const Field2D rhs{makeField<Field2D>(
      [](const Ind2D& i) { return static_cast<BoutReal>(10 + i.x() - i.y()); })};

  const auto expr = if_else(true, lhs, rhs);

  static_assert(std::is_same_v<std::decay_t<decltype(expr)>,
                               BinaryExpr<Field2D, Field2D, Field2D, bout::op::IfElse>>);
  EXPECT_TRUE(IsFieldEqual(expr, lhs));
  EXPECT_TRUE(IsFieldEqual(if_else(false, lhs, rhs), rhs));
}

TEST_F(IfElseTest, Field3DMixesField2DAndField3D) {
  const Field2D lhs{makeField<Field2D>(
      [](const Ind2D& i) { return static_cast<BoutReal>(i.x() + 2 * i.y()); })};
  const Field3D rhs{makeField<Field3D>(
      [](const Ind3D& i) { return static_cast<BoutReal>(100 + i.x() + i.y() + i.z()); })};

  const auto expr = if_else(true, lhs, rhs);
  const Field3D expected{lhs};

  static_assert(std::is_same_v<std::decay_t<decltype(expr)>,
                               BinaryExpr<Field3D, Field2D, Field3D, bout::op::IfElse>>);
  EXPECT_TRUE(IsFieldEqual(expr, expected));
  EXPECT_TRUE(IsFieldEqual(if_else(false, lhs, rhs), rhs));
}

TEST_F(IfElseTest, IfElseZeroKeepsExpressionWhenConditionTrue) {
  const Field3D field{makeField<Field3D>(
      [](const Ind3D& i) { return static_cast<BoutReal>(1 + i.x() + i.y() + i.z()); })};
  const auto source = 2.0 * field + 1.0;

  EXPECT_TRUE(IsFieldEqual(if_else_zero(true, source), source));
  EXPECT_TRUE(IsFieldEqual(if_else_zero(false, source), 0.0));
}

TEST_F(IfElseTest, InactiveBranchIsNotEvaluatedThroughMaskedArithmetic) {
  const Field2D lhs{makeField<Field2D>(
      [](const Ind2D& i) { return static_cast<BoutReal>(1 + i.x() + i.y()); })};
  const Field2D rhs{filledFrom(lhs, BoutNaN)};

  EXPECT_TRUE(IsFieldEqual(if_else(true, lhs, rhs), lhs));
}
