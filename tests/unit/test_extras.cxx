#include "test_extras.hxx"
#include "fake_mesh_fixture.hxx"
#include "bout/bout_types.hxx"
#include "bout/field2d.hxx"
#include "bout/field3d.hxx"
#include "bout/region.hxx"
#include "gtest/gtest.h"

#include <cmath>
#include <string>
#include <vector>

::testing::AssertionResult IsSubString(const std::string& str,
                                       const std::string& substring) {
  if (str.find(substring) != std::string::npos) {
    return ::testing::AssertionSuccess();
  } else {
    return ::testing::AssertionFailure()
           << '"' << substring << "\" not found in \"" << str << '"';
  }
}

void fillField(Field3D& f, std::vector<std::vector<std::vector<BoutReal>>> values) {
  f.allocate();
  Ind3D i{0};
  for (auto& x : values) {
    for (auto& y : x) {
      for (auto& z : y) {
        f[i] = z;
        ++i;
      }
    }
  }
}

void fillField(Field2D& f, std::vector<std::vector<BoutReal>> values) {
  f.allocate();
  Ind2D i{0};
  for (auto& x : values) {
    for (auto& y : x) {
      f[i] = y;
      ++i;
    }
  }
}

using TestExtrasFieldExpr = FakeMeshFixture;

TEST_F(TestExtrasFieldExpr, IsFieldEqualHandlesBinaryExprOnEitherSide) {
  const Field2D field{1.0};
  const Field2D expected{3.0};

  EXPECT_TRUE(IsFieldEqual(field + 2.0, expected));
  EXPECT_TRUE(IsFieldEqual(expected, field + 2.0));
}

TEST_F(TestExtrasFieldExpr, BinaryExprCanBeIndexedWithRegionIndex) {
  const Field3D lhs{
      makeField<Field3D>([](const Ind3D& i) { return static_cast<BoutReal>(i.x()); })};
  const Field3D rhs{
      makeField<Field3D>([](const Ind3D& i) { return static_cast<BoutReal>(i.y()); })};

  const auto expr = lhs + 2.0 * rhs;
  Field3D result{emptyFrom(lhs)};

  BOUT_FOR_SERIAL(i, result.getRegion("RGN_ALL")) { result[i] = expr[i]; }

  BOUT_FOR(i, result.getRegion("RGN_ALL")) {
    EXPECT_DOUBLE_EQ(result[i], lhs[i] + 2.0 * rhs[i]);
  }
}
