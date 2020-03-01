#include "gtest/gtest.h"

#include "bout/sys/variant.hxx"

#include <string>
#include <typeinfo>

// Check that we can compare values which can and can't be stored in the variant
TEST(VariantTest, IsEqualTest) {
  bout::utils::variant<bool, int> var;

  var = true;

  EXPECT_TRUE(bout::utils::variantEqualTo(var, true));
  EXPECT_FALSE(bout::utils::variantEqualTo(var, false));
  EXPECT_FALSE(bout::utils::variantEqualTo(var, 2));
  EXPECT_FALSE(bout::utils::variantEqualTo(var, "test"));

  var = 2;

  EXPECT_FALSE(bout::utils::variantEqualTo(var, true));
  EXPECT_FALSE(bout::utils::variantEqualTo(var, false));
  EXPECT_TRUE(bout::utils::variantEqualTo(var, 2));
  EXPECT_FALSE(bout::utils::variantEqualTo(var, 3));
  EXPECT_FALSE(bout::utils::variantEqualTo(var, "test"));
}

// Check casting to allowed types
TEST(VariantTest, StaticCastTest) {
  using VarType = bout::utils::variant<bool, int>;
  VarType var;

  var = 3;

  // Note: Extra brackets needed to avoid passing macro 3 arguments
  EXPECT_EQ((bout::utils::variantStaticCastOrThrow<VarType, int>(var)), 3);
  EXPECT_DOUBLE_EQ((bout::utils::variantStaticCastOrThrow<VarType, BoutReal>(var)), 3.0);
  EXPECT_THROW((bout::utils::variantStaticCastOrThrow<VarType, std::string>(var)),
               std::bad_cast);
}

// Check conversion of variant to std::string
TEST(VariantTest, ToStringTest) {
  bout::utils::variant<bool, int> var;

  var = 42;

  EXPECT_EQ(bout::utils::variantToString(var), "42");

  var = true;

  EXPECT_EQ(bout::utils::variantToString(var), "true");
}
