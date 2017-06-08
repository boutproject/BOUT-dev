#include "gtest/gtest.h"
#include "options.hxx"

#include <string>

TEST(OptionsTest, SetGetInt) {
  Options options;
  options.set("int_key", 42, "code");

  ASSERT_TRUE(options.isSet("int_key"));

  int value;
  options.get("int_key", value, 99, false);

  EXPECT_EQ(value, 42);
}

TEST(OptionsTest, DefaultValueInt) {
  Options options;

  int value;
  options.get("int_key", value, 99, false);

  EXPECT_EQ(value, 99);
}

TEST(OptionsTest, SetGetReal) {
  Options options;
  options.set("real_key", 6.7e8, "code");

  ASSERT_TRUE(options.isSet("real_key"));

  BoutReal value;
  options.get("real_key", value, -78.0, false);

  EXPECT_DOUBLE_EQ(value, 6.7e8);
}

TEST(OptionsTest, DefaultValueReal) {
  Options options;

  BoutReal value;
  options.get("real_key", value, -78.0, false);

  EXPECT_DOUBLE_EQ(value, -78.0);
}

TEST(OptionsTest, SetGetBool) {
  Options options;
  options.set("bool_key", true, "code");

  ASSERT_TRUE(options.isSet("bool_key"));

  bool value;
  options.get("bool_key", value, false, false);

  EXPECT_EQ(value, true);
}

TEST(OptionsTest, GetBoolFromString) {
  Options options;
  options.set("bool_key", "true", "code");
  options.set("bool_key2", "yes", "code");

  ASSERT_TRUE(options.isSet("bool_key"));

  bool value;
  options.get("bool_key", value, false, false);

  EXPECT_EQ(value, true);

  bool value2;
  options.get("bool_key", value2, false, false);

  EXPECT_EQ(value2, true);
}

TEST(OptionsTest, DefaultValueBool) {
  Options options;

  bool value;
  options.get("bool_key", value, false, false);

  EXPECT_EQ(value, false);
}

TEST(OptionsTest, SetGetString) {
  Options options;
  options.set("string_key", "abcdef", "code");

  ASSERT_TRUE(options.isSet("string_key"));

  std::string value;
  options.get("string_key", value, "ghijkl", false);

  EXPECT_EQ(value, "abcdef");
}

TEST(OptionsTest, DefaultValueString) {
  Options options;

  std::string value;
  options.get("string_key", value, "ghijkl", false);

  EXPECT_EQ(value, "ghijkl");
}
