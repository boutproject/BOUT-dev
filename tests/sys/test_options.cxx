#include "gtest/gtest.h"
#include "options.hxx"

#include <string>

TEST(OptionsTest, IsSet) {
  Options options;
  options.set("int_key", 42, "code");

  ASSERT_TRUE(options.isSet("int_key"));
}

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

TEST(OptionsTest, SingletonTest) {
  Options *root = Options::getRoot();
  Options *second = Options::getRoot();

  EXPECT_EQ(root, second);
}

::testing::AssertionResult IsSubString(const std::string &str,
                                       const std::string &substring) {
  if (str.find(substring) != std::string::npos) {
    return ::testing::AssertionSuccess();
  } else {
    return ::testing::AssertionFailure() << substring << " not found in " << str;
  }
}

TEST(OptionsTest, CheckUsed) {
  // stdout redirection code from https://stackoverflow.com/a/4043813/2043465

  std::stringstream buffer;
  // Save cout's buffer here
  std::streambuf *sbuf = std::cout.rdbuf();
  // Redirect cout to our stringstream buffer or any other ostream
  std::cout.rdbuf(buffer.rdbuf());

  Options options;
  options.set("key1", "a", "code");
  options.set("key2", "b", "code");
  options.set("key3", "c", "code");
  options.set("key4", "d", "code");

  options.printUnused();

  // Check keys are all in buffer
  EXPECT_TRUE(IsSubString(buffer.str(), "key1"));
  EXPECT_TRUE(IsSubString(buffer.str(), "key2"));
  EXPECT_TRUE(IsSubString(buffer.str(), "key3"));
  EXPECT_TRUE(IsSubString(buffer.str(), "key4"));

  // Clear buffer
  buffer.str("");

  std::string value;
  options.get("key1", value, "--", false);
  options.get("key2", value, "--", false);

  options.printUnused();

  // Check only key3, key4 are in buffer
  EXPECT_FALSE(IsSubString(buffer.str(), "key1"));
  EXPECT_FALSE(IsSubString(buffer.str(), "key2"));
  EXPECT_TRUE(IsSubString(buffer.str(), "key3"));
  EXPECT_TRUE(IsSubString(buffer.str(), "key4"));

  // Clear buffer
  buffer.str("");

  options.get("key3", value, "--", false);
  options.get("key4", value, "--", false);

  options.printUnused();

  EXPECT_TRUE(IsSubString(buffer.str(), "All options used"));

  // When done redirect cout to its old self
  std::cout.rdbuf(sbuf);
}

TEST(OptionsTest, GetEmptySection) {
  Options options;
  Options *new_section = options.getSection("");

  EXPECT_EQ(new_section, &options);
}

TEST(OptionsTest, MakeNewSection) {
  Options options;
  Options *new_section = options.getSection("section1");

  EXPECT_NE(new_section, &options);
  EXPECT_EQ(new_section->getParent(), &options);
  EXPECT_EQ(new_section->str(), "section1");
}

TEST(OptionsTest, GetExistingSection) {
  Options options;
  Options *new_section = options.getSection("section1");
  Options *old_section = options.getSection("section1");

  EXPECT_EQ(new_section, old_section);
}

TEST(OptionsTest, CheckCaseSensitivity) {
  Options options;
  Options *new_section = options.getSection("section1");
  Options *old_section = options.getSection("SECTION1");

  EXPECT_EQ(new_section, old_section);
}

TEST(OptionsTest, GetCorrectSection) {
  Options options;
  Options *section1 = options.getSection("section1");
  options.getSection("section2");

  Options *old_section = options.getSection("section1");

  EXPECT_EQ(section1, old_section);
}

TEST(OptionsTest, MakeNestedSection) {
  Options options;
  Options *section1 = options.getSection("section1");
  Options *section2 = section1->getSection("section2");

  EXPECT_NE(section2, section1);
  EXPECT_EQ(section2->getParent(), section1);
  EXPECT_EQ(section2->str(), "section1:section2");
}
