#include "gtest/gtest.h"
#include "test_extras.hxx"

#include "options.hxx"
#include "output.hxx"
#include <boutexception.hxx>

#include <string>

class OptionsTest : public ::testing::Test {
public:
  OptionsTest() {
    output_info.disable();
  }

  ~OptionsTest() {
    output_info.enable();
  }
};

TEST_F(OptionsTest, IsSet) {
  Options options;
  options.set("int_key", 42, "code");

  ASSERT_TRUE(options.isSet("int_key"));
}

TEST_F(OptionsTest, SetGetInt) {
  Options options;
  options.set("int_key", 42, "code");

  ASSERT_TRUE(options.isSet("int_key"));

  int value;
  options.get("int_key", value, 99, false);

  EXPECT_EQ(value, 42);
}

TEST_F(OptionsTest, SetGetIntFromReal) {
  Options options;
  options.set("int_key", 42.00001, "code");

  ASSERT_TRUE(options.isSet("int_key"));

  int value;
  options.get("int_key", value, 99, false);

  EXPECT_EQ(value, 42);

  options.set("int_key2", 12.5, "code");
  EXPECT_THROW(options.get("int_key2", value, 99, false), BoutException);
  // Note we expect to get the rounded value despite the throw as we
  // pass by reference and modify the passed variable in options.get.
  EXPECT_EQ(value, 13);
}

TEST_F(OptionsTest, DefaultValueInt) {
  Options options;

  int value;
  options.get("int_key", value, 99, false);

  EXPECT_EQ(value, 99);
}

TEST_F(OptionsTest, InconsistentDefaultValueInt) {
  Options options;

  int value;
  options.get("int_key", value, 99, false);
  EXPECT_THROW(options.get("int_key", value, 98, false), BoutException);

  EXPECT_EQ(value, 99);
}

TEST_F(OptionsTest, SetGetReal) {
  Options options;
  options.set("real_key", 6.7e8, "code");

  ASSERT_TRUE(options.isSet("real_key"));

  BoutReal value;
  options.get("real_key", value, -78.0, false);

  EXPECT_DOUBLE_EQ(value, 6.7e8);
}

TEST_F(OptionsTest, SetGetDouble) {
  Options options;
  options.set("real_key", 0.7853981633974483, "code");

  ASSERT_TRUE(options.isSet("real_key"));

  BoutReal value;
  options.get("real_key", value, -78.0, false);

  EXPECT_DOUBLE_EQ(value, 0.7853981633974483);
}

TEST_F(OptionsTest, SetGetNegativeDouble) {
  Options options;
  options.set("real_key", -0.7853981633974483, "code");

  ASSERT_TRUE(options.isSet("real_key"));

  BoutReal value;
  options.get("real_key", value, -78.0, false);

  EXPECT_DOUBLE_EQ(value, -0.7853981633974483);
}

TEST_F(OptionsTest, DefaultValueReal) {
  Options options;

  BoutReal value;
  options.get("real_key", value, -78.0, false);

  EXPECT_DOUBLE_EQ(value, -78.0);
}

TEST_F(OptionsTest, InconsistentDefaultValueReal) {
  Options options;

  BoutReal value;
  options.get("real_key", value, -78.0, false);
  EXPECT_THROW(options.get("real_key", value, -68.0, false), BoutException);

  EXPECT_EQ(value, -78.0);
}

TEST_F(OptionsTest, GetBool) {
  Options options;
  bool value;
  options.get("bool_key", value, true, false);
  EXPECT_EQ(value, true);
}

TEST_F(OptionsTest, SetGetBool) {
  Options options;
  options.set("bool_key", true, "code");

  ASSERT_TRUE(options.isSet("bool_key"));

  bool value;
  options.get("bool_key", value, false, false);

  EXPECT_EQ(value, true);
}

TEST_F(OptionsTest, SetGetBoolFalse) {
  Options options;
  options.set("bool_key", false, "code");

  ASSERT_TRUE(options.isSet("bool_key"));

  bool value;
  options.get("bool_key", value, true, false);

  EXPECT_EQ(value, false);
}

TEST_F(OptionsTest, GetBoolFromString) {
  Options options;
  options.set("bool_key", "true", "code");
  options.set("bool_key2", "yes", "code");

  ASSERT_TRUE(options.isSet("bool_key"));

  bool value;
  options.get("bool_key", value, false, false);

  EXPECT_EQ(value, true);

  bool value2;
  options.get("bool_key2", value2, false, false);

  EXPECT_EQ(value2, true);

  bool value3;
  // Note we only test the first character so "not_a_bool" is treated as
  // a bool that is false.
  options.set("bool_key3", "A_bool_starts_with_T_or_N_or_Y_or_F_or_1_or_0", "code");
  EXPECT_THROW(options.get("bool_key3", value3, false, false), BoutException);
  // Surprise true
  options.set("bool_key3", "yes_this_is_a_bool", "code");
  EXPECT_NO_THROW(options.get("bool_key3", value3, false, false));
  EXPECT_EQ(value3, true);
}

TEST_F(OptionsTest, DefaultValueBool) {
  Options options;

  bool value;
  options.get("bool_key", value, false, false);

  EXPECT_EQ(value, false);
}

TEST_F(OptionsTest, SetGetString) {
  Options options;
  options.set("string_key", "abcdef", "code");

  ASSERT_TRUE(options.isSet("string_key"));

  std::string value;
  options.get("string_key", value, "ghijkl", false);

  EXPECT_EQ(value, "abcdef");
}

TEST_F(OptionsTest, DefaultValueString) {
  Options options;

  std::string value;
  options.get("string_key", value, "ghijkl", false);

  EXPECT_EQ(value, "ghijkl");
}

TEST_F(OptionsTest, InconsistentDefaultValueString) {
  Options options;

  std::string value;
  options.get("string_key", value, "ghijkl", false);

  EXPECT_EQ(value, "ghijkl");

  EXPECT_THROW(options.get("string_key", value, "_ghijkl", false), BoutException);

  EXPECT_EQ(value, "ghijkl");
}

TEST_F(OptionsTest, SingletonTest) {
  Options *root = Options::getRoot();
  Options *second = Options::getRoot();

  EXPECT_EQ(root, second);
}

TEST_F(OptionsTest, CheckUsed) {
  // stdout redirection code from https://stackoverflow.com/a/4043813/2043465

  // Need output_info enabled, as Options::printUnused writes to it
  output_info.enable();

  std::stringstream buffer;
  // Save cout's buffer here
  std::streambuf *sbuf = std::cout.rdbuf();
  // Redirect cout to our stringstream buffer or any other ostream
  std::cout.rdbuf(buffer.rdbuf());

  Options options;
  Options *section1 = options.getSection("section1");
  options.set("key1", "a", "code");
  section1->set("key2", "b", "code");
  options.set("key3", "c", "code");
  section1->set("key4", "d", "code");

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
  section1->get("key2", value, "--", false);

  // Clear buffer
  buffer.str("");
  options.printUnused();

  // Check only key3, key4 are in buffer
  EXPECT_FALSE(IsSubString(buffer.str(), "key1"));
  EXPECT_FALSE(IsSubString(buffer.str(), "section1:key2"));
  EXPECT_TRUE(IsSubString(buffer.str(), "key3"));
  EXPECT_TRUE(IsSubString(buffer.str(), "section1:key4"));

  // Clear buffer
  buffer.str("");

  options.get("key3", value, "--", false);
  section1->get("key4", value, "--", false);

  options.printUnused();

  EXPECT_TRUE(IsSubString(buffer.str(), "All options used"));

  // When done redirect cout to its old self
  std::cout.rdbuf(sbuf);
}

TEST_F(OptionsTest, GetEmptySection) {
  Options options;
  Options *new_section = options.getSection("");

  EXPECT_EQ(new_section, &options);
}

TEST_F(OptionsTest, MakeNewSection) {
  Options options;
  Options *new_section = options.getSection("section1");

  EXPECT_NE(new_section, &options);
  EXPECT_EQ(new_section->getParent(), &options);
  EXPECT_EQ(new_section->str(), "section1");
}

TEST_F(OptionsTest, GetExistingSection) {
  Options options;
  Options *new_section = options.getSection("section1");
  Options *old_section = options.getSection("section1");

  EXPECT_EQ(new_section, old_section);
}

TEST_F(OptionsTest, CheckCaseSensitivity) {
  Options options;
  Options *new_section = options.getSection("section1");
  Options *old_section = options.getSection("SECTION1");

  EXPECT_EQ(new_section, old_section);
}

TEST_F(OptionsTest, GetCorrectSection) {
  Options options;
  Options *section1 = options.getSection("section1");
  options.getSection("section2");

  Options *old_section = options.getSection("section1");

  EXPECT_EQ(section1, old_section);
}

TEST_F(OptionsTest, MakeNestedSection) {
  Options options;
  Options *section1 = options.getSection("section1");
  Options *section2 = section1->getSection("section2");

  EXPECT_NE(section2, section1);
  EXPECT_EQ(section2->getParent(), section1);
  EXPECT_EQ(section2->str(), "section1:section2");
}
