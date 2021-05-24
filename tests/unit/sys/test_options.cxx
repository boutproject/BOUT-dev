#include "gtest/gtest.h"
#include "test_extras.hxx"

#include "options.hxx"
#include "output.hxx"
#include <boutexception.hxx>

#include <string>

class OptionsTest : public FakeMeshFixture {
public:
  virtual ~OptionsTest() = default;
  WithQuietOutput quiet_info{output_info};
  WithQuietOutput quiet_warn{output_warn};
  WithQuietOutput quiet_progress{output_progress};
};

TEST_F(OptionsTest, IsSet) {
  Options options;

  ASSERT_FALSE(options.isSet("int_key"));
  
  options.set("int_key", 42, "code");

  ASSERT_TRUE(options.isSet("int_key"));
}

TEST_F(OptionsTest, IsSetDefault) {
  Options options;
  int value;
  ASSERT_FALSE(options.isSet("default_value"));
  options.get("default_value", value, 42);
  ASSERT_FALSE(options.isSet("default_value"));
}

TEST_F(OptionsTest, IsSection) {
  Options options;

  // make sure options is initialized as a section
  options["testkey"] = 1.;

  ASSERT_TRUE(options.isSection());
  ASSERT_FALSE(options["testkey"].isSection());
  ASSERT_TRUE(options.isSection(""));
  ASSERT_FALSE(options.isSection("subsection"));

  options["subsection"]["testkey"] = 1.;

  ASSERT_TRUE(options.isSection("subsection"));
}

TEST_F(OptionsTest, CompoundName) {
  Options options;

  // make sure options is initialized as a section
  options["compoundkey"] = 321.;

  ASSERT_TRUE(options.isSection());
  ASSERT_FALSE(options["compoundkey"].isSection());
  ASSERT_TRUE(options.isSection(""));
  ASSERT_FALSE(options.isSection("compoundsubsection"));

  options["compoundsubsection:compoundkey"] = 321.;

  ASSERT_TRUE(options.isSection("compoundsubsection"));

  BoutReal value = options["compoundsubsection"]["compoundkey"];
  EXPECT_EQ(value, 321.);

  BoutReal value2 = options["compoundsubsection:compoundkey"];
  EXPECT_EQ(value2, 321.);
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
  
  // value is not changed
  EXPECT_EQ(value, 42);
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

TEST_F(OptionsTest, DefaultValueOptions) {
  Options options, default_options;

  default_options.set("int_key", 99);

  int value = options["int_key"].withDefault(default_options["int_key"]).as<int>();

  EXPECT_EQ(value, 99);
}

TEST_F(OptionsTest, InconsistentDefaultValueOptions) {
  Options options, default_options;

  default_options.set("int_key", 99);

  EXPECT_EQ(options["int_key"].withDefault(42), 42);

  int value = 0;
  EXPECT_THROW(
      value = options["int_key"].withDefault(default_options["int_key"]).as<int>(),
      BoutException);

  EXPECT_EQ(value, 0);
}

TEST_F(OptionsTest, OverrideDefaultValueOptions) {
  Options options;

  options["override_key"].overrideDefault("override_value");

  std::string value = options["override_key"].withDefault("default_value");

  EXPECT_EQ(value, "override_value");
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
  EXPECT_FALSE(IsSubString(buffer.str(), "key2"));
  EXPECT_TRUE(IsSubString(buffer.str(), "key3"));
  EXPECT_TRUE(IsSubString(buffer.str(), "key4"));

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

TEST_F(OptionsTest, SetSameOptionTwice) {
  Options options;
  options.set("key", "value", "code");
  EXPECT_THROW(options.set("key", "new value", "code"),BoutException);

  options.set("key", "value", "code");
  EXPECT_NO_THROW(options.forceSet("key", "new value", "code"));
  EXPECT_NO_THROW(options.set("key", "value", "code",true));
}

/// New interface


TEST_F(OptionsTest, NewIsSet) {
  Options options;

  ASSERT_FALSE(options["int_key"].isSet());

  options["int_key"].assign(42, "code");

  ASSERT_TRUE(options["int_key"].isSet());
}

TEST_F(OptionsTest, NewSubSection) {
  Options options;
  
  options["sub-section"]["int_key"].assign(42, "code");
  
  ASSERT_FALSE(options["int_key"].isSet());
  ASSERT_TRUE(options["sub-section"]["int_key"].isSet());
  
  int value = options["sub-section"]["int_key"].withDefault(99);
  EXPECT_EQ(value, 42);
}

TEST_F(OptionsTest, NewIsSetDefault) {
  Options options;
  ASSERT_FALSE(options.isSet());
  int value = options.withDefault(42);
  ASSERT_EQ(value, 42);
  ASSERT_FALSE(options.isSet());
}

TEST_F(OptionsTest, NewSetGetInt) {
  Options options;
  options.assign(42, "code");

  ASSERT_TRUE(options.isSet());

  int value = options.withDefault(99);

  EXPECT_EQ(value, 42);
}

TEST_F(OptionsTest, NewSetGetIntFromReal) {
  Options options;
  options["key1"] = 42.00001;

  ASSERT_TRUE(options["key1"].isSet());

  int value = options["key1"].withDefault(99);

  EXPECT_EQ(value, 42);

  options["key2"] = 12.5;
  EXPECT_THROW(options["key2"].as<int>(), BoutException);
}

TEST_F(OptionsTest, NewDefaultValueInt) {
  Options options;

  int value = options.withDefault(99);
  EXPECT_EQ(value, 99);
}

TEST_F(OptionsTest, WithDefaultString) {
  Options options;

  std::string value = options.withDefault("hello");
  EXPECT_EQ(value, "hello");
}

TEST_F(OptionsTest, WithDefaultStringCaseSensitive) {
  Options options;

  std::string value = options.withDefault("Hello");
  EXPECT_NE(value, "hello");
}

TEST_F(OptionsTest, OptionsMacroPointer) {
  Options options;

  options["val"] = 42;

  int val{0};
  OPTION(&options, val, 3);
  EXPECT_EQ(val, 42);
}

TEST_F(OptionsTest, OptionsMacroConstPointer) {
  Options options;

  options["val"] = 42;

  int val{0};
  OPTION(const_cast<const Options*>(&options), val, 3);
  EXPECT_EQ(val, 42);
}

TEST_F(OptionsTest, OptionsMacroReference) {
  Options options;

  options["val"] = 42;

  int val{0};
  OPTION(options, val, 3);
  EXPECT_EQ(val, 42);
}

TEST_F(OptionsTest, OptionsMacroConstReference) {
  Options options;

  options["val"] = 42;

  int val{0};
  OPTION(const_cast<const Options&>(options), val, 3);
  EXPECT_EQ(val, 42);
}

/// Copy constructor copies value
TEST_F(OptionsTest, CopyOption) {
  Options option1;

  option1 = 42;

  Options option2(option1);

  EXPECT_EQ(option2.as<int>(), 42);
}

/// Copy constructor makes independent copy
TEST_F(OptionsTest, CopyOptionDistinct) {
  Options option1;
  option1 = 42;

  Options option2(option1);

  option1.force(23);
  
  EXPECT_EQ(option1.as<int>(), 23);
  EXPECT_EQ(option2.as<int>(), 42);
}

/// Copies of sections get values
TEST_F(OptionsTest, CopySection) {
  Options option1;

  option1["key"] = 42;   // option1 now a section

  Options option2(option1);

  EXPECT_EQ(option2["key"].as<int>(), 42);
}

/// The parent should be updated when copied
TEST_F(OptionsTest, CopySectionParent) {
  Options option1;

  option1["key"] = 42;

  Options option2(option1);
  
  EXPECT_TRUE( &option2["key"].parent() == &option2 );
}

TEST_F(OptionsTest, AssignOption) {
  Options option1, option2;

  option1 = 42;
  
  option2 = option1;

  EXPECT_EQ(option2.as<int>(), 42);
}

TEST_F(OptionsTest, AssignSection) {
  Options option1, option2;

  option1["key"] = 42;
  
  option2 = option1;

  EXPECT_EQ(option2["key"].as<int>(), 42);
}

TEST_F(OptionsTest, AssignSectionReplace) {
  Options option1, option2;

  option1["key"] = 42;
  option2["key"] = 23;
  
  option2 = option1;

  EXPECT_EQ(option2["key"].as<int>(), 42);
}

TEST_F(OptionsTest, AssignSectionParent) {
  Options option1, option2;

  option1["key"] = 42;
  
  option2 = option1;
  
  EXPECT_TRUE( &option2["key"].parent() == &option2 );
}

TEST_F(OptionsTest, AssignSubSection) {
  Options option1, option2;

  option1["key1"] = 42;
  
  option2["key2"] = option1;

  EXPECT_EQ(option2["key2"]["key1"].as<int>(), 42);
}

TEST_F(OptionsTest, AssignSubSectionParent) {
  Options option1, option2;

  option1["key1"] = 42;
  
  option2["key2"] = option1;

  EXPECT_EQ(&option2["key2"].parent(), &option2);
  EXPECT_EQ(&option2["key2"]["key1"].parent(), &option2["key2"]);
}

TEST_F(OptionsTest, HasAttribute) {
  Options option;

  EXPECT_FALSE(option.hasAttribute("not here"));

  option.attributes["here"] = true;
  EXPECT_TRUE(option.hasAttribute("here"));
}

TEST_F(OptionsTest, AttributeMissingBool) {
  Options option;

  bool a = option.attributes["test"];
  EXPECT_EQ(a, false);
}

TEST_F(OptionsTest, AttributeMissingInt) {
  Options option;

  int a = option.attributes["test"];
  EXPECT_EQ(a, 0);
}

TEST_F(OptionsTest, AttributeMissingBoutReal) {
  Options option;

  BoutReal a = option.attributes["test"];
  EXPECT_DOUBLE_EQ(a, 0.0);
}

TEST_F(OptionsTest, AttributeMissingString) {
  Options option;

  EXPECT_THROW(option.attributes["test"].as<std::string>(), std::bad_cast);
}

TEST_F(OptionsTest, AttributeStoreBool) {
  Options option;
  option.attributes["test"] = true;

  EXPECT_TRUE(option.attributes["test"].as<bool>());

  option.attributes["test"] = false;
  EXPECT_FALSE(option.attributes["test"].as<bool>());
}

TEST_F(OptionsTest, AttributeStoreBoolCaseSensitive) {
  Options option;
  option.attributes["Test"] = true;

  EXPECT_FALSE(option.attributes["test"].as<bool>());
  EXPECT_TRUE(option.attributes["Test"].as<bool>());
}

TEST_F(OptionsTest, AttributeStoreInt) {
  Options option;
  option.attributes["test"] = 42;

  int value = option.attributes["test"];
  EXPECT_EQ(value, 42);
}

TEST_F(OptionsTest, AttributeStoreIntCaseSensitive) {
  Options option;
  option.attributes["Test"] = 42;

  int value = option.attributes["tEst"];
  EXPECT_NE(value, 42);
}

TEST_F(OptionsTest, AttributeStoreBoutReal) {
  Options option;
  option.attributes["test"] = 3.1415;

  BoutReal value = option.attributes["test"];
  EXPECT_DOUBLE_EQ(value, 3.1415);
}

TEST_F(OptionsTest, AttributeStoreBoutRealCaseSensitive) {
  Options option;
  option.attributes["Test"] = 3.1415;

  BoutReal value = option.attributes["tEst"];
  EXPECT_DOUBLE_EQ(value, 0.);
}

TEST_F(OptionsTest, AttributeStoreConstChars) {
  Options option;
  option.attributes["test"] = "hello";

  std::string test = option.attributes["test"];
  EXPECT_EQ(test, "hello");
}

TEST_F(OptionsTest, AttributeStoreConstCharsCaseSensitive) {
  Options option;
  option.attributes["Test"] = "HeLlO";

  std::string test = option.attributes["Test"];
  EXPECT_EQ(test, "HeLlO");
  EXPECT_NE(test, "hello");
}

TEST_F(OptionsTest, AttributeTimeDimension) {
  Options option;

  option = 3;
  EXPECT_EQ(option.as<int>(), 3);
  
  option.attributes["time_dimension"] = "t";

  option = 4;

  EXPECT_EQ(option.as<int>(), 4);
}

TEST_F(OptionsTest, EqualityBool) {
  Options option;

  option = true;

  EXPECT_TRUE(option == true);
  EXPECT_FALSE(option == false);

  option.force(false);

  EXPECT_TRUE(option == false);
  EXPECT_FALSE(option == true);
}

TEST_F(OptionsTest, EqualityInt) {
  Options option;

  option = 3;

  EXPECT_TRUE(option == 3);
  EXPECT_FALSE(option == 4);
}

TEST_F(OptionsTest, EqualityString) {
  Options option;

  option = "hello";

  EXPECT_TRUE(option == "hello");
  EXPECT_FALSE(option == "goodbye");
}

TEST_F(OptionsTest, EqualityStringCaseSensitive) {
  Options option;

  option = "HeLlO";

  EXPECT_TRUE(option == "HeLlO");
  EXPECT_FALSE(option == "hello");
  EXPECT_FALSE(option == "goodbye");
}

TEST_F(OptionsTest, ComparisonInt) {
  Options option;

  option = 3;

  EXPECT_TRUE(option < 4);
  EXPECT_FALSE(option < 3);
}

TEST_F(OptionsTest, ComparisonString) {
  Options option;

  option = "bbb";

  EXPECT_TRUE(option < "ccc");
  EXPECT_FALSE(option < "aaa");
}

TEST_F(OptionsTest, WithDefaultIntThrow) {
  // If given an integer as default, will try to cast to int

  Options option;
  option = "4.32";
  
  EXPECT_THROW(option.withDefault(0), BoutException);
}

TEST_F(OptionsTest, TypeAttributeBool) {
  Options option;
  option = "true";

  // Getting into bool using withDefault should modify the "type" attribute
  bool value = option.withDefault(false);

  EXPECT_TRUE(value);
  EXPECT_EQ(option.attributes["type"].as<std::string>(), "bool");
}

TEST_F(OptionsTest, AsNoTypeAttribute) {
  Options option;
  option = "true";

  // as is const so doesn't set the type attribute
  bool value = option.as<bool>();

  EXPECT_TRUE(value);
  EXPECT_EQ(option.attributes.count("type"), 0);
}

TEST_F(OptionsTest, TypeAttributeInt) {
  Options option;
  option = "42";

  // Casting to int should modify the "type" attribute
  int value = option.withDefault<int>(-1);

  EXPECT_EQ(value, 42);
  EXPECT_EQ(option.attributes["type"].as<std::string>(), "int");
}

TEST_F(OptionsTest, TypeAttributeField2D) {
  Options option;
  option = "42";

  // Casting to Field2D should modify the "type" attribute
  Field2D value = option.withDefault<Field2D>(Field2D(-1, bout::globals::mesh));

  EXPECT_EQ(value(0,0), 42);
  EXPECT_EQ(option.attributes["type"].as<std::string>(), "Field2D");
}

TEST_F(OptionsTest, TypeAttributeField3D) {
  Options option;
  option = "42";

  // Casting to Field3D should modify the "type" attribute
  Field3D value = option.withDefault<Field3D>(Field3D(-1, bout::globals::mesh));

  EXPECT_EQ(value(0,0,0), 42);
  EXPECT_EQ(option.attributes["type"].as<std::string>(), "Field3D");
}

TEST_F(OptionsTest, TypeAttributeFieldPerp) {
  Options option;
  option = "36";

  // Casting to FieldPerp should modify the "type" attribute
  FieldPerp value = option.withDefault<FieldPerp>(FieldPerp(-1, bout::globals::mesh));

  EXPECT_EQ(value(0,0,0), 36);
  EXPECT_EQ(option.attributes["type"].as<std::string>(), "FieldPerp");
}

TEST_F(OptionsTest, DocString) {
  Options option;

  option.doc("test string");

  EXPECT_EQ(option.attributes["doc"].as<std::string>(), "test string");
}

TEST_F(OptionsTest, DocStringAssignTo) {
  Options option;
  
  option.doc("test string") = 42;

  EXPECT_EQ(option.attributes["doc"].as<std::string>(), "test string");
  EXPECT_EQ(option.as<int>(), 42);
}

TEST_F(OptionsTest, DocStringAssignFrom) {
  Options option;
  option = 42;
  
  int value = option.doc("test string");

  EXPECT_EQ(option.attributes["doc"].as<std::string>(), "test string");
  EXPECT_EQ(value, 42);
}

TEST_F(OptionsTest, DocStringWithDefault) {
  Options option;
  option = 42;

  int value = option.doc("some value").withDefault(2);

  EXPECT_EQ(value, 42);
  EXPECT_EQ(option.attributes["doc"].as<std::string>(), "some value"); 
}

TEST_F(OptionsTest, DocStringNotCopied) {
  Options option;
  option = 32;

  Options option2 = option;

  int value = option2.doc("test value");
  
  EXPECT_EQ(value, 32);
  EXPECT_EQ(option2.attributes["doc"].as<std::string>(), "test value");
  EXPECT_EQ(option.attributes.count("doc"), 0);
}

TEST_F(OptionsTest, InitializeInt) {
  Options option {3};
  EXPECT_EQ(option.as<int>(), 3);
}

TEST_F(OptionsTest, InitialiseTree) {
  Options option {{"section1", {{"value1", 42},
                                {"value2", "hello"}}},
                  {"section2", {{"subsection1", {{"value3", true},
                                                 {"value4", 3.2}}},
                                {"value5", 3}}}};
  
  EXPECT_EQ(option["section1"]["value1"].as<int>(), 42);
  EXPECT_EQ(option["section1"]["value2"].as<std::string>(), "hello");
  EXPECT_EQ(option["section2"]["subsection1"]["value3"].as<bool>(), true);
  EXPECT_DOUBLE_EQ(option["section2"]["subsection1"]["value4"].as<BoutReal>(), 3.2);
  EXPECT_EQ(option["section2"]["value5"].as<int>(), 3);
}

TEST_F(OptionsTest, ToString) {
  Options option{
      {"section1", {{"value1", 42}, {"value2", "hello"}}},
      {"section2", {{"subsection1", {{"value3", true}, {"value4", 3.2}}}, {"value5", 3}}},
      {"section3", {{"subsection2", {{"value6", 12}}}}}};

  // It's plausible this test is fragile if the internal storage
  // changes the order -- at time of writing (Jan 2020) it's
  // lexographical rather than insertion order
  std::string expected = R"(
[section1]
value1 = 42
value2 = hello

[section2]
value5 = 3

[section2:subsection1]
value3 = true
value4 = 3.2

[section3:subsection2]
value6 = 12
)";

  EXPECT_EQ(toString(option), expected);
}

TEST_F(OptionsTest, GetUnused) {
  Options option{{"section1", {{"value1", 42}, {"value2", "hello"}}},
                 {"section2",
                  {{"subsection1", {{"value3", true}, {"value4", 3.2}}}, {"value5", 3}}}};

  // This shouldn't count as unused
  option["section2"]["value5"].attributes["source"] = "Output";

  MAYBE_UNUSED(auto value1) = option["section1"]["value1"].as<int>();
  MAYBE_UNUSED(auto value3) = option["section2"]["subsection1"]["value3"].as<bool>();

  Options expected_unused{
      {"section1", {{"value2", "hello"}}},
      {"section2", {{"subsection1", {{"value4", 3.2}}}}}};

  EXPECT_EQ(option.getUnused(), expected_unused);

  MAYBE_UNUSED(auto value2) = option["section1"]["value2"].as<std::string>();
  MAYBE_UNUSED(auto value4) = option["section2"]["subsection1"]["value4"].as<double>();
  MAYBE_UNUSED(auto value5) = option["section2"]["value5"].as<int>();

  Options expected_empty{};

  EXPECT_EQ(option.getUnused(), expected_empty);
}

TEST_F(OptionsTest, SetConditionallyUsed) {
  Options option{{"section1", {{"value1", 42}, {"value2", "hello"}}},
                 {"section2",
                  {{"subsection1", {{"value3", true}, {"value4", 3.2}}}, {"value5", 3}}}};

  option.setConditionallyUsed();

  Options expected_empty{};

  EXPECT_EQ(option.getUnused(), expected_empty);
}

TEST_F(OptionsTest, FuzzyFind) {
  Options option{{"value1", 21},
                 {"section1", {{"value1", 42}, {"value2", "hello"}, {"not this", 1}}},
                 {"section2",
                  {{"subsection1", {{"value3", true}, {"value4", 3.2}}}, {"value_5", 3}}}};

  auto fuzzy_matches = option.fuzzyFind("value1");
  EXPECT_EQ(fuzzy_matches.size(), 6);
  auto first_match = fuzzy_matches.begin();
  EXPECT_EQ(first_match->match.str(), "value1");
  EXPECT_EQ(first_match->distance, 0);
  auto second_match = ++first_match;
  EXPECT_EQ(second_match->match.str(), "section1:value1");
  EXPECT_EQ(second_match->distance, 1);

  auto fuzzy_section_matches = option.fuzzyFind("section1:subsection2:value__3");
  EXPECT_EQ(fuzzy_section_matches.size(), 1);
  auto first_section_match = fuzzy_section_matches.begin();
  EXPECT_EQ(first_section_match->match.str(), "section2:subsection1:value3");
  EXPECT_EQ(first_section_match->distance, 4);

  auto fuzzy_CAPS_matches = option.fuzzyFind("section2:VALUE_5");
  EXPECT_EQ(fuzzy_CAPS_matches.size(), 1);
  auto first_CAPS_match = fuzzy_CAPS_matches.begin();
  EXPECT_EQ(first_CAPS_match->match.str(), "section2:value_5");
  EXPECT_EQ(first_CAPS_match->distance, 1);
}

TEST_F(OptionsTest, GetFlattenedKeys) {
  Options option{
      {"value1", 21},
      {"section1", {{"value1", 42}, {"value2", "hello"}}},
      {"section2",
       {{"subsection1", {{"value3", true}, {"value4", 3.2}}}, {"value_5", 3}}}};

  auto flat_keys = option.getFlattenedKeys();

  std::vector<std::string> expected_keys{"value1",
                                         "section1:value1",
                                         "section1:value2",
                                         "section2:subsection1:value3",
                                         "section2:subsection1:value4",
                                         "section2:value_5"};

  std::sort(flat_keys.begin(), flat_keys.end());
  std::sort(expected_keys.begin(), expected_keys.end());

  EXPECT_EQ(flat_keys, expected_keys);
}

TEST_F(OptionsTest, CheckForUnusedOptions) {
  Options option{{"section1", {{"value1", 42}, {"value2", "hello"}}},
                 {"section2",
                  {{"subsection1", {{"value3", true}, {"value4", 3.2}}}, {"value5", 3}}}};

  // This shouldn't count as unused
  option["section2"]["value5"].attributes["source"] = "Output";

  MAYBE_UNUSED(auto value1) = option["section1"]["value1"].as<int>();
  MAYBE_UNUSED(auto value3) = option["section2"]["subsection1"]["value3"].as<bool>();

  EXPECT_THROW(bout::checkForUnusedOptions(option, "data", "BOUT.inp"), BoutException);
}

TEST_F(OptionsTest, CheckForUnusedOptionsGlobalRoot) {
  Options::root()["unused"] = 42;

  EXPECT_THROW(bout::checkForUnusedOptions(), BoutException);

  Options::root()["input"]["error_on_unused_options"] = false;
  EXPECT_NO_THROW(bout::checkForUnusedOptions());
}

class BoolTrueTestParametrized : public OptionsTest,
                                public ::testing::WithParamInterface<std::string> {};

TEST_P(BoolTrueTestParametrized, BoolTrueFromString) {
  std::string testval = GetParam();
  Options options;
  options["bool_key"] = testval;
  ASSERT_TRUE(options.isSet("bool_key"));
  ASSERT_TRUE(options["bool_key"].as<bool>());
}

INSTANTIATE_TEST_CASE_P(
    BoolTrueTests,
    BoolTrueTestParametrized,
    ::testing::Values(
      "y", "Y", "yes", "Yes", "yeS", "t", "true", "T", "True", "tRuE", "1"
    )
);

class BoolFalseTestParametrized : public OptionsTest,
                                  public ::testing::WithParamInterface<std::string> {};

TEST_P(BoolFalseTestParametrized, BoolFalseFromString) {
  std::string testval = GetParam();
  Options options;
  options["bool_key"] = testval;
  ASSERT_TRUE(options.isSet("bool_key"));
  ASSERT_FALSE(options["bool_key"].as<bool>());
}

INSTANTIATE_TEST_CASE_P(
    BoolFalseTests,
    BoolFalseTestParametrized,
    ::testing::Values(
      "n", "N", "no", "No", "nO", "f", "false", "F", "False", "fAlSe", "0"
    )
);

class BoolInvalidTestParametrized : public OptionsTest,
                                  public ::testing::WithParamInterface<std::string> {};

TEST_P(BoolInvalidTestParametrized, BoolInvalidFromString) {
  std::string testval = GetParam();
  Options options;
  options["bool_key"] = testval;
  ASSERT_TRUE(options.isSet("bool_key"));
  EXPECT_THROW(options["bool_key"].as<bool>(), BoutException);
}

INSTANTIATE_TEST_CASE_P(
    BoolInvalidTests,
    BoolInvalidTestParametrized,
    ::testing::Values(
      "a", "B", "yellow", "Yogi", "test", "truelong", "Tim", "2", "not", "No bool",
      "nOno", "falsebuttoolong", "-1"
    )
);
