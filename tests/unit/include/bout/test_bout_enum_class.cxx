#include <gtest/gtest.h>

#include "bout/bout_enum_class.hxx"
#include "bout/boutexception.hxx"
#include "bout/options.hxx"
#include "bout/output.hxx"

BOUT_ENUM_CLASS(TestEnum, foo, bar);

TEST(BoutEnumClass, toString) {
  EXPECT_EQ(toString(TestEnum::foo), "foo");
  EXPECT_EQ(toString(TestEnum::bar), "bar");
}

TEST(BoutEnumClass, fromString) {
  EXPECT_EQ(TestEnumFromString("foo"), TestEnum::foo);
  EXPECT_EQ(TestEnumFromString("bar"), TestEnum::bar);
  EXPECT_THROW(TestEnumFromString("expect_fail"), BoutException);
}

TEST(BoutEnumClass, options) {
  WithQuietOutput quiet_info{output_info};

  Options options;

  auto opt1 = options["opt"].withDefault(TestEnum::foo);
  EXPECT_EQ(opt1, TestEnum::foo);
  EXPECT_NE(opt1, TestEnum::bar);

  options["opt"] = "bar";

  auto opt2 = options["opt"].as<TestEnum>();
  EXPECT_EQ(opt2, TestEnum::bar);
  EXPECT_NE(opt2, TestEnum::foo);

  options["optfail"] = "expect_fail";

  EXPECT_THROW(options["optfail"].as<TestEnum>(), BoutException);
}

TEST(BoutEnumClass, ostream) {
  auto sstream = std::stringstream();

  sstream << TestEnum::foo << TestEnum::bar;

  EXPECT_EQ(sstream.str(), "foobar");
}
