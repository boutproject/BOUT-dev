#include <gtest/gtest.h>

#include "bout/bout_enum_class.hxx"
#include "bout/boutexception.hxx"
#include "bout/options.hxx"
#include "bout/output.hxx"

BOUT_ENUM_CLASS(TestEnum, foo, bar);
BOUT_ENUM_CLASS_NS(bout, TestEnum, foo, bar);
BOUT_ENUM_CLASS_NS(bout::test, TestEnum, foo, bar);

TEST(BoutEnumClassNS, toString) {
  EXPECT_EQ(toString(bout::TestEnum::foo), "foo");
  EXPECT_EQ(toString(bout::TestEnum::bar), "bar");
}

TEST(BoutEnumClass, toString) {
  EXPECT_EQ(toString(TestEnum::foo), "foo");
  EXPECT_EQ(toString(TestEnum::bar), "bar");
}

TEST(BoutEnumClassNS, fromString) {
  EXPECT_EQ(bout::TestEnumFromString("foo"), bout::TestEnum::foo);
  EXPECT_EQ(bout::TestEnumFromString("bar"), bout::TestEnum::bar);
  EXPECT_THROW(bout::TestEnumFromString("expect_fail"), BoutException);
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

TEST(BoutEnumClassNS, options) {
  WithQuietOutput quiet_info{output_info};

  Options options;

  auto opt1 = options["opt"].withDefault(bout::TestEnum::foo);
  EXPECT_EQ(opt1, bout::TestEnum::foo);
  EXPECT_NE(opt1, bout::TestEnum::bar);

  options["opt"] = "bar";

  auto opt2 = options["opt"].as<bout::TestEnum>();
  EXPECT_EQ(opt2, bout::TestEnum::bar);
  EXPECT_NE(opt2, bout::TestEnum::foo);

  options["optfail"] = "expect_fail";

  EXPECT_THROW(options["optfail"].as<bout::TestEnum>(), BoutException);
}

TEST(BoutEnumClassNS2, options) {
  WithQuietOutput quiet_info{output_info};

  Options options;

  auto opt1 = options["opt"].withDefault(bout::test::TestEnum::foo);
  EXPECT_EQ(opt1, bout::test::TestEnum::foo);
  EXPECT_NE(opt1, bout::test::TestEnum::bar);

  options["opt"] = "bar";

  auto opt2 = options["opt"].as<bout::test::TestEnum>();
  EXPECT_EQ(opt2, bout::test::TestEnum::bar);
  EXPECT_NE(opt2, bout::test::TestEnum::foo);

  options["optfail"] = "expect_fail";

  EXPECT_THROW(options["optfail"].as<bout::test::TestEnum>(), BoutException);
}

TEST(BoutEnumClass, ostream) {
  auto sstream = std::stringstream();

  sstream << TestEnum::foo << TestEnum::bar;

  EXPECT_EQ(sstream.str(), "foobar");
}

TEST(BoutEnumClassNS, ostream) {
  auto sstream = std::stringstream();

  sstream << bout::TestEnum::foo << bout::TestEnum::bar;

  EXPECT_EQ(sstream.str(), "foobar");
}
