#include "gtest/gtest.h"
#include "utils.hxx"

#include <string>

TEST(MatrixTest, CreateAndFree) {
  BoutReal **test_matrix = matrix<BoutReal>(5, 10);

  EXPECT_NE(nullptr, test_matrix);

  free_matrix(test_matrix);
}

TEST(NumberUtilitiesTest, SquareInt) {
  EXPECT_EQ(4, SQ(2));
  EXPECT_EQ(4, SQ(-2));
}

TEST(NumberUtilitiesTest, SquareReal) {
  EXPECT_DOUBLE_EQ(4.0, SQ(BoutReal(2.0)));
  EXPECT_DOUBLE_EQ(4.0, SQ(BoutReal(-2.0)));
}

TEST(NumberUtilitiesTest, Round) {
  EXPECT_EQ(3, ROUND(3.142));
  EXPECT_EQ(5, ROUND(4.566));
  EXPECT_EQ(-42, ROUND(-41.789));
  EXPECT_EQ(-42, ROUND(-42.123));
}

TEST(NumberUtilitiesTest, BoutMaxInt) {
  EXPECT_EQ(37, BOUTMAX(37));
  EXPECT_EQ(37, BOUTMAX(37, 12));
  EXPECT_EQ(37, BOUTMAX(-89, 37, 12));
}

TEST(NumberUtilitiesTest, BoutMaxReal) {
  EXPECT_DOUBLE_EQ(37, BOUTMAX(BoutReal(37)));
  EXPECT_DOUBLE_EQ(37, BOUTMAX(BoutReal(37), BoutReal(12)));
  EXPECT_DOUBLE_EQ(37, BOUTMAX(BoutReal(-89), BoutReal(37), BoutReal(12)));
}

TEST(NumberUtilitiesTest, BoutMinInt) {
  EXPECT_EQ(37, BOUTMIN(37));
  EXPECT_EQ(12, BOUTMIN(37, 12));
  EXPECT_EQ(-89, BOUTMIN(-89, 37, 12));
}

TEST(NumberUtilitiesTest, BoutMinReal) {
  EXPECT_DOUBLE_EQ(37, BOUTMIN(BoutReal(37)));
  EXPECT_DOUBLE_EQ(12, BOUTMIN(BoutReal(37), BoutReal(12)));
  EXPECT_DOUBLE_EQ(-89, BOUTMIN(BoutReal(-89), BoutReal(37), BoutReal(12)));
}

TEST(NumberUtilitiesTest, IsPow2) {
  EXPECT_TRUE(is_pow2(512));
  EXPECT_FALSE(is_pow2(887));
}

TEST(NumberUtilitiesTest, SignInt) {
  EXPECT_EQ(1, SIGN(5));
  EXPECT_EQ(1, SIGN(0));
  EXPECT_EQ(-1, SIGN(-5));
}

TEST(NumberUtilitiesTest, SignBoutReal) {
  EXPECT_EQ(1, SIGN(BoutReal(5.0)));
  EXPECT_EQ(1, SIGN(BoutReal(0.0)));
  EXPECT_EQ(-1, SIGN(BoutReal(-5.0)));
}

TEST(NumberUtilitiesTest, MinModInt) {
  EXPECT_EQ(5, MINMOD(5, 10));
  EXPECT_EQ(0, MINMOD(5, -10));
  EXPECT_EQ(0, MINMOD(-5, 10));
  EXPECT_EQ(-5, MINMOD(-5, -10));
  EXPECT_EQ(5, MINMOD(10, 5));
  EXPECT_EQ(0, MINMOD(10, -5));
  EXPECT_EQ(0, MINMOD(-10, 5));
  EXPECT_EQ(-5, MINMOD(-10, -5));
}

TEST(StringUtilitiesTest, CopyString) {
  const char hello[] = "Hello, world";
  char* copy = copy_string(hello);

  EXPECT_STREQ(hello, copy);

  free(copy);
}

TEST(StringUtilitiesTest, LowerCaseString) {
  std::string upper = "UPPERCASE";
  std::string lower = "uppercase";

  EXPECT_EQ(lower, lowercase(upper));
}

TEST(StringUtilitiesTest, LowerCaseStringWithSingleQuote) {
  std::string upper = "UPPER'CASE'";
  std::string lower = "upper'CASE'";

  EXPECT_EQ(lower, lowercasequote(upper));
}

TEST(StringUtilitiesTest, LowerCaseStringWithDoubleQuote) {
  std::string upper = "UPPER\"CASE\"";
  std::string lower = "upper\"CASE\"";

  EXPECT_EQ(lower, lowercasequote(upper));
}

TEST(StringUtilitiesTest, StringToReal) {
  std::string number_string = "0.3142e1";
  BoutReal number_real = 3.142e0;

  EXPECT_DOUBLE_EQ(number_real, stringToReal(number_string));
}

TEST(StringUtilitiesTest, StringToRealFail) {
  std::string number_string = "Not a number";

  EXPECT_THROW(stringToReal(number_string), BoutException);
}

TEST(StringUtilitiesTest, StringToInt) {
  std::string number_string = "42";
  int number_int = 42;

  EXPECT_DOUBLE_EQ(number_int, stringToInt(number_string));
}

TEST(StringUtilitiesTest, StringToIntFail) {
  std::string number_string = "Not a number";

  EXPECT_THROW(stringToInt(number_string), BoutException);
}

TEST(StringUtilitiesTest, IntToString) {
  int number_int = 42;
  std::string number_string = "42";

  EXPECT_EQ(number_string, toString(number_int));
}

TEST(StringUtilitiesTest, RealToString) {
  BoutReal number_real = 3.142e8;
  std::string number_string = "3.142e+08";

  EXPECT_EQ(number_string, toString(number_real));
}

TEST(StringUtilitiesTest, StringSplit) {
  std::string string_list = "a,b,c,d";
  std::list<std::string> split_list = {"a", "b", "c", "d"};

  EXPECT_EQ(split_list, strsplit(string_list, ','));
}

TEST(StringUtilitiesTest, StringTrim) {
  std::string input = "    space    ";

  EXPECT_EQ("space", trim(input));
}

TEST(StringUtilitiesTest, StringTrimRight) {
  std::string input = "    space    ";

  EXPECT_EQ("    space", trimRight(input));
}

TEST(StringUtilitiesTest, StringTrimLeft) {
  std::string input = "    space    ";

  EXPECT_EQ("space    ", trimLeft(input));
}

TEST(StringUtilitiesTest, StringTrimComments) {
  std::string input = "space  # comment";

  EXPECT_EQ("space  ", trimComments(input, "#"));
}
