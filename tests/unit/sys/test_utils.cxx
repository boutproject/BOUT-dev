#include "gtest/gtest.h"
#include "utils.hxx"

#include <string>

TEST(MatrixTest, DefaultShape) {
  Matrix<int> matrix;

  int shape0, shape1;
  std::tie(shape0, shape1) = matrix.shape();
  EXPECT_EQ(shape0, 0);
  EXPECT_EQ(shape1, 0);
}

TEST(MatrixTest, CreateGivenSize) {
  Matrix<int> matrix(3, 5);

  int shape0, shape1;
  std::tie(shape0, shape1) = matrix.shape();
  EXPECT_EQ(shape0, 3);
  EXPECT_EQ(shape1, 5);
}

TEST(MatrixTest, Reallocate) {
  Matrix<int> matrix{};

  ASSERT_TRUE(matrix.empty());

  matrix.reallocate(3, 5);

  int shape0, shape1;
  std::tie(shape0, shape1) = matrix.shape();

  EXPECT_EQ(shape0, 3);
  EXPECT_EQ(shape1, 5);

  matrix.reallocate(5, 3);
  std::tie(shape0, shape1) = matrix.shape();

  EXPECT_EQ(shape0, 5);
  EXPECT_EQ(shape1, 3);
}

TEST(MatrixTest, Empty) {
  Matrix<int> matrix;
  EXPECT_TRUE(matrix.empty());

  Matrix<int> matrix2(3, 5);
  EXPECT_FALSE(matrix2.empty());
}

TEST(MatrixTest, CopyConstuctor) {
  Matrix<int> matrix(3, 5);
  matrix = 0;
  Matrix<int> matrix2(matrix);

  int shape0, shape1;
  std::tie(shape0, shape1) = matrix2.shape();
  EXPECT_EQ(shape0, 3);
  EXPECT_EQ(shape1, 5);

  // Now check that matrix and matrix2 are unique
  for (const auto i : matrix2) {
    EXPECT_EQ(i, 0);
  }

  matrix = 2;

  for (const auto i : matrix2) {
    EXPECT_EQ(i, 0);
  }
}

TEST(MatrixTest, CopyAssignment) {
  Matrix<int> matrix(3, 5);
  matrix = 0;

  Matrix<int> matrix2;
  ASSERT_TRUE(matrix2.empty());

  matrix2 = matrix;

  int shape0, shape1;
  std::tie(shape0, shape1) = matrix2.shape();
  EXPECT_EQ(shape0, 3);
  EXPECT_EQ(shape1, 5);

  // Now check that matrix and matrix2 are unique
  for (const auto i : matrix2) {
    EXPECT_EQ(i, 0);
  }

  matrix = 2;

  for (const auto i : matrix2) {
    EXPECT_EQ(i, 0);
  }
}

TEST(MatrixTest, Iterator) {
  Matrix<int> matrix(3, 5);
  int count = 0;

  for (auto iter = std::begin(matrix); iter < std::end(matrix); ++iter) {
    ++count;
  }

  EXPECT_EQ(count, 15);
}

TEST(MatrixTest, ConstIterator) {
  const Matrix<int> matrix(3, 5);
  int count = 0;

  for (auto iter = std::begin(matrix); iter < std::end(matrix); ++iter) {
    ++count;
  }

  EXPECT_EQ(count, 15);
}

TEST(MatrixTest, AssignmentValue) {
  Matrix<int> matrix(3, 5);
  int count = 0;

  matrix = 2;

  for (auto iter = std::begin(matrix); iter < std::end(matrix); ++iter) {
    count += *iter;
  }

  EXPECT_EQ(count, 30);
}

TEST(MatrixTest, Indexing) {
  Matrix<int> matrix(3, 5);
  matrix = 3;

  EXPECT_EQ(matrix(1, 1), 3);

  matrix(1, 1) = 4;

  EXPECT_EQ(matrix(1, 1), 4);
}

#if CHECK > 1
TEST(MatrixTest, OutOfBoundsIndexing) {
  Matrix<int> matrix(3, 5);
  matrix = 3;

  EXPECT_THROW(matrix(-1, 0), BoutException);
  EXPECT_THROW(matrix(0, -1), BoutException);
  EXPECT_THROW(matrix(10, 0), BoutException);
  EXPECT_THROW(matrix(0, 10), BoutException);
}
#endif

TEST(MatrixTest, ConstIndexing) {
  Matrix<int> matrix(3, 5);
  matrix = 3;
  const Matrix<int> matrix2(matrix);

  EXPECT_EQ(matrix2(1, 1), 3);
}

TEST(MatrixTest, GetData) {
  Matrix<int> matrix(3, 5);
  matrix = 3;

  auto data = matrix.getData();

  EXPECT_TRUE(std::all_of(std::begin(data), std::end(data), [](int a) { return a == 3; }));

  data[0] = 4;

  EXPECT_EQ(matrix(0, 0), 4);
}

TEST(MatrixTest, ConstGetData) {
  Matrix<int> matrix(3, 5);
  matrix = 3;

  const auto data = matrix.getData();

  EXPECT_TRUE(std::all_of(std::begin(data), std::end(data), [](int a) { return a == 3; }));
  EXPECT_TRUE(std::all_of(std::begin(matrix), std::end(matrix), [](int a) { return a == 3; }));
}

TEST(TensorTest, DefaultShape) {
  Tensor<int> tensor;

  int shape0, shape1, shape2;
  std::tie(shape0, shape1, shape2) = tensor.shape();
  EXPECT_EQ(shape0, 0);
  EXPECT_EQ(shape1, 0);
  EXPECT_EQ(shape2, 0);
}

TEST(TensorTest, CreateGivenSize) {
  Tensor<int> tensor(3, 5, 7);

  int shape0, shape1, shape2;
  std::tie(shape0, shape1, shape2) = tensor.shape();
  EXPECT_EQ(shape0, 3);
  EXPECT_EQ(shape1, 5);
  EXPECT_EQ(shape2, 7);
}

TEST(TensorTest, Reallocate) {
  Tensor<int> tensor{};

  ASSERT_TRUE(tensor.empty());

  tensor.reallocate(3, 5, 7);

  int shape0, shape1, shape2;
  std::tie(shape0, shape1, shape2) = tensor.shape();

  EXPECT_EQ(shape0, 3);
  EXPECT_EQ(shape1, 5);
  EXPECT_EQ(shape2, 7);

  tensor.reallocate(5, 7, 3);
  std::tie(shape0, shape1, shape2) = tensor.shape();

  EXPECT_EQ(shape0, 5);
  EXPECT_EQ(shape1, 7);
  EXPECT_EQ(shape2, 3);
}

TEST(TensorTest, Empty) {
  Tensor<int> tensor;
  EXPECT_TRUE(tensor.empty());

  Tensor<int> tensor2(3, 5, 7);
  EXPECT_FALSE(tensor2.empty());
}

TEST(TensorTest, CopyConstuctor) {
  Tensor<int> tensor(3, 5, 7);
  tensor = 0;

  Tensor<int> tensor2(tensor);

  int shape0, shape1, shape2;
  std::tie(shape0, shape1, shape2) = tensor2.shape();
  EXPECT_EQ(shape0, 3);
  EXPECT_EQ(shape1, 5);
  EXPECT_EQ(shape2, 7);

  // Now check that matrix and matrix2 are unique
  for (const auto i : tensor2) {
    EXPECT_EQ(i, 0);
  }

  tensor = 2;

  for (const auto i : tensor2) {
    EXPECT_EQ(i, 0);
  }
}

TEST(TensorTest, CopyAssignment) {
  Tensor<int> tensor(3, 5, 7);
  tensor = 0;

  Tensor<int> tensor2;

  ASSERT_TRUE(tensor2.empty());

  tensor2 = tensor;

  int shape0, shape1, shape2;
  std::tie(shape0, shape1, shape2) = tensor2.shape();
  EXPECT_EQ(shape0, 3);
  EXPECT_EQ(shape1, 5);
  EXPECT_EQ(shape2, 7);

  // Now check that matrix and matrix2 are unique
  for (const auto i : tensor2) {
    EXPECT_EQ(i, 0);
  }

  tensor = 2;

  for (const auto i : tensor2) {
    EXPECT_EQ(i, 0);
  }
}

TEST(TensorTest, Iterator) {
  Tensor<int> tensor(3, 5, 7);
  int count = 0;

  for (auto iter = std::begin(tensor); iter < std::end(tensor); ++iter) {
    ++count;
  }

  EXPECT_EQ(count, 105);
}

TEST(TensorTest, ConstIterator) {
  const Tensor<int> tensor(3, 5, 7);
  int count = 0;

  for (auto iter = std::begin(tensor); iter < std::end(tensor); ++iter) {
    ++count;
  }

  EXPECT_EQ(count, 105);
}

TEST(TensorTest, AssignmentValue) {
  Tensor<int> tensor(3, 5, 7);
  int count = 0;

  tensor = 2;

  for (auto iter = std::begin(tensor); iter < std::end(tensor); ++iter) {
    count += *iter;
  }

  EXPECT_EQ(count, 210);
}

TEST(TensorTest, Indexing) {
  Tensor<int> tensor(3, 5, 7);
  tensor = 3;

  EXPECT_EQ(tensor(1, 1, 1), 3);

  tensor(1, 1, 1) = 4;

  EXPECT_EQ(tensor(1, 1, 1), 4);
}

#if CHECK > 1
TEST(TensorTest, OutOfBoundsIndexing) {
  Tensor<int> tensor(3, 5, 7);
  tensor = 3;

  EXPECT_THROW(tensor(-1, 0, 0), BoutException);
  EXPECT_THROW(tensor(0, -1, 0), BoutException);
  EXPECT_THROW(tensor(0, 0, -1), BoutException);
  EXPECT_THROW(tensor(10, 0, 0), BoutException);
  EXPECT_THROW(tensor(0, 10, 0), BoutException);
  EXPECT_THROW(tensor(0, 0, 10), BoutException);
}
#endif

TEST(TensorTest, ConstIndexing) {
  Tensor<int> tensor(3, 5, 7);
  tensor = 3;
  const Tensor<int> tensor2(tensor);

  EXPECT_EQ(tensor2(1, 1, 1), 3);
}

TEST(TensorTest, GetData) {
  Tensor<int> tensor(3, 5, 7);
  tensor = 3;

  auto data = tensor.getData();

  EXPECT_TRUE(std::all_of(std::begin(data), std::end(data), [](int a) { return a == 3; }));

  data[0] = 4;

  EXPECT_EQ(tensor(0, 0, 0), 4);
}

TEST(TensorTest, ConstGetData) {
  Tensor<int> tensor(3, 5, 7);
  tensor = 3;

  const auto data = tensor.getData();

  EXPECT_TRUE(std::all_of(std::begin(data), std::end(data), [](int a) { return a == 3; }));
  EXPECT_TRUE(std::all_of(std::begin(tensor), std::end(tensor), [](int a) { return a == 3; }));
}

TEST(Invert3x3Test, Identity) {
  Matrix<BoutReal> input(3, 3);
  input = 0;
  for (int i = 0; i < 3; i++) {
    input(i, i) = 1.0;
  }
  auto expected = input;
  invert3x3(input);

  for (int j = 0; j < 3; j++) {
    for (int i = 0; i < 3; i++) {
      EXPECT_EQ(input(i, j), expected(i, j));
    }
  }
}

TEST(Invert3x3Test, InvertTwice) {
  std::vector<BoutReal> rawDataMat = {0.05567105, 0.92458227, 0.19954631,
                                      0.28581972, 0.54009039, 0.13234403,
                                      0.8841194,  0.161224,   0.74853209};
  std::vector<BoutReal> rawDataInv = {-2.48021781, 4.27410022,  -0.09449605,
                                      0.6278449,   0.87275842,  -0.32168092,
                                      2.79424897,  -5.23628123, 1.51684677};

  Matrix<BoutReal> input(3, 3);
  Matrix<BoutReal> expected(3, 3);

  int counter = 0;
  for (int j = 0; j < 3; j++) {
    for (int i = 0; i < 3; i++) {
      input(i, j) = rawDataMat[counter];
      expected(i, j) = rawDataInv[counter];
      counter++;
    }
  }

  // Invert twice to check if we get back to where we started
  invert3x3(input);

  for (int j = 0; j < 3; j++) {
    for (int i = 0; i < 3; i++) {
      // Note we only check to single tolerance here
      EXPECT_FLOAT_EQ(input(i, j), expected(i, j));
    }
  }
}

TEST(Invert3x3Test, Singular) {
  Matrix<BoutReal> input(3, 3);
  input = 0;
  EXPECT_THROW(invert3x3(input), BoutException);
}

TEST(Invert3x3Test, BadCondition) {
  Matrix<BoutReal> input(3, 3);

  // Default small
  input = 0.;
  input(0, 0) = 1.0e-16;
  input(1, 1) = 1.0;
  input(2, 2) = 1.0;
  EXPECT_THROW(invert3x3(input), BoutException);

  // Default small -- not quite bad enough condition
  input = 0.;
  input(0, 0) = 1.0e-12;
  input(1, 1) = 1.0;
  input(2, 2) = 1.0;
  EXPECT_NO_THROW(invert3x3(input));

  // Non-default small
  input = 0.;
  input(0, 0) = 1.0e-12;
  input(1, 1) = 1.0;
  input(2, 2) = 1.0;
  EXPECT_THROW(invert3x3(input, 1.0e-10), BoutException);

  // Non-default small
  input = 0.;
  input(0, 0) = 1.0e-12;
  input(1, 1) = 1.0; input(2, 2) = 1.0;
  EXPECT_NO_THROW(invert3x3(input, -1.0e-10));
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

#if CHECK > 0
TEST(NumberUtilitiesTest, CheckDataGood) {
  EXPECT_NO_THROW(checkData(5.0));
}

TEST(NumberUtilitiesTest, CheckDataBad) {
  EXPECT_THROW(checkData(nan("")), BoutException);
}
#else
TEST(NumberUtilitiesTest, CheckDataGoodDisabled) {
  EXPECT_NO_THROW(checkData(5.0));
}

TEST(NumberUtilitiesTest, CheckDataBadDisabled) {
  EXPECT_NO_THROW(checkData(nan("")));
}
#endif

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

TEST(StringUtilitiesTest, BoolToString) {
  std::string true_string = "true";
  std::string false_string = "false";

  EXPECT_EQ(true_string, toString(true));
  EXPECT_EQ(false_string, toString(false));
}

TEST(StringUtilitiesTest, ConstCharToString) {
  EXPECT_EQ(std::string("hello"), toString("hello"));
}

TEST(StringUtilitiesTest, StringToString) {
  std::string test_string = "dlkjl872kj";

  EXPECT_EQ( test_string, toString(test_string) );
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

namespace {
using function_typedef = int(*)(char, int, double);
int function_pointer(double, char) {return 0;};
template <typename T>
void function_template(T) {}
} // namespace

TEST(FunctionTraitsTest, ResultType) {
  using bout::utils::function_traits;
  static_assert(std::is_same<function_traits<function_typedef>::result_type, int>::value,
                "Wrong result_type for function_traits of a typedef");
  static_assert(
      std::is_same<function_traits<decltype(&function_pointer)>::result_type, int>::value,
      "Wrong result_type for function_traits of a function pointer");
  static_assert(
      std::is_same<function_traits<decltype(&function_template<int>)>::result_type,
                   void>::value,
      "Wrong result_type for function_traits of a template function");

  // Use foo to suppress warning
  function_pointer(0, 0);
}

TEST(FunctionTraitsTest, NumberOfArgs) {
  using bout::utils::function_traits;
  static_assert(function_traits<function_typedef>::nargs == 3,
                "Wrong number of arguments for function_traits of a typedef>");
  static_assert(function_traits<decltype(&function_pointer)>::nargs == 2,
                "Wrong number of arguments for function_traits of a function pointer");
  static_assert(function_traits<decltype(&function_template<int>)>::nargs == 1,
                "Wrong number of arguments for function_traits of a template function");
}

TEST(FunctionTraitsTest, FirstArg) {
  using bout::utils::function_traits;
  static_assert(
      std::is_same<function_traits<function_typedef>::arg<0>::type, char>::value,
      "Wrong first argument type for function_traits of a typedef");
  static_assert(std::is_same<function_traits<decltype(&function_pointer)>::arg<0>::type,
                             double>::value,
                "Wrong first argument type for function_traits of a function pointer");
  static_assert(
      std::is_same<function_traits<decltype(&function_template<int>)>::arg<0>::type,
                   int>::value,
      "Wrong first argument type for function_traits of a template function");

  static_assert(std::is_same<function_traits<function_typedef>::arg_t<0>, char>::value,
                "Wrong first argument type for function_traits of a typedef using arg_t");
  static_assert(
      std::is_same<function_traits<decltype(&function_pointer)>::arg_t<0>, double>::value,
      "Wrong first argument type for function_traits of a function pointer using arg_t");
  static_assert(
      std::is_same<function_traits<decltype(&function_template<int>)>::arg_t<0>,
                   int>::value,
      "Wrong first argument type for function_traits of a template function using arg_t");
}

TEST(FunctionTraitsTest, SecondArg) {
  using bout::utils::function_traits;
  static_assert(std::is_same<function_traits<function_typedef>::arg<1>::type, int>::value,
                "Wrong second argument type for function_traits of a typedef");
  static_assert(
      std::is_same<function_traits<function_typedef>::arg_t<1>, int>::value,
      "Wrong second argument type for function_traits of a typedef using arg_t");
}
