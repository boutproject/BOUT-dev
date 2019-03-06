#include "gtest/gtest.h"

#include "boutcomm.hxx"
#include "cyclic_reduction.hxx"
#include "test_extras.hxx"
#include "bout/array.hxx"

#include <algorithm>
#include <vector>

namespace bout {
namespace testing {
constexpr int reduction_size{5};
constexpr BoutReal CyclicReduceTolerance{1.e-14};
} // namespace testing
} // namespace bout

Array<BoutReal> makeArrayFromVector(const std::vector<BoutReal>& values) {
  using std::begin;
  using std::end;

  Array<BoutReal> array{static_cast<Array<BoutReal>::size_type>(values.size())};
  std::copy(begin(values), end(values), begin(array));
  return array;
}

Matrix<BoutReal> makeMatrixFromVector(const std::vector<std::vector<BoutReal>>& values) {
  using std::begin;
  using std::end;

  Matrix<BoutReal> matrix{static_cast<Matrix<BoutReal>::size_type>(values.size()),
                          static_cast<Matrix<BoutReal>::size_type>(values[0].size())};
  auto start = begin(matrix);
  for (const auto& sub_values : values) {
    start = std::copy(begin(sub_values), end(sub_values), start);
  }
  return matrix;
}

TEST(CyclicReduction, SerialSolveSingleArray) {
  CyclicReduce<BoutReal> reduce{BoutComm::get(), bout::testing::reduction_size};

  auto a = makeArrayFromVector({0., 1., 1., 1., 1.});
  auto b = makeArrayFromVector({5., 4., 3., 2., 1.});
  auto c = makeArrayFromVector({2., 2., 2., 2., 0.});

  reduce.setCoefs(a, b, c);

  auto rhs = makeArrayFromVector({0., 1., 2., 2., 3.});
  Array<BoutReal> x{bout::testing::reduction_size};

  reduce.solve(rhs, x);

  EXPECT_NEAR(x[0], -1., bout::testing::CyclicReduceTolerance);
  EXPECT_NEAR(x[1], 2.5, bout::testing::CyclicReduceTolerance);
  EXPECT_NEAR(x[2], -4., bout::testing::CyclicReduceTolerance);
  EXPECT_NEAR(x[3], 5.75, bout::testing::CyclicReduceTolerance);
  EXPECT_NEAR(x[4], -2.75, bout::testing::CyclicReduceTolerance);
}

TEST(CyclicReduction, SerialSolveSingleMatrix) {
  CyclicReduce<BoutReal> reduce{BoutComm::get(), bout::testing::reduction_size};

  auto a = makeMatrixFromVector({{0., 1., 1., 1., 1.}});
  auto b = makeMatrixFromVector({{5., 4., 3., 2., 1.}});
  auto c = makeMatrixFromVector({{2., 2., 2., 2., 0.}});

  reduce.setCoefs(a, b, c);

  auto rhs = makeMatrixFromVector({{0., 1., 2., 2., 3.}});
  Matrix<BoutReal> x{1, bout::testing::reduction_size};

  reduce.solve(rhs, x);

  EXPECT_NEAR(x(0, 0), -1., bout::testing::CyclicReduceTolerance);
  EXPECT_NEAR(x(0, 1), 2.5, bout::testing::CyclicReduceTolerance);
  EXPECT_NEAR(x(0, 2), -4., bout::testing::CyclicReduceTolerance);
  EXPECT_NEAR(x(0, 3), 5.75, bout::testing::CyclicReduceTolerance);
  EXPECT_NEAR(x(0, 4), -2.75, bout::testing::CyclicReduceTolerance);
}

TEST(CyclicReduction, SerialSolveDoubleMatrix) {
  CyclicReduce<BoutReal> reduce{BoutComm::get(), bout::testing::reduction_size};

  auto a = makeMatrixFromVector({{0., 1., 1., 1., 1.}, {0., -2., -2., -2., -2.}});
  auto b = makeMatrixFromVector({{5., 4., 3., 2., 1.}, {1., 1., 1., 1., 1.}});
  auto c = makeMatrixFromVector({{2., 2., 2., 2., 0.}, {2., 2., 2., 2., 0.}});

  reduce.setCoefs(a, b, c);

  auto rhs = makeMatrixFromVector({{0., 1., 2., 2., 3.}, {5., 4., 5., 4., 5.}});
  Matrix<BoutReal> x{2, bout::testing::reduction_size};

  reduce.solve(rhs, x);

  EXPECT_NEAR(x(0, 0), -1., bout::testing::CyclicReduceTolerance);
  EXPECT_NEAR(x(0, 1), 2.5, bout::testing::CyclicReduceTolerance);
  EXPECT_NEAR(x(0, 2), -4., bout::testing::CyclicReduceTolerance);
  EXPECT_NEAR(x(0, 3), 5.75, bout::testing::CyclicReduceTolerance);
  EXPECT_NEAR(x(0, 4), -2.75, bout::testing::CyclicReduceTolerance);
  EXPECT_NEAR(x(1, 0), 3.4, bout::testing::CyclicReduceTolerance);
  EXPECT_NEAR(x(1, 1), 0.8, bout::testing::CyclicReduceTolerance);
  EXPECT_NEAR(x(1, 2), 5., bout::testing::CyclicReduceTolerance);
  EXPECT_NEAR(x(1, 3), 0.8, bout::testing::CyclicReduceTolerance);
  EXPECT_NEAR(x(1, 4), 6.6, bout::testing::CyclicReduceTolerance);
}
