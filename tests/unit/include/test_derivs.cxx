#include "gtest/gtest.h"

#include "bout_types.hxx"
#include "fft.hxx"
#include "field3d.hxx"
#include "test_extras.hxx"
#include "bout/constants.hxx"
#include "bout/deriv_store.hxx"
#include "bout/paralleltransform.hxx"

#include <algorithm>
#include <string>
#include <tuple>
#include <vector>

namespace {
// These are merely sanity checks, so don't expect them to agree very well!
constexpr BoutReal derivatives_tolerance{1.e-2};
}

class DerivativesTest
    : public ::testing::TestWithParam<std::pair<DIRECTION, std::string>> {
public:
  DerivativesTest() : input{mesh}, first_order_expected{mesh} {

    // Make sure fft functions are quiet by setting fft_measure to false
    bout::fft::fft_init(false);

    using Index = Field3D::ind_type;

    // Pointer to index method that converts to index-space
    using DirectionFunction = int (Index::*)() const;
    DirectionFunction dir;

    // Number of guard cells in (x, y). The "other" direction will
    // have none
    int x_guards{0};
    int y_guards{0};

    // This must be a balance between getting any kind of accuracy and
    // each derivative running in ~1ms or less
    constexpr int grid_size{64};
    const BoutReal box_length{TWOPI / grid_size};

    switch (std::get<0>(GetParam())) {
    case DIRECTION::X:
      nx = grid_size;
      dir = &Index::x;
      x_guards = 2;
      region = RGN_NOX;
      break;
    case DIRECTION::Y:
      ny = grid_size;
      dir = &Index::y;
      y_guards = 2;
      region = RGN_NOY;
      break;
    case DIRECTION::Z:
      nz = grid_size;
      dir = &Index::z;
      region = RGN_ALL;
      break;
    default:
      throw BoutException("bad direction");
    }

    if (mesh != nullptr) {
      delete mesh;
      mesh = nullptr;
    }

    mesh = new FakeMesh(nx, ny, nz);

    mesh->xstart = x_guards;
    mesh->xend = nx - (x_guards + 1);
    mesh->ystart = y_guards;
    mesh->yend = ny - (y_guards + 1);

    output_info.disable();
    mesh->createDefaultRegions();
    output_info.enable();

    Field3D input_{mesh};
    input_.allocate();

    for (auto i : input_) {
      input_[i] = std::sin((i.*dir)() * box_length);
    }

    input = input_;

    // We need the parallel slices for the y-direction
    ParallelTransformIdentity identity{};
    identity.calcYUpDown(input);

    Field3D first_order_expected_{mesh};
    first_order_expected_.allocate();

    for (auto i : first_order_expected_) {
      first_order_expected_[i] = std::cos((i.*dir)() * box_length) * box_length;
    }

    first_order_expected = first_order_expected_;

    Field3D second_order_expected_{mesh};
    second_order_expected_.allocate();

    for (auto i : second_order_expected_) {
      second_order_expected_[i] = -std::sin((i.*dir)() * box_length) * pow(box_length, 2);
    }

    second_order_expected = second_order_expected_;

    DerivativeStore<Field3D>::getInstance().initialise(Options::getRoot());
  };

  Field3D input;
  Field3D first_order_expected;
  Field3D second_order_expected;

  int nx{3};
  int ny{3};
  int nz{2};

  REGION region;
};

// Get all the available methods for this direction and turn it from a
// collection of strings to a collection of pairs of the direction and
// strings so that the test fixture knows which direction we're using
auto getMethodsForDirection(DERIV derivative_order, DIRECTION direction)
    -> std::vector<std::pair<DIRECTION, std::string>> {

  auto available_methods = DerivativeStore<Field3D>::getInstance().getAvailableMethods(
      derivative_order, direction);

  // Method names paired with the current direction
  std::vector<std::pair<DIRECTION, std::string>> methods{};

  std::transform(std::begin(available_methods), std::end(available_methods),
                 std::back_inserter(methods), [&direction](std::string method) {
                   return std::make_pair(direction, method);
                 });

  return methods;
};

// Returns the method out of the direction/method pair for printing
// the test name
auto methodDirectionPairToString(
    const ::testing::TestParamInfo<std::pair<DIRECTION, std::string>>& param)
    -> std::string {
  return std::get<1>(param.param);
}

using FirstDerivatives = DerivativesTest;

// Instantiate the test for X, Y, Z for first derivatives
INSTANTIATE_TEST_CASE_P(X, FirstDerivatives,
                        ::testing::ValuesIn(getMethodsForDirection(DERIV::Standard,
                                                                   DIRECTION::X)),
                        methodDirectionPairToString);

INSTANTIATE_TEST_CASE_P(Y, FirstDerivatives,
                        ::testing::ValuesIn(getMethodsForDirection(DERIV::Standard,
                                                                   DIRECTION::Y)),
                        methodDirectionPairToString);

INSTANTIATE_TEST_CASE_P(Z, FirstDerivatives,
                        ::testing::ValuesIn(getMethodsForDirection(DERIV::Standard,
                                                                   DIRECTION::Z)),
                        methodDirectionPairToString);

// The actual first derivative test!
TEST_P(FirstDerivatives, FirstOrder) {
  auto derivative = DerivativeStore<Field3D>::getInstance().getStandardDerivative(
      std::get<1>(GetParam()), std::get<0>(GetParam()));

  Field3D result{mesh};
  result.allocate();
  derivative(input, result, region);

  EXPECT_TRUE(IsField3DEqualField3D(result, first_order_expected, "RGN_NOBNDRY",
                                    derivatives_tolerance));
}

using SecondDerivatives = DerivativesTest;

// Instantiate the test for X, Y, Z for second derivatives
INSTANTIATE_TEST_CASE_P(X, SecondDerivatives,
                        ::testing::ValuesIn(getMethodsForDirection(DERIV::StandardSecond,
                                                                   DIRECTION::X)),
                        methodDirectionPairToString);

INSTANTIATE_TEST_CASE_P(Y, SecondDerivatives,
                        ::testing::ValuesIn(getMethodsForDirection(DERIV::StandardSecond,
                                                                   DIRECTION::Y)),
                        methodDirectionPairToString);

INSTANTIATE_TEST_CASE_P(Z, SecondDerivatives,
                        ::testing::ValuesIn(getMethodsForDirection(DERIV::StandardSecond,
                                                                   DIRECTION::Z)),
                        methodDirectionPairToString);

// The actual second derivative test!
TEST_P(SecondDerivatives, SecondOrder) {
  auto derivative = DerivativeStore<Field3D>::getInstance().getStandard2ndDerivative(
      std::get<1>(GetParam()), std::get<0>(GetParam()));

  Field3D result{mesh};
  result.allocate();
  derivative(input, result, region);

  EXPECT_TRUE(IsField3DEqualField3D(result, second_order_expected, "RGN_NOBNDRY",
                                    derivatives_tolerance));
}
