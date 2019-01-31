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

// Some basic sanity checks for the derivative kernels. Checks the
// derivatives of sin(R) where R = {X, Y, Z} for each R
// individually. To make this as fast as possible, we use only a
// couple of points in the non-tested directions -- not just one
// though, as this allows us to check we're not introducing spurious
// variation in the other directions.
//
// This is one of the more complicated uses of googletest! We need to
// all combinations of methods and directions. Unfortunately, Z has
// one more method than X and Y -- FFT. This means we can't just use
// the provided `Combine` to produce the Cartesian product of
// directions and methods as the latter depends on the
// former. Instead, we instantiate the tests separately for each
// direction and test all the methods for that direction in that
// instantiation.

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

    // Set all the variables for this direction
    // In C++14 this can be the more explicit std::get<DIRECTION>()
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

    // Make the input and expected output fields
    // Weird `(i.*dir)()` syntax here in order to call the direction method
    // C++17 makes this nicer with std::invoke
    input = makeField<Field3D>([&](Index& i) { return std::sin((i.*dir)() * box_length); }, mesh);

    first_order_expected = makeField<Field3D>(
        [&](Index& i) { return std::cos((i.*dir)() * box_length) * box_length; }, mesh);

    second_order_expected = makeField<Field3D>(
        [&](Index& i) { return -std::sin((i.*dir)() * box_length) * pow(box_length, 2); },
        mesh);

    // We need the parallel slices for the y-direction
    ParallelTransformIdentity identity{};
    identity.calcYUpDown(input);

    // FIXME: remove when defaults are set in the DerivativeStore ctor
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

// Use an alias to distinguish the first and second derivative tests
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

// Use an alias to distinguish the first and second derivative tests
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
