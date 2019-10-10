#include "gtest/gtest.h"

#include "bout_types.hxx"
#include "fft.hxx"
#include "field3d.hxx"
#include "test_extras.hxx"
#include "bout/constants.hxx"
#include "bout/deriv_store.hxx"
#include "bout/index_derivs_interface.hxx"
#include "bout/paralleltransform.hxx"

#include <algorithm>
#include <string>
#include <tuple>
#include <vector>

// The unit tests use the global mesh
using namespace bout::globals;

// Some basic sanity checks for the derivative kernels. Checks the
// derivatives of sin(R) where R = {X, Y, Z} for each R
// individually. To make this as fast as possible, we use only a
// couple of points in the non-tested directions -- not just one
// though, as this allows us to check we're not introducing spurious
// variation in the other directions.
//
// This is one of the more complicated uses of googletest! We need to
// test all combinations of methods, directions and (standard)
// derivative types. Unfortunately, Z has one more method than X and Y
// (FFT), and the different orders have different sets of
// methods. This means we can't just use the provided `Combine` to
// produce the Cartesian product of methods, directions, and types as
// the first depends on the second two. Instead, we instantiate the
// tests separately for each direction and order and test all the
// methods for that direction in that instantiation.

namespace {
// These are merely sanity checks, so don't expect them to agree very
// well! This has to be sufficiently loose for the least accurate
// method to pass, or we need to also match up tolerances to methods
constexpr BoutReal derivatives_tolerance{5.e-3};
}

class DerivativesTest
    : public ::testing::TestWithParam<std::tuple<DIRECTION, DERIV, std::string>> {
public:
  DerivativesTest() : input{mesh}, expected{mesh} {
    WithQuietOutput quiet_info{output_info};
    WithQuietOutput quiet_warn{output_warn};

    using Index = Field3D::ind_type;

    // Pointer to index method that converts single-index to 3-index space
    using DirectionFunction = int (Index::*)() const;
    DirectionFunction dir;

    // Number of guard cells in (x, y). The "other" direction will
    // have none
    int x_guards{0};
    int y_guards{0};

    // Grid sizes
    int nx{3};
    int ny{3};
    int nz{2};

    // This must be a balance between getting any kind of accuracy and
    // each derivative running in ~1ms or less
    constexpr int grid_size{128};
    const BoutReal box_length{TWOPI / grid_size};

    // Set all the variables for this direction
    // In C++14 this can be the more explicit std::get<DIRECTION>()
    switch (std::get<0>(GetParam())) {
    case DIRECTION::X:
      nx = grid_size;
      dir = &Index::x;
      x_guards = 2;
      region = "RGN_NOX";
      break;
    case DIRECTION::Y:
      ny = grid_size;
      dir = &Index::y;
      y_guards = 2;
      region = "RGN_NOY";
      break;
    case DIRECTION::Z:
      nz = grid_size;
      dir = &Index::z;
      region = "RGN_ALL";
      break;
    default:
      throw BoutException("bad direction");
    }

    mesh = new FakeMesh(nx, ny, nz);
    static_cast<FakeMesh*>(mesh)->setCoordinates(nullptr);

    mesh->xstart = x_guards;
    mesh->xend = nx - (x_guards + 1);
    mesh->ystart = y_guards;
    mesh->yend = ny - (y_guards + 1);

    mesh->createDefaultRegions();

    // Make the input and expected output fields
    // Weird `(i.*dir)()` syntax here in order to call the direction method
    // C++17 makes this nicer with std::invoke
    input = makeField<Field3D>([&](Index& i) { return std::sin((i.*dir)() * box_length); }, mesh);

    // Make the velocity field
    velocity = makeField<Field3D>([&](Index& UNUSED(i)) { return 2.0; }, mesh);

    // Get the expected result for this order of derivative
    // Again, could be nicer in C++17 with std::get<DERIV>(GetParam())
    switch (std::get<1>(GetParam())) {
    case DERIV::Standard:
      expected = makeField<Field3D>(
        [&](Index& i) { return std::cos((i.*dir)() * box_length) * box_length; }, mesh);
      break;
    case DERIV::StandardSecond:
      expected = makeField<Field3D>(
        [&](Index& i) { return -std::sin((i.*dir)() * box_length) * pow(box_length, 2); },
        mesh);
      break;
    case DERIV::StandardFourth:
      expected = makeField<Field3D>(
        [&](Index& i) { return std::sin((i.*dir)() * box_length) * pow(box_length, 4); },
        mesh);
      break;
    // For now advection derivatives (upwind, flux) can have the same expected
    // result as the velocity field is constant
    case DERIV::Upwind:
    case DERIV::Flux:
      expected = makeField<Field3D>(
          [&](Index& i) { return 2.0 * std::cos((i.*dir)() * box_length) * box_length; },
          mesh);
      break;
    default:
      throw BoutException("Sorry, don't we test that type of derivative yet!");
    }

    // We need the parallel slices for the y-direction
    ParallelTransformIdentity identity{*mesh, std::vector<BoutReal>()};
    identity.calcParallelSlices(input);
    identity.calcParallelSlices(velocity);
  };

  virtual ~DerivativesTest() {
    delete mesh;
    mesh = nullptr;
  }

  Field3D input, velocity;
  Field3D expected;

  // Region not including the guard cells in current direction
  std::string region;
};

using DerivativesTestAdvection = DerivativesTest;

// Get all the available methods for this direction and turn it from a
// collection of strings to a collection of tuples of the direction,
// order, and strings so that the test fixture knows which direction
// we're using
auto getMethodsForDirection(DERIV derivative_order, DIRECTION direction)
    -> std::vector<std::tuple<DIRECTION, DERIV, std::string>> {

  auto available_methods = DerivativeStore<Field3D>::getInstance().getAvailableMethods(
      derivative_order, direction);

  // Method names together with the current direction and derivative type
  std::vector<std::tuple<DIRECTION, DERIV, std::string>> methods{};

  std::transform(std::begin(available_methods), std::end(available_methods),
                 std::back_inserter(methods), [&](std::string method) {
                   return std::make_tuple(direction, derivative_order, method);
                 });

  return methods;
};

// Returns the method out of the direction/order/method tuple for
// printing the test name
auto methodDirectionTupleToString(
    const ::testing::TestParamInfo<std::tuple<DIRECTION, DERIV, std::string>>& param)
    -> std::string {
  return std::get<2>(param.param);
}

// Instantiate the test for X, Y, Z for first derivatives
INSTANTIATE_TEST_SUITE_P(FirstX, DerivativesTest,
                        ::testing::ValuesIn(getMethodsForDirection(DERIV::Standard,
                                                                   DIRECTION::X)),
                        methodDirectionTupleToString);

INSTANTIATE_TEST_SUITE_P(FirstY, DerivativesTest,
                        ::testing::ValuesIn(getMethodsForDirection(DERIV::Standard,
                                                                   DIRECTION::Y)),
                        methodDirectionTupleToString);

INSTANTIATE_TEST_SUITE_P(FirstZ, DerivativesTest,
                        ::testing::ValuesIn(getMethodsForDirection(DERIV::Standard,
                                                                   DIRECTION::Z)),
                        methodDirectionTupleToString);

// Instantiate the test for X, Y, Z for second derivatives
INSTANTIATE_TEST_SUITE_P(SecondX, DerivativesTest,
                        ::testing::ValuesIn(getMethodsForDirection(DERIV::StandardSecond,
                                                                   DIRECTION::X)),
                        methodDirectionTupleToString);

INSTANTIATE_TEST_SUITE_P(SecondY, DerivativesTest,
                        ::testing::ValuesIn(getMethodsForDirection(DERIV::StandardSecond,
                                                                   DIRECTION::Y)),
                        methodDirectionTupleToString);

INSTANTIATE_TEST_SUITE_P(SecondZ, DerivativesTest,
                        ::testing::ValuesIn(getMethodsForDirection(DERIV::StandardSecond,
                                                                   DIRECTION::Z)),
                        methodDirectionTupleToString);

// Instantiate the test for X, Y, Z for fourth derivatives
INSTANTIATE_TEST_SUITE_P(FourthX, DerivativesTest,
                        ::testing::ValuesIn(getMethodsForDirection(DERIV::StandardFourth,
                                                                   DIRECTION::X)),
                        methodDirectionTupleToString);

INSTANTIATE_TEST_SUITE_P(FourthY, DerivativesTest,
                        ::testing::ValuesIn(getMethodsForDirection(DERIV::StandardFourth,
                                                                   DIRECTION::Y)),
                        methodDirectionTupleToString);

INSTANTIATE_TEST_SUITE_P(FourthZ, DerivativesTest,
                        ::testing::ValuesIn(getMethodsForDirection(DERIV::StandardFourth,
                                                                   DIRECTION::Z)),
                        methodDirectionTupleToString);

// Instantiate the test for X, Y, Z for upwind derivatives
INSTANTIATE_TEST_SUITE_P(UpwindX, DerivativesTestAdvection,
                        ::testing::ValuesIn(getMethodsForDirection(DERIV::Upwind,
                                                                   DIRECTION::X)),
                        methodDirectionTupleToString);

INSTANTIATE_TEST_SUITE_P(UpwindY, DerivativesTestAdvection,
                        ::testing::ValuesIn(getMethodsForDirection(DERIV::Upwind,
                                                                   DIRECTION::Y)),
                        methodDirectionTupleToString);

INSTANTIATE_TEST_SUITE_P(UpwindZ, DerivativesTestAdvection,
                        ::testing::ValuesIn(getMethodsForDirection(DERIV::Upwind,
                                                                   DIRECTION::Z)),
                        methodDirectionTupleToString);

// Instantiate the test for X, Y, Z for flux derivatives
INSTANTIATE_TEST_SUITE_P(FluxX, DerivativesTestAdvection,
                        ::testing::ValuesIn(getMethodsForDirection(DERIV::Flux,
                                                                   DIRECTION::X)),
                        methodDirectionTupleToString);

INSTANTIATE_TEST_SUITE_P(FluxY, DerivativesTestAdvection,
                        ::testing::ValuesIn(getMethodsForDirection(DERIV::Flux,
                                                                   DIRECTION::Y)),
                        methodDirectionTupleToString);

INSTANTIATE_TEST_SUITE_P(FluxZ, DerivativesTestAdvection,
                        ::testing::ValuesIn(getMethodsForDirection(DERIV::Flux,
                                                                   DIRECTION::Z)),
                        methodDirectionTupleToString);

// All standard derivatives have the same signature, so we can use a
// single test, just instantiate it for each direction/order combination
TEST_P(DerivativesTest, Sanity) {
  auto derivative = DerivativeStore<Field3D>::getInstance().getStandardDerivative(
      std::get<2>(GetParam()), std::get<0>(GetParam()), STAGGER::None,
      std::get<1>(GetParam()));

  Field3D result{mesh};
  result.allocate();
  derivative(input, result, region);

  EXPECT_TRUE(IsFieldEqual(result, expected, "RGN_NOBNDRY", derivatives_tolerance));
}

// All advection (upwind/flux) derivatives have the same signature, so we can use a
// single test, just instantiate it for each direction/order combination
TEST_P(DerivativesTestAdvection, Sanity) {
  auto derivative = DerivativeStore<Field3D>::getInstance().getFlowDerivative(
      std::get<2>(GetParam()), std::get<0>(GetParam()), STAGGER::None,
      std::get<1>(GetParam()));

  Field3D result{mesh};
  result.allocate();
  derivative(velocity, input, result, region);

  EXPECT_TRUE(
      IsFieldEqual(result, expected, "RGN_NOBNDRY", derivatives_tolerance));
}

/////////////////////////////////////////////////////////////////////
// The following tests are essentially identical to the above, expect
// that we test the derivatives through the actual (index) derivative
// interface. This makes things a little more awkward to do completely
// generically, so let's not bother

using FirstDerivativesInterfaceTest = DerivativesTest;

INSTANTIATE_TEST_SUITE_P(X, FirstDerivativesInterfaceTest,
                        ::testing::ValuesIn(getMethodsForDirection(DERIV::Standard,
                                                                   DIRECTION::X)),
                        methodDirectionTupleToString);

INSTANTIATE_TEST_SUITE_P(FirstY, FirstDerivativesInterfaceTest,
                        ::testing::ValuesIn(getMethodsForDirection(DERIV::Standard,
                                                                   DIRECTION::Y)),
                        methodDirectionTupleToString);

INSTANTIATE_TEST_SUITE_P(FirstZ, FirstDerivativesInterfaceTest,
                        ::testing::ValuesIn(getMethodsForDirection(DERIV::Standard,
                                                                   DIRECTION::Z)),
                        methodDirectionTupleToString);

TEST_P(FirstDerivativesInterfaceTest, Sanity) {
  Field3D result;
  switch (std::get<0>(GetParam())) {
    case DIRECTION::X:
      result = bout::derivatives::index::DDX(input);
      break;
    case DIRECTION::Y:
      result = bout::derivatives::index::DDY(input);
      break;
    case DIRECTION::Z:
      result = bout::derivatives::index::DDZ(input);
      break;
  default:
    break;
  }

  EXPECT_TRUE(IsFieldEqual(result, expected, "RGN_NOBNDRY", derivatives_tolerance));
}

using SecondDerivativesInterfaceTest = DerivativesTest;

INSTANTIATE_TEST_SUITE_P(X, SecondDerivativesInterfaceTest,
                        ::testing::ValuesIn(getMethodsForDirection(DERIV::StandardSecond,
                                                                   DIRECTION::X)),
                        methodDirectionTupleToString);

INSTANTIATE_TEST_SUITE_P(Y, SecondDerivativesInterfaceTest,
                        ::testing::ValuesIn(getMethodsForDirection(DERIV::StandardSecond,
                                                                   DIRECTION::Y)),
                        methodDirectionTupleToString);

INSTANTIATE_TEST_SUITE_P(Z, SecondDerivativesInterfaceTest,
                        ::testing::ValuesIn(getMethodsForDirection(DERIV::StandardSecond,
                                                                   DIRECTION::Z)),
                        methodDirectionTupleToString);

TEST_P(SecondDerivativesInterfaceTest, Sanity) {
  Field3D result;
  switch (std::get<0>(GetParam())) {
    case DIRECTION::X:
      result = bout::derivatives::index::D2DX2(input);
      break;
    case DIRECTION::Y:
      result = bout::derivatives::index::D2DY2(input);
      break;
    case DIRECTION::Z:
      result = bout::derivatives::index::D2DZ2(input);
      break;
  default:
    break;
  }

  EXPECT_TRUE(IsFieldEqual(result, expected, "RGN_NOBNDRY", derivatives_tolerance));
}

using FourthDerivativesInterfaceTest = DerivativesTest;

INSTANTIATE_TEST_SUITE_P(X, FourthDerivativesInterfaceTest,
                        ::testing::ValuesIn(getMethodsForDirection(DERIV::StandardFourth,
                                                                   DIRECTION::X)),
                        methodDirectionTupleToString);

INSTANTIATE_TEST_SUITE_P(Y, FourthDerivativesInterfaceTest,
                        ::testing::ValuesIn(getMethodsForDirection(DERIV::StandardFourth,
                                                                   DIRECTION::Y)),
                        methodDirectionTupleToString);

INSTANTIATE_TEST_SUITE_P(Z, FourthDerivativesInterfaceTest,
                        ::testing::ValuesIn(getMethodsForDirection(DERIV::StandardFourth,
                                                                   DIRECTION::Z)),
                        methodDirectionTupleToString);

TEST_P(FourthDerivativesInterfaceTest, Sanity) {
  Field3D result;
  switch (std::get<0>(GetParam())) {
    case DIRECTION::X:
      result = bout::derivatives::index::D4DX4(input);
      break;
    case DIRECTION::Y:
      result = bout::derivatives::index::D4DY4(input);
      break;
    case DIRECTION::Z:
      result = bout::derivatives::index::D4DZ4(input);
      break;
  default:
    break;
  }

  EXPECT_TRUE(IsFieldEqual(result, expected, "RGN_NOBNDRY", derivatives_tolerance));
}


using UpwindDerivativesInterfaceTest = DerivativesTest;

// Instantiate the test for X, Y, Z for upwind derivatives
INSTANTIATE_TEST_SUITE_P(X, UpwindDerivativesInterfaceTest,
                        ::testing::ValuesIn(getMethodsForDirection(DERIV::Upwind,
                                                                   DIRECTION::X)),
                        methodDirectionTupleToString);

INSTANTIATE_TEST_SUITE_P(Y, UpwindDerivativesInterfaceTest,
                        ::testing::ValuesIn(getMethodsForDirection(DERIV::Upwind,
                                                                   DIRECTION::Y)),
                        methodDirectionTupleToString);

INSTANTIATE_TEST_SUITE_P(Z, UpwindDerivativesInterfaceTest,
                        ::testing::ValuesIn(getMethodsForDirection(DERIV::Upwind,
                                                                   DIRECTION::Z)),
                        methodDirectionTupleToString);

TEST_P(UpwindDerivativesInterfaceTest, Sanity) {
  Field3D result;

  switch (std::get<0>(GetParam())) {
  case DIRECTION::X:
    result = bout::derivatives::index::VDDX(velocity, input);
    break;
  case DIRECTION::Y:
    result = bout::derivatives::index::VDDY(velocity, input);
    break;
  case DIRECTION::Z:
    result = bout::derivatives::index::VDDZ(velocity, input);
    break;
  default:
    break;
  }

  EXPECT_TRUE(IsFieldEqual(result, expected, "RGN_NOBNDRY", derivatives_tolerance));
}

using FluxDerivativesInterfaceTest = DerivativesTest;

// Instantiate the test for X, Y, Z for flux derivatives
INSTANTIATE_TEST_SUITE_P(X, FluxDerivativesInterfaceTest,
                        ::testing::ValuesIn(getMethodsForDirection(DERIV::Flux,
                                                                   DIRECTION::X)),
                        methodDirectionTupleToString);

INSTANTIATE_TEST_SUITE_P(Y, FluxDerivativesInterfaceTest,
                        ::testing::ValuesIn(getMethodsForDirection(DERIV::Flux,
                                                                   DIRECTION::Y)),
                        methodDirectionTupleToString);

INSTANTIATE_TEST_SUITE_P(Z, FluxDerivativesInterfaceTest,
                        ::testing::ValuesIn(getMethodsForDirection(DERIV::Flux,
                                                                   DIRECTION::Z)),
                        methodDirectionTupleToString);

TEST_P(FluxDerivativesInterfaceTest, Sanity) {
  Field3D result;

  switch (std::get<0>(GetParam())) {
  case DIRECTION::X:
    result = bout::derivatives::index::FDDX(velocity, input);
    break;
  case DIRECTION::Y:
    result = bout::derivatives::index::FDDY(velocity, input);
    break;
  case DIRECTION::Z:
    result = bout::derivatives::index::FDDZ(velocity, input);
    break;
  default:
    break;
  }

  EXPECT_TRUE(IsFieldEqual(result, expected, "RGN_NOBNDRY", derivatives_tolerance));
}
