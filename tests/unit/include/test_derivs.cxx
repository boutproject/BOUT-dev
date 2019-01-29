#include "gtest/gtest.h"

#include "bout/constants.hxx"
#include "bout/deriv_store.hxx"
#include "bout_types.hxx"
#include "field3d.hxx"
#include "test_extras.hxx"

#include <string>

// We don't just inherit from FakeMeshTest as we want to increase nx
// to get reasonable agreement
class DerivativesTest : public ::testing::TestWithParam<std::string> {
public:
  DerivativesTest() : input{mesh}, expected{mesh} {
    if (mesh != nullptr) {
      delete mesh;
      mesh = nullptr;
    }
    mesh = new FakeMesh(nx, ny, nz);
    mesh->xstart = 2;
    mesh->xend = nx - 3;

    mesh->ystart = 0;
    mesh->yend = ny - 1;

    output_info.disable();
    mesh->createDefaultRegions();
    output_info.enable();

    Field3D input_{mesh};
    input_.allocate();

    for (auto i : input_) {
      input_[i] = std::sin(i.x() * TWOPI / (nx - 1));
    }

    input = input_;

    Field3D expected_{mesh};
    expected_.allocate();

    for (auto i : expected_) {
      expected_[i] = std::cos(i.x() * TWOPI / (nx - 1)) * TWOPI / (nx - 1);
    }

    expected = expected_;

    DerivativeStore<Field3D>::getInstance().initialise(Options::getRoot());
  };

  Field3D input;
  Field3D expected;

  static constexpr int nx = 64;
  static constexpr int ny = 3;
  static constexpr int nz = 2;
};

INSTANTIATE_TEST_CASE_P(SanityCheck, DerivativesTest,
                        ::testing::ValuesIn(DerivativeStore<Field3D>::getInstance()
                                                .getAvailableMethods(DERIV::Standard,
                                                                     DIRECTION::X)));

TEST_P(DerivativesTest, Basic) {
  auto derivative = DerivativeStore<Field3D>::getInstance().getStandardDerivative(
      GetParam(), DIRECTION::X);

  Field3D result{mesh};
  result.allocate();
  derivative(input, result, RGN_NOX);

  EXPECT_TRUE(IsField3DEqualField3D(result, expected, "RGN_NOBNDRY", 1.e-2));
}
