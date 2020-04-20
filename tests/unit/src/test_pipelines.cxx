#include "gtest/gtest.h"


#include "test_extras.hxx"

#include "bout/pipelines.hxx"
#include "difops.hxx"
#include "derivs.hxx"

#include <random>

/// Global mesh
namespace bout{
namespace globals{
extern Mesh *mesh;
} // namespace globals
} // namespace bout

// The unit tests use the global mesh
using namespace bout::globals;

// Reuse the "standard" fixture for FakeMesh
using PipelinesTest = FakeMeshFixture;

using namespace bout::experimental;

TEST_F(PipelinesTest, Identity) {
  Field3D field {1.0};

  // Create a generator, then convert back to a field
  Field3D result = field | toField();

  ASSERT_TRUE(IsFieldEqual(result, field));
}

TEST_F(PipelinesTest, IdentityChain) {
  Field3D field {1.0};

  // Create a generator, convert into a Field, and repeat
  Field3D result = field | toField() | toField();

  ASSERT_TRUE(IsFieldEqual(result, field));
}

TEST_F(PipelinesTest, AddConstant) {
  Field3D field {1.0};

  // Create a generator, add and collect
  Field3D result = field | add(2.0) | toField();

  ASSERT_TRUE(IsFieldEqual(result, 3.0));
}

TEST_F(PipelinesTest, AddConstantToConstant) {
  // Create a generator from a constant, add and convert to field
  Field3D result = 1.0 | add(2.0) | toField();

  ASSERT_TRUE(IsFieldEqual(result, 3.0));
}
TEST_F(PipelinesTest, MakeFieldRegion) {
  // Only iterate over the core region (no boundaries)
  Field3D result = 1.0 | toField(mesh->getRegion3D("RGN_NOBNDRY"));

  ASSERT_TRUE(IsFieldEqual(result, 1.0, "RGN_NOBNDRY"));
  ASSERT_FALSE(IsFieldEqual(result, 1.0));
}

TEST_F(PipelinesTest, SubtractConstant) {
  // Create a generator from a constant, add and convert to field
  Field3D result = 3.0 | subtract(2.0) | toField();

  ASSERT_TRUE(IsFieldEqual(result, 1.0));
}

TEST_F(PipelinesTest, MultiplyConstant) {
  // Create a generator from a constant, add and convert to field
  Field3D result = 3.0 | multiply(2.0) | toField();

  ASSERT_TRUE(IsFieldEqual(result, 6.0));
}

TEST_F(PipelinesTest, DivideConstant) {
  // Create a generator from a constant, add and convert to field
  Field3D result = 3.0 | divide(2.0) | toField();

  ASSERT_TRUE(IsFieldEqual(result, 1.5));
}

TEST_F(PipelinesTest, Arithmetic) {
  // Create a generator from a constant, add and convert to field
  Field3D result =
    3.0
    | divide(2.0
             | multiply(1.5)
             | subtract(4.0)
             | multiply(2.0))
    | add(5.0 | divide(2.0))
    | toField();
  
  ASSERT_TRUE(IsFieldEqual(result, 3./((2 * 1.5 - 4.)*2) + 5./2));
}

TEST_F(PipelinesTest, bracket) {
  // Fill a field with random numbers
  std::uniform_real_distribution<double> unif(-1.0, 1.0);
  std::default_random_engine re;
  
  Field3D input; input.allocate();
  BOUT_FOR(i, input.getRegion("RGN_ALL")) {
    input[i] = unif(re);
  }
  // Check that the field is not zero
  EXPECT_FALSE(IsFieldEqual(input, 0.0, "RGN_NOBNDRY"));

  Field3D input2; input2.allocate();
  BOUT_FOR(i, input.getRegion("RGN_ALL")) {
    input2[i] = unif(re);
  }
  ASSERT_FALSE(IsFieldEqual(input2, 0.0, "RGN_NOBNDRY"));
  
  // Differentiate whole field
  Field3D difops = bracket(input, input2, BRACKET_ARAKAWA);

  input.fast_coords = input.getCoordinates(); // Set coordinates for single-index ops
  input2.fast_coords = input2.getCoordinates();
  Field3D pipes = bp::bracket(input, input2) | toField(input.getRegion("RGN_NOBNDRY"));

  // Check the answer is the same
  ASSERT_TRUE(IsFieldEqual(difops, pipes, "RGN_NOBNDRY"));
}

TEST_F(PipelinesTest, DDZ) {
  // Fill a field with random numbers
  std::uniform_real_distribution<double> unif(-1.0, 1.0);
  std::default_random_engine re;
  
  Field3D input; input.allocate();
  BOUT_FOR(i, input.getRegion("RGN_ALL")) {
    input[i] = unif(re);
  }
  // Check that the field is not zero
  EXPECT_FALSE(IsFieldEqual(input, 0.0, "RGN_NOBNDRY"));

  // Differentiate whole field
  Field3D difops = DDZ(input);

  input.fast_coords = input.getCoordinates(); // Set coordinates for single-index ops
  Field3D pipes = bp::DDZ(input) | toField(input.getRegion("RGN_NOBNDRY"));

  // Check the answer is the same
  ASSERT_TRUE(IsFieldEqual(difops, pipes, "RGN_NOBNDRY"));
}
