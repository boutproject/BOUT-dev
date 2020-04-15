#include "gtest/gtest.h"

#include "test_extras.hxx"

#include "bout/single_index_ops.hxx"
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
using SingleIndexOpsTest = FakeMeshFixture;

TEST_F(SingleIndexOpsTest, DDX) {
  // Fill a field with random numbers
  std::uniform_real_distribution<double> unif(-1.0, 1.0);
  std::default_random_engine re;
  
  Field3D input; input.allocate();
  BOUT_FOR(i, input.getRegion("RGN_ALL")) {
    input[i] = unif(re);
  }

  // Check that the field is not zero
  ASSERT_FALSE(IsFieldEqual(input, 0.0, "RGN_NOBNDRY"));
  
  // Differentiate whole field
  Field3D difops = DDX(input);

  // Differentiate using index operations
  Field3D indexops; indexops.allocate();
  BOUT_FOR(i, input.getRegion("RGN_NOBNDRY")) {
    indexops[i] = DDX(input, i);
  }

  // Check the answer is the same
  ASSERT_TRUE(IsFieldEqual(difops, indexops, "RGN_NOBNDRY"));
}

TEST_F(SingleIndexOpsTest, DDY) {
  // Fill a field with random numbers
  std::uniform_real_distribution<double> unif(-1.0, 1.0);
  std::default_random_engine re;
  
  Field3D input; input.allocate();
  BOUT_FOR(i, input.getRegion("RGN_ALL")) {
    input[i] = unif(re);
  }

  // Check that the field is not zero
  ASSERT_FALSE(IsFieldEqual(input, 0.0, "RGN_NOBNDRY"));

  // need the yup/ydown fields
  mesh->communicate(input);
  
  // Differentiate whole field
  Field3D difops = DDY(input);

  // Differentiate using index operations
  Field3D indexops; indexops.allocate();
  BOUT_FOR(i, input.getRegion("RGN_NOBNDRY")) {
    indexops[i] = DDY(input, i);
  }

  // Check the answer is the same
  ASSERT_TRUE(IsFieldEqual(difops, indexops, "RGN_NOBNDRY"));
}

TEST_F(SingleIndexOpsTest, DDZ) {
  // Fill a field with random numbers
  std::uniform_real_distribution<double> unif(-1.0, 1.0);
  std::default_random_engine re;
  
  Field3D input; input.allocate();
  BOUT_FOR(i, input.getRegion("RGN_ALL")) {
    input[i] = unif(re);
  }

  // Check that the field is not zero
  ASSERT_FALSE(IsFieldEqual(input, 0.0, "RGN_NOBNDRY"));
  
  // Differentiate whole field
  Field3D difops = DDZ(input);

  // Differentiate using index operations
  Field3D indexops; indexops.allocate();
  BOUT_FOR(i, input.getRegion("RGN_NOBNDRY")) {
    indexops[i] = DDZ(input, i);
  }

  // Check the answer is the same
  ASSERT_TRUE(IsFieldEqual(difops, indexops, "RGN_NOBNDRY"));
}

TEST_F(SingleIndexOpsTest, bracket) {
  // Fill a field with random numbers
  std::uniform_real_distribution<double> unif(-1.0, 1.0);
  std::default_random_engine re;
  
  Field3D input; input.allocate();
  BOUT_FOR(i, input.getRegion("RGN_ALL")) {
    input[i] = unif(re);
  }
  // Check that the field is not zero
  ASSERT_FALSE(IsFieldEqual(input, 0.0, "RGN_NOBNDRY"));
  
  Field3D input2; input2.allocate();
  BOUT_FOR(i, input.getRegion("RGN_ALL")) {
    input2[i] = unif(re);
  }
  ASSERT_FALSE(IsFieldEqual(input2, 0.0, "RGN_NOBNDRY"));

  // Differentiate whole field
  Field3D difops = bracket(input, input2, BRACKET_ARAKAWA);

  // Differentiate using index operations
  Field3D indexops; indexops.allocate();
  BOUT_FOR(i, input.getRegion("RGN_NOBNDRY")) {
    indexops[i] = bracket(input, input2, i);
  }

  // Check the answer is the same
  ASSERT_TRUE(IsFieldEqual(difops, indexops, "RGN_NOBNDRY"));
}
