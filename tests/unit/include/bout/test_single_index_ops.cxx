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

/// Creates a Field filled with random, uniformly distributed numbers
///
/// Usage
///   std::default_random_engine re;
///   auto field = random_field<Field3D>(re);
///
template <class FieldType>
FieldType random_field(std::default_random_engine& re) {
  std::uniform_real_distribution<double> unif(-1.0, 1.0);
  // Create a field of requested type
  FieldType result;
  result.allocate();
  // Fill with random numbers
  BOUT_FOR(i, result.getRegion("RGN_ALL")) { result[i] = unif(re); }
  return result;
}

TEST_F(SingleIndexOpsTest, DDX) {
  // Fill a field with random numbers
  std::default_random_engine re;

  auto input = random_field<Field3D>(re);

  // Check that the field is not zero
  ASSERT_FALSE(IsFieldEqual(input, 0.0, "RGN_NOBNDRY"));
  
  // Differentiate whole field
  Field3D difops = DDX(input);

  // Differentiate using index operations
  Field3D indexops; indexops.allocate();
  auto input_acc = FieldAccessor<>(input);
  BOUT_FOR(i, input.getRegion("RGN_NOBNDRY")) {
    indexops[i] = DDX(input_acc, i);
  }

  // Check the answer is the same
  ASSERT_TRUE(IsFieldEqual(difops, indexops, "RGN_NOBNDRY"));
}

TEST_F(SingleIndexOpsTest, DDY) {
  // Fill a field with random numbers
  std::default_random_engine re;

  auto input = random_field<Field3D>(re);

  // Check that the field is not zero
  ASSERT_FALSE(IsFieldEqual(input, 0.0, "RGN_NOBNDRY"));

  // need the yup/ydown fields
  mesh->communicate(input);
  
  // Differentiate whole field
  Field3D difops = DDY(input);

  // Differentiate using index operations
  Field3D indexops; indexops.allocate();
  auto input_acc = FieldAccessor<>(input);
  BOUT_FOR(i, input.getRegion("RGN_NOBNDRY")) {
    indexops[i] = DDY(input_acc, i);
  }

  // Check the answer is the same
  ASSERT_TRUE(IsFieldEqual(difops, indexops, "RGN_NOBNDRY"));
}

TEST_F(SingleIndexOpsTest, DDZ) {
  // Fill a field with random numbers
  std::default_random_engine re;

  auto input = random_field<Field3D>(re);

  // Check that the field is not zero
  ASSERT_FALSE(IsFieldEqual(input, 0.0, "RGN_NOBNDRY"));
  
  // Differentiate whole field
  Field3D difops = DDZ(input);

  // Differentiate using index operations
  Field3D indexops; indexops.allocate();
  auto input_acc = FieldAccessor<>(input);
  BOUT_FOR(i, input.getRegion("RGN_NOBNDRY")) {
    indexops[i] = DDZ(input_acc, i);
  }

  // Check the answer is the same
  ASSERT_TRUE(IsFieldEqual(difops, indexops, "RGN_NOBNDRY"));
}

TEST_F(SingleIndexOpsTest, bracket) {
  // Fill a field with random numbers
  std::default_random_engine re;

  auto input = random_field<Field3D>(re);
  // Check that the field is not zero
  ASSERT_FALSE(IsFieldEqual(input, 0.0, "RGN_NOBNDRY"));

  auto input2 = random_field<Field3D>(re);
  ASSERT_FALSE(IsFieldEqual(input2, 0.0, "RGN_NOBNDRY"));

  // Differentiate whole field
  Field3D difops = bracket(input, input2, BRACKET_ARAKAWA);

  // Differentiate using index operations
  Field3D indexops; indexops.allocate();
  auto input_acc = FieldAccessor<>(input);
  auto input2_acc = FieldAccessor<>(input2);
  BOUT_FOR(i, input.getRegion("RGN_NOBNDRY")) {
    indexops[i] = bracket(input_acc, input2_acc, i);
  }

  // Check the answer is the same
  ASSERT_TRUE(IsFieldEqual(difops, indexops, "RGN_NOBNDRY"));
}

TEST_F(SingleIndexOpsTest, Delp2_3D) {
  // Fill a field with random numbers
  std::default_random_engine re;

  auto input = random_field<Field3D>(re);

  // Check that the field is not zero
  ASSERT_FALSE(IsFieldEqual(input, 0.0, "RGN_NOBNDRY"));

  // Differentiate whole field
  // Note: Uses non-FFT version for comparison to the single index operator
  Field3D difops = Delp2(input, CELL_DEFAULT, false);

  // Differentiate using index operations
  Field3D indexops; indexops.allocate();
  auto input_acc = FieldAccessor<>(input);
  BOUT_FOR(i, input.getRegion("RGN_NOBNDRY")) {
    indexops[i] = Delp2(input_acc, i);
  }

  // Check the answer is the same
  ASSERT_TRUE(IsFieldEqual(difops, indexops, "RGN_NOBNDRY"));
}

TEST_F(SingleIndexOpsTest, b0xGrad_dot_Grad_3D_2D) {
  std::default_random_engine re;

  auto input1 = random_field<Field3D>(re);
  auto input2 = random_field<Field2D>(re);

  ASSERT_FALSE(IsFieldEqual(input1, 0.0, "RGN_NOBNDRY"));
  ASSERT_FALSE(IsFieldEqual(input2, 0.0, "RGN_NOBNDRY"));

  // Need parallel derivatives of input1 and input2
  // input2 is Field2D, so always has parallel slices
  input1.calcParallelSlices();

  // Set the numerical methods
  // Note: Important setting is Field3D upwind (U1 -> C2)
  Options diff_options {{"ddx", {{"first", "C2"},
                                 {"upwind", "C2"}}},
                        {"ddy", {{"first", "C2"},
                                 {"upwind", "C2"}}},
                        {"ddz", {{"first", "C2"},
                                 {"upwind", "C2"}}}};
  DerivativeStore<Field2D>::getInstance().initialise(&diff_options);
  DerivativeStore<Field3D>::getInstance().initialise(&diff_options);

  // Differentiate whole field
  Field3D difops = b0xGrad_dot_Grad(input1, input2);

  // Differentiate using index operations
  Field3D indexops; indexops.allocate();
  auto input1_acc = FieldAccessor<>(input1);
  auto input2_acc = Field2DAccessor<>(input2);
  BOUT_FOR(i, input1.getRegion("RGN_NOBNDRY")) {
    indexops[i] = b0xGrad_dot_Grad(input1_acc, input2_acc, i);
  }

  // Check the answer is the same
  ASSERT_TRUE(IsFieldEqual(difops, indexops, "RGN_NOBNDRY"));
}

TEST_F(SingleIndexOpsTest, D2DY2) {
  std::default_random_engine re;

  auto input = random_field<Field3D>(re);

  ASSERT_FALSE(IsFieldEqual(input, 0.0, "RGN_NOBNDRY"));

  // Need parallel derivatives of input
  input.calcParallelSlices();

  // Differentiate whole field
  Field3D difops = D2DY2(input);

  // Differentiate using index operations
  Field3D indexops; indexops.allocate();
  auto input_acc = FieldAccessor<>(input);
  BOUT_FOR(i, input.getRegion("RGN_NOBNDRY")) {
    indexops[i] = D2DY2(input_acc, i);
  }

  // Check the answer is the same
  ASSERT_TRUE(IsFieldEqual(difops, indexops, "RGN_NOBNDRY"));
}

TEST_F(SingleIndexOpsTest, Grad_par) {
  std::default_random_engine re;

  auto input = random_field<Field3D>(re);

  ASSERT_FALSE(IsFieldEqual(input, 0.0, "RGN_NOBNDRY"));

  // Need parallel derivatives of input
  input.calcParallelSlices();

  // Differentiate whole field
  Field3D difops = Grad_par(input);

  // Differentiate using index operations
  Field3D indexops; indexops.allocate();
  auto input_acc = FieldAccessor<>(input);
  BOUT_FOR(i, input.getRegion("RGN_NOBNDRY")) {
    indexops[i] = Grad_par(input_acc, i);
  }

  // Check the answer is the same
  ASSERT_TRUE(IsFieldEqual(difops, indexops, "RGN_NOBNDRY"));
}

TEST_F(SingleIndexOpsTest, Div_par) {
  std::default_random_engine re;

  auto input = random_field<Field3D>(re);

  ASSERT_FALSE(IsFieldEqual(input, 0.0, "RGN_NOBNDRY"));

  // Need parallel derivatives of input
  input.calcParallelSlices();

  // Differentiate whole field
  Field3D difops = Div_par(input);

  // Differentiate using index operations
  Field3D indexops; indexops.allocate();
  auto input_acc = FieldAccessor<>(input);
  BOUT_FOR(i, input.getRegion("RGN_NOBNDRY")) {
    indexops[i] = Div_par(input_acc, i);
  }

  // Check the answer is the same
  ASSERT_TRUE(IsFieldEqual(difops, indexops, "RGN_NOBNDRY"));
}
