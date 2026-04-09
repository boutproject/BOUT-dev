#include "bout/build_defines.hxx"

#if BOUT_HAS_PETSC

#include "gtest/gtest.h"

#include "bout/array.hxx"
#include "bout/bout_types.hxx"
#include "bout/output.hxx"
#include "bout/output_bout_types.hxx"
#include "bout/petsc_operators.hxx"
#include "bout/region.hxx"

#include "fake_mesh_fixture.hxx"
#include "test_extras.hxx"

#include <memory>

// Reuse the "standard" fixture for FakeMesh
using PetscMappingTest = FakeMeshFixture;

TEST_F(PetscMappingTest, create_region_empty) {
  const Field3D cell_number{-1.0}; // No cells >= 0

  auto rgn = PetscCellMapping::create_region(cell_number);

  ASSERT_EQ(0, rgn.size());
}

TEST_F(PetscMappingTest, create_region) {
  Field3D cell_number{-1.0}; // No cells >= 0

  // Note: 1 boundary cell in X and Y
  cell_number(0, 0, 0) = 0; // Corner
  cell_number(1, 1, 0) = 1; // In domain
  cell_number(0, 1, 1) = 2; // Xin boundary

  auto rgn = PetscCellMapping::create_region(cell_number);
  ASSERT_EQ(1, rgn.size());

  const Ind3D expected_ind = cell_number.indexAt(1, 1, 0);
  ASSERT_TRUE(expected_ind == *rgn.begin());
}

TEST_F(PetscMappingTest, create_region_xin) {
  Field3D cell_number{-1.0}; // No cells >= 0

  // Note: 1 boundary cell in X and Y
  cell_number(0, 0, 0) = 0; // Corner
  cell_number(1, 1, 0) = 1; // In domain
  cell_number(0, 1, 1) = 2; // Xin boundary

  auto rgn = PetscCellMapping::create_region_xin(cell_number);
  ASSERT_EQ(1, rgn.size());

  const Ind3D expected_ind = cell_number.indexAt(0, 1, 1);
  output.write("Expecting {} == {}\n", expected_ind, *rgn.begin());
  ASSERT_TRUE(expected_ind == *rgn.begin());
}

TEST_F(PetscMappingTest, mapping) {
  Field3D cell_number{-1.0}; // No cells >= 0

  // Note: 1 boundary cell in X and Y
  cell_number(0, 0, 0) = 2; // Corner
  cell_number(1, 1, 0) = 0; // In domain
  cell_number(0, 1, 1) = 1; // Xin boundary

  const Field3D forward_cell_number{-1.0};
  const Field3D backward_cell_number{-1.0};

  const PetscCellMapping mapping(cell_number, forward_cell_number, backward_cell_number,
                                 2);

  // Two cells: one evolving and one in xin
  ASSERT_EQ(2, mapping.localSize());
  ASSERT_EQ(2, mapping.globalSize());
}

using PetscOperatorTest = FakeMeshFixture;

TEST_F(PetscOperatorTest, identity) {
  Field3D cell_number{-1.0};
  const Field3D forward_cell_number{-1.0};
  const Field3D backward_cell_number{-1.0};

  // Three cells active
  cell_number(1, 1, 0) = 0;
  cell_number(1, 2, 0) = 1;
  cell_number(1, 2, 1) = 2;

  auto mapping = std::make_shared<const PetscCellMapping>(
      cell_number, forward_cell_number, backward_cell_number, 3);
  ASSERT_EQ(3, mapping->localSize());

  const auto rows = Array<int>::fromValues({0, 1, 2, 3});
  const auto cols = Array<int>::fromValues({0, 1, 2});
  const auto weights = Array<BoutReal>::fromValues({1.0, 1.0, 1.0});

  const PetscCellOperator identity(mapping, mapping, rows, cols, weights);

  Field3D value{0.0};
  value(1, 1, 0) = 10;
  value(1, 2, 0) = 21;
  value(1, 2, 1) = 32;
  value.splitParallelSlices();
  value.yup() = -1.0;
  value.ydown() = -1.0;

  const Field3D result = identity(value);

  EXPECT_EQ(10, result(1, 1, 0));
  EXPECT_EQ(21, result(1, 2, 0));
  EXPECT_EQ(32, result(1, 2, 1));
}

TEST_F(PetscOperatorTest, identity_yupdown) {
  Field3D cell_number{-1.0};
  Field3D forward_cell_number{-1.0};
  Field3D backward_cell_number{-1.0};

  // Three cells active
  cell_number(1, 1, 0) = 0;
  cell_number(1, 2, 0) = 1;
  cell_number(1, 2, 1) = 2;

  forward_cell_number(1, 2, 0) = 3;
  backward_cell_number(1, 1, 2) = 4;

  auto mapping = std::make_shared<const PetscCellMapping>(
      cell_number, forward_cell_number, backward_cell_number, 5);
  ASSERT_EQ(5, mapping->localSize());

  const auto rows = Array<int>::fromValues({0, 1, 2, 3});
  const auto cols = Array<int>::fromValues({0, 1, 2});
  const auto weights = Array<BoutReal>::fromValues({1.0, 1.0, 1.0});

  const PetscCellOperator identity(mapping, mapping, rows, cols, weights);

  Field3D value{0.0};
  value(1, 1, 0) = 10;
  value(1, 2, 0) = 21;
  value(1, 2, 1) = 32;
  value.splitParallelSlices();
  value.yup() = -1.0;
  value.ydown() = -1.0;

  const Field3D result = identity(value);

  EXPECT_EQ(10, result(1, 1, 0));
  EXPECT_EQ(21, result(1, 2, 0));
  EXPECT_EQ(32, result(1, 2, 1));
}

#endif // BOUT_HAS_PETSC
