#include "gtest/gtest.h"

#include "bout/output.hxx"
#include "bout/output_bout_types.hxx"
#include "bout/petsc_operators.hxx"
#include "bout/region.hxx"

#include "fake_mesh_fixture.hxx"
#include "test_extras.hxx"

// Reuse the "standard" fixture for FakeMesh
using PetscMappingTest = FakeMeshFixture;

TEST_F(PetscMappingTest, create_region_empty) {
  const Field3D cell_number{-1.0}; // No cells >= 0

  auto rgn = PetscMapping::create_region(cell_number);

  ASSERT_EQ(0, rgn.size());
}

TEST_F(PetscMappingTest, create_region) {
  Field3D cell_number{-1.0}; // No cells >= 0

  // Note: 1 boundary cell in X and Y
  cell_number(0, 0, 0) = 0; // Corner
  cell_number(1, 1, 0) = 1; // In domain
  cell_number(0, 1, 1) = 2; // Xin boundary

  auto rgn = PetscMapping::create_region(cell_number);
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

  auto rgn = PetscMapping::create_region_xin(cell_number);
  ASSERT_EQ(1, rgn.size());

  const Ind3D expected_ind = cell_number.indexAt(0, 1, 1);
  output.write("Expecting {} == {}\n", expected_ind, *rgn.begin());
  ASSERT_TRUE(expected_ind == *rgn.begin());
}

TEST_F(PetscMappingTest, mapping) {
  Field3D cell_number{-1.0}; // No cells >= 0

  // Note: 1 boundary cell in X and Y
  cell_number(0, 0, 0) = 0; // Corner
  cell_number(1, 1, 0) = 1; // In domain
  cell_number(0, 1, 1) = 2; // Xin boundary

  const Field3D forward_cell_number{-1.0};
  const Field3D backward_cell_number{-1.0};

  const PetscMapping mapping(cell_number, forward_cell_number, backward_cell_number);

  // Two cells: one evolving and one in xin
  ASSERT_EQ(2, mapping.size());
}
