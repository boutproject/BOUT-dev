#include "gtest/gtest.h"

#include "bout/output_bout_types.hxx"
#include "bout/petsc_operators.hxx"

#include "fake_mesh_fixture.hxx"
#include "test_extras.hxx"

namespace {
Ind3D indexAt(const Field3D& f, int x, int y, int z) {
  int ny = f.getNy();
  int nz = f.getNz();
  return Ind3D{(((x * ny) + y) * nz) + z, ny, nz};
}
} // namespace

using bout::globals::mesh;

// Reuse the "standard" fixture for FakeMesh
using PetscOperatorsTest = FakeMeshFixture;

TEST_F(PetscOperatorsTest, create_region_empty) {
  Field3D cell_number{-1.0}; // No cells >= 0

  auto rgn = PetscOperators::create_region(cell_number);

  ASSERT_EQ(0, rgn.size());
}

TEST_F(PetscOperatorsTest, create_region) {
  Field3D cell_number{-1.0}; // No cells >= 0

  // Note: 1 boundary cell in X and Y
  cell_number(0, 0, 0) = 0; // Corner
  cell_number(1, 1, 0) = 1; // In domain
  cell_number(0, 1, 1) = 2; // Xin boundary

  auto rgn = PetscOperators::create_region(cell_number);
  ASSERT_EQ(1, rgn.size());

  Ind3D expected_ind = indexAt(cell_number, 1, 1, 0);
  output.write("{} : {}\n", expected_ind, *rgn.begin());
  ASSERT_TRUE(expected_ind == *rgn.begin());
}
