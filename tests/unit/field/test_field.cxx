#include "gtest/gtest.h"

#include "bout/constants.hxx"
#include "bout/mesh.hxx"
#include "boutexception.hxx"
#include "field.hxx"
#include "output.hxx"
#include "test_extras.hxx"

/// Global mesh
namespace bout{
namespace globals{
extern Mesh *mesh;
} // namespace globals
} // namespace bout

// The unit tests use the global mesh
using namespace bout::globals;

// Reuse the "standard" fixture for FakeMesh
using FieldTest = FakeMeshFixture;


TEST_F(FieldTest, GetGlobalMesh) {
  Field field;

  auto localmesh = field.getMesh();

  EXPECT_EQ(localmesh, mesh);
}

TEST_F(FieldTest, GetLocalMesh) {
  FakeMesh myMesh{nx + 1, ny + 2, nz + 3};
  Field field(&myMesh, CELL_CENTRE, DIRECTION::X, DIRECTION::Y, DIRECTION::Z);

  auto localmesh = field.getMesh();

  EXPECT_EQ(localmesh, &myMesh);
}

TEST_F(FieldTest, GetGridSizes) {
  Field field;

  EXPECT_EQ(field.getNx(), nx);
  EXPECT_EQ(field.getNy(), ny);
  EXPECT_EQ(field.getNz(), nz);
}
