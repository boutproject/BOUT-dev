#include "gtest/gtest.h"

#include "bout/field_accessor.hxx"

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
using FieldAccessorTest = FakeMeshFixture;

TEST_F(FieldAccessorTest, CreationDefault) {
  Field3D var = 0.0;
  auto var_acc = FieldAccessor<>(var);
  EXPECT_EQ(var_acc.data, &var(0,0,0));
}

TEST_F(FieldAccessorTest, CreationUnallocated) {
  Field3D var;
  // var not allocated
  EXPECT_THROW(auto var_acc = FieldAccessor<>(var), BoutException);
}

TEST_F(FieldAccessorTest, CreateYlow) {
  Field3D field(mesh_staggered);
  field = 0.0;
  
  field.getMesh()->StaggerGrids = true;
  
  field.setLocation(CELL_YLOW);

  auto field_acc = FieldAccessor<CELL_YLOW>(field);
  EXPECT_EQ(field_acc.data, &field(0,0,0));
}

TEST_F(FieldAccessorTest, CreateYlowWrongLocation) {
  Field3D field(mesh_staggered);
  field = 0.0;

  field.getMesh()->StaggerGrids = true;
  field.setLocation(CELL_CENTRE);

  // Trying to create a YLOW field access from a CENTRE field
  EXPECT_THROW(auto field_acc = FieldAccessor<CELL_YLOW>(field), BoutException);
}

