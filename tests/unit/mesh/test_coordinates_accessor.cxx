#include "gtest/gtest.h"

#include "bout/build_config.hxx"
#include "bout/coordinates_accessor.hxx"

#include "test_extras.hxx"

/// Global mesh
namespace bout {
namespace globals {
extern Mesh* mesh;
}
} // namespace bout

using bout::globals::mesh;

using CoordinatesAccessorTest = FakeMeshFixture;
using FieldMetric = Coordinates::FieldMetric;

TEST_F(CoordinatesAccessorTest, CreateAccessor) {
  CoordinatesAccessor acc(mesh->getCoordinates());

  // Basic sanity checks
  EXPECT_TRUE(acc.mesh_nz == mesh->LocalNz);
  EXPECT_TRUE(acc.data != nullptr);
  EXPECT_FLOAT_EQ(mesh->getCoordinates()->dx(0,0,0), acc.dx(0));
}

TEST_F(CoordinatesAccessorTest, CreateTwoAccessors) {
  CoordinatesAccessor acc(mesh->getCoordinates());

  // Create another accessor from the same Coordinates
  CoordinatesAccessor acc2(mesh->getCoordinates());

  // Should share the same data pointer
  EXPECT_EQ(acc.data, acc2.data);
}
