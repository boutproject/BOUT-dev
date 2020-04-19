#include "gtest/gtest.h"


#include "test_extras.hxx"

#include "bout/pipelines.hxx"

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
