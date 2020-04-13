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

  // Create a generator, then collect it
  Field3D result = field | collect();

  ASSERT_TRUE(IsFieldEqual(result, field));
}

TEST_F(PipelinesTest, IdentityChain) {
  Field3D field {1.0};

  // Create a generator, collect into a Field, and repeat
  Field3D result = field | collect() | collect();

  ASSERT_TRUE(IsFieldEqual(result, field));
}

TEST_F(PipelinesTest, AddConstant) {
  Field3D field {1.0};

  // Create a generator, add and collect
  Field3D result = field | add(2.0) | collect();

  ASSERT_TRUE(IsFieldEqual(result, 3.0));
}
