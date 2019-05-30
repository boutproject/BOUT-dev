#include "gtest/gtest.h"

#include "../src/mesh/impls/bout/boutmesh.hxx"
#include "bout/mesh.hxx"
#include "output.hxx"
#include "unused.hxx"

#include "test_extras.hxx"

TEST(BoutMeshTest, NullOptionsCheck) {
  // Temporarily turn off outputs to make test quiet
  output_info.disable();
  output_warn.disable();
  EXPECT_NO_THROW(BoutMesh mesh(new FakeGridDataSource, nullptr));
  output_info.enable();
  output_warn.enable();
}
