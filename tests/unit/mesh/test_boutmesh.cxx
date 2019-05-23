#include "gtest/gtest.h"

#include "../src/mesh/impls/bout/boutmesh.hxx"
#include "bout/mesh.hxx"
#include "output.hxx"
#include "unused.hxx"
#include "bout/griddata.hxx"

#include "test_extras.hxx"

TEST(BoutMeshTest, NullOptionsCheck) {
  // Temporarily turn off outputs to make test quiet
  output_info.disable();
  output_warn.disable();
  EXPECT_NO_THROW(BoutMesh mesh(new FakeGridDataSource, nullptr));
  output_info.enable();
  output_warn.enable();
// Not a great test as it's not specific to the thing we want to test,
// and also takes a whopping ~300ms!
TEST(BoutMeshTest, SingleCoreDecomposition) {
  WithQuietOutput info{output_info};
  WithQuietOutput warn{output_warn};
  WithQuietOutput progress{output_progress};

  Options options{};
  options["ny"] = 1;
  options["nx"] = 4;
  options["nz"] = 1;
  options["MXG"] = 1;
  options["MYG"] = 0;

  BoutMesh mesh{new GridFromOptions{&options}, &options};
  EXPECT_NO_THROW(mesh.load());
}
