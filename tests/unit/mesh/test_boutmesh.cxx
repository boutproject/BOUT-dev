#include "gtest/gtest.h"

#include "../src/mesh/impls/bout/boutmesh.hxx"
#include "options.hxx"
#include "output.hxx"
#include "bout/griddata.hxx"

#include "test_extras.hxx"

TEST(BoutMeshTest, NullOptionsCheck) {
  WithQuietOutput info{output_info};
  WithQuietOutput warn{output_warn};

  EXPECT_NO_THROW(BoutMesh mesh(new FakeGridDataSource, nullptr));
}

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
