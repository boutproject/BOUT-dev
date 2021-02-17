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

TEST(BoutMeshTest, SingleCoreYDecomposition) {
  const int total_processors = 1;
  const int num_y_processors = 1;
  const int ny = 1;
  const int num_local_y_points = ny / num_y_processors;
  const int num_y_guards = 1;
  const int jyseps1_1 = -1;
  const int jyseps2_1 = ny / 2;
  const int jyseps1_2 = ny / 2;
  const int jyseps2_2 = ny - 1;
  const int ny_inner = ny / 2;

  auto result = bout::checkBoutMeshYDecomposition(
      total_processors, num_y_processors, num_local_y_points, ny, num_y_guards, jyseps1_1,
      jyseps2_1, jyseps1_2, jyseps2_2, ny_inner);

  EXPECT_TRUE(result.success);
}

TEST(BoutMeshTest, BadTwoCoreYDecomposition) {
  const int total_processors = 2;
  const int num_y_processors = 2;
  const int ny = 1;
  const int num_local_y_points = ny / num_y_processors;
  const int num_y_guards = 1;
  const int jyseps1_1 = -1;
  const int jyseps2_1 = ny / 2;
  const int jyseps1_2 = ny / 2;
  const int jyseps2_2 = ny - 1;
  const int ny_inner = ny / 2;

  auto result = bout::checkBoutMeshYDecomposition(
      total_processors, num_y_processors, num_local_y_points, ny, num_y_guards, jyseps1_1,
      jyseps2_1, jyseps1_2, jyseps2_2, ny_inner);

  EXPECT_FALSE(result.success);
  EXPECT_FALSE(result.reason.empty());
}

TEST(BoutMeshTest, TwoCoreYDecomposition) {
  const int total_processors = 2;
  const int num_y_processors = 2;
  const int ny = 2;
  const int num_local_y_points = ny / num_y_processors;
  const int num_y_guards = 1;
  const int jyseps1_1 = -1;
  const int jyseps2_1 = ny / 2;
  const int jyseps1_2 = ny / 2;
  const int jyseps2_2 = ny - 1;
  const int ny_inner = ny / 2;

  auto result = bout::checkBoutMeshYDecomposition(
      total_processors, num_y_processors, num_local_y_points, ny, num_y_guards, jyseps1_1,
      jyseps2_1, jyseps1_2, jyseps2_2, ny_inner);

  EXPECT_TRUE(result.success);
}
