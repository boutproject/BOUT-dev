#include "gtest/gtest.h"

#include "bout/mesh.hxx"
#include "boutexception.hxx"
#include "test_extras.hxx"

/// Test fixture to ensure a different mesh object for each test
///
/// We only test the non-virtual methods in Mesh in these tests
class MeshTest : public ::testing::Test {
public:
  const int nx = 3;
  const int ny = 5;
  const int nz = 7;

  /// Not called just "mesh" to avoid any ambiguity with the global name "mesh"
  FakeMesh testmesh{nx, ny, nz};
};

TEST_F(MeshTest, MakeSingleIndexRegion) {
  auto region_all = testmesh.makeSingleIndexRegion(
      0, testmesh.LocalNx - 1, 0, testmesh.LocalNy - 1, 0, testmesh.LocalNz - 1);
  auto total_size = nx * ny * nz;

  EXPECT_EQ(region_all.size(), total_size);
}

TEST_F(MeshTest, AddRegion) {
  auto region_all = testmesh.makeSingleIndexRegion(
      0, testmesh.LocalNx - 1, 0, testmesh.LocalNy - 1, 0, testmesh.LocalNz - 1);

  int starting_number_of_regions = testmesh.getRegionMap().size();
  EXPECT_EQ(starting_number_of_regions, 0);
  testmesh.addRegion("test_add_region", region_all);
  int ending_number_of_regions = testmesh.getRegionMap().size();
  EXPECT_EQ(ending_number_of_regions, 1);
}

TEST_F(MeshTest, GetRegion) {
  auto region_all = testmesh.makeSingleIndexRegion(
      0, testmesh.LocalNx - 1, 0, testmesh.LocalNy - 1, 0, testmesh.LocalNz - 1);

  testmesh.addRegion("test_get_region", region_all);

  auto returned_region = testmesh.getRegion("test_get_region");

  EXPECT_EQ(region_all.size(), returned_region.size());
}

TEST_F(MeshTest, GetNonExistantRegion) {
  EXPECT_THROW(testmesh.getRegion("this will throw"), BoutException);
}

TEST_F(MeshTest, CreateDefaultRegions) {
  int starting_number_of_regions = testmesh.getRegionMap().size();
  EXPECT_EQ(starting_number_of_regions, 0);

  testmesh.createDefaultRegions();

  int ending_number_of_regions = testmesh.getRegionMap().size();
  EXPECT_EQ(ending_number_of_regions, 4);

  auto region_map = testmesh.getRegionMap();
  EXPECT_EQ(region_map.count("RGN_ALL"), 1);
  EXPECT_EQ(region_map.count("RGN_NOBNDRY"), 1);
  EXPECT_EQ(region_map.count("RGN_NOX"), 1);
  EXPECT_EQ(region_map.count("RGN_NOY"), 1);
}
