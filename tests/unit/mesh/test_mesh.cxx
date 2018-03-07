#include "gtest/gtest.h"

#include "bout/mesh.hxx"
#include "boutexception.hxx"

#include "test_extras.hxx"

/// Test fixture to make sure the global mesh is our fake one
class MeshTest : public ::testing::Test {
public:
  MeshTest() : localmesh(nx, ny, nz) {}
  static const int nx = 3;
  static const int ny = 5;
  static const int nz = 7;
  FakeMesh localmesh;
};

TEST_F(MeshTest, CreateDefaultRegions) {
  EXPECT_NO_THROW(localmesh.createDefaultRegions());
  EXPECT_THROW(localmesh.createDefaultRegions(), BoutException);
}

TEST_F(MeshTest, GetRegionFromMesh) {
  localmesh.createDefaultRegions();
  EXPECT_NO_THROW(localmesh.getRegion("RGN_ALL"));
  EXPECT_NO_THROW(localmesh.getRegion("RGN_NOBNDRY"));
  EXPECT_NO_THROW(localmesh.getRegion("RGN_NOX"));
  EXPECT_NO_THROW(localmesh.getRegion("RGN_NOY"));
  EXPECT_THROW(localmesh.getRegion("SOME_MADE_UP_REGION_NAME"), BoutException);
}

TEST_F(MeshTest, GetRegion3DFromMesh) {
  localmesh.createDefaultRegions();
  EXPECT_NO_THROW(localmesh.getRegion3D("RGN_ALL"));
  EXPECT_NO_THROW(localmesh.getRegion3D("RGN_NOBNDRY"));
  EXPECT_THROW(localmesh.getRegion3D("SOME_MADE_UP_REGION_NAME"), BoutException);
}

TEST_F(MeshTest, GetRegion2DFromMesh) {
  localmesh.createDefaultRegions();
  EXPECT_NO_THROW(localmesh.getRegion2D("RGN_ALL"));
  EXPECT_NO_THROW(localmesh.getRegion2D("RGN_NOBNDRY"));
  EXPECT_THROW(localmesh.getRegion2D("SOME_MADE_UP_REGION_NAME"), BoutException);
}

TEST_F(MeshTest, AddRegionToMesh) {
  Region<Ind3D> junk(0, 0, 0, 0, 0, 0, 1, 1);
  EXPECT_NO_THROW(localmesh.addRegion("RGN_JUNK", junk));
  EXPECT_THROW(localmesh.addRegion("RGN_JUNK", junk), BoutException);
  Region<Ind2D> junk2D(0, 0, 0, 0, 0, 0, 1, 1);
  EXPECT_NO_THROW(localmesh.addRegion("RGN_JUNK", junk2D));
}

TEST_F(MeshTest, AddRegion3DToMesh) {
  Region<Ind3D> junk(0, 0, 0, 0, 0, 0, 1, 1);
  EXPECT_NO_THROW(localmesh.addRegion3D("RGN_JUNK_3D", junk));
  EXPECT_THROW(localmesh.addRegion3D("RGN_JUNK_3D", junk), BoutException);
}

TEST_F(MeshTest, AddRegion2DToMesh) {
  Region<Ind2D> junk(0, 0, 0, 0, 0, 0, 1, 1);
  EXPECT_NO_THROW(localmesh.addRegion2D("RGN_JUNK_2D", junk));
  EXPECT_THROW(localmesh.addRegion2D("RGN_JUNK_2D", junk), BoutException);
}
