#include "gtest/gtest.h"

#include "bout/mesh.hxx"
#include "bout/region.hxx"
#include "boutexception.hxx"
#include "output.hxx"

#include "test_extras.hxx"

/// Test fixture to make sure the global mesh is our fake one
class MeshTest : public ::testing::Test {
protected:
  static void SetUpTestCase() { output_info.disable(); }

  static void TearDownTestCase() { output_info.enable(); }

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

TEST_F(MeshTest, Ind2DTo3D) {
  Ind2D index2d_0(0);
  Ind2D index2d_7(7);
  Ind2D index2d_14(14);

  EXPECT_EQ(localmesh.ind2Dto3D(index2d_0, 0), Ind3D(0));
  EXPECT_EQ(localmesh.ind2Dto3D(index2d_0, 1), Ind3D(1));

  EXPECT_EQ(localmesh.ind2Dto3D(index2d_7, 0), Ind3D(49));
  EXPECT_EQ(localmesh.ind2Dto3D(index2d_7, 1), Ind3D(50));

  EXPECT_EQ(localmesh.ind2Dto3D(index2d_14, 0), Ind3D(98));
  EXPECT_EQ(localmesh.ind2Dto3D(index2d_14, 1), Ind3D(99));
}

TEST_F(MeshTest, Ind3DTo2D) {
  Ind3D index3d_0(0);
  Ind3D index3d_49(49);
  Ind3D index3d_98(98);

  EXPECT_EQ(localmesh.ind3Dto2D(index3d_0), Ind2D(0));
  EXPECT_EQ(localmesh.ind3Dto2D(index3d_49), Ind2D(7));
  EXPECT_EQ(localmesh.ind3Dto2D(index3d_98), Ind2D(14));
}

TEST_F(MeshTest, MapInd3DTo2D) {
  Ind3D index3d_0(0);
  Ind3D index3d_49(49);
  Ind3D index3d_98(98);

  EXPECT_EQ(localmesh.ind3Dto2D(index3d_0), 0);
  EXPECT_EQ(localmesh.ind3Dto2D(index3d_49), 7);
  EXPECT_EQ(localmesh.ind3Dto2D(index3d_98), 14);
}
