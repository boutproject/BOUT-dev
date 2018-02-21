#include "gtest/gtest.h"

#include "bout/constants.hxx"
#include "bout/mesh.hxx"
#include "bout/boutexception.hxx"
#include "bout/output.hxx"
#include "test_extras.hxx"
#include "bout/unused.hxx"

#include <cmath>
#include <set>
#include <vector>

/// Global mesh
extern Mesh *mesh;

/// Test fixture to make sure the global mesh is our fake one
class MeshTest : public ::testing::Test {
protected:
  static void SetUpTestCase() {
    // Delete any existing mesh
    if (mesh != nullptr) {
      delete mesh;
      mesh = nullptr;
    }
    mesh = new FakeMesh(nx, ny, nz);
    mesh->createDefaultRegions();
  }

  static void TearDownTestCase() {
    delete mesh;
    mesh = nullptr;
  }

public:
  static const int nx;
  static const int ny;
  static const int nz;
};

const int MeshTest::nx = 3;
const int MeshTest::ny = 5;
const int MeshTest::nz = 7;

TEST(MeshTestNoFixture, createDefaultRegions) {
  FakeMesh localmesh(MeshTest::nx, MeshTest::ny, MeshTest::nz);
  EXPECT_NO_THROW(localmesh.createDefaultRegions());
  EXPECT_THROW(localmesh.createDefaultRegions(), BoutException);
}

TEST_F(MeshTest, getRegionFromMesh) {
  EXPECT_NO_THROW(mesh->getRegion("RGN_ALL"));
  EXPECT_NO_THROW(mesh->getRegion("RGN_NOBNDRY"));
  EXPECT_NO_THROW(mesh->getRegion("RGN_NOX"));
  EXPECT_NO_THROW(mesh->getRegion("RGN_NOY"));
  EXPECT_THROW(mesh->getRegion("SOME_MADE_UP_REGION_NAME"), BoutException);
}

TEST_F(MeshTest, getRegion3DFromMesh) {
  EXPECT_NO_THROW(mesh->getRegion3D("RGN_ALL"));
  EXPECT_NO_THROW(mesh->getRegion3D("RGN_NOBNDRY"));
  EXPECT_THROW(mesh->getRegion3D("SOME_MADE_UP_REGION_NAME"), BoutException);
}

TEST_F(MeshTest, getRegion2DFromMesh) {
  EXPECT_NO_THROW(mesh->getRegion2D("RGN_ALL"));
  EXPECT_NO_THROW(mesh->getRegion2D("RGN_NOBNDRY"));
  EXPECT_THROW(mesh->getRegion2D("SOME_MADE_UP_REGION_NAME"), BoutException);
}

TEST_F(MeshTest, addRegionToMesh) {
  Region<Ind3D> junk(0, 0, 0, 0, 0, 0, 1, 1);
  EXPECT_NO_THROW(mesh->addRegion("RGN_JUNK", junk));
  EXPECT_THROW(mesh->addRegion("RGN_JUNK", junk), BoutException);
  Region<Ind2D> junk2D(0, 0, 0, 0, 0, 0, 1, 1);
  EXPECT_NO_THROW(mesh->addRegion("RGN_JUNK", junk2D));
}

TEST_F(MeshTest, addRegion3DToMesh) {
  Region<Ind3D> junk(0, 0, 0, 0, 0, 0, 1, 1);
  EXPECT_NO_THROW(mesh->addRegion3D("RGN_JUNK_3D", junk));
  EXPECT_THROW(mesh->addRegion3D("RGN_JUNK_3D", junk), BoutException);
}

TEST_F(MeshTest, addRegion2DToMesh) {
  Region<Ind2D> junk(0, 0, 0, 0, 0, 0, 1, 1);
  EXPECT_NO_THROW(mesh->addRegion2D("RGN_JUNK_2D", junk));
  EXPECT_THROW(mesh->addRegion2D("RGN_JUNK_2D", junk), BoutException);
}
