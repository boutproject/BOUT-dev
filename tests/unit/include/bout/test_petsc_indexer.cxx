#include <set>

#include "gtest/gtest.h"
#include "test_extras.hxx"

#include "bout/petsc_interface.hxx"
#include "bout/region.hxx"

#ifdef BOUT_HAS_PETSC

/// Global mesh
namespace bout{
namespace globals{
extern Mesh *mesh;
} // namespace globals
} // namespace bout

// The unit tests use the global mesh
using namespace bout::globals;

using IndexerTest = FakeMeshFixture;

TEST_F(IndexerTest, TestIndexerConstructor){
  FakeMesh mesh2(1, 1, 1);
  IndexerPtr i1 = GlobalIndexer::getInstance(mesh);
  IndexerPtr i2 = GlobalIndexer::getInstance(mesh);
  IndexerPtr i3 = GlobalIndexer::getInstance(&mesh2);
  IndexerPtr i4 = GlobalIndexer::getInstance(&mesh2);

  // Check only one instance for the global mesh
  EXPECT_EQ(&i1, &i2);
  // Check all other instances are unique
  EXPECT_NE(&i1, &i3);
  EXPECT_NE(&i3, &i4);
}


class SingleProcIndexerTest : public testing::TestWithParam<Mesh*> {
};

class FakeNonperiodicMesh : public FakeMesh {
public:
  FakeNonperiodicMesh(int nx, int ny, int nz): FakeMesh(nx, ny, nz) {
    periodicX = false;
  }
  bool periodicY(int UNUSED(jx)) const { return false; }
};

class FakeXPeriodicMesh : public FakeMesh {
public:
  FakeXPeriodicMesh(int nx, int ny, int nz): FakeMesh(nx, ny, nz) {
    periodicX = true;
  }
  bool periodicY(int UNUSED(jx)) const { return false; }
};

using FakeYPeriodicMesh = FakeMesh;

class FakePeriodicMesh : public FakeMesh {
public:
  FakePeriodicMesh(int nx, int ny, int nz): FakeMesh(nx, ny, nz) {
    periodicX = true;
  }
  bool periodicY(int UNUSED(jx)) const { return true; }
};

FakeNonperiodicMesh mesh_nonperiod(FakeMeshFixture::nx, FakeMeshFixture::ny,
				   FakeMeshFixture::nz);
FakeXPeriodicMesh mesh_Xperiod(FakeMeshFixture::nx, FakeMeshFixture::ny,
				   FakeMeshFixture::nz);
FakeYPeriodicMesh mesh_Yperiod(FakeMeshFixture::nx, FakeMeshFixture::ny,
				   FakeMeshFixture::nz);
FakePeriodicMesh mesh_XYperiod(FakeMeshFixture::nx, FakeMeshFixture::ny,
				   FakeMeshFixture::nz);

INSTANTIATE_TEST_SUITE_P(UniqueIndices, SingleProcIndexerTest,
			 testing::Values(&mesh_nonperiod, &mesh_Xperiod,
					 &mesh_Yperiod, &mesh_XYperiod));

TEST_P(SingleProcIndexerTest, TestConvertIndex2D) {
  Mesh* localmesh = GetParam();
  IndexerPtr index = GlobalIndexer::getInstance(localmesh);
  set <int, greater <int>> returnedIndices;
  localmesh->addRegion2D("RGN_XGUARDS", mask(localmesh->getRegion2D("RGN_ALL"),
                                        localmesh->getRegion2D("RGN_NOX")));
  localmesh->addRegion2D("RGN_YGUARDS", mask(localmesh->getRegion2D("RGN_ALL"),
                                        localmesh->getRegion2D("RGN_NOY")));

  // Check each of the interior global indices is unique
  BOUT_FOR(i, localmesh->getRegion2D("RGN_NOBNDRY")) {
    int global = index->getGlobal(i);
    EXPECT_GE(global, 0);
    EXPECT_EQ(returnedIndices.count(global), 0);
    returnedIndices.insert(global);
  }

  // If periodic in X, check indices of X guard cells are not unique;
  // otherwise check that they are
  BOUT_FOR(i, localmesh->getRegion2D("RGN_XGUARDS")) {
    int global = index->getGlobal(i);
    EXPECT_GE(global, 0);
    if (localmesh->periodicX) {
      EXPECT_NE(returnedIndices.count(global), 0);
    } else {
      EXPECT_EQ(returnedIndices.count(global), 0);
      returnedIndices.insert(global);
    }
  }

  // If periodic in Y, check indices of Y guard cells are not unique;
  // otherwise check that they are
  BOUT_FOR(i, localmesh->getRegion2D("RGN_YGUARDS")) {
    int global = index->getGlobal(i);
    EXPECT_GE(global, 0);
    if (localmesh->periodicY(i.y())) {
      EXPECT_NE(returnedIndices.count(global), 0);
    } else {
      EXPECT_EQ(returnedIndices.count(global), 0);
      returnedIndices.insert(global);
    }
  }

  ASSERT_LT(*returnedIndices.rbegin(), localmesh->LocalNx*localmesh->LocalNy);
}


TEST_P(SingleProcIndexerTest, TestConvertIndex3D) {
  Mesh* localmesh = GetParam();
  IndexerPtr index = GlobalIndexer::getInstance(localmesh);
  set <int, greater <int>> returnedIndices;
  localmesh->addRegion3D("RGN_XGUARDS", mask(localmesh->getRegion3D("RGN_ALL"),
                                        localmesh->getRegion3D("RGN_NOX")));
  localmesh->addRegion3D("RGN_YGUARDS", mask(localmesh->getRegion3D("RGN_ALL"),
                                        localmesh->getRegion3D("RGN_NOY")));

  // Check each of the interior global indices is unique
  BOUT_FOR(i, localmesh->getRegion3D("RGN_NOBNDRY")) {
    int global = index->getGlobal(i);
    EXPECT_EQ(returnedIndices.count(global), 0);
    returnedIndices.insert(global);
  }

  // If periodic in X, check indices of X guard cells are not unique;
  // otherwise check that they are
  BOUT_FOR(i, localmesh->getRegion3D("RGN_XGUARDS")) {
    int global = index->getGlobal(i);
    if (localmesh->periodicX) {
      EXPECT_NE(returnedIndices.count(global), 0);
    } else {
      EXPECT_EQ(returnedIndices.count(global), 0);
      returnedIndices.insert(global);
    }
  }

  // If periodic in Y, check indices of Y guard cells are not unique;
  // otherwise check that they are
  BOUT_FOR(i, localmesh->getRegion3D("RGN_YGUARDS")) {
    int global = index->getGlobal(i);
    if (localmesh->periodicY(i.y())) {
      EXPECT_NE(returnedIndices.count(global), 0);
    } else {
      EXPECT_EQ(returnedIndices.count(global), 0);
      returnedIndices.insert(global);
    }
  }
  ASSERT_LT(*returnedIndices.rbegin(),
	    localmesh->LocalNx*localmesh->LocalNy*localmesh->LocalNz);
}


TEST_P(SingleProcIndexerTest, TestConvertIndexPerp) {
  Mesh* localmesh = GetParam();
  IndexerPtr index = GlobalIndexer::getInstance(localmesh);
  set <int, greater <int>> returnedIndices;
  localmesh->addRegionPerp("RGN_XGUARDS", mask(localmesh->getRegionPerp("RGN_ALL"),
					       localmesh->getRegionPerp("RGN_NOX")));

  // Check each of the interior global indices is unique
  BOUT_FOR(i, localmesh->getRegionPerp("RGN_NOBNDRY")) {
    int global = index->getGlobal(i);
    EXPECT_EQ(returnedIndices.count(global), 0);
    returnedIndices.insert(global);
  }

  // If periodic in X, check indices of X guard cells are not unique;
  // otherwise check that they are
  BOUT_FOR(i, localmesh->getRegionPerp("RGN_XGUARDS")) {
    int global = index->getGlobal(i);
    if (localmesh->periodicX) {
      EXPECT_NE(returnedIndices.count(global), 0);
    } else {
      EXPECT_EQ(returnedIndices.count(global), 0);
      returnedIndices.insert(global);
    }
  }

  ASSERT_LT(*returnedIndices.rbegin(), localmesh->LocalNx*localmesh->LocalNz);
}


// TODO: Add tests for communicating indices across processors, somehow...

#endif // BOUT_HAS_PETSC
