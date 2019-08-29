#include <set>
#include <vector>
#include <tuple>

#include "gtest/gtest.h"
#include "test_extras.hxx"
#include "fake_parallel_mesh.hxx"

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
class FakeParallelIndexer : public GlobalIndexer {
public:
  FakeParallelIndexer(Mesh* localmesh) : GlobalIndexer(localmesh) {}

private:
  virtual void registerFieldForTest(FieldData& f) {
    auto* meshPtr = static_cast<FakeParallelMesh*>(fieldmesh);
    if (meshPtr) {
      int idnum = f.is3D() ? 0 : 1;
      meshPtr->registerField(f, idnum);
    }
  }
  virtual void registerFieldsForTest(FieldPerp& f) {
    auto* meshPtr = static_cast<FakeParallelMesh*>(fieldmesh);
    if (meshPtr) {
      meshPtr->registerField(f, 2);
    }
  }
};

class ParallelIndexerTest : public ::testing::TestWithParam<std::tuple<int, int, int, int> > {
public:
  ParallelIndexerTest() {
    nxpe = std::get<0>(GetParam());
    nype = std::get<1>(GetParam());
    pe_xind = std::get<2>(GetParam());
    pe_yind = std::get<3>(GetParam());
    meshes = createFakeProcessors(nx, ny, nz, nxpe, nype);
    xstart = meshes[0].xstart;
    xend = meshes[0].xend;
    ystart = meshes[0].ystart;
    yend = meshes[0].yend;
    int i = 0;
    for (auto &ind : indexers) {
      ind = std::make_unique<FakeParallelIndexer>(&meshes[i]);
      ind->initialiseTest();
      i++;
    }
    i = 0;
    for (auto &ind : indexers) {
      ind->initialise();
      i++;
    }
    localmesh = &meshes[nype*pe_yind + pe_xind];
  }

  static constexpr int nx = 7;
  static constexpr int ny = 5;
  static constexpr int nz = 3;

  int nxpe, nype;
  int pe_xind, pe_yind;
  int xstart, xend, ystart, yend;
  std::vector<FakeParallelMesh> meshes;
  std::vector<IndexerPtr> indexers;
  Mesh* localmesh;

  IndexerPtr getIndexer(int xind, int yind) {
    return indexers[nype*yind + xind];
  }
};

std::vector<std::tuple<int, int, int, int> > makeParallelTestCases(int nxpe, int nype) {
  std::vector<std::tuple<int, int, int, int> > cases;
  for (int i = 0; i < nxpe; i++) {
    for (int j = 0; j < nype; j++) {
      cases.push_back(std::make_tuple(nxpe, nype, i, j));
    }
  }
  return cases;
}

INSTANTIATE_TEST_SUITE_P(1by3, ParallelIndexerTest,
			 testing::ValuesIn(makeParallelTestCases(1, 3)));
INSTANTIATE_TEST_SUITE_P(3by1, ParallelIndexerTest,
			 testing::ValuesIn(makeParallelTestCases(3, 1)));
INSTANTIATE_TEST_SUITE_P(3by3, ParallelIndexerTest,
			 testing::ValuesIn(makeParallelTestCases(3, 3)));

TEST_P(ParallelIndexerTest, TestConvertIndex3D) {
  IndexerPtr index = getIndexer(pe_xind, pe_yind);
  localmesh->addRegion3D("RGN_XGUARDS", mask(localmesh->getRegion3D("RGN_ALL"),
                                        localmesh->getRegion3D("RGN_NOX")));
  localmesh->addRegion3D("RGN_YGUARDS", mask(localmesh->getRegion3D("RGN_ALL"),
                                        localmesh->getRegion3D("RGN_NOY")));

  // Test xIn boundary
  if (pe_xind > 0) {
    IndexerPtr xInIndex = getIndexer(pe_xind - 1, pe_yind);
    BOUT_FOR(i, localmesh->getRegion3D("RGN_XGUARDS")) {
      if (i.x() < xstart) {
        int global = index->getGlobal(i);
	int otherGlobal = xInIndex->getGlobal(i.xp(xend - xstart));
	ASSERT_EQ(global, otherGlobal);
      }
    }
  }
  
  // Test xOut boundary
  if (pe_xind < nxpe - 1) {
    IndexerPtr xOutIndex = getIndexer(pe_xind + 1, pe_yind);
    BOUT_FOR(i, localmesh->getRegion3D("RGN_XGUARDS")) {
      if (i.x() >= xend) {
        int global = index->getGlobal(i);
	int otherGlobal = xOutIndex->getGlobal(i.xm(xend - xstart));
	ASSERT_EQ(global, otherGlobal);
      }
    }
  }
  
  // Test yDown boundary
  if (pe_yind > 0) {
    IndexerPtr yDownIndex = getIndexer(pe_xind, pe_yind - 1);
    BOUT_FOR(i, localmesh->getRegion3D("RGN_YGUARDS")) {
      if (i.y() < ystart) {
        int global = index->getGlobal(i);
	int otherGlobal = yDownIndex->getGlobal(i.yp(yend - ystart));
	ASSERT_EQ(global, otherGlobal);
      }
    }
    
  }
  
  // Test yUp boundary
  if (pe_yind < nype - 1) {
    IndexerPtr yUpIndex = getIndexer(pe_xind, pe_yind + 1);
    BOUT_FOR(i, localmesh->getRegion3D("RGN_YGUARDS")) {
      if (i.y() >= yend) {
        int global = index->getGlobal(i);
	int otherGlobal = yUpIndex->getGlobal(i.ym(yend - ystart));
	ASSERT_EQ(global, otherGlobal);
      }
    }
    
  }
}

TEST_P(ParallelIndexerTest, TestConvertIndex2D) {
  IndexerPtr index = getIndexer(pe_xind, pe_yind);
  localmesh->addRegion2D("RGN_XGUARDS", mask(localmesh->getRegion2D("RGN_ALL"),
                                        localmesh->getRegion2D("RGN_NOX")));
  localmesh->addRegion2D("RGN_YGUARDS", mask(localmesh->getRegion2D("RGN_ALL"),
                                        localmesh->getRegion2D("RGN_NOY")));

  // Test xIn boundary
  if (pe_xind > 0) {
    IndexerPtr xInIndex = getIndexer(pe_xind - 1, pe_yind);
    BOUT_FOR(i, localmesh->getRegion2D("RGN_XGUARDS")) {
      if (i.x() < xstart) {
        int global = index->getGlobal(i);
	int otherGlobal = xInIndex->getGlobal(i.xp(xend - xstart));
	ASSERT_EQ(global, otherGlobal);
      }
    }
  }
  
  // Test xOut boundary
  if (pe_xind < nxpe - 1) {
    IndexerPtr xOutIndex = getIndexer(pe_xind + 1, pe_yind);
    BOUT_FOR(i, localmesh->getRegion2D("RGN_XGUARDS")) {
      if (i.x() >= xend) {
        int global = index->getGlobal(i);
	int otherGlobal = xOutIndex->getGlobal(i.xm(xend - xstart));
	ASSERT_EQ(global, otherGlobal);
      }
    }
  }
  
  // Test yDown boundary
  if (pe_yind > 0) {
    IndexerPtr yDownIndex = getIndexer(pe_xind, pe_yind - 1);
    BOUT_FOR(i, localmesh->getRegion2D("RGN_YGUARDS")) {
      if (i.y() < ystart) {
        int global = index->getGlobal(i);
	int otherGlobal = yDownIndex->getGlobal(i.yp(yend - ystart));
	ASSERT_EQ(global, otherGlobal);
      }
    }
    
  }
  
  // Test yUp boundary
  if (pe_yind < nype - 1) {
    IndexerPtr yUpIndex = getIndexer(pe_xind, pe_yind + 1);
    BOUT_FOR(i, localmesh->getRegion2D("RGN_YGUARDS")) {
      if (i.y() >= yend) {
        int global = index->getGlobal(i);
	int otherGlobal = yUpIndex->getGlobal(i.ym(yend - ystart));
	ASSERT_EQ(global, otherGlobal);
      }
    }
  }
}


TEST_P(ParallelIndexerTest, TestConvertIndexPerp) {
  IndexerPtr index = getIndexer(pe_xind, pe_yind);
  localmesh->addRegionPerp("RGN_XGUARDS", mask(localmesh->getRegionPerp("RGN_ALL"),
                                        localmesh->getRegionPerp("RGN_NOX")));

  // Test xIn boundary
  if (pe_xind > 0) {
    IndexerPtr xInIndex = getIndexer(pe_xind - 1, pe_yind);
    BOUT_FOR(i, localmesh->getRegionPerp("RGN_XGUARDS")) {
      if (i.x() < xstart) {
        int global = index->getGlobal(i);
	int otherGlobal = xInIndex->getGlobal(i.xp(xend - xstart));
	ASSERT_EQ(global, otherGlobal);
      }
    }
  }
  
  // Test xOut boundary
  if (pe_xind < nxpe - 1) {
    IndexerPtr xOutIndex = getIndexer(pe_xind + 1, pe_yind);
    BOUT_FOR(i, localmesh->getRegionPerp("RGN_XGUARDS")) {
      if (i.x() >= xend) {
        int global = index->getGlobal(i);
	int otherGlobal = xOutIndex->getGlobal(i.xm(xend - xstart));
	ASSERT_EQ(global, otherGlobal);
      }
    }
  }
}

#endif // BOUT_HAS_PETSC
