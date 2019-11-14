#include <set>
#include <tuple>
#include <vector>

#include "fake_parallel_mesh.hxx"
#include "test_extras.hxx"
#include "gtest/gtest.h"

#include "bout/petsc_interface.hxx"
#include "bout/region.hxx"

#ifdef BOUT_HAS_PETSC

/// Global mesh
namespace bout {
namespace globals {
extern Mesh* mesh;
} // namespace globals
} // namespace bout

// The unit tests use the global mesh
using namespace bout::globals;

using IndexerTest = FakeMeshFixture;

TEST_F(IndexerTest, TestIndexerConstructor) {
  FakeMesh mesh2(2, 2, 2);
  mesh2.createDefaultRegions();
  mesh2.setCoordinates(nullptr);
  test_coords = std::make_shared<Coordinates>(
      bout::globals::mesh, Field2D{1.0}, Field2D{1.0}, BoutReal{1.0}, Field2D{1.0},
      Field2D{0.0}, Field2D{1.0}, Field2D{1.0}, Field2D{1.0}, Field2D{0.0}, Field2D{0.0},
      Field2D{0.0}, Field2D{1.0}, Field2D{1.0}, Field2D{1.0}, Field2D{0.0}, Field2D{0.0},
      Field2D{0.0}, Field2D{0.0}, Field2D{0.0}, false);
  mesh2.setCoordinates(test_coords);
  // May need a ParallelTransform to create fields, because create3D calls
  // fromFieldAligned
  test_coords->setParallelTransform(
      bout::utils::make_unique<ParallelTransformIdentity>(*bout::globals::mesh));
  mesh2.createBoundaryRegions();

  IndexerPtr i1 = GlobalIndexer::getInstance(mesh);
  IndexerPtr i2 = GlobalIndexer::getInstance(mesh);
  IndexerPtr i3 = GlobalIndexer::getInstance(&mesh2);
  IndexerPtr i4 = GlobalIndexer::getInstance(&mesh2);

  // Check only one instance for the global mesh
  EXPECT_EQ(i1.get(), i2.get());
  // Check all other instances are unique
  EXPECT_NE(i1.get(), i3.get());
  EXPECT_NE(i3.get(), i4.get());

  GlobalIndexer::recreateGlobalInstance();
}

TEST_F(IndexerTest, TestConvertIndex3D) {
  Mesh* localmesh = bout::globals::mesh;
  IndexerPtr index = GlobalIndexer::getInstance(localmesh);
  std::set<int, std::greater<int>> returnedIndices;

  // Check each of the interior global indices is unique
  BOUT_FOR(i, localmesh->getRegion3D("RGN_NOBNDRY")) {
    int global = index->getGlobal(i);
    EXPECT_GE(global, 0);
    EXPECT_EQ(returnedIndices.count(global), 0);
    BOUT_OMP(critical) returnedIndices.insert(global);
  }

  // Check indices of X guard cells are unique
  BOUT_FOR(i, localmesh->getRegion3D("RGN_XGUARDS")) {
    int global = index->getGlobal(i);
    EXPECT_GE(global, 0);
    EXPECT_EQ(returnedIndices.count(global), 0);
    BOUT_OMP(critical) returnedIndices.insert(global);
  }

  // Check indices of Y guard cells are unique
  BOUT_FOR(i, localmesh->getRegion3D("RGN_YGUARDS")) {
    int global = index->getGlobal(i);
    EXPECT_GE(global, 0);
    EXPECT_EQ(returnedIndices.count(global), 0);
    BOUT_OMP(critical) returnedIndices.insert(global);
  }
  ASSERT_LT(*returnedIndices.rbegin(),
            localmesh->LocalNx * localmesh->LocalNy * localmesh->LocalNz);

  GlobalIndexer::recreateGlobalInstance();
}

TEST_F(IndexerTest, TestConvertIndex2D) {
  Mesh* localmesh = bout::globals::mesh;
  IndexerPtr index = GlobalIndexer::getInstance(localmesh);
  std::set<int, std::greater<int>> returnedIndices;

  // Check each of the interior global indices is unique
  BOUT_FOR(i, localmesh->getRegion2D("RGN_NOBNDRY")) {
    int global = index->getGlobal(i);
    EXPECT_GE(global, 0);
    EXPECT_EQ(returnedIndices.count(global), 0);
    BOUT_OMP(critical) returnedIndices.insert(global);
  }

  // If periodic in X, check indices of X guard cells are not unique;
  // otherwise check that they are
  BOUT_FOR(i, localmesh->getRegion2D("RGN_XGUARDS")) {
    int global = index->getGlobal(i);
    EXPECT_GE(global, 0);
    EXPECT_EQ(returnedIndices.count(global), 0);
    BOUT_OMP(critical) returnedIndices.insert(global);
  }

  // If periodic in Y, check indices of Y guard cells are not unique;
  // otherwise check that they are
  BOUT_FOR(i, localmesh->getRegion2D("RGN_YGUARDS")) {
    int global = index->getGlobal(i);
    EXPECT_GE(global, 0);
    EXPECT_EQ(returnedIndices.count(global), 0);
    BOUT_OMP(critical) returnedIndices.insert(global);
  }

  ASSERT_LT(*returnedIndices.rbegin(), localmesh->LocalNx * localmesh->LocalNy);

  GlobalIndexer::recreateGlobalInstance();
}

TEST_F(IndexerTest, TestConvertIndexPerp) {
  Mesh* localmesh = bout::globals::mesh;
  IndexerPtr index = GlobalIndexer::getInstance(localmesh);
  std::set<int, std::greater<int>> returnedIndices;

  // Check each of the interior global indices is unique
  BOUT_FOR(i, localmesh->getRegionPerp("RGN_NOBNDRY")) {
    int global = index->getGlobal(i);
    EXPECT_EQ(returnedIndices.count(global), 0);
    BOUT_OMP(critical) returnedIndices.insert(global);
  }

  // If periodic in X, check indices of X guard cells are not unique;
  // otherwise check that they are
  BOUT_FOR(i, localmesh->getRegionPerp("RGN_XGUARDS")) {
    int global = index->getGlobal(i);
    EXPECT_EQ(returnedIndices.count(global), 0);
    BOUT_OMP(critical) returnedIndices.insert(global);
  }

  ASSERT_LT(*returnedIndices.rbegin(), localmesh->LocalNx * localmesh->LocalNz);

  GlobalIndexer::recreateGlobalInstance();
}

// TODO: Add tests for communicating indices across processors, somehow...
class FakeParallelIndexer : public GlobalIndexer {
public:
  FakeParallelIndexer(Mesh* localmesh) : GlobalIndexer(localmesh) {}
  void initialiseTest() {
    registerFieldForTest(getIndices3D());
    registerFieldForTest(getIndices2D());
    registerFieldForTest(getIndicesPerp());
  }

private:
  virtual void registerFieldForTest(FieldData& f) {
    auto* meshPtr = static_cast<FakeParallelMesh*>(getMesh());
    if (meshPtr) {
      int idnum = f.is3D() ? 0 : 1;
      meshPtr->registerField(f, idnum);
    }
  }
  virtual void registerFieldForTest(FieldPerp& f) {
    auto* meshPtr = static_cast<FakeParallelMesh*>(getMesh());
    if (meshPtr) {
      meshPtr->registerField(f, 2);
    }
  }
};

class ParallelIndexerTest
    : public ::testing::TestWithParam<std::tuple<int, int, int, int>> {
public:
  WithQuietOutput info{output_info}, warn{output_warn}, all{output};

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
    for (i = 0; i < nxpe * nype; i++) {
      auto ind = std::make_shared<FakeParallelIndexer>(&meshes[i]);
      ind->initialiseTest();
      indexers.push_back(ind);
    }
    //    i = 0;
    //    for (auto &ind : indexers) {
    //      ind->initialise();
    //      i++;
    //    }
    indexers[pe_yind + pe_xind * nype]->initialise();
    localmesh = &meshes[pe_yind + pe_xind * nype];
  }

  virtual ~ParallelIndexerTest() {
    bout::globals::mesh = nullptr;
    GlobalIndexer::recreateGlobalInstance();
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
    return indexers[(yind + nype) % nype + xind * nype];
  }
};

std::vector<std::tuple<int, int, int, int>> makeParallelTestCases(int nxpe, int nype) {
  std::vector<std::tuple<int, int, int, int>> cases;
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

  // Test xIn boundary
  BOUT_FOR(i, localmesh->getRegion3D("RGN_XGUARDS")) {
    if (i.x() < xstart && i.y() >= ystart && i.y() <= yend) {
      int global = index->getGlobal(i);
      if (pe_xind > 0) {
        int otherGlobal =
            getIndexer(pe_xind - 1, pe_yind)->getGlobal(i.xp(xend - xstart + 1));
        EXPECT_EQ(global, otherGlobal);
      } else {
        EXPECT_NE(global, -1);
      }
    }
  }

  // Test xOut boundary
  BOUT_FOR(i, localmesh->getRegion3D("RGN_XGUARDS")) {
    if (i.x() >= xend && i.y() >= ystart && i.y() <= yend) {
      int global = index->getGlobal(i);
      if (pe_xind < nxpe - 1) {
        int otherGlobal =
            getIndexer(pe_xind + 1, pe_yind)->getGlobal(i.xm(xend - xstart + 1));
        EXPECT_EQ(global, otherGlobal);
      } else {
        EXPECT_NE(global, -1);
      }
    }
  }

  // Test yDown boundary
  BOUT_FOR(i, localmesh->getRegion3D("RGN_YGUARDS")) {
    if (i.y() < ystart && i.x() >= xstart && i.x() <= xend) {
      int global = index->getGlobal(i);
      int otherGlobal =
          getIndexer(pe_xind, pe_yind - 1)->getGlobal(i.yp(yend - ystart + 1));
      EXPECT_EQ(global, otherGlobal);
    }
  }

  // Test yUp boundary
  BOUT_FOR(i, localmesh->getRegion3D("RGN_YGUARDS")) {
    if (i.y() >= yend && i.x() >= xstart && i.x() <= xend) {
      int global = index->getGlobal(i);
      int otherGlobal =
          getIndexer(pe_xind, pe_yind + 1)->getGlobal(i.ym(yend - ystart + 1));
      EXPECT_EQ(global, otherGlobal);
    }
  }
}

TEST_P(ParallelIndexerTest, TestConvertIndex2D) {
  IndexerPtr index = getIndexer(pe_xind, pe_yind);

  // Test xIn boundary
  BOUT_FOR(i, localmesh->getRegion2D("RGN_XGUARDS")) {
    if (i.x() < xstart && i.y() >= ystart && i.y() <= yend) {
      int global = index->getGlobal(i);
      if (pe_xind > 0) {
        int otherGlobal =
            getIndexer(pe_xind - 1, pe_yind)->getGlobal(i.xp(xend - xstart + 1));
        EXPECT_EQ(global, otherGlobal);
      } else {
        EXPECT_NE(global, -1);
      }
    }
  }

  // Test xOut boundary
  BOUT_FOR(i, localmesh->getRegion2D("RGN_XGUARDS")) {
    if (i.x() >= xend && i.y() >= ystart && i.y() <= yend) {
      int global = index->getGlobal(i);
      if (pe_xind < nxpe - 1) {
        int otherGlobal =
            getIndexer(pe_xind + 1, pe_yind)->getGlobal(i.xm(xend - xstart + 1));
        EXPECT_EQ(global, otherGlobal);
      } else {
        EXPECT_NE(global, -1);
      }
    }
  }

  // Test yDown boundary
  BOUT_FOR(i, localmesh->getRegion2D("RGN_YGUARDS")) {
    if (i.y() < ystart && i.x() >= xstart && i.x() <= xend) {
      int global = index->getGlobal(i);
      int otherGlobal =
          getIndexer(pe_xind, pe_yind - 1)->getGlobal(i.yp(yend - ystart + 1));
      EXPECT_EQ(global, otherGlobal);
    }
  }

  // Test yUp boundary
  BOUT_FOR(i, localmesh->getRegion2D("RGN_YGUARDS")) {
    if (i.y() >= yend && i.x() >= xstart && i.x() <= xend) {
      int global = index->getGlobal(i);
      int otherGlobal =
          getIndexer(pe_xind, pe_yind + 1)->getGlobal(i.ym(yend - ystart + 1));
      EXPECT_EQ(global, otherGlobal);
    }
  }
}

TEST_P(ParallelIndexerTest, TestConvertIndexPerp) {
  IndexerPtr index = getIndexer(pe_xind, pe_yind);

  // Test xIn boundary
  BOUT_FOR(i, localmesh->getRegionPerp("RGN_XGUARDS")) {
    if (i.x() < xstart) {
      int global = index->getGlobal(i);
      if (pe_xind > 0) {
        int otherGlobal =
            getIndexer(pe_xind - 1, pe_yind)->getGlobal(i.xp(xend - xstart + 1));
        EXPECT_EQ(global, otherGlobal);
      } else {
        EXPECT_NE(global, -1);
      }
    }
  }

  // Test xOut boundary
  BOUT_FOR(i, localmesh->getRegionPerp("RGN_XGUARDS")) {
    if (i.x() >= xend) {
      int global = index->getGlobal(i);
      if (pe_xind < nxpe - 1) {
        int otherGlobal =
            getIndexer(pe_xind + 1, pe_yind)->getGlobal(i.xm(xend - xstart + 1));
        EXPECT_EQ(global, otherGlobal);
      } else {
        EXPECT_NE(global, -1);
      }
    }
  }
}

#endif // BOUT_HAS_PETSC
