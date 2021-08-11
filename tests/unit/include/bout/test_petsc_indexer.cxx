#include "bout/build_config.hxx"

#include <set>
#include <tuple>
#include <vector>

#include "fake_parallel_mesh.hxx"
#include "test_extras.hxx"
#include "gtest/gtest.h"

#include "bout/petsc_interface.hxx"
#include "bout/region.hxx"

#if BOUT_HAS_PETSC

/// Global mesh
namespace bout {
namespace globals {
extern Mesh* mesh;
} // namespace globals
} // namespace bout

// The unit tests use the global mesh
using namespace bout::globals;

template <typename T>
class IndexerTest : public FakeMeshFixture {
public:
  using ind_type = typename T::ind_type;
  GlobalIndexer<T> globalSquareIndexer, globalStarIndexer, globalDefaultIndexer;
  int guardx, guardy, nx, ny, nz;
  FakeMesh mesh2;
  GlobalIndexer<T> localIndexer;

  IndexerTest()
      : FakeMeshFixture(),
        globalSquareIndexer(bout::globals::mesh,
                            squareStencil<ind_type>(bout::globals::mesh)),
        globalStarIndexer(bout::globals::mesh,
                          starStencil<ind_type>(bout::globals::mesh)),
        globalDefaultIndexer(bout::globals::mesh), guardx(bout::globals::mesh->getNXPE()),
        guardy(std::is_same<T, FieldPerp>::value ? 0 : bout::globals::mesh->getNYPE()),
        nx(bout::globals::mesh->LocalNx - 2 * guardx),
        ny(std::is_same<T, FieldPerp>::value ? 1
                                             : bout::globals::mesh->LocalNy - 2 * guardy),
        nz(std::is_same<T, Field2D>::value ? 1 : bout::globals::mesh->LocalNz),
        mesh2(2, 2, 2) {
    mesh2.createDefaultRegions();
    mesh2.setCoordinates(nullptr);
    test_coords = std::make_shared<Coordinates>(
        bout::globals::mesh, Field2D{1.0}, Field2D{1.0}, BoutReal{1.0}, Field2D{1.0},
        Field2D{0.0}, Field2D{1.0}, Field2D{1.0}, Field2D{1.0}, Field2D{0.0},
        Field2D{0.0}, Field2D{0.0}, Field2D{1.0}, Field2D{1.0}, Field2D{1.0},
        Field2D{0.0}, Field2D{0.0}, Field2D{0.0}, Field2D{0.0}, Field2D{0.0});
    // No call to Coordinates::geometry() needed here
    mesh2.setCoordinates(test_coords);
    // May need a ParallelTransform to create fields, because create3D calls
    // fromFieldAligned
    test_coords->setParallelTransform(
        bout::utils::make_unique<ParallelTransformIdentity>(*bout::globals::mesh));
    mesh2.createBoundaryRegions();
    localIndexer = GlobalIndexer<T>(&mesh2);
  }
};

using FieldTypes = ::testing::Types<Field2D, Field3D, FieldPerp>;
TYPED_TEST_SUITE(IndexerTest, FieldTypes);

TYPED_TEST(IndexerTest, TestGetMesh) {
  EXPECT_EQ(this->globalSquareIndexer.getMesh(), bout::globals::mesh);
  EXPECT_EQ(this->globalStarIndexer.getMesh(), bout::globals::mesh);
  EXPECT_EQ(this->globalDefaultIndexer.getMesh(), bout::globals::mesh);
  EXPECT_EQ(this->localIndexer.getMesh(), &(this->mesh2));
}

TYPED_TEST(IndexerTest, TestConvertIndex) {
  TypeParam f(bout::globals::mesh);
  std::set<int, std::greater<int>> indicesGlobalSquare, indicesGlobalStar,
      indicesGlobalDefault;

  // Check each of the interior global indices is unique
  BOUT_FOR(i, f.getRegion("RGN_NOBNDRY")) {
    int global = this->globalSquareIndexer.getGlobal(i);
    EXPECT_GE(global, 0);
    BOUT_OMP(critical) EXPECT_TRUE(indicesGlobalSquare.insert(global).second);
    global = this->globalStarIndexer.getGlobal(i);
    EXPECT_GE(global, 0);
    BOUT_OMP(critical) EXPECT_TRUE(indicesGlobalStar.insert(global).second);
    global = this->globalDefaultIndexer.getGlobal(i);
    EXPECT_GE(global, 0);
    BOUT_OMP(critical) EXPECT_TRUE(indicesGlobalDefault.insert(global).second);
  }

  // Check indices of X guard cells are unique
  BOUT_FOR(i, f.getRegion("RGN_XGUARDS")) {
    int global = this->globalSquareIndexer.getGlobal(i);
    EXPECT_GE(global, 0);
    BOUT_OMP(critical) EXPECT_TRUE(indicesGlobalSquare.insert(global).second);
    global = this->globalStarIndexer.getGlobal(i);
    EXPECT_GE(global, 0);
    BOUT_OMP(critical) EXPECT_TRUE(indicesGlobalStar.insert(global).second);
    EXPECT_LT(this->globalDefaultIndexer.getGlobal(i), 0);
  }

  // Check indices of Y guard cells are unique
  if (!std::is_same<TypeParam, FieldPerp>::value) {
    BOUT_FOR(i, f.getRegion("RGN_YGUARDS")) {
      int global = this->globalSquareIndexer.getGlobal(i);
      EXPECT_GE(global, 0);
      BOUT_OMP(critical) EXPECT_TRUE(indicesGlobalSquare.insert(global).second);
      global = this->globalStarIndexer.getGlobal(i);
      EXPECT_GE(global, 0);
      BOUT_OMP(critical) EXPECT_TRUE(indicesGlobalStar.insert(global).second);
      EXPECT_LT(this->globalDefaultIndexer.getGlobal(i), 0);
    }
  }

  ASSERT_LT(*indicesGlobalSquare.rbegin(),
            this->globalSquareIndexer.getMesh()->LocalNx
                * this->globalSquareIndexer.getMesh()->LocalNy
                * this->globalSquareIndexer.getMesh()->LocalNz);
  ASSERT_LT(*indicesGlobalStar.rbegin(),
            (this->globalSquareIndexer.getMesh()->LocalNx
                 * this->globalSquareIndexer.getMesh()->LocalNy
             - 4)
                * this->globalSquareIndexer.getMesh()->LocalNz);
  ASSERT_LT(*indicesGlobalStar.rbegin(), this->nx * this->ny * this->nz);
}

TYPED_TEST(IndexerTest, TestIsLocal) {
  BOUT_FOR(i, this->globalSquareIndexer.getRegionAll()) {
    EXPECT_TRUE(this->globalSquareIndexer.isLocal(i));
  }
  BOUT_FOR(i, this->globalStarIndexer.getRegionAll()) {
    EXPECT_TRUE(this->globalStarIndexer.isLocal(i));
  }
  BOUT_FOR(i, this->globalDefaultIndexer.getRegionAll()) {
    EXPECT_TRUE(this->globalDefaultIndexer.isLocal(i));
  }
  BOUT_FOR(i, this->localIndexer.getRegionAll()) {
    EXPECT_TRUE(this->localIndexer.isLocal(i));
  }
}

TYPED_TEST(IndexerTest, TestGetGlobalStart) {
  EXPECT_EQ(this->globalSquareIndexer.getGlobalStart(), 0);
  EXPECT_EQ(this->globalStarIndexer.getGlobalStart(), 0);
  EXPECT_EQ(this->globalDefaultIndexer.getGlobalStart(), 0);
  EXPECT_EQ(this->localIndexer.getGlobalStart(), 0);
}

TYPED_TEST(IndexerTest, TestGetRegionAll) {
  EXPECT_EQ(size(this->globalSquareIndexer.getRegionAll()),
            size(this->globalSquareIndexer.getRegionNobndry()
                 + this->globalSquareIndexer.getRegionBndry()));
  EXPECT_EQ(size(this->globalStarIndexer.getRegionAll()),
            size(this->globalStarIndexer.getRegionNobndry()
                 + this->globalStarIndexer.getRegionBndry()));
  EXPECT_EQ(size(this->globalDefaultIndexer.getRegionAll()),
            size(this->globalDefaultIndexer.getRegionNobndry()
                 + this->globalDefaultIndexer.getRegionBndry()));
  EXPECT_EQ(
      size(this->localIndexer.getRegionAll()),
      size(this->localIndexer.getRegionNobndry() + this->localIndexer.getRegionBndry()));
}

TYPED_TEST(IndexerTest, TestGetRegionNobndry) {
  Region<typename TypeParam::ind_type> rgn;
  rgn = this->globalSquareIndexer.getRegionNobndry();
  EXPECT_EQ(rgn.asUnique().size(), this->nx * this->ny * this->nz);
  BOUT_FOR(i, rgn) {
    EXPECT_GE(i.x(), this->globalSquareIndexer.getMesh()->xstart);
    EXPECT_LE(i.x(), this->globalSquareIndexer.getMesh()->xend);
    if (!std::is_same<TypeParam, FieldPerp>::value) {
      EXPECT_GE(i.y(), this->globalSquareIndexer.getMesh()->ystart);
      EXPECT_LE(i.y(), this->globalSquareIndexer.getMesh()->yend);
    }
  }
  rgn = this->globalStarIndexer.getRegionNobndry();
  EXPECT_EQ(rgn.asUnique().size(), this->nx * this->ny * this->nz);
  BOUT_FOR(i, rgn) {
    EXPECT_GE(i.x(), this->globalStarIndexer.getMesh()->xstart);
    EXPECT_LE(i.x(), this->globalStarIndexer.getMesh()->xend);
    if (!std::is_same<TypeParam, FieldPerp>::value) {
      EXPECT_GE(i.y(), this->globalStarIndexer.getMesh()->ystart);
      EXPECT_LE(i.y(), this->globalStarIndexer.getMesh()->yend);
    }
  }
  rgn = this->globalDefaultIndexer.getRegionNobndry();
  EXPECT_EQ(rgn.asUnique().size(), this->nx * this->ny * this->nz);
  BOUT_FOR(i, rgn) {
    EXPECT_GE(i.x(), this->globalDefaultIndexer.getMesh()->xstart);
    EXPECT_LE(i.x(), this->globalDefaultIndexer.getMesh()->xend);
    if (!std::is_same<TypeParam, FieldPerp>::value) {
      EXPECT_GE(i.y(), this->globalDefaultIndexer.getMesh()->ystart);
      EXPECT_LE(i.y(), this->globalDefaultIndexer.getMesh()->yend);
    }
  }
  rgn = this->localIndexer.getRegionNobndry();
  EXPECT_EQ(rgn.asUnique().size(), 0);
  BOUT_FOR(i, rgn) {
    EXPECT_GE(i.x(), this->localIndexer.getMesh()->xstart);
    EXPECT_LE(i.x(), this->localIndexer.getMesh()->xend);
    if (!std::is_same<TypeParam, FieldPerp>::value) {
      EXPECT_GE(i.y(), this->localIndexer.getMesh()->ystart);
      EXPECT_LE(i.y(), this->localIndexer.getMesh()->yend);
    }
  }
}

#if 0 // fails compilation under nvcc -- need to reconcile later
TYPED_TEST(IndexerTest, TestGetRegionBndry) {
  Region<typename TypeParam::ind_type> bounds;
  for (auto& indexer : {this->globalSquareIndexer, this->globalStarIndexer,
                        this->globalDefaultIndexer, this->localIndexer}) {
    bounds = indexer.getRegionLowerY() + indexer.getRegionUpperY()
             + indexer.getRegionInnerX() + indexer.getRegionOuterX();
    bounds.sort();
    EXPECT_EQ(indexer.getRegionBndry().getIndices(), bounds.getIndices());
  }
}
#endif

TYPED_TEST(IndexerTest, TestGetRegionLowerY) {
  Region<typename TypeParam::ind_type> rgn;
  if (std::is_same<TypeParam, FieldPerp>::value) {
    rgn = this->globalSquareIndexer.getRegionLowerY();
    EXPECT_EQ(rgn.asUnique().size(), 0);
    rgn = this->globalStarIndexer.getRegionLowerY();
    EXPECT_EQ(rgn.size(), 0);
  } else {
    rgn = this->globalSquareIndexer.getRegionLowerY();
    EXPECT_EQ(rgn.asUnique().size(), (this->nx + 2 * this->guardx) * this->nz);
    BOUT_FOR(i, rgn) { EXPECT_LT(i.y(), this->globalSquareIndexer.getMesh()->ystart); }
    rgn = this->globalStarIndexer.getRegionLowerY();
    EXPECT_EQ(rgn.asUnique().size(), this->nx * this->nz);
    BOUT_FOR(i, rgn) { EXPECT_LT(i.y(), this->globalStarIndexer.getMesh()->ystart); }
  }
  rgn = this->globalDefaultIndexer.getRegionLowerY();
  EXPECT_EQ(rgn.asUnique().size(), 0);
  rgn = this->localIndexer.getRegionLowerY();
  EXPECT_EQ(rgn.size(), 0);
}

TYPED_TEST(IndexerTest, TestGetRegionUpperY) {
  Region<typename TypeParam::ind_type> rgn;
  if (std::is_same<TypeParam, FieldPerp>::value) {
    rgn = this->globalSquareIndexer.getRegionUpperY();
    EXPECT_EQ(rgn.asUnique().size(), 0);
    rgn = this->globalStarIndexer.getRegionUpperY();
    EXPECT_EQ(rgn.size(), 0);
  } else {
    rgn = this->globalSquareIndexer.getRegionUpperY();
    EXPECT_EQ(rgn.asUnique().size(), (this->nx + 2 * this->guardx) * this->nz);
    BOUT_FOR(i, rgn) { EXPECT_GT(i.y(), this->globalSquareIndexer.getMesh()->yend); }
    rgn = this->globalStarIndexer.getRegionUpperY();
    EXPECT_EQ(rgn.asUnique().size(), this->nx * this->nz);
    BOUT_FOR(i, rgn) { EXPECT_GT(i.y(), this->globalStarIndexer.getMesh()->yend); }
  }
  rgn = this->globalDefaultIndexer.getRegionUpperY();
  EXPECT_EQ(rgn.asUnique().size(), 0);
  rgn = this->localIndexer.getRegionUpperY();
  EXPECT_EQ(rgn.asUnique().size(), 0);
}

TYPED_TEST(IndexerTest, TestGetRegionInnerX) {
  Region<typename TypeParam::ind_type> rgn;
  rgn = this->globalSquareIndexer.getRegionInnerX();
  EXPECT_EQ(rgn.asUnique().size(), this->ny * this->nz);
  BOUT_FOR(i, rgn) { EXPECT_LT(i.x(), this->globalSquareIndexer.getMesh()->xstart); }
  rgn = this->globalStarIndexer.getRegionInnerX();
  EXPECT_EQ(rgn.asUnique().size(), this->ny * this->nz);
  BOUT_FOR(i, rgn) { EXPECT_LT(i.x(), this->globalStarIndexer.getMesh()->xstart); }
  rgn = this->globalDefaultIndexer.getRegionInnerX();
  EXPECT_EQ(rgn.asUnique().size(), 0);
  rgn = this->localIndexer.getRegionInnerX();
  EXPECT_EQ(rgn.asUnique().size(), 0);
}

TYPED_TEST(IndexerTest, TestGetRegionOuterX) {
  Region<typename TypeParam::ind_type> rgn;
  rgn = this->globalSquareIndexer.getRegionOuterX();
  EXPECT_EQ(rgn.asUnique().size(), this->ny * this->nz);
  BOUT_FOR(i, rgn) { EXPECT_GT(i.x(), this->globalSquareIndexer.getMesh()->xend); }
  rgn = this->globalStarIndexer.getRegionOuterX();
  EXPECT_EQ(rgn.asUnique().size(), this->ny * this->nz);
  BOUT_FOR(i, rgn) { EXPECT_GT(i.x(), this->globalStarIndexer.getMesh()->xend); }
  rgn = this->globalDefaultIndexer.getRegionOuterX();
  EXPECT_EQ(rgn.asUnique().size(), 0);
  rgn = this->localIndexer.getRegionOuterX();
  EXPECT_EQ(rgn.asUnique().size(), 0);
}

TYPED_TEST(IndexerTest, TestSparsityPatternAvailable) {
  EXPECT_TRUE(this->globalSquareIndexer.sparsityPatternAvailable());
  EXPECT_TRUE(this->globalStarIndexer.sparsityPatternAvailable());
  EXPECT_FALSE(this->globalDefaultIndexer.sparsityPatternAvailable());
  EXPECT_FALSE(this->localIndexer.sparsityPatternAvailable());
}

TYPED_TEST(IndexerTest, TestGetNumDiagonal) {
  const int squareStencilInteriorSize = std::is_same<TypeParam, Field3D>::value ? 19 : 9,
            starStencilInteriorSize = std::is_same<TypeParam, Field3D>::value ? 7 : 5,
            boundsSize = 1;
  int numBoundaryCells = 0;
  for (int i : this->globalSquareIndexer.getNumDiagonal()) {
    if (i != boundsSize) {
      EXPECT_EQ(i, squareStencilInteriorSize);
    } else {
      numBoundaryCells++;
    }
  }
  EXPECT_EQ(numBoundaryCells, this->globalSquareIndexer.getRegionBndry().size());
  numBoundaryCells = 0;
  for (int i : this->globalStarIndexer.getNumDiagonal()) {
    if (i != boundsSize) {
      EXPECT_EQ(i, starStencilInteriorSize);
    } else {
      numBoundaryCells++;
    }
  }
  EXPECT_EQ(numBoundaryCells, this->globalStarIndexer.getRegionBndry().size());
}

TYPED_TEST(IndexerTest, TestGetNumOffDiagonal) {
  for (int i : this->globalSquareIndexer.getNumOffDiagonal()) {
    EXPECT_EQ(i, 0);
  }
  for (int i : this->globalStarIndexer.getNumOffDiagonal()) {
    EXPECT_EQ(i, 0);
  }
}

TYPED_TEST(IndexerTest, TestSize) {
  EXPECT_EQ(this->globalSquareIndexer.size(),
            (this->nx + 2 * this->guardx) * (this->ny + 2 * this->guardy) * this->nz);
  EXPECT_EQ(this->globalStarIndexer.size(), this->nx * this->ny * this->nz
                                                + 2 * this->nz * this->ny * this->guardx
                                                + 2 * this->nz * this->nx * this->guardy);
  EXPECT_EQ(this->globalDefaultIndexer.size(), this->nx * this->ny * this->nz);
  EXPECT_EQ(this->localIndexer.size(), 0);
}

template <class T>
class FakeParallelIndexer : public GlobalIndexer<T> {
public:
  FakeParallelIndexer(Mesh* localmesh, OperatorStencil<typename T::ind_type> stencil)
      : GlobalIndexer<T>(localmesh, stencil, false) {}
  void initialiseTest() { registerFieldForTest(this->getIndices()); }

private:
  void registerFieldForTest(T& f) override {
    auto* meshPtr = static_cast<FakeParallelMesh*>(this->getMesh());
    if (meshPtr) {
      meshPtr->registerField(f, 0);
    }
  }
};

// Need to parameterise this over different field types, I guess

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
    localmesh = &meshes[pe_yind + pe_xind * nype];

    bout::globals::mpi = &localmesh->getMpi();
  }

  virtual ~ParallelIndexerTest() {
    bout::globals::mesh = nullptr;
    bout::globals::mpi = nullptr;
  }

  static const int nx;
  static const int ny;
  static const int nz;

  int nxpe, nype;
  int pe_xind, pe_yind;
  int xstart, xend, ystart, yend;
  std::vector<FakeParallelMesh> meshes;
  Mesh* localmesh;

  template <class T>
  std::vector<IndexerPtr<T>> makeIndexers() {
    using ind_type = typename T::ind_type;
    std::vector<IndexerPtr<T>> result;
    OperatorStencil<ind_type> stencil = squareStencil<ind_type>(bout::globals::mesh);
    for (int i = 0; i < nxpe * nype; i++) {
      auto ind = std::make_shared<FakeParallelIndexer<T>>(&meshes[i], stencil);
      ind->initialiseTest();
      result.push_back(ind);
    }
    result[pe_yind + pe_xind * nype]->initialise();
    return result;
  }

  template <class T>
  IndexerPtr<T> getIndexer(std::vector<IndexerPtr<T>> indexers, int xind, int yind) {
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

const int ParallelIndexerTest::nx = 7;
const int ParallelIndexerTest::ny = 5;
const int ParallelIndexerTest::nz = 3;

INSTANTIATE_TEST_SUITE_P(1by3, ParallelIndexerTest,
                         testing::ValuesIn(makeParallelTestCases(1, 3)));
INSTANTIATE_TEST_SUITE_P(3by1, ParallelIndexerTest,
                         testing::ValuesIn(makeParallelTestCases(3, 1)));
INSTANTIATE_TEST_SUITE_P(3by3, ParallelIndexerTest,
                         testing::ValuesIn(makeParallelTestCases(3, 3)));

TEST_P(ParallelIndexerTest, TestConvertIndex3D) {
  std::vector<IndexerPtr<Field3D>> indexers = this->makeIndexers<Field3D>();
  IndexerPtr<Field3D> index = this->getIndexer(indexers, this->pe_xind, this->pe_yind);

  // Test x boundaries
  BOUT_FOR(i, this->localmesh->getRegion3D("RGN_XGUARDS")) {
    if (i.x() < this->xstart && i.y() >= this->ystart && i.y() <= this->yend) {
      int global = index->getGlobal(i);
      if (this->pe_xind > 0) {
        int otherGlobal =
            this->getIndexer<Field3D>(indexers, this->pe_xind - 1, this->pe_yind)
                ->getGlobal(i.xp(this->xend - this->xstart + 1));
        EXPECT_EQ(global, otherGlobal);
      } else {
        EXPECT_NE(global, -1);
      }
    }
    if (i.x() >= this->xend && i.y() >= this->ystart && i.y() <= this->yend) {
      int global = index->getGlobal(i);
      if (this->pe_xind < this->nxpe - 1) {
        int otherGlobal =
            this->getIndexer<Field3D>(indexers, this->pe_xind + 1, this->pe_yind)
                ->getGlobal(i.xm(this->xend - this->xstart + 1));
        EXPECT_EQ(global, otherGlobal);
      } else {
        EXPECT_NE(global, -1);
      }
    }
  }

  // Test y boundaries
  BOUT_FOR(i, this->localmesh->getRegion3D("RGN_YGUARDS")) {
    if (i.y() < this->ystart) {
      
      if (i.x() == this->xstart) {
	const int global = index->getGlobal(i.xm()),
	  otherXind = (this->pe_xind == 0) ? this->pe_xind : this->pe_xind - 1,
	  otherYind = this->pe_yind - 1;
	const Ind3D otherInd = i.offset((this->pe_xind == 0) ? -1 :
					this->xend - this->xstart,
					this->yend - this->ystart + 1, 0);
	const int otherGlobal =
	  this->getIndexer<Field3D>(indexers, otherXind, otherYind)->getGlobal(otherInd);
	EXPECT_NE(global, -1);
	EXPECT_EQ(global, otherGlobal);
      }
      
      int global = index->getGlobal(i);
      int otherGlobal =
          this->getIndexer<Field3D>(indexers, this->pe_xind, this->pe_yind - 1)
              ->getGlobal(i.yp(this->yend - this->ystart + 1));
      EXPECT_EQ(global, otherGlobal);

      if (i.x() == this->xend) {
	const int global = index->getGlobal(i.xp()),
	  otherXind = (this->pe_xind == this->nxpe - 1) ? this->pe_xind : this->pe_xind + 1,
	  otherYind = this->pe_yind - 1;
	const Ind3D otherInd = i.offset((this->pe_xind == this->nxpe - 1) ? 1 :
					this->xstart - this->xend,
					this->yend - this->ystart + 1, 0);
	const int otherGlobal =
	  this->getIndexer<Field3D>(indexers, otherXind, otherYind)->getGlobal(otherInd);
	EXPECT_NE(global, -1);
	EXPECT_EQ(global, otherGlobal);
      }
      
    } else if (i.y() >= this->yend && i.x() >= this->xstart && i.x() <= this->xend) {

      if (i.x() == this->xstart) {
	const int global = index->getGlobal(i.xm()),
	  otherXind = (this->pe_xind == 0) ? this->pe_xind : this->pe_xind - 1,
	  otherYind = this->pe_yind + 1;
	const Ind3D otherInd = i.offset((this->pe_xind == 0) ? -1 :
					this->xend - this->xstart,
					this->ystart - this->yend - 1, 0);
	const int otherGlobal =
	  this->getIndexer<Field3D>(indexers, otherXind, otherYind)->getGlobal(otherInd);
	EXPECT_NE(global, -1);
	EXPECT_EQ(global, otherGlobal);
      }

      int global = index->getGlobal(i);
      int otherGlobal =
          this->getIndexer<Field3D>(indexers, this->pe_xind, this->pe_yind + 1)
              ->getGlobal(i.ym(this->yend - this->ystart + 1));
      EXPECT_EQ(global, otherGlobal);

      if (i.x() == this->xend) {
	const int global = index->getGlobal(i.xp()),
	  otherXind = (this->pe_xind == this->nxpe - 1) ? this->pe_xind : this->pe_xind + 1,
	  otherYind = this->pe_yind + 1;
	const Ind3D otherInd = i.offset((this->pe_xind == this->nxpe - 1) ? 1 :
					this->xstart - this->xend,
					this->ystart - this->yend - 1, 0);
	const int otherGlobal =
	  this->getIndexer<Field3D>(indexers, otherXind, otherYind)->getGlobal(otherInd);
	EXPECT_NE(global, -1);
	EXPECT_EQ(global, otherGlobal);
      }

    }
  }
}

TEST_P(ParallelIndexerTest, TestConvertIndex2D) {
  std::vector<IndexerPtr<Field2D>> indexers = this->makeIndexers<Field2D>();
  IndexerPtr<Field2D> index = this->getIndexer(indexers, this->pe_xind, this->pe_yind);

  // Test x boundaries
  BOUT_FOR(i, this->localmesh->getRegion2D("RGN_XGUARDS")) {
    if (i.x() < this->xstart) {
      int global = index->getGlobal(i);
      if (this->pe_xind > 0) {
        int otherGlobal =
            this->getIndexer<Field2D>(indexers, this->pe_xind - 1, this->pe_yind)
                ->getGlobal(i.xp(this->xend - this->xstart + 1));
        EXPECT_EQ(global, otherGlobal);
      } else {
        EXPECT_NE(global, -1);
      }
    }
    if (i.x() >= this->xend) {
      int global = index->getGlobal(i);
      if (this->pe_xind < this->nxpe - 1) {
        int otherGlobal =
            this->getIndexer<Field2D>(indexers, this->pe_xind + 1, this->pe_yind)
                ->getGlobal(i.xm(this->xend - this->xstart + 1));
        EXPECT_EQ(global, otherGlobal);
      } else {
        EXPECT_NE(global, -1);
      }
    }
  }

  // Test y boundaries
  BOUT_FOR(i, this->localmesh->getRegion2D("RGN_YGUARDS")) {
    if (i.y() < this->ystart && i.x() >= this->xstart && i.x() <= this->xend) {
      int global = index->getGlobal(i);
      int otherGlobal =
          this->getIndexer<Field2D>(indexers, this->pe_xind, this->pe_yind - 1)
              ->getGlobal(i.yp(this->yend - this->ystart + 1));
      EXPECT_EQ(global, otherGlobal);

      if (i.x() == this->xend) {
	const int global = index->getGlobal(i.xp()),
	  otherXind = (this->pe_xind == this->nxpe - 1) ? this->pe_xind : this->pe_xind + 1,
	  otherYind = this->pe_yind - 1;
	const Ind2D otherInd = i.offset((this->pe_xind == this->nxpe - 1) ? 1 :
					this->xstart - this->xend,
					this->yend - this->ystart + 1, 0);
	const int otherGlobal =
	  this->getIndexer<Field2D>(indexers, otherXind, otherYind)->getGlobal(otherInd);
	EXPECT_NE(global, -1);
	EXPECT_EQ(global, otherGlobal);
      }
      
    } else if (i.y() >= this->yend) {

      if (i.x() == this->xstart) {
	const int global = index->getGlobal(i.xm()),
	  otherXind = (this->pe_xind == 0) ? this->pe_xind : this->pe_xind - 1,
	  otherYind = this->pe_yind + 1;
	const Ind2D otherInd = i.offset((this->pe_xind == 0) ? -1 :
					this->xend - this->xstart,
					this->ystart - this->yend - 1, 0);
	const int otherGlobal =
	  this->getIndexer<Field2D>(indexers, otherXind, otherYind)->getGlobal(otherInd);
	EXPECT_NE(global, -1);
	EXPECT_EQ(global, otherGlobal);
      }

  int global = index->getGlobal(i);
      int otherGlobal =
          this->getIndexer<Field2D>(indexers, this->pe_xind, this->pe_yind + 1)
              ->getGlobal(i.ym(this->yend - this->ystart + 1));
      EXPECT_EQ(global, otherGlobal);

      if (i.x() == this->xend) {
	const int global = index->getGlobal(i.xp()),
	  otherXind = (this->pe_xind == this->nxpe - 1) ? this->pe_xind : this->pe_xind + 1,
	  otherYind = this->pe_yind + 1;
	const Ind2D otherInd = i.offset((this->pe_xind == this->nxpe - 1) ? 1 :
					this->xstart - this->xend,
					this->ystart - this->yend - 1, 0);
	const int otherGlobal =
	  this->getIndexer<Field2D>(indexers, otherXind, otherYind)->getGlobal(otherInd);
	EXPECT_NE(global, -1);
	EXPECT_EQ(global, otherGlobal);
      }
      
    }
  }
}

TEST_P(ParallelIndexerTest, TestConvertIndexPerp) {
  std::vector<IndexerPtr<FieldPerp>> indexers = this->makeIndexers<FieldPerp>();
  IndexerPtr<FieldPerp> index = this->getIndexer(indexers, this->pe_xind, this->pe_yind);

  // Test x boundaries
  BOUT_FOR(i, this->localmesh->getRegionPerp("RGN_XGUARDS")) {
    if (i.x() < this->xstart) {
      int global = index->getGlobal(i);
      if (this->pe_xind > 0) {
        int otherGlobal =
            this->getIndexer<FieldPerp>(indexers, this->pe_xind - 1, this->pe_yind)
                ->getGlobal(i.xp(this->xend - this->xstart + 1));
        EXPECT_EQ(global, otherGlobal);
      } else {
        EXPECT_NE(global, -1);
      }
    }
    if (i.x() >= this->xend) {
      int global = index->getGlobal(i);
      if (this->pe_xind < this->nxpe - 1) {
        int otherGlobal =
            this->getIndexer<FieldPerp>(indexers, this->pe_xind + 1, this->pe_yind)
                ->getGlobal(i.xm(this->xend - this->xstart + 1));
        EXPECT_EQ(global, otherGlobal);
      } else {
        EXPECT_NE(global, -1);
      }
    }
  }
}

TEST_P(ParallelIndexerTest, TestRegions3D) {
  std::vector<IndexerPtr<Field3D>> indexers = this->makeIndexers<Field3D>();
  IndexerPtr<Field3D> index = this->getIndexer(indexers, this->pe_xind, this->pe_yind);
  Region<Ind3D> rgn;

  rgn = index->getRegionAll();
  EXPECT_EQ(size(index->getRegionAll()),
            size(index->getRegionNobndry() + index->getRegionBndry()));

  rgn = index->getRegionNobndry();
  EXPECT_EQ(rgn.size(), (this->nx - 2) * this->ny * this->nz);
  BOUT_FOR(i, rgn) {
    EXPECT_GE(i.x(), index->getMesh()->xstart);
    EXPECT_LE(i.x(), index->getMesh()->xend);
    EXPECT_GE(i.y(), index->getMesh()->ystart);
    EXPECT_LE(i.y(), index->getMesh()->yend);
  }

  Region<Ind3D> bounds = index->getRegionLowerY() + index->getRegionUpperY()
                         + index->getRegionInnerX() + index->getRegionOuterX();
  bounds.sort();
  EXPECT_EQ(index->getRegionBndry().getIndices(), bounds.getIndices());

  rgn = index->getRegionLowerY();
  EXPECT_EQ(rgn.size(), 0);
  rgn = index->getRegionUpperY();
  EXPECT_EQ(rgn.size(), 0);

  rgn = index->getRegionInnerX();
  if (this->pe_xind == 0) {
    EXPECT_EQ(rgn.size(), this->ny * this->nz);
    BOUT_FOR(i, rgn) {
      EXPECT_EQ(i.x(), index->getMesh()->xstart - 1);
      EXPECT_GE(i.y(), index->getMesh()->ystart);
      EXPECT_LE(i.y(), index->getMesh()->yend);
    }
  } else {
    EXPECT_EQ(rgn.size(), 0);
  }
  rgn = index->getRegionOuterX();
  if (this->pe_xind == this->nxpe - 1) {
    EXPECT_EQ(rgn.size(), this->ny * this->nz);
    BOUT_FOR(i, rgn) {
      EXPECT_EQ(i.x(), index->getMesh()->xend + 1);
      EXPECT_GE(i.y(), index->getMesh()->ystart);
      EXPECT_LE(i.y(), index->getMesh()->yend);
    }
  } else {
    EXPECT_EQ(rgn.size(), 0);
  }
}

TEST_P(ParallelIndexerTest, TestRegions2D) {
  std::vector<IndexerPtr<Field2D>> indexers = this->makeIndexers<Field2D>();
  IndexerPtr<Field2D> index = this->getIndexer(indexers, this->pe_xind, this->pe_yind);
  Region<Ind2D> rgn;

  rgn = index->getRegionAll();
  EXPECT_EQ(size(index->getRegionAll()),
            size(index->getRegionNobndry() + index->getRegionBndry()));

  rgn = index->getRegionNobndry();
  EXPECT_EQ(rgn.size(), (this->nx - 2) * this->ny);
  BOUT_FOR(i, rgn) {
    EXPECT_GE(i.x(), index->getMesh()->xstart);
    EXPECT_LE(i.x(), index->getMesh()->xend);
    EXPECT_GE(i.y(), index->getMesh()->ystart);
    EXPECT_LE(i.y(), index->getMesh()->yend);
  }

  Region<Ind2D> bounds = index->getRegionLowerY() + index->getRegionUpperY()
                         + index->getRegionInnerX() + index->getRegionOuterX();
  bounds.sort();
  EXPECT_EQ(index->getRegionBndry().getIndices(), bounds.getIndices());

  rgn = index->getRegionLowerY();
  EXPECT_EQ(rgn.size(), 0);
  rgn = index->getRegionUpperY();
  EXPECT_EQ(rgn.size(), 0);

  rgn = index->getRegionInnerX();
  if (this->pe_xind == 0) {
    EXPECT_EQ(rgn.size(), this->ny);
    BOUT_FOR(i, rgn) {
      EXPECT_EQ(i.x(), index->getMesh()->xstart - 1);
      EXPECT_GE(i.y(), index->getMesh()->ystart);
      EXPECT_LE(i.y(), index->getMesh()->yend);
    }
  } else {
    EXPECT_EQ(rgn.size(), 0);
  }
  rgn = index->getRegionOuterX();
  if (this->pe_xind == this->nxpe - 1) {
    EXPECT_EQ(rgn.size(), this->ny);
    BOUT_FOR(i, rgn) {
      EXPECT_EQ(i.x(), index->getMesh()->xend + 1);
      EXPECT_GE(i.y(), index->getMesh()->ystart);
      EXPECT_LE(i.y(), index->getMesh()->yend);
    }
  } else {
    EXPECT_EQ(rgn.size(), 0);
  }
}

TEST_P(ParallelIndexerTest, TestRegionsPerp) {
  std::vector<IndexerPtr<FieldPerp>> indexers = this->makeIndexers<FieldPerp>();
  IndexerPtr<FieldPerp> index = this->getIndexer(indexers, this->pe_xind, this->pe_yind);
  Region<IndPerp> rgn;

  rgn = index->getRegionAll();
  EXPECT_EQ(size(index->getRegionAll()),
            size(index->getRegionNobndry() + index->getRegionBndry()));

  rgn = index->getRegionNobndry();
  EXPECT_EQ(rgn.size(), (this->nx - 2) * this->nz);
  BOUT_FOR(i, rgn) {
    EXPECT_GE(i.x(), index->getMesh()->xstart);
    EXPECT_LE(i.x(), index->getMesh()->xend);
  }

  Region<IndPerp> bounds = index->getRegionLowerY() + index->getRegionUpperY()
                           + index->getRegionInnerX() + index->getRegionOuterX();
  bounds.sort();
  EXPECT_EQ(index->getRegionBndry().getIndices(), bounds.getIndices());

  rgn = index->getRegionLowerY();
  EXPECT_EQ(rgn.size(), 0);
  rgn = index->getRegionUpperY();
  EXPECT_EQ(rgn.size(), 0);

  rgn = index->getRegionInnerX();
  if (this->pe_xind == 0) {
    EXPECT_EQ(rgn.size(), this->nz);
    BOUT_FOR(i, rgn) { EXPECT_EQ(i.x(), index->getMesh()->xstart - 1); }
  } else {
    EXPECT_EQ(rgn.size(), 0);
  }
  rgn = index->getRegionOuterX();
  if (this->pe_xind == this->nxpe - 1) {
    EXPECT_EQ(rgn.size(), this->nz);
    BOUT_FOR(i, rgn) { EXPECT_EQ(i.x(), index->getMesh()->xend + 1); }
  } else {
    EXPECT_EQ(rgn.size(), 0);
  }
}

#endif // BOUT_HAS_PETSC
