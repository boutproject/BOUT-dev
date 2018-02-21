#include "gtest/gtest.h"

#include "bout/constants.hxx"
#include "bout/mesh.hxx"
#include "bout/region.hxx"
#include "boutexception.hxx"
#include "output.hxx"
#include "test_extras.hxx"
#include "unused.hxx"

#include <cmath>
#include <set>
#include <vector>

/// Global mesh
extern Mesh *mesh;

/// Test fixture to make sure the global mesh is our fake one
class RegionTest : public ::testing::Test {
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

const int RegionTest::nx = 3;
const int RegionTest::ny = 5;
const int RegionTest::nz = 7;

TEST_F(RegionTest, maxBlockSize) { EXPECT_TRUE(MAXREGIONBLOCKSIZE > 0); }

TEST_F(RegionTest, regionFromRange) {
  Region<Ind3D> region(0, mesh->LocalNx - 1, 0, mesh->LocalNy - 1, 0, mesh->LocalNz - 1,
                       mesh->LocalNy, mesh->LocalNz);

  int nmesh = RegionTest::nx * RegionTest::ny * RegionTest::nz;
  auto regionIndices = region.getIndices();

  // Check index vector is expected length
  EXPECT_EQ(regionIndices.size(), nmesh);

  for (int i = 0; i < nmesh; i++) {
    EXPECT_EQ(regionIndices[i].ind, i);
  }
}

TEST_F(RegionTest, regionFromIndices) {
  std::vector<std::pair<int, int>> blocksIn = {
      {0, 3}, {5, 7}, {9, 10}, {12, 12}, {14, 20}};
  std::vector<Ind3D> indicesIn;
  int maxContiguousSizeUsed = 0;
  for (auto &block : blocksIn) {
    int currBlockSize = 1 + block.second - block.first;
    maxContiguousSizeUsed =
        currBlockSize > maxContiguousSizeUsed ? currBlockSize : maxContiguousSizeUsed;
    for (int i = block.first; i <= block.second; i++) {
      indicesIn.push_back(i);
    }
  }

  // If our region block size is too small then this test won't work
  if (maxContiguousSizeUsed > MAXREGIONBLOCKSIZE) {
    return;
  }

  Region<Ind3D> region(indicesIn);

  auto regionIndices = region.getIndices();

  // Check index vector is expected length
  EXPECT_EQ(regionIndices.size(), indicesIn.size());

  for (unsigned int i = 0; i < indicesIn.size(); i++) {
    EXPECT_EQ(regionIndices[i], indicesIn[i]);
  }

  auto regionBlocks = region.getBlocks();

  // Check expected number of blocks
  EXPECT_EQ(regionBlocks.size(), blocksIn.size());

  for (unsigned int i = 0; i < blocksIn.size(); i++) {
    EXPECT_EQ(regionBlocks[i].first.ind, blocksIn[i].first);
    // The actual block second is exclusive, blocksIn.second is inclusive
    EXPECT_EQ(regionBlocks[i].second.ind - 1, blocksIn[i].second);
  }
}

TEST_F(RegionTest, numberOfBlocks) {
  Region<Ind3D> region(0, mesh->LocalNx - 1, 0, mesh->LocalNy - 1, 0, mesh->LocalNz - 1,
                       mesh->LocalNy, mesh->LocalNz);

  auto blocks = region.getBlocks();
  int nmesh = RegionTest::nx * RegionTest::ny * RegionTest::nz;
  int nblocks = blocks.size();

  // Calculate expected number of blocks, assuming region is entirely
  // contiguous, as is the case here.
  int nexpected = nmesh / MAXREGIONBLOCKSIZE; // Integer division required
  // Above is only correct when MAXREGIONBLOCKSIZE is a factor of nmesh,
  // so check if there needs to be an extra block for remainder
  if ((1.0 * nmesh) / MAXREGIONBLOCKSIZE > 0) {
    nexpected++;
  }

  EXPECT_EQ(nexpected, nblocks);
}

TEST_F(RegionTest, contiguousBlockSize) {
  Region<Ind3D> region(0, mesh->LocalNx - 1, 0, mesh->LocalNy - 1, 0, mesh->LocalNz - 1,
                       mesh->LocalNy, mesh->LocalNz);

  const int nmesh = RegionTest::nx * RegionTest::ny * RegionTest::nz;

  auto blocks = region.getBlocks();
  const int nblocks = blocks.size();
  auto firstBlock = blocks[0];
  auto lastBlock = blocks[nblocks - 1];

  // Calculate expected size of blocks, assuming region is entirely
  // contiguous, as is the case here.
  const int expectedFirstBlockSize =
      nmesh >= MAXREGIONBLOCKSIZE ? MAXREGIONBLOCKSIZE : nmesh;
  EXPECT_EQ(firstBlock.second.ind - firstBlock.first.ind, expectedFirstBlockSize);
  const int expectedLastBlockSize = nmesh % MAXREGIONBLOCKSIZE;
  if (expectedLastBlockSize != 0) {
    EXPECT_EQ(lastBlock.second.ind - lastBlock.first.ind, expectedLastBlockSize);
  } else {
    // If no remainder then expect same block size as in first block
    EXPECT_EQ(lastBlock.second.ind - lastBlock.first.ind, expectedFirstBlockSize);
  }
}

TEST_F(RegionTest, defaultRegions) {
  const int nmesh = RegionTest::nx * RegionTest::ny * RegionTest::nz;
  EXPECT_EQ(mesh->getRegion("RGN_ALL").getIndices().size(), nmesh);

  const int nnobndry =
      RegionTest::nz * (mesh->yend - mesh->ystart + 1) * (mesh->xend - mesh->xstart + 1);
  EXPECT_EQ(mesh->getRegion("RGN_NOBNDRY").getIndices().size(), nnobndry);

  const int nnox = RegionTest::nz * RegionTest::ny * (mesh->xend - mesh->xstart + 1);
  EXPECT_EQ(mesh->getRegion("RGN_NOX").getIndices().size(), nnox);

  const int nnoy = RegionTest::nz * RegionTest::nx * (mesh->yend - mesh->ystart + 1);
  EXPECT_EQ(mesh->getRegion("RGN_NOY").getIndices().size(), nnoy);
}

TEST_F(RegionTest, regionLoopAll) {
  auto region = mesh->getRegion3D("RGN_ALL");

  Field3D a = 0.0;
  BLOCK_REGION_LOOP(region, i, a[i] = 1.0;);

  for (const auto &i : a.region(RGN_ALL)) {
    EXPECT_EQ(a[i], 1.0);
  }
}

TEST_F(RegionTest, regionLoopNoBndry) {
  auto region = mesh->getRegion3D("RGN_NOBNDRY");

  Field3D a = 0.0;
  BLOCK_REGION_LOOP(region, i, a[i] = 1.0;);

  const int nmesh = RegionTest::nx * RegionTest::ny * RegionTest::nz;
  const int ninner =
      (mesh->LocalNz * (1 + mesh->xend - mesh->xstart) * (1 + mesh->yend - mesh->ystart));
  int numExpectNotMatching = nmesh - ninner;

  int numNotMatching = 0;
  int numMatching = 0;
  for (const auto &i : a.region(RGN_ALL)) {
    if (a[i] != 1.0) {
      numNotMatching++;
    } else {
      numMatching++;
    }
  }
  EXPECT_EQ(numNotMatching, numExpectNotMatching);
  EXPECT_EQ(numMatching, ninner);
}

TEST_F(RegionTest, regionLoopAllSerial) {
  auto region = mesh->getRegion3D("RGN_ALL");

  Field3D a = 0.0;
  BLOCK_REGION_LOOP_SERIAL(region, i, a[i] = 1.0;);

  for (const auto &i : a.region(RGN_ALL)) {
    EXPECT_EQ(a[i], 1.0);
  }
}

TEST_F(RegionTest, regionLoopNoBndrySerial) {
  auto region = mesh->getRegion3D("RGN_NOBNDRY");

  Field3D a = 0.0;
  BLOCK_REGION_LOOP_SERIAL(region, i, a[i] = 1.0;);

  const int nmesh = RegionTest::nx * RegionTest::ny * RegionTest::nz;
  const int ninner =
      (mesh->LocalNz * (1 + mesh->xend - mesh->xstart) * (1 + mesh->yend - mesh->ystart));
  int numExpectNotMatching = nmesh - ninner;

  int numNotMatching = 0;
  int numMatching = 0;
  for (const auto &i : a.region(RGN_ALL)) {
    if (a[i] != 1.0) {
      numNotMatching++;
    } else {
      numMatching++;
    }
  }
  EXPECT_EQ(numNotMatching, numExpectNotMatching);
  EXPECT_EQ(numMatching, ninner);
}
