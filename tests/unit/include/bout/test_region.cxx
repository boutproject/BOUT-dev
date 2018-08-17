#include "gtest/gtest.h"

#include "bout/constants.hxx"
#include "bout/mesh.hxx"
#include "bout/region.hxx"
#include "boutexception.hxx"
#include "output.hxx"
#include "test_extras.hxx"
#include "unused.hxx"

#include <algorithm>
#include <list>
#include <vector>
#include <sstream>

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
    output_info.disable();
    mesh->createDefaultRegions();
    output_info.enable();
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
    EXPECT_EQ(regionIndices[i], i);
  }

  // Single point
  Region<Ind3D> region2(0, 0, 0, 0, 0, 0, 1, 1);
  EXPECT_EQ(region2.getIndices().size(), 1);

  // Invalid range results in empty region
  { Region<Ind3D> region3(0, -1, 0, 0, 0, 0, 1, 1);
    EXPECT_EQ(region3.size(), 0);}
  { Region<Ind3D> region3(0, 0, 1, 0, 0, 0, 1, 1);
    EXPECT_EQ(region3.size(), 0);}
  { Region<Ind3D> region3(0, 0, 0, 0, 20, 10, 1, 1);
    EXPECT_EQ(region3.size(), 0);}

  // Invalid size throws if CHECK >= 1
#if CHECK >= 1
  EXPECT_THROW(Region<Ind3D> region3(0, 0, 0, 0, 0, 0, 0, 1), BoutException);
  EXPECT_THROW(Region<Ind3D> region3(0, 0, 0, 0, 0, 0, 1, 0), BoutException);
#else
  EXPECT_NO_THROW(Region<Ind3D> region3(0, 0, 0, 0, 0, 0, 0, 1));
  EXPECT_NO_THROW(Region<Ind3D> region3(0, 0, 0, 0, 0, 0, 1, 0));
#endif
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
    EXPECT_EQ(regionBlocks[i].first, blocksIn[i].first);
    // The actual block second is exclusive, blocksIn.second is inclusive
    EXPECT_EQ(regionBlocks[i].second - 1, blocksIn[i].second);
  }
}

TEST_F(RegionTest, regionFromBlocks) {
  Region<Ind3D> region(0, mesh->LocalNx - 1, 0, mesh->LocalNy - 1, 0, mesh->LocalNz - 1,
                       mesh->LocalNy, mesh->LocalNz);

  const int nmesh = RegionTest::nx * RegionTest::ny * RegionTest::nz;

  auto blocks = region.getBlocks();

  Region<Ind3D> region2(blocks);
  auto regionIndices = region2.getIndices();

  EXPECT_EQ(regionIndices.size(), nmesh);
  for (int i = 0; i < nmesh; i++) {
    EXPECT_EQ(regionIndices[i], i);
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
  EXPECT_EQ(firstBlock.second - firstBlock.first, expectedFirstBlockSize);
  const int expectedLastBlockSize = nmesh % MAXREGIONBLOCKSIZE;
  if (expectedLastBlockSize != 0) {
    EXPECT_EQ(lastBlock.second - lastBlock.first, expectedLastBlockSize);
  } else {
    // If no remainder then expect same block size as in first block
    EXPECT_EQ(lastBlock.second - lastBlock.first, expectedFirstBlockSize);
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
  BLOCK_REGION_LOOP(i, region) { a[i] = 1.0; }

  for (const auto &i : a.region(RGN_ALL)) {
    EXPECT_EQ(a[i], 1.0);
  }
}

TEST_F(RegionTest, regionLoopNoBndry) {
  auto region = mesh->getRegion3D("RGN_NOBNDRY");

  Field3D a = 0.0;
  BLOCK_REGION_LOOP(i, region) { a[i] = 1.0; }

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
  BLOCK_REGION_LOOP_SERIAL(i, region) { a[i] = 1.0; }

  for (const auto &i : a.region(RGN_ALL)) {
    EXPECT_EQ(a[i], 1.0);
  }
}

TEST_F(RegionTest, regionLoopNoBndrySerial) {
  auto region = mesh->getRegion3D("RGN_NOBNDRY");

  Field3D a = 0.0;
  BLOCK_REGION_LOOP_SERIAL(i, region) { a[i] = 1.0; }

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

TEST_F(RegionTest, regionAsSorted) {
  // Contiguous blocks to insert
  std::vector<std::pair<int, int>> blocksIn = {
      {0, 3}, {5, 7}, {9, 10}, {12, 12}, {14, 20}};
  Region<Ind3D>::RegionIndices indicesIn;

  for (auto &block : blocksIn) {
    for (int i = block.first; i <= block.second; i++) {
      indicesIn.push_back(i);
    }
  }

  // This is the sorted region and indices
  Region<Ind3D> regionSortedIn(indicesIn);
  Region<Ind3D>::RegionIndices regionIndicesSortedIn = regionSortedIn.getIndices();

  // Now shuffle the order and create a new region
  std::random_shuffle(std::begin(indicesIn), std::end(indicesIn));
  Region<Ind3D> regionShuffledIn(indicesIn);
  Region<Ind3D>::RegionIndices regionIndicesShuffledIn = regionShuffledIn.getIndices();
  // Should we check the shuffle has actually changed the index order?
  // If it hasn't the test isn't really testing sorted

  // Check index vectors haven't changed length
  EXPECT_EQ(regionIndicesShuffledIn.size(), regionIndicesSortedIn.size());

  // Now get a sorted verion of the region
  Region<Ind3D> regionAsSorted = regionShuffledIn.asSorted();
  Region<Ind3D>::RegionIndices regionIndicesAsSorted = regionAsSorted.getIndices();

  for (unsigned int i = 0; i < regionIndicesSortedIn.size(); i++) {
    EXPECT_EQ(regionIndicesAsSorted[i], regionIndicesSortedIn[i]);
  }
}

TEST_F(RegionTest, regionAsUnique) {
  // Values to insert
  std::vector<int> rawIndicesIn1 = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9,
                                    8, 7, 6, 5, 4, 3, 2, 1, 0};

  Region<Ind3D>::RegionIndices indicesIn1;
  for (auto i : rawIndicesIn1) {
    indicesIn1.push_back(i);
  }

  // This is the sorted region and indices
  Region<Ind3D> regionIn1(indicesIn1);
  Region<Ind3D>::RegionIndices regionIndicesIn1 = regionIn1.getIndices();

  // Now get a unique version of the region
  Region<Ind3D> regionUnique1 = regionIn1.asUnique();
  Region<Ind3D>::RegionIndices regionIndicesUnique1 = regionUnique1.getIndices();

  EXPECT_EQ(regionIndicesUnique1.size(), 10);
  for (unsigned int i = 0; i < regionIndicesUnique1.size(); i++) {
    EXPECT_EQ(regionIndicesUnique1[i], i);
  }

  std::vector<int> rawIndicesIn2 = {0, 0, 0, 2, 2, 2, 5, 4, 4, 5, 10, 12, 11, 9};

  Region<Ind3D>::RegionIndices indicesIn2;
  for (auto i : rawIndicesIn2) {
    indicesIn2.push_back(i);
  }

  // This is the sorted region and indices
  Region<Ind3D> regionIn2(indicesIn2);
  Region<Ind3D>::RegionIndices regionIndicesIn2 = regionIn2.getIndices();

  // Now get a unique version of the region
  Region<Ind3D> regionUnique2 = regionIn2.asUnique();
  Region<Ind3D>::RegionIndices regionIndicesUnique2 = regionUnique2.getIndices();

  EXPECT_EQ(regionIndicesUnique2.size(), 8);

  std::vector<int> rawIndicesIn3 = {9, 8, 7, 6, 5, 4, 3, 2, 1, 0};

  Region<Ind3D>::RegionIndices indicesIn3;
  for (auto i : rawIndicesIn3) {
    indicesIn3.push_back(i);
  }

  // This is the sorted region and indices
  Region<Ind3D> regionIn3(indicesIn3);
  Region<Ind3D>::RegionIndices regionIndicesIn3 = regionIn3.getIndices();

  // Now get a unique version of the region
  Region<Ind3D> regionUnique3 = regionIn3.asUnique();
  Region<Ind3D>::RegionIndices regionIndicesUnique3 = regionUnique3.getIndices();

  EXPECT_EQ(regionIndicesUnique3.size(), 10);
  for (unsigned int i = 0; i < regionIndicesUnique3.size(); i++) {
    EXPECT_EQ(regionIndicesUnique3[i], i);
  }
}

TEST_F(RegionTest, regionSetIndices) {
  // Values to insert
  std::vector<int> rawIndicesIn1 = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};

  std::vector<int> rawIndicesIn2 = {10, 11, 12};

  Region<Ind3D>::RegionIndices indicesIn1;
  for (auto i : rawIndicesIn1) {
    indicesIn1.push_back(i);
  }

  Region<Ind3D>::RegionIndices indicesIn2;
  for (auto i : rawIndicesIn2) {
    indicesIn2.push_back(i);
  }

  // This is the sorted region and indices
  Region<Ind3D> regionIn1(indicesIn1);
  Region<Ind3D>::RegionIndices regionIndicesIn1 = regionIn1.getIndices();
  EXPECT_EQ(regionIndicesIn1.size(), indicesIn1.size());

  for (unsigned int i = 0; i < regionIndicesIn1.size(); i++) {
    EXPECT_EQ(regionIndicesIn1[i], indicesIn1[i]);
  }

  regionIn1.setIndices(indicesIn2);
  Region<Ind3D>::RegionIndices regionIndicesIn2 = regionIn1.getIndices();
  EXPECT_EQ(regionIndicesIn2.size(), indicesIn2.size());

  for (unsigned int i = 0; i < regionIndicesIn2.size(); i++) {
    EXPECT_EQ(regionIndicesIn2[i], indicesIn2[i]);
  }
}

TEST_F(RegionTest, regionSetBlocks) {
  const int nmesh = RegionTest::nx * RegionTest::ny * RegionTest::nz;

  Region<Ind3D> region(0, mesh->LocalNx - 1, 0, mesh->LocalNy - 1, 0, mesh->LocalNz - 1,
                       mesh->LocalNy, mesh->LocalNz);
  auto blocks = region.getBlocks();
  auto indices = region.getIndices();

  EXPECT_EQ(indices.size(), nmesh);

  Region<Ind3D> region2(0, 0, 0, 0, 0, 0, 1, 1);
  auto blocks2 = region2.getBlocks();
  auto indices2 = region2.getIndices();

  EXPECT_EQ(indices2.size(), 1);

  region2.setBlocks(blocks);
  auto indices3 = region2.getIndices();

  EXPECT_EQ(indices3.size(), nmesh);
  for (int i = 0; i < nmesh; i++) {
    EXPECT_EQ(indices3[i], i);
  }
}

TEST_F(RegionTest, regionSortInPlace) {
  // Values to insert
  std::vector<int> rawIndicesBwd = {9, 8, 7, 6, 5, 4, 3, 2, 1, 0};
  std::vector<int> rawIndicesFwd = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};

  Region<Ind3D>::RegionIndices indicesIn1;
  for (auto i : rawIndicesBwd) {
    indicesIn1.push_back(i);
  }

  Region<Ind3D> region(indicesIn1);

  // Check initial order
  auto regionIndices = region.getIndices();
  for (unsigned int i = 0; i < indicesIn1.size(); i++) {
    EXPECT_EQ(regionIndices[i], rawIndicesBwd[i]);
  }

  // Sort in place
  region.sort();

  // Check new order
  regionIndices = region.getIndices();
  for (unsigned int i = 0; i < indicesIn1.size(); i++) {
    EXPECT_EQ(regionIndices[i], rawIndicesFwd[i]);
  }
}

TEST_F(RegionTest, regionFriendSort) {
  // Values to insert
  std::vector<int> rawIndicesBwd = {9, 8, 7, 6, 5, 4, 3, 2, 1, 0};
  std::vector<int> rawIndicesFwd = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};

  Region<Ind3D>::RegionIndices indicesIn1;
  for (auto i : rawIndicesBwd) {
    indicesIn1.push_back(i);
  }

  Region<Ind3D> region(indicesIn1);

  // Check initial order
  auto regionIndices = region.getIndices();
  for (unsigned int i = 0; i < indicesIn1.size(); i++) {
    EXPECT_EQ(regionIndices[i], rawIndicesBwd[i]);
  }

  // Sort with friend
  auto region2 = sort(region);

  // Check new order
  regionIndices = region2.getIndices();
  for (unsigned int i = 0; i < indicesIn1.size(); i++) {
    EXPECT_EQ(regionIndices[i], rawIndicesFwd[i]);
  }
}

TEST_F(RegionTest, regionUniqueInPlace) {
  // Values to insert
  std::vector<int> rawIndicesIn1 = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9,
                                    8, 7, 6, 5, 4, 3, 2, 1, 0};

  Region<Ind3D>::RegionIndices indicesIn1;
  for (auto i : rawIndicesIn1) {
    indicesIn1.push_back(i);
  }

  Region<Ind3D> region(indicesIn1);
  auto regionIndices = region.getIndices();
  EXPECT_EQ(regionIndices.size(), 19);

  for (unsigned int i = 0; i < regionIndices.size(); i++) {
    EXPECT_EQ(regionIndices[i], rawIndicesIn1[i]);
  }

  // Make unique in place
  region.unique();

  auto regionIndices2 = region.getIndices();
  EXPECT_EQ(regionIndices2.size(), 10);

  for (unsigned int i = 0; i < regionIndices2.size(); i++) {
    EXPECT_EQ(regionIndices2[i], i);
  }
}

TEST_F(RegionTest, regionFriendUnique) {
  // Values to insert
  std::vector<int> rawIndicesIn1 = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9,
                                    8, 7, 6, 5, 4, 3, 2, 1, 0};

  Region<Ind3D>::RegionIndices indicesIn1;
  for (auto i : rawIndicesIn1) {
    indicesIn1.push_back(i);
  }

  Region<Ind3D> region(indicesIn1);
  auto regionIndices = region.getIndices();
  EXPECT_EQ(regionIndices.size(), 19);

  for (unsigned int i = 0; i < regionIndices.size(); i++) {
    EXPECT_EQ(regionIndices[i], rawIndicesIn1[i]);
  }

  // Make unique in place
  auto region2 = unique(region);

  auto regionIndices2 = region2.getIndices();
  EXPECT_EQ(regionIndices2.size(), 10);

  for (unsigned int i = 0; i < regionIndices2.size(); i++) {
    EXPECT_EQ(regionIndices2[i], i);
  }
}

TEST_F(RegionTest, regionMask) {
  // Values to insert
  std::vector<int> rawIndicesIn = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
  Region<Ind3D>::RegionIndices indicesIn;
  for (auto i : rawIndicesIn) {
    indicesIn.push_back(i);
  }

  std::vector<int> rawIndicesMask1 = {1, 3, 5, 7, 9};
  Region<Ind3D>::RegionIndices indicesMask1;
  for (auto i : rawIndicesMask1) {
    indicesMask1.push_back(i);
  }

  std::vector<int> rawIndicesMask2 = {11, 13, 15, 17, 19};
  Region<Ind3D>::RegionIndices indicesMask2;
  for (auto i : rawIndicesMask2) {
    indicesMask2.push_back(i);
  }

  // Create base region and two masks
  Region<Ind3D> regionIn(indicesIn);
  Region<Ind3D> mask1(indicesMask1);
  Region<Ind3D> mask2(indicesMask2);

  EXPECT_EQ(regionIn.getIndices().size(), indicesIn.size());
  EXPECT_EQ(mask1.getIndices().size(), indicesMask1.size());

  auto masked1Indices = regionIn.mask(mask1).getIndices();
  EXPECT_EQ(masked1Indices.size(), indicesIn.size() - indicesMask1.size());

  // Check values
  for (unsigned int i = 0; i < masked1Indices.size(); i++) {
    EXPECT_EQ(masked1Indices[i] % 2, 0);
  }

  // Check size of other regions not changed
  EXPECT_EQ(mask1.getIndices().size(), indicesMask1.size());

  // Reset region
  regionIn.setIndices(indicesIn);

  // Check values of mask and region haven't changed
  auto regionIndices = regionIn.getIndices();
  for (unsigned int i = 0; i < regionIndices.size(); i++) {
    EXPECT_EQ(regionIndices[i], indicesIn[i]);
  }
  auto mask1Indices = mask1.getIndices();
  for (unsigned int i = 0; i < mask1Indices.size(); i++) {
    EXPECT_EQ(mask1Indices[i], indicesMask1[i]);
  }

  auto masked2Indices = regionIn.mask(mask2).getIndices();
  EXPECT_EQ(masked2Indices.size(), indicesIn.size());

  // Check size of other regions not changed
  EXPECT_EQ(mask2.getIndices().size(), indicesMask2.size());
}

TEST_F(RegionTest, regionFriendMask) {
  // Values to insert
  std::vector<int> rawIndicesIn = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
  Region<Ind3D>::RegionIndices indicesIn;
  for (auto i : rawIndicesIn) {
    indicesIn.push_back(i);
  }

  std::vector<int> rawIndicesMask1 = {1, 3, 5, 7, 9};
  Region<Ind3D>::RegionIndices indicesMask1;
  for (auto i : rawIndicesMask1) {
    indicesMask1.push_back(i);
  }

  std::vector<int> rawIndicesMask2 = {11, 13, 15, 17, 19};
  Region<Ind3D>::RegionIndices indicesMask2;
  for (auto i : rawIndicesMask2) {
    indicesMask2.push_back(i);
  }

  // Create base region and two masks
  Region<Ind3D> regionIn(indicesIn);
  Region<Ind3D> mask1(indicesMask1);
  Region<Ind3D> mask2(indicesMask2);

  EXPECT_EQ(regionIn.getIndices().size(), indicesIn.size());
  EXPECT_EQ(mask1.getIndices().size(), indicesMask1.size());

  auto masked1 = mask(regionIn, mask1);
  auto masked1Indices = masked1.getIndices();
  EXPECT_EQ(masked1Indices.size(), indicesIn.size() - indicesMask1.size());

  // Check values
  for (unsigned int i = 0; i < masked1Indices.size(); i++) {
    EXPECT_EQ(masked1Indices[i] % 2, 0);
  }

  // Check size of other regions not changed
  EXPECT_EQ(regionIn.getIndices().size(), indicesIn.size());
  EXPECT_EQ(mask1.getIndices().size(), indicesMask1.size());

  // Check values of mask and region haven't changed
  auto regionIndices = regionIn.getIndices();
  for (unsigned int i = 0; i < regionIndices.size(); i++) {
    EXPECT_EQ(regionIndices[i], indicesIn[i]);
  }
  auto mask1Indices = mask1.getIndices();
  for (unsigned int i = 0; i < mask1Indices.size(); i++) {
    EXPECT_EQ(mask1Indices[i], indicesMask1[i]);
  }

  auto masked2 = mask(regionIn, mask2);
  auto masked2Indices = masked2.getIndices();
  EXPECT_EQ(masked2Indices.size(), indicesIn.size());

  // Check size of other regions not changed
  EXPECT_EQ(regionIn.getIndices().size(), indicesIn.size());
  EXPECT_EQ(mask2.getIndices().size(), indicesMask2.size());
}

TEST_F(RegionTest, regionOperatorAdd) {
  // Values to insert
  std::vector<int> rawIndicesIn1 = {0, 1, 2, 3, 4};
  Region<Ind3D>::RegionIndices indicesIn1;
  for (auto i : rawIndicesIn1) {
    indicesIn1.push_back(i);
  }

  std::vector<int> rawIndicesIn2 = {5, 6, 7, 8, 9};
  Region<Ind3D>::RegionIndices indicesIn2;
  for (auto i : rawIndicesIn2) {
    indicesIn2.push_back(i);
  }

  // Create base regions
  Region<Ind3D> region1(indicesIn1);
  Region<Ind3D> region2(indicesIn2);

  // Check size of other regions not changed
  EXPECT_EQ(region1.getIndices().size(), indicesIn1.size());
  EXPECT_EQ(region2.getIndices().size(), indicesIn2.size());

  auto region3 = region1 + region2;
  auto indices3 = region3.getIndices();
  EXPECT_EQ(indices3.size(), indicesIn1.size() + indicesIn2.size());

  // Check values
  for (unsigned int i = 0; i < indices3.size(); i++) {
    EXPECT_EQ(indices3[i], i);
  }

  auto region4 = region1 + region2 + region2;
  auto indices4 = region4.getIndices();
  EXPECT_EQ(indices4.size(), indicesIn1.size() + 2 * indicesIn2.size());
  EXPECT_EQ(region1.getIndices().size(), indicesIn1.size());
  EXPECT_EQ(region2.getIndices().size(), indicesIn2.size());
}

TEST_F(RegionTest, regionOperatorAccumulate) {
  // Values to insert
  std::vector<int> rawIndicesIn1 = {0, 1, 2, 3, 4};
  Region<Ind3D>::RegionIndices indicesIn1;
  for (auto i : rawIndicesIn1) {
    indicesIn1.push_back(i);
  }

  std::vector<int> rawIndicesIn2 = {5, 6, 7, 8, 9};
  Region<Ind3D>::RegionIndices indicesIn2;
  for (auto i : rawIndicesIn2) {
    indicesIn2.push_back(i);
  }

  // Create base regions
  Region<Ind3D> region1(indicesIn1);
  Region<Ind3D> region2(indicesIn2);

  // Check size of other regions not changed
  EXPECT_EQ(region1.getIndices().size(), indicesIn1.size());
  EXPECT_EQ(region2.getIndices().size(), indicesIn2.size());

  auto region3 = region1;
  region3 += region2;
  auto indices3 = region3.getIndices();
  EXPECT_EQ(indices3.size(), indicesIn1.size() + indicesIn2.size());

  // Check values
  for (unsigned int i = 0; i < indices3.size(); i++) {
    EXPECT_EQ(indices3[i], i);
  }

  auto region4 = region1;
  region4 += region2 + region2;
  auto indices4 = region4.getIndices();
  EXPECT_EQ(indices4.size(), indicesIn1.size() + 2 * indicesIn2.size());
  EXPECT_EQ(region1.getIndices().size(), indicesIn1.size());
  EXPECT_EQ(region2.getIndices().size(), indicesIn2.size());
}

TEST_F(RegionTest, regionOffset) {
  // Values to insert
  std::vector<int> rawIndicesIn1 = {0, 1, 2, 3, 4};
  Region<Ind3D>::RegionIndices indicesIn1;
  for (auto i : rawIndicesIn1) {
    indicesIn1.push_back(i);
  }

  // Create base regions
  Region<Ind3D> region1(indicesIn1);
  auto inputInd = region1.getIndices();

  // Check size of region
  EXPECT_EQ(inputInd.size(), indicesIn1.size());

  // Check values
  for (unsigned int i = 0; i < inputInd.size(); i++) {
    EXPECT_EQ(inputInd[i], rawIndicesIn1[i]);
  }

  // Now offset
  region1.offset(2);
  inputInd = region1.getIndices();

  // Check size of region
  EXPECT_EQ(inputInd.size(), indicesIn1.size());

  // Check values
  for (unsigned int i = 0; i < inputInd.size(); i++) {
    EXPECT_EQ(inputInd[i], rawIndicesIn1[i] + 2);
  }

  region1.offset(-5);
  inputInd = region1.getIndices();

  // Check size of region
  EXPECT_EQ(inputInd.size(), indicesIn1.size());

  // Check values
  for (unsigned int i = 0; i < inputInd.size(); i++) {
    EXPECT_EQ(inputInd[i], rawIndicesIn1[i] - 3);
  }

  // Reset
  region1.setIndices(indicesIn1);
  region1.offset(0);
  inputInd = region1.getIndices();

  // Check size of region
  EXPECT_EQ(inputInd.size(), indicesIn1.size());

  // Check values
  for (unsigned int i = 0; i < inputInd.size(); i++) {
    EXPECT_EQ(inputInd[i], rawIndicesIn1[i]);
  }
}

TEST_F(RegionTest, regionFriendOffset) {
  // Values to insert
  std::vector<int> rawIndicesIn1 = {0, 1, 2, 3, 4};
  Region<Ind3D>::RegionIndices indicesIn1;
  for (auto i : rawIndicesIn1) {
    indicesIn1.push_back(i);
  }

  // Create base regions
  Region<Ind3D> region1(indicesIn1);
  auto inputInd = region1.getIndices();

  // Check size of region
  EXPECT_EQ(inputInd.size(), indicesIn1.size());

  // Check values
  for (unsigned int i = 0; i < inputInd.size(); i++) {
    EXPECT_EQ(inputInd[i], rawIndicesIn1[i]);
  }

  // Now offset
  auto region2 = offset(region1, 2);
  inputInd = region2.getIndices();

  // Check size of region
  EXPECT_EQ(inputInd.size(), indicesIn1.size());

  // Check values
  for (unsigned int i = 0; i < inputInd.size(); i++) {
    EXPECT_EQ(inputInd[i], rawIndicesIn1[i] + 2);
  }

  auto region3 = offset(region1, -5);
  inputInd = region3.getIndices();

  // Check size of region
  EXPECT_EQ(inputInd.size(), indicesIn1.size());

  // Check values
  for (unsigned int i = 0; i < inputInd.size(); i++) {
    EXPECT_EQ(inputInd[i], rawIndicesIn1[i] - 5);
  }

  // Reset
  auto region4 = offset(region1, 0);
  inputInd = region4.getIndices();

  // Check size of region
  EXPECT_EQ(inputInd.size(), indicesIn1.size());

  // Check values
  for (unsigned int i = 0; i < inputInd.size(); i++) {
    EXPECT_EQ(inputInd[i], rawIndicesIn1[i]);
  }

  inputInd = region1.getIndices();

  // Check size of region
  EXPECT_EQ(inputInd.size(), indicesIn1.size());

  // Check values
  for (unsigned int i = 0; i < inputInd.size(); i++) {
    EXPECT_EQ(inputInd[i], rawIndicesIn1[i]);
  }
}

TEST_F(RegionTest, regionPeriodicShift) {
  const int nmesh = RegionTest::nx * RegionTest::ny * RegionTest::nz;
  auto regionReference = mesh->getRegion("RGN_ALL");

  auto regionShiftInZ = mesh->getRegion("RGN_ALL");
  auto regionShiftInZNeg = mesh->getRegion("RGN_ALL");
  auto regionShiftInY = mesh->getRegion("RGN_ALL");
  auto regionShiftInX = mesh->getRegion("RGN_ALL");

  regionShiftInZ.periodicShift(1, RegionTest::nz);
  regionShiftInZNeg.periodicShift(-1, RegionTest::nz);
  regionShiftInY.periodicShift(RegionTest::nz, RegionTest::nz * RegionTest::ny);
  regionShiftInX.periodicShift(RegionTest::nz * RegionTest::ny, nmesh);

  auto indicesReg = regionReference.getIndices();
  auto shiftZReg = regionShiftInZ.getIndices();
  auto shiftZRegNeg = regionShiftInZNeg.getIndices();
  auto shiftYReg = regionShiftInY.getIndices();
  auto shiftXReg = regionShiftInX.getIndices();

  // Check size
  EXPECT_EQ(indicesReg.size(), nmesh);
  EXPECT_EQ(shiftZReg.size(), nmesh);
  EXPECT_EQ(shiftZRegNeg.size(), nmesh);
  EXPECT_EQ(shiftYReg.size(), nmesh);
  EXPECT_EQ(shiftXReg.size(), nmesh);

  // Specific values -- first point
  EXPECT_EQ(indicesReg[0], 0);
  EXPECT_EQ(shiftZReg[0], 1);
  EXPECT_EQ(shiftZRegNeg[0], RegionTest::nz - 1);
  EXPECT_EQ(shiftYReg[0], RegionTest::nz);
  EXPECT_EQ(shiftXReg[0], RegionTest::nz * RegionTest::ny);

  // Specific values -- last point (nmesh -1 - period) + shift
  // Note this test could probably be made more robust
  EXPECT_EQ(indicesReg[nmesh - 1], nmesh - 1);
  EXPECT_EQ(shiftZReg[nmesh - 1], (nmesh - RegionTest::nz));
  EXPECT_EQ(shiftZRegNeg[nmesh - 1], (nmesh - 1 - 1));
  EXPECT_EQ(shiftYReg[nmesh - 1],
            (nmesh - RegionTest::nz * RegionTest::ny) + RegionTest::nz - 1);
  EXPECT_EQ(shiftXReg[nmesh - 1], RegionTest::nz * RegionTest::ny - 1);

  // Value range
  for (unsigned int i = 0; i < nmesh; i++) {
    EXPECT_TRUE(indicesReg[i] < nmesh);
    EXPECT_TRUE(indicesReg[i] >= 0);
    EXPECT_TRUE(shiftZReg[i] < nmesh);
    EXPECT_TRUE(shiftZReg[i] >= 0);
    EXPECT_TRUE(shiftZRegNeg[i] < nmesh);
    EXPECT_TRUE(shiftZRegNeg[i] >= 0);
    EXPECT_TRUE(shiftYReg[i] < nmesh);
    EXPECT_TRUE(shiftYReg[i] >= 0);
    EXPECT_TRUE(shiftXReg[i] < nmesh);
    EXPECT_TRUE(shiftXReg[i] >= 0);
  }
}

TEST_F(RegionTest, regionGetStatsHomogenous) {
  int numBlocks = 10;
  int maxBlockSize = 64;
  int minBlockSize = 64;
  int smallBlockSize = static_cast<int>(maxBlockSize * 0.5) - 1;
  smallBlockSize = smallBlockSize < 1 ? 1 : smallBlockSize;
  int numMaxBlocks = numBlocks;
  int numMinBlocks = 0;
  int numSmallBlocks = 0;
  BoutReal maxImbalance =
      static_cast<BoutReal>(maxBlockSize) / static_cast<BoutReal>(minBlockSize);

  Region<Ind3D>::ContiguousBlocks blocks;

  for (int i = 0; i < numMaxBlocks; i++) {
    blocks.push_back(Region<Ind3D>::ContiguousBlock(0, maxBlockSize));
  }
  for (int i = 0; i < numMinBlocks; i++) {
    blocks.push_back(Region<Ind3D>::ContiguousBlock(0, minBlockSize));
  }
  for (int i = 0; i < numSmallBlocks; i++) {
    blocks.push_back(Region<Ind3D>::ContiguousBlock(0, smallBlockSize));
  }

  Region<Ind3D> region(blocks);
  auto stats = region.getStats();

  EXPECT_EQ(stats.numBlocks, numBlocks);
  EXPECT_EQ(stats.minBlockSize, minBlockSize);
  EXPECT_EQ(stats.numMinBlocks, numMaxBlocks); // As maxBlockSize == minBlockSize
  EXPECT_EQ(stats.maxBlockSize, maxBlockSize);
  EXPECT_EQ(stats.numMaxBlocks, numMaxBlocks);
  EXPECT_EQ(stats.numSmallBlocks, numSmallBlocks);
  EXPECT_EQ(stats.maxImbalance, maxImbalance);

  std::ostringstream strRepresentation;
  strRepresentation << stats;
  std::ostringstream expectedStrRepresentation;
  expectedStrRepresentation << "Total blocks : "<< numBlocks;
  expectedStrRepresentation << ", " << "min(count)/max(count) :";
  expectedStrRepresentation << " " << minBlockSize << " (" << numMaxBlocks << ")/";
  expectedStrRepresentation << " " << maxBlockSize << " (" << numMaxBlocks << ")";
  expectedStrRepresentation << ", " << "Max imbalance : " << maxImbalance;
  expectedStrRepresentation << ", " << "Small block count : " << numSmallBlocks;
  EXPECT_EQ(strRepresentation.str(), expectedStrRepresentation.str());

}

TEST_F(RegionTest, regionGetStatsHeterogenous) {
  int numBlocks = 20;
  int maxBlockSize = 64;
  int minBlockSize = 4;
  int smallBlockSize = static_cast<int>(maxBlockSize * 0.5) - 1;
  smallBlockSize = smallBlockSize < 1 ? 1 : smallBlockSize;
  int numMaxBlocks = 10;
  int numMinBlocks = 6;
  int numExtraSmallBlocks = 4;
  BoutReal maxImbalance =
      static_cast<BoutReal>(maxBlockSize) / static_cast<BoutReal>(minBlockSize);

  Region<Ind3D>::ContiguousBlocks blocks;

  for (int i = 0; i < numMaxBlocks; i++) {
    blocks.push_back(Region<Ind3D>::ContiguousBlock(0, maxBlockSize));
  }
  for (int i = 0; i < numMinBlocks; i++) {
    blocks.push_back(Region<Ind3D>::ContiguousBlock(0, minBlockSize));
  }
  for (int i = 0; i < numExtraSmallBlocks; i++) {
    blocks.push_back(Region<Ind3D>::ContiguousBlock(0, smallBlockSize));
  }

  Region<Ind3D> region(blocks);
  auto stats = region.getStats();

  EXPECT_EQ(stats.numBlocks, numBlocks);
  EXPECT_EQ(stats.minBlockSize, minBlockSize);
  EXPECT_EQ(stats.numMinBlocks, numMinBlocks);
  EXPECT_EQ(stats.maxBlockSize, maxBlockSize);
  EXPECT_EQ(stats.numMaxBlocks, numMaxBlocks);
  EXPECT_EQ(stats.numSmallBlocks, numExtraSmallBlocks + numMinBlocks);
  EXPECT_EQ(stats.maxImbalance, maxImbalance);
}

TEST_F(RegionTest, regionGetStatsEmpty) {
  int numBlocks = 0;
  int maxBlockSize = 0;
  int minBlockSize = 0;
  int numMaxBlocks = 0;
  int numMinBlocks = 0;
  int numExtraSmallBlocks = 0;
  BoutReal maxImbalance = 0.0;

  Region<Ind3D>::ContiguousBlocks blocks;

  Region<Ind3D> region(blocks);
  auto stats = region.getStats();

  EXPECT_EQ(stats.numBlocks, numBlocks);
  EXPECT_EQ(stats.minBlockSize, minBlockSize);
  EXPECT_EQ(stats.numMinBlocks, numMinBlocks);
  EXPECT_EQ(stats.maxBlockSize, maxBlockSize);
  EXPECT_EQ(stats.numMaxBlocks, numMaxBlocks);
  EXPECT_EQ(stats.numSmallBlocks, numExtraSmallBlocks + numMinBlocks);
  EXPECT_EQ(stats.maxImbalance, maxImbalance);

  std::ostringstream strRepresentation;
  strRepresentation << stats;
  EXPECT_EQ(strRepresentation.str(), "Empty");
}

template <typename T> class RegionIndexTest : public ::testing::Test {
public:
  typedef std::list<T> List;
  static T shared_;
  T value_;
};

typedef ::testing::Types<Ind2D, Ind3D> RegionIndexTypes;
TYPED_TEST_CASE(RegionIndexTest, RegionIndexTypes);

TYPED_TEST(RegionIndexTest, MemberSize) {
  Region<TypeParam> region(0, 1, 0, 2, 0, 4, 3, 5);

  int nmesh = 2 * 3 * 5;

  EXPECT_EQ(region.size(), nmesh);
}

TYPED_TEST(RegionIndexTest, NonMemberSize) {
  Region<TypeParam> region(0, 1, 0, 2, 0, 4, 3, 5);

  int nmesh = 2 * 3 * 5;

  EXPECT_EQ(size(region), nmesh);
}

TYPED_TEST(RegionIndexTest, Begin) {
  typename Region<TypeParam>::RegionIndices region{0, 2, 4, 6, 8, 10, 12, 14, 16};
  Region<TypeParam> range(region);

  auto iter = range.begin();
  EXPECT_EQ(*iter, 0);
}

// Dereferencing an end() iterator is an error, so we need to test
// end() works a little indirectly. If the addition and less-than
// tests fail, this one is suspect!
TYPED_TEST(RegionIndexTest, End) {
  typename Region<TypeParam>::RegionIndices region{0, 2, 4, 6, 8, 10, 12, 14, 16};
  Region<TypeParam> range(region);

  auto iter = range.begin() + 8;
  auto iter_end = range.end();
  EXPECT_EQ(*iter, region[8]);
  EXPECT_EQ(*iter, 16);
  // iter_end is one-past the last element of region
  EXPECT_TRUE(iter < iter_end);
}

TYPED_TEST(RegionIndexTest, PrefixIncrement) {
  typename Region<TypeParam>::RegionIndices region{0, 2, 4, 6, 8, 10, 12, 14, 16};
  Region<TypeParam> range(region);

  auto iter = range.begin();
  ++iter;

  EXPECT_EQ(*iter, region[1]);
  EXPECT_EQ(*iter, 2);
}

TYPED_TEST(RegionIndexTest, PostfixIncrement) {
  typename Region<TypeParam>::RegionIndices region{0, 2, 4, 6, 8, 10, 12, 14, 16};
  Region<TypeParam> range(region);

  auto iter = range.begin();
  auto iter2 = iter++;

  EXPECT_EQ(*iter, region[1]);
  EXPECT_EQ(*iter, 2);
  EXPECT_EQ(*iter2, region[0]);
  EXPECT_EQ(*iter2, 0);
}

TYPED_TEST(RegionIndexTest, PrefixDecrement) {
  typename Region<TypeParam>::RegionIndices region{0, 2, 4, 6, 8, 10, 12, 14, 16};
  Region<TypeParam> range(region);

  auto iter = range.end();
  --iter;

  EXPECT_EQ(*iter, region[8]);
  EXPECT_EQ(*iter, 16);
}

TYPED_TEST(RegionIndexTest, PostfixDecrement) {
  typename Region<TypeParam>::RegionIndices region{0, 2, 4, 6, 8, 10, 12, 14, 16};
  Region<TypeParam> range(region);

  // end() is one-past-the-last element, so we need to decrement it in
  // order to be able to dereference it
  auto iter = --(range.end());
  auto iter2 = iter--;

  EXPECT_EQ(*iter, region[7]);
  EXPECT_EQ(*iter, 14);
  EXPECT_EQ(*iter2, region[8]);
  EXPECT_EQ(*iter2, 16);
}

TYPED_TEST(RegionIndexTest, NotEquals) {
  typename Region<TypeParam>::RegionIndices region{0, 2, 4, 6, 8, 10, 12, 14, 16};
  Region<TypeParam> range(region);

  auto iter = range.begin();
  auto iter2 = iter;
  ++iter2;

  EXPECT_TRUE(iter != iter2);
  ++iter;
  EXPECT_FALSE(iter != iter2);
}

TYPED_TEST(RegionIndexTest, Equals) {
  typename Region<TypeParam>::RegionIndices region{0, 2, 4, 6, 8, 10, 12, 14, 16};
  Region<TypeParam> range(region);

  auto iter = range.begin();
  auto iter2 = iter;
  ++iter2;
  ++iter;

  EXPECT_TRUE(iter == iter2);
  ++iter;
  EXPECT_FALSE(iter == iter2);
}

TYPED_TEST(RegionIndexTest, LessThan) {
  typename Region<TypeParam>::RegionIndices region{0, 2, 4, 6, 8, 10, 12, 14, 16};
  Region<TypeParam> range(region);

  auto iter = range.begin();
  auto iter2 = iter;
  ++iter2;

  EXPECT_TRUE(iter < iter2);
  ++iter;
  EXPECT_FALSE(iter < iter2);
  ++iter;
  EXPECT_FALSE(iter < iter2);
}

TYPED_TEST(RegionIndexTest, MoreThan) {
  typename Region<TypeParam>::RegionIndices region{0, 2, 4, 6, 8, 10, 12, 14, 16};
  Region<TypeParam> range(region);

  auto iter = range.begin();
  auto iter2 = iter;
  ++iter2;

  EXPECT_TRUE(iter2 > iter);
  ++iter;
  EXPECT_FALSE(iter2 > iter);
  ++iter;
  EXPECT_FALSE(iter2 > iter);
}

TYPED_TEST(RegionIndexTest, LessThanOrEqualTo) {
  typename Region<TypeParam>::RegionIndices region{0, 2, 4, 6, 8, 10, 12, 14, 16};
  Region<TypeParam> range(region);

  auto iter = range.begin();
  auto iter2 = iter;
  ++iter2;

  EXPECT_TRUE(iter <= iter2);
  ++iter;
  EXPECT_TRUE(iter <= iter2);
  ++iter;
  EXPECT_FALSE(iter <= iter2);
}

TYPED_TEST(RegionIndexTest, MoreThanOrEqualTo) {
  typename Region<TypeParam>::RegionIndices region{0, 2, 4, 6, 8, 10, 12, 14, 16};
  Region<TypeParam> range(region);

  auto iter = range.begin();
  auto iter2 = iter;
  ++iter2;

  EXPECT_TRUE(iter2 >= iter);
  ++iter;
  EXPECT_TRUE(iter2 >= iter);
  ++iter;
  EXPECT_FALSE(iter2 >= iter);
}

TYPED_TEST(RegionIndexTest, PlusEqualsInt) {
  typename Region<TypeParam>::RegionIndices region{0, 2, 4, 6, 8, 10, 12, 14, 16};
  Region<TypeParam> range(region);

  auto iter = range.begin();

  iter += 2;
  EXPECT_EQ(*iter, region[2]);
  EXPECT_EQ(*iter, 4);
}

TYPED_TEST(RegionIndexTest, PlusInt) {
  typename Region<TypeParam>::RegionIndices region{0, 2, 4, 6, 8, 10, 12, 14, 16};
  Region<TypeParam> range(region);

  auto iter = range.begin();

  auto iter2 = iter + 3;
  EXPECT_EQ(*iter2, region[3]);
  EXPECT_EQ(*iter2, 6);
}

TYPED_TEST(RegionIndexTest, IntPlus) {
  typename Region<TypeParam>::RegionIndices region{0, 2, 4, 6, 8, 10, 12, 14, 16};
  Region<TypeParam> range(region);

  auto iter = range.begin();

  auto iter2 = 5 + iter;
  EXPECT_EQ(*iter2, region[5]);
  EXPECT_EQ(*iter2, 10);
}

TYPED_TEST(RegionIndexTest, MinusEqualsInt) {
  typename Region<TypeParam>::RegionIndices region{0, 2, 4, 6, 8, 10, 12, 14, 16};
  Region<TypeParam> range(region);

  auto iter = range.end();

  iter -= 2;
  EXPECT_EQ(*iter, region[7]);
  EXPECT_EQ(*iter, 14);
}

TYPED_TEST(RegionIndexTest, MinusInt) {
  typename Region<TypeParam>::RegionIndices region{0, 2, 4, 6, 8, 10, 12, 14, 16};
  Region<TypeParam> range(region);

  auto iter = range.end();

  auto iter2 = iter - 3;
  EXPECT_EQ(*iter2, region[6]);
  EXPECT_EQ(*iter2, 12);
}

TYPED_TEST(RegionIndexTest, MinusIterator) {
  typename Region<TypeParam>::RegionIndices region{0, 2, 4, 6, 8, 10, 12, 14, 16};
  Region<TypeParam> range(region);

  auto start = range.begin();
  auto end = range.end();

  EXPECT_EQ(end - start, region.size());
}

TYPED_TEST(RegionIndexTest, IndexInt) {
  typename Region<TypeParam>::RegionIndices region{0, 2, 4, 6, 8, 10, 12, 14, 16};
  Region<TypeParam> range(region);

  auto iter = range.begin();

  EXPECT_EQ(iter[4], region[4]);
  EXPECT_EQ(iter[4], 8);
}

TYPED_TEST(RegionIndexTest, Iteration) {
  typename Region<TypeParam>::RegionIndices region{0, 2, 4, 6, 8, 10, 12};
  Region<TypeParam> range(region);
  typename Region<TypeParam>::RegionIndices region2;

  int count = 0;
  for (auto iter2 = range.begin(); iter2 != range.end(); ++iter2) {
    ++count;
    region2.push_back(*iter2);
  }

  EXPECT_EQ(count, region.size());
  EXPECT_EQ(region2, region);
}

TYPED_TEST(RegionIndexTest, RangeBasedForLoop) {
  typename Region<TypeParam>::RegionIndices region{0, 2, 4, 6, 8, 10, 12};
  Region<TypeParam> range(region);
  typename Region<TypeParam>::RegionIndices region2;

  int count = 0;
  for (const auto &iter : range) {
    ++count;
    region2.push_back(iter);
  }

  EXPECT_EQ(count, region.size());
  EXPECT_EQ(region2, region);
}

/////////////////////////////////////////////////////////
// Type-parameterised tests for Ind2D, Ind3D

template <typename T>
class FieldIndexTest : public ::testing::Test {
 public:
  typedef std::list<T> List;
  static T shared_;
  T value_;
};

typedef ::testing::Types<Ind2D, Ind3D> FieldIndexTypes;
TYPED_TEST_CASE(FieldIndexTest, FieldIndexTypes);

TYPED_TEST(FieldIndexTest, Constructor) {
  TypeParam index(1);
  EXPECT_EQ(index, 1);
}

TYPED_TEST(FieldIndexTest, CopyConstructor) {
  TypeParam index(2);
  TypeParam index2(index);
  EXPECT_EQ(index2, 2);
}

TYPED_TEST(FieldIndexTest, Assignment) {
  TypeParam index;
  index = 3;
  EXPECT_EQ(index, 3);
}

TYPED_TEST(FieldIndexTest, CopyAssignment) {
  TypeParam index(4);
  TypeParam index2;
  index2 = index;
  EXPECT_EQ(index2, 4);
}

TYPED_TEST(FieldIndexTest, PreIncrement) {
  TypeParam index(5);
  TypeParam index2;
  index2 = ++index;
  EXPECT_EQ(index2, 6);
  EXPECT_EQ(index, 6);
}

TYPED_TEST(FieldIndexTest, PostIncrement) {
  TypeParam index(7);
  TypeParam index2;
  index2 = index++;
  EXPECT_EQ(index2, 7);
  EXPECT_EQ(index, 8);
}

TYPED_TEST(FieldIndexTest, PreDecrement) {
  TypeParam index(9);
  TypeParam index2;
  index2 = --index;
  EXPECT_EQ(index2, 8);
  EXPECT_EQ(index, 8);
}

TYPED_TEST(FieldIndexTest, PostDecrement) {
  TypeParam index(10);
  TypeParam index2;
  index2 = index--;
  EXPECT_EQ(index2, 10);
  EXPECT_EQ(index, 9);
}

TYPED_TEST(FieldIndexTest, Equality) {
  TypeParam index(11);
  TypeParam index2(11);
  EXPECT_EQ(index, index2);
}

TYPED_TEST(FieldIndexTest, EqualityInt) {
  TypeParam index(11);
  EXPECT_EQ(index, 11);
}

TYPED_TEST(FieldIndexTest, Inequality) {
  TypeParam index(12);
  TypeParam index2(13);
  EXPECT_NE(index, index2);
}

TYPED_TEST(FieldIndexTest, InequalityInt) {
  TypeParam index(12);
  EXPECT_NE(index, 13);
}

TYPED_TEST(FieldIndexTest, LessThan) {
  TypeParam index(14);
  TypeParam index2(15);
  EXPECT_TRUE(index < index2);
}

TYPED_TEST(FieldIndexTest, LessThanEqualTo) {
  TypeParam index(16);
  TypeParam index2(17);
  TypeParam index3(16);
  EXPECT_TRUE(index <= index2);
  EXPECT_TRUE(index <= index3);
}

TYPED_TEST(FieldIndexTest, GreaterThan) {
  TypeParam index(19);
  TypeParam index2(18);
  EXPECT_TRUE(index > index2);
}

TYPED_TEST(FieldIndexTest, GreaterThanEqualTo) {
  TypeParam index(21);
  TypeParam index2(20);
  TypeParam index3(21);
  EXPECT_TRUE(index >= index2);
  EXPECT_TRUE(index >= index3);
}

TYPED_TEST(FieldIndexTest, Addition) {
  TypeParam index(22);
  TypeParam index2(3);
  TypeParam index3;
  index3 = index + index2;
  EXPECT_EQ(index, 22);
  EXPECT_EQ(index2, 3);
  EXPECT_EQ(index3, 25);
}

TYPED_TEST(FieldIndexTest, AdditionInt) {
  TypeParam index(22);
  TypeParam index2, index3;
  index2 = index + 1;
  index3 = 1 + index;
  EXPECT_EQ(index2, 23);
  EXPECT_EQ(index3, 23);
}

TYPED_TEST(FieldIndexTest, Subtraction) {
  TypeParam index(23);
  TypeParam index2(3);
  TypeParam index3;
  index3 = index - index2;
  EXPECT_EQ(index, 23);
  EXPECT_EQ(index2, 3);
  EXPECT_EQ(index3, 20);
}

TYPED_TEST(FieldIndexTest, SubtractionInt) {
  TypeParam index(24);
  TypeParam index2;
  index2 = index - 1;
  EXPECT_EQ(index2, 23);
}

TYPED_TEST(FieldIndexTest, InPlaceAddition) {
  TypeParam index(25);
  index += 1;
  EXPECT_EQ(index, 26);
}

TYPED_TEST(FieldIndexTest, InPlaceSubtraction) {
  TypeParam index(27);
  index -= 1;
  EXPECT_EQ(index, 26);
}

TYPED_TEST(FieldIndexTest, Modulus) {
  TypeParam index(28);
  TypeParam index2;
  index2 = index % 27;
  EXPECT_EQ(index, 28);
  EXPECT_EQ(index2, 1);
}
