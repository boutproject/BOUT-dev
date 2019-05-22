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
#include <random>
#include <sstream>
#include <type_traits>
#include <vector>

/// Global mesh
namespace bout{
namespace globals{
extern Mesh *mesh;
} // namespace globals
} // namespace bout

// The unit tests use the global mesh
using namespace bout::globals;

/// Test fixture to make sure the global mesh is our fake one
using RegionTest = FakeMeshFixture;

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
      indicesIn.emplace_back(i);
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
    EXPECT_EQ((regionBlocks[i].second - 1).ind, blocksIn[i].second);
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
    EXPECT_EQ(regionIndices[i].ind, i);
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
  const auto &region = mesh->getRegion3D("RGN_ALL");

  // Need to use a Field3D as a jig as OpenMP complicates things here
  Field3D a{0.};
  BOUT_FOR(i, region) {
    a[i] = 1.0;
  }

  for (int i = 0; i < mesh->LocalNx; ++i) {
    for (int j = 0; j < mesh->LocalNy; ++j) {
      for (int k = 0; k < mesh->LocalNz; ++k) {
        EXPECT_DOUBLE_EQ(a(i, j, k), 1.0);
      }
    }
  }
}

TEST_F(RegionTest, regionLoopNoBndry) {
  const auto &region = mesh->getRegion3D("RGN_NOBNDRY");

  Field3D a{0.};
  BOUT_FOR(i, region) {
    a[i] = 1.0;
  }

  const int nmesh = RegionTest::nx * RegionTest::ny * RegionTest::nz;
  const int ninner =
      (mesh->LocalNz * (1 + mesh->xend - mesh->xstart) * (1 + mesh->yend - mesh->ystart));
  int numExpectNotMatching = nmesh - ninner;

  int numNotMatching = 0;
  int numMatching = 0;
  for (int i = 0; i < mesh->LocalNx; ++i) {
    for (int j = 0; j < mesh->LocalNy; ++j) {
      for (int k = 0; k < mesh->LocalNz; ++k) {
        if (a(i, j, k) != 1.0) {
          numNotMatching++;
        } else {
          numMatching++;
        }
      }
    }
  }
  EXPECT_EQ(numNotMatching, numExpectNotMatching);
  EXPECT_EQ(numMatching, ninner);
}

TEST_F(RegionTest, regionLoopAllSerial) {
  const auto &region = mesh->getRegion3D("RGN_ALL");

  int count = 0;
  BOUT_FOR_SERIAL(i, region) {
    ++count;
  }

  const int nmesh = RegionTest::nx * RegionTest::ny * RegionTest::nz;

  EXPECT_EQ(count, nmesh);
}

TEST_F(RegionTest, regionLoopNoBndrySerial) {
  const auto &region = mesh->getRegion3D("RGN_NOBNDRY");

  int count = 0;
  BOUT_FOR_SERIAL(i, region) {
    ++count;
  }

  const int ninner =
      (mesh->LocalNz * (1 + mesh->xend - mesh->xstart) * (1 + mesh->yend - mesh->ystart));

  EXPECT_EQ(count, ninner);
}

TEST_F(RegionTest, regionLoopAllSection) {
  const auto &region = mesh->getRegion3D("RGN_ALL");

  int count = 0;
  BOUT_OMP(parallel) {
    BOUT_FOR_OMP(i, region, for reduction(+:count)) {
      ++count;
    }
  }

  const int nmesh = RegionTest::nx * RegionTest::ny * RegionTest::nz;

  EXPECT_EQ(count, nmesh);
}

TEST_F(RegionTest, regionLoopNoBndrySection) {
  const auto &region = mesh->getRegion3D("RGN_NOBNDRY");

  int count = 0;
  BOUT_OMP(parallel) {
    BOUT_FOR_OMP(i, region, for reduction(+:count)) {
      ++count;
    }
  }

  const int ninner =
      (mesh->LocalNz * (1 + mesh->xend - mesh->xstart) * (1 + mesh->yend - mesh->ystart));

  EXPECT_EQ(count, ninner);
}

TEST_F(RegionTest, regionLoopAllInner) {
  const auto &region = mesh->getRegion3D("RGN_ALL");

  Field3D a{0.};
  BOUT_OMP(parallel) {
    BOUT_FOR_INNER(i, region) {
      a[i] = 1.0;
    }
  }

  for (int i = 0; i < mesh->LocalNx; ++i) {
    for (int j = 0; j < mesh->LocalNy; ++j) {
      for (int k = 0; k < mesh->LocalNz; ++k) {
        EXPECT_DOUBLE_EQ(a(i, j, k), 1.0);
      }
    }
  }
}

TEST_F(RegionTest, regionLoopNoBndryInner) {
  const auto &region = mesh->getRegion3D("RGN_NOBNDRY");

  Field3D a{0.};
  BOUT_OMP(parallel) {
    BOUT_FOR_INNER(i, region) {
      a[i] = 1.0;
    }
  }

  const int nmesh = RegionTest::nx * RegionTest::ny * RegionTest::nz;
  const int ninner =
      (mesh->LocalNz * (1 + mesh->xend - mesh->xstart) * (1 + mesh->yend - mesh->ystart));
  int numExpectNotMatching = nmesh - ninner;

  int numNotMatching = 0;
  int numMatching = 0;
  for (int i = 0; i < mesh->LocalNx; ++i) {
    for (int j = 0; j < mesh->LocalNy; ++j) {
      for (int k = 0; k < mesh->LocalNz; ++k) {
        if (a(i, j, k) != 1.0) {
          numNotMatching++;
        } else {
          numMatching++;
        }
      }
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
      indicesIn.push_back(Ind3D{i});
    }
  }

  // This is the sorted region and indices
  Region<Ind3D> regionSortedIn(indicesIn);
  Region<Ind3D>::RegionIndices regionIndicesSortedIn = regionSortedIn.getIndices();

  // Now shuffle the order and create a new region
  std::shuffle(std::begin(indicesIn), std::end(indicesIn), std::mt19937());
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
    indicesIn1.push_back(Ind3D{i});
  }

  // This is the sorted region and indices
  Region<Ind3D> regionIn1(indicesIn1);
  Region<Ind3D>::RegionIndices regionIndicesIn1 = regionIn1.getIndices();

  // Now get a unique version of the region
  Region<Ind3D> regionUnique1 = regionIn1.asUnique();
  Region<Ind3D>::RegionIndices regionIndicesUnique1 = regionUnique1.getIndices();

  EXPECT_EQ(regionIndicesUnique1.size(), 10);
  for (unsigned int i = 0; i < regionIndicesUnique1.size(); i++) {
    EXPECT_EQ(regionIndicesUnique1[i].ind, i);
  }

  std::vector<int> rawIndicesIn2 = {0, 0, 0, 2, 2, 2, 5, 4, 4, 5, 10, 12, 11, 9};

  Region<Ind3D>::RegionIndices indicesIn2;
  for (auto i : rawIndicesIn2) {
    indicesIn2.push_back(Ind3D{i});
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
    indicesIn3.push_back(Ind3D{i});
  }

  // This is the sorted region and indices
  Region<Ind3D> regionIn3(indicesIn3);
  Region<Ind3D>::RegionIndices regionIndicesIn3 = regionIn3.getIndices();

  // Now get a unique version of the region
  Region<Ind3D> regionUnique3 = regionIn3.asUnique();
  Region<Ind3D>::RegionIndices regionIndicesUnique3 = regionUnique3.getIndices();

  EXPECT_EQ(regionIndicesUnique3.size(), 10);
  for (unsigned int i = 0; i < regionIndicesUnique3.size(); i++) {
    EXPECT_EQ(regionIndicesUnique3[i].ind, i);
  }
}

TEST_F(RegionTest, regionSetIndices) {
  // Values to insert
  std::vector<int> rawIndicesIn1 = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};

  std::vector<int> rawIndicesIn2 = {10, 11, 12};

  Region<Ind3D>::RegionIndices indicesIn1;
  for (auto i : rawIndicesIn1) {
    indicesIn1.push_back(Ind3D{i});
  }

  Region<Ind3D>::RegionIndices indicesIn2;
  for (auto i : rawIndicesIn2) {
    indicesIn2.push_back(Ind3D{i});
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
    EXPECT_EQ(indices3[i].ind, i);
  }
}

TEST_F(RegionTest, regionSortInPlace) {
  // Values to insert
  std::vector<int> rawIndicesBwd = {9, 8, 7, 6, 5, 4, 3, 2, 1, 0};
  std::vector<int> rawIndicesFwd = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};

  Region<Ind3D>::RegionIndices indicesIn1;
  for (auto i : rawIndicesBwd) {
    indicesIn1.push_back(Ind3D{i});
  }

  Region<Ind3D> region(indicesIn1);

  // Check initial order
  auto regionIndices = region.getIndices();
  for (unsigned int i = 0; i < indicesIn1.size(); i++) {
    EXPECT_EQ(regionIndices[i].ind, rawIndicesBwd[i]);
  }

  // Sort in place
  region.sort();

  // Check new order
  regionIndices = region.getIndices();
  for (unsigned int i = 0; i < indicesIn1.size(); i++) {
    EXPECT_EQ(regionIndices[i].ind, rawIndicesFwd[i]);
  }
}

TEST_F(RegionTest, regionFriendSort) {
  // Values to insert
  std::vector<int> rawIndicesBwd = {9, 8, 7, 6, 5, 4, 3, 2, 1, 0};
  std::vector<int> rawIndicesFwd = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};

  Region<Ind3D>::RegionIndices indicesIn1;
  for (auto i : rawIndicesBwd) {
    indicesIn1.push_back(Ind3D{i});
  }

  Region<Ind3D> region(indicesIn1);

  // Check initial order
  auto regionIndices = region.getIndices();
  for (unsigned int i = 0; i < indicesIn1.size(); i++) {
    EXPECT_EQ(regionIndices[i].ind, rawIndicesBwd[i]);
  }

  // Sort with friend
  auto region2 = sort(region);

  // Check new order
  regionIndices = region2.getIndices();
  for (unsigned int i = 0; i < indicesIn1.size(); i++) {
    EXPECT_EQ(regionIndices[i].ind, rawIndicesFwd[i]);
  }
}

TEST_F(RegionTest, regionUniqueInPlace) {
  // Values to insert
  std::vector<int> rawIndicesIn1 = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9,
                                    8, 7, 6, 5, 4, 3, 2, 1, 0};

  Region<Ind3D>::RegionIndices indicesIn1;
  for (auto i : rawIndicesIn1) {
    indicesIn1.push_back(Ind3D{i});
  }

  Region<Ind3D> region(indicesIn1);
  auto regionIndices = region.getIndices();
  EXPECT_EQ(regionIndices.size(), 19);

  for (unsigned int i = 0; i < regionIndices.size(); i++) {
    EXPECT_EQ(regionIndices[i].ind, rawIndicesIn1[i]);
  }

  // Make unique in place
  region.unique();

  auto regionIndices2 = region.getIndices();
  EXPECT_EQ(regionIndices2.size(), 10);

  for (unsigned int i = 0; i < regionIndices2.size(); i++) {
    EXPECT_EQ(regionIndices2[i].ind, i);
  }
}

TEST_F(RegionTest, regionFriendUnique) {
  // Values to insert
  std::vector<int> rawIndicesIn1 = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9,
                                    8, 7, 6, 5, 4, 3, 2, 1, 0};

  Region<Ind3D>::RegionIndices indicesIn1;
  for (auto i : rawIndicesIn1) {
    indicesIn1.push_back(Ind3D{i});
  }

  Region<Ind3D> region(indicesIn1);
  auto regionIndices = region.getIndices();
  EXPECT_EQ(regionIndices.size(), 19);

  for (unsigned int i = 0; i < regionIndices.size(); i++) {
    EXPECT_EQ(regionIndices[i].ind, rawIndicesIn1[i]);
  }

  // Make unique in place
  auto region2 = unique(region);

  auto regionIndices2 = region2.getIndices();
  EXPECT_EQ(regionIndices2.size(), 10);

  for (unsigned int i = 0; i < regionIndices2.size(); i++) {
    EXPECT_EQ(regionIndices2[i].ind, i);
  }
}

TEST_F(RegionTest, regionMask) {
  // Values to insert
  std::vector<int> rawIndicesIn = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
  Region<Ind3D>::RegionIndices indicesIn;
  for (auto i : rawIndicesIn) {
    indicesIn.push_back(Ind3D{i});
  }

  std::vector<int> rawIndicesMask1 = {1, 3, 5, 7, 9};
  Region<Ind3D>::RegionIndices indicesMask1;
  for (auto i : rawIndicesMask1) {
    indicesMask1.push_back(Ind3D{i});
  }

  std::vector<int> rawIndicesMask2 = {11, 13, 15, 17, 19};
  Region<Ind3D>::RegionIndices indicesMask2;
  for (auto i : rawIndicesMask2) {
    indicesMask2.push_back(Ind3D{i});
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
  for (auto& masked1Index : masked1Indices) {
    EXPECT_EQ((masked1Index % 2).ind, 0);
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
    indicesIn.push_back(Ind3D{i});
  }

  std::vector<int> rawIndicesMask1 = {1, 3, 5, 7, 9};
  Region<Ind3D>::RegionIndices indicesMask1;
  for (auto i : rawIndicesMask1) {
    indicesMask1.push_back(Ind3D{i});
  }

  std::vector<int> rawIndicesMask2 = {11, 13, 15, 17, 19};
  Region<Ind3D>::RegionIndices indicesMask2;
  for (auto i : rawIndicesMask2) {
    indicesMask2.push_back(Ind3D{i});
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
  for (auto& masked1Index : masked1Indices) {
    EXPECT_EQ((masked1Index % 2).ind, 0);
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
    indicesIn1.push_back(Ind3D{i});
  }

  std::vector<int> rawIndicesIn2 = {5, 6, 7, 8, 9};
  Region<Ind3D>::RegionIndices indicesIn2;
  for (auto i : rawIndicesIn2) {
    indicesIn2.push_back(Ind3D{i});
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
    EXPECT_EQ(indices3[i].ind, i);
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
    indicesIn1.push_back(Ind3D{i});
  }

  std::vector<int> rawIndicesIn2 = {5, 6, 7, 8, 9};
  Region<Ind3D>::RegionIndices indicesIn2;
  for (auto i : rawIndicesIn2) {
    indicesIn2.push_back(Ind3D{i});
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
    EXPECT_EQ(indices3[i].ind, i);
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
    indicesIn1.push_back(Ind3D{i});
  }

  // Create base regions
  Region<Ind3D> region1(indicesIn1);
  auto inputInd = region1.getIndices();

  // Check size of region
  EXPECT_EQ(inputInd.size(), indicesIn1.size());

  // Check values
  for (unsigned int i = 0; i < inputInd.size(); i++) {
    EXPECT_EQ(inputInd[i].ind, rawIndicesIn1[i]);
  }

  // Now offset
  region1.offset(2);
  inputInd = region1.getIndices();

  // Check size of region
  EXPECT_EQ(inputInd.size(), indicesIn1.size());

  // Check values
  for (unsigned int i = 0; i < inputInd.size(); i++) {
    EXPECT_EQ(inputInd[i].ind, rawIndicesIn1[i] + 2);
  }

  region1.offset(-5);
  inputInd = region1.getIndices();

  // Check size of region
  EXPECT_EQ(inputInd.size(), indicesIn1.size());

  // Check values
  for (unsigned int i = 0; i < inputInd.size(); i++) {
    EXPECT_EQ(inputInd[i].ind, rawIndicesIn1[i] - 3);
  }

  // Reset
  region1.setIndices(indicesIn1);
  region1.offset(0);
  inputInd = region1.getIndices();

  // Check size of region
  EXPECT_EQ(inputInd.size(), indicesIn1.size());

  // Check values
  for (unsigned int i = 0; i < inputInd.size(); i++) {
    EXPECT_EQ(inputInd[i].ind, rawIndicesIn1[i]);
  }
}

TEST_F(RegionTest, regionFriendOffset) {
  // Values to insert
  std::vector<int> rawIndicesIn1 = {0, 1, 2, 3, 4};
  Region<Ind3D>::RegionIndices indicesIn1;
  for (auto i : rawIndicesIn1) {
    indicesIn1.push_back(Ind3D{i});
  }

  // Create base regions
  Region<Ind3D> region1(indicesIn1);
  auto inputInd = region1.getIndices();

  // Check size of region
  EXPECT_EQ(inputInd.size(), indicesIn1.size());

  // Check values
  for (unsigned int i = 0; i < inputInd.size(); i++) {
    EXPECT_EQ(inputInd[i].ind, rawIndicesIn1[i]);
  }

  // Now offset
  auto region2 = offset(region1, 2);
  inputInd = region2.getIndices();

  // Check size of region
  EXPECT_EQ(inputInd.size(), indicesIn1.size());

  // Check values
  for (unsigned int i = 0; i < inputInd.size(); i++) {
    EXPECT_EQ(inputInd[i].ind, rawIndicesIn1[i] + 2);
  }

  auto region3 = offset(region1, -5);
  inputInd = region3.getIndices();

  // Check size of region
  EXPECT_EQ(inputInd.size(), indicesIn1.size());

  // Check values
  for (unsigned int i = 0; i < inputInd.size(); i++) {
    EXPECT_EQ(inputInd[i].ind, rawIndicesIn1[i] - 5);
  }

  // Reset
  auto region4 = offset(region1, 0);
  inputInd = region4.getIndices();

  // Check size of region
  EXPECT_EQ(inputInd.size(), indicesIn1.size());

  // Check values
  for (unsigned int i = 0; i < inputInd.size(); i++) {
    EXPECT_EQ(inputInd[i].ind, rawIndicesIn1[i]);
  }

  inputInd = region1.getIndices();

  // Check size of region
  EXPECT_EQ(inputInd.size(), indicesIn1.size());

  // Check values
  for (unsigned int i = 0; i < inputInd.size(); i++) {
    EXPECT_EQ(inputInd[i].ind, rawIndicesIn1[i]);
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
  EXPECT_EQ(indicesReg[0].ind, 0);
  EXPECT_EQ(shiftZReg[0].ind, 1);
  EXPECT_EQ(shiftZRegNeg[0].ind, RegionTest::nz - 1);
  EXPECT_EQ(shiftYReg[0].ind, RegionTest::nz);
  EXPECT_EQ(shiftXReg[0].ind, RegionTest::nz * RegionTest::ny);

  // Specific values -- last point (nmesh -1 - period) + shift
  // Note this test could probably be made more robust
  EXPECT_EQ(indicesReg[nmesh - 1].ind, nmesh - 1);
  EXPECT_EQ(shiftZReg[nmesh - 1].ind, (nmesh - RegionTest::nz));
  EXPECT_EQ(shiftZRegNeg[nmesh - 1].ind, (nmesh - 1 - 1));
  EXPECT_EQ(shiftYReg[nmesh - 1].ind,
            (nmesh - RegionTest::nz * RegionTest::ny) + RegionTest::nz - 1);
  EXPECT_EQ(shiftXReg[nmesh - 1].ind, RegionTest::nz * RegionTest::ny - 1);

  // Value range
  for (unsigned int i = 0; i < nmesh; i++) {
    EXPECT_TRUE(indicesReg[i].ind < nmesh);
    EXPECT_TRUE(indicesReg[i].ind >= 0);
    EXPECT_TRUE(shiftZReg[i].ind < nmesh);
    EXPECT_TRUE(shiftZReg[i].ind >= 0);
    EXPECT_TRUE(shiftZRegNeg[i].ind < nmesh);
    EXPECT_TRUE(shiftZRegNeg[i].ind >= 0);
    EXPECT_TRUE(shiftYReg[i].ind < nmesh);
    EXPECT_TRUE(shiftYReg[i].ind >= 0);
    EXPECT_TRUE(shiftXReg[i].ind < nmesh);
    EXPECT_TRUE(shiftXReg[i].ind >= 0);
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
    blocks.push_back(Region<Ind3D>::ContiguousBlock(Ind3D{0}, Ind3D{maxBlockSize}));
  }
  for (int i = 0; i < numMinBlocks; i++) {
    blocks.push_back(Region<Ind3D>::ContiguousBlock(Ind3D{0}, Ind3D{minBlockSize}));
  }
  for (int i = 0; i < numSmallBlocks; i++) {
    blocks.push_back(Region<Ind3D>::ContiguousBlock(Ind3D{0}, Ind3D{smallBlockSize}));
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
    blocks.push_back(Region<Ind3D>::ContiguousBlock(Ind3D{0}, Ind3D{maxBlockSize}));
  }
  for (int i = 0; i < numMinBlocks; i++) {
    blocks.push_back(Region<Ind3D>::ContiguousBlock(Ind3D{0}, Ind3D{minBlockSize}));
  }
  for (int i = 0; i < numExtraSmallBlocks; i++) {
    blocks.push_back(Region<Ind3D>::ContiguousBlock(Ind3D{0}, Ind3D{smallBlockSize}));
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

TEST(RegionIndexConversionTest, Ind3DtoInd2D) {
  // This could just be:
  //     EXPECT_FALSE(std::is_convertible<Ind3D, Ind2D>::value());
  // but requires C++14
  bool convert = std::is_convertible<Ind3D, Ind2D>::value;
  EXPECT_FALSE(convert);
}

TEST(RegionIndexConversionTest, Ind2DtoInd3D) {
  bool convert = std::is_convertible<Ind2D, Ind3D>::value;
  EXPECT_FALSE(convert);
}

TEST(RegionIndexConversionTest, Ind2Dtoint) {
  bool convert = std::is_convertible<Ind2D, int>::value;
  EXPECT_FALSE(convert);
}

TEST(RegionIndexConversionTest, Ind3Dtoint) {
  bool convert = std::is_convertible<Ind3D, int>::value;
  EXPECT_FALSE(convert);
}

TEST(RegionIndexConversionTest, inttoInd2D) {
  bool convert = std::is_convertible<int, Ind2D>::value;
  EXPECT_FALSE(convert);
}

TEST(RegionIndexConversionTest, inttoInd3D) {
  bool convert = std::is_convertible<int, Ind3D>::value;
  EXPECT_FALSE(convert);
}

TEST(RegionIndex2DTest, MemberSize) {
  const int nx = 2, ny = 3, nz = 1;
  Region<Ind2D> region(0, nx - 1, 0, ny - 1, 0, nz - 1, ny, nz);

  int nmesh = nx * ny * nz;

  EXPECT_EQ(region.size(), nmesh);
}

TEST(RegionIndex2DTest, NonMemberSize) {
  const int nx = 2, ny = 3, nz = 1;
  Region<Ind2D> region(0, nx - 1, 0, ny - 1, 0, nz - 1, ny, nz);

  int nmesh = nx * ny * nz;

  EXPECT_EQ(size(region), nmesh);
}

TEST(RegionIndex3DTest, MemberSize) {
  const int nx = 2, ny = 3, nz = 5;
  Region<Ind3D> region(0, nx - 1, 0, ny - 1, 0, nz - 1, ny, nz);

  int nmesh = nx * ny * nz;

  EXPECT_EQ(region.size(), nmesh);
}

TEST(RegionIndex3DTest, NonMemberSize) {
  const int nx = 2, ny = 3, nz = 5;
  Region<Ind3D> region(0, nx - 1, 0, ny - 1, 0, nz - 1, ny, nz);

  int nmesh = nx * ny * nz;

  EXPECT_EQ(size(region), nmesh);
}

template <typename T>
class RegionIndexTest : public ::testing::Test {
public:
  ~RegionIndexTest() override = default;
};

using RegionIndexTypes = ::testing::Types<Ind2D, Ind3D, IndPerp>;
TYPED_TEST_SUITE(RegionIndexTest, RegionIndexTypes);

TYPED_TEST(RegionIndexTest, Begin) {
  typename Region<TypeParam>::RegionIndices region{
      TypeParam{0},  TypeParam{2},  TypeParam{4},  TypeParam{6}, TypeParam{8},
      TypeParam{10}, TypeParam{12}, TypeParam{14}, TypeParam{16}};
  Region<TypeParam> range(region);

  auto iter = range.begin();
  EXPECT_EQ(iter->ind, 0);
}

// Dereferencing an end() iterator is an error, so we need to test
// end() works a little indirectly. If the addition and less-than
// tests fail, this one is suspect!
TYPED_TEST(RegionIndexTest, End) {
  typename Region<TypeParam>::RegionIndices region{
      TypeParam{0},  TypeParam{2},  TypeParam{4},  TypeParam{6}, TypeParam{8},
      TypeParam{10}, TypeParam{12}, TypeParam{14}, TypeParam{16}};
  Region<TypeParam> range(region);

  auto iter = range.begin() + 8;
  auto iter_end = range.end();
  EXPECT_EQ(*iter, region[8]);
  EXPECT_EQ(iter->ind, 16);
  // iter_end is one-past the last element of region
  EXPECT_TRUE(iter < iter_end);
}

TYPED_TEST(RegionIndexTest, PrefixIncrement) {
  typename Region<TypeParam>::RegionIndices region{
      TypeParam{0},  TypeParam{2},  TypeParam{4},  TypeParam{6}, TypeParam{8},
      TypeParam{10}, TypeParam{12}, TypeParam{14}, TypeParam{16}};
  Region<TypeParam> range(region);

  auto iter = range.begin();
  ++iter;

  EXPECT_EQ(*iter, region[1]);
  EXPECT_EQ(iter->ind, 2);
}

TYPED_TEST(RegionIndexTest, PostfixIncrement) {
  typename Region<TypeParam>::RegionIndices region{
      TypeParam{0},  TypeParam{2},  TypeParam{4},  TypeParam{6}, TypeParam{8},
      TypeParam{10}, TypeParam{12}, TypeParam{14}, TypeParam{16}};
  Region<TypeParam> range(region);

  auto iter = range.begin();
  auto iter2 = iter++;

  EXPECT_EQ(*iter, region[1]);
  EXPECT_EQ(iter->ind, 2);
  EXPECT_EQ(*iter2, region[0]);
  EXPECT_EQ(iter2->ind, 0);
}

TYPED_TEST(RegionIndexTest, PrefixDecrement) {
  typename Region<TypeParam>::RegionIndices region{
      TypeParam{0},  TypeParam{2},  TypeParam{4},  TypeParam{6}, TypeParam{8},
      TypeParam{10}, TypeParam{12}, TypeParam{14}, TypeParam{16}};
  Region<TypeParam> range(region);

  auto iter = range.end();
  --iter;

  EXPECT_EQ(*iter, region[8]);
  EXPECT_EQ(iter->ind, 16);
}

TYPED_TEST(RegionIndexTest, PostfixDecrement) {
  typename Region<TypeParam>::RegionIndices region{
      TypeParam{0},  TypeParam{2},  TypeParam{4},  TypeParam{6}, TypeParam{8},
      TypeParam{10}, TypeParam{12}, TypeParam{14}, TypeParam{16}};
  Region<TypeParam> range(region);

  // end() is one-past-the-last element, so we need to decrement it in
  // order to be able to dereference it
  auto iter = --(range.end());
  auto iter2 = iter--;

  EXPECT_EQ(*iter, region[7]);
  EXPECT_EQ(iter->ind, 14);
  EXPECT_EQ(*iter2, region[8]);
  EXPECT_EQ(iter2->ind, 16);
}

TYPED_TEST(RegionIndexTest, NotEquals) {
  typename Region<TypeParam>::RegionIndices region{
      TypeParam{0},  TypeParam{2},  TypeParam{4},  TypeParam{6}, TypeParam{8},
      TypeParam{10}, TypeParam{12}, TypeParam{14}, TypeParam{16}};
  Region<TypeParam> range(region);

  auto iter = range.begin();
  auto iter2 = iter;
  ++iter2;

  EXPECT_TRUE(iter != iter2);
  ++iter;
  EXPECT_FALSE(iter != iter2);
}

TYPED_TEST(RegionIndexTest, Equals) {
  typename Region<TypeParam>::RegionIndices region{
      TypeParam{0},  TypeParam{2},  TypeParam{4},  TypeParam{6}, TypeParam{8},
      TypeParam{10}, TypeParam{12}, TypeParam{14}, TypeParam{16}};
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
  typename Region<TypeParam>::RegionIndices region{
      TypeParam{0},  TypeParam{2},  TypeParam{4},  TypeParam{6}, TypeParam{8},
      TypeParam{10}, TypeParam{12}, TypeParam{14}, TypeParam{16}};
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
  typename Region<TypeParam>::RegionIndices region{
      TypeParam{0},  TypeParam{2},  TypeParam{4},  TypeParam{6}, TypeParam{8},
      TypeParam{10}, TypeParam{12}, TypeParam{14}, TypeParam{16}};
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
  typename Region<TypeParam>::RegionIndices region{
      TypeParam{0},  TypeParam{2},  TypeParam{4},  TypeParam{6}, TypeParam{8},
      TypeParam{10}, TypeParam{12}, TypeParam{14}, TypeParam{16}};
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
  typename Region<TypeParam>::RegionIndices region{
      TypeParam{0},  TypeParam{2},  TypeParam{4},  TypeParam{6}, TypeParam{8},
      TypeParam{10}, TypeParam{12}, TypeParam{14}, TypeParam{16}};
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
  typename Region<TypeParam>::RegionIndices region{
      TypeParam{0},  TypeParam{2},  TypeParam{4},  TypeParam{6}, TypeParam{8},
      TypeParam{10}, TypeParam{12}, TypeParam{14}, TypeParam{16}};
  Region<TypeParam> range(region);

  auto iter = range.begin();

  iter += 2;
  EXPECT_EQ(*iter, region[2]);
  EXPECT_EQ(iter->ind, 4);
}

TYPED_TEST(RegionIndexTest, PlusInt) {
  typename Region<TypeParam>::RegionIndices region{
      TypeParam{0},  TypeParam{2},  TypeParam{4},  TypeParam{6}, TypeParam{8},
      TypeParam{10}, TypeParam{12}, TypeParam{14}, TypeParam{16}};
  Region<TypeParam> range(region);

  auto iter = range.begin();

  auto iter2 = iter + 3;
  EXPECT_EQ(*iter2, region[3]);
  EXPECT_EQ(iter2->ind, 6);
}

TYPED_TEST(RegionIndexTest, IntPlus) {
  typename Region<TypeParam>::RegionIndices region{
      TypeParam{0},  TypeParam{2},  TypeParam{4},  TypeParam{6}, TypeParam{8},
      TypeParam{10}, TypeParam{12}, TypeParam{14}, TypeParam{16}};
  Region<TypeParam> range(region);

  auto iter = range.begin();

  auto iter2 = 5 + iter;
  EXPECT_EQ(*iter2, region[5]);
  EXPECT_EQ(iter2->ind, 10);
}

TYPED_TEST(RegionIndexTest, MinusEqualsInt) {
  typename Region<TypeParam>::RegionIndices region{
      TypeParam{0},  TypeParam{2},  TypeParam{4},  TypeParam{6}, TypeParam{8},
      TypeParam{10}, TypeParam{12}, TypeParam{14}, TypeParam{16}};
  Region<TypeParam> range(region);

  auto iter = range.end();

  iter -= 2;
  EXPECT_EQ(*iter, region[7]);
  EXPECT_EQ(iter->ind, 14);
}

TYPED_TEST(RegionIndexTest, MinusInt) {
  typename Region<TypeParam>::RegionIndices region{
      TypeParam{0},  TypeParam{2},  TypeParam{4},  TypeParam{6}, TypeParam{8},
      TypeParam{10}, TypeParam{12}, TypeParam{14}, TypeParam{16}};
  Region<TypeParam> range(region);

  auto iter = range.end();

  auto iter2 = iter - 3;
  EXPECT_EQ(*iter2, region[6]);
  EXPECT_EQ(iter2->ind, 12);
}

TYPED_TEST(RegionIndexTest, MinusIterator) {
  typename Region<TypeParam>::RegionIndices region{
      TypeParam{0},  TypeParam{2},  TypeParam{4},  TypeParam{6}, TypeParam{8},
      TypeParam{10}, TypeParam{12}, TypeParam{14}, TypeParam{16}};
  Region<TypeParam> range(region);

  auto start = range.begin();
  auto end = range.end();

  EXPECT_EQ(end - start, region.size());
}

TYPED_TEST(RegionIndexTest, IndexInt) {
  typename Region<TypeParam>::RegionIndices region{
      TypeParam{0},  TypeParam{2},  TypeParam{4},  TypeParam{6}, TypeParam{8},
      TypeParam{10}, TypeParam{12}, TypeParam{14}, TypeParam{16}};
  Region<TypeParam> range(region);

  auto iter = range.begin();

  EXPECT_EQ(iter[4], region[4]);
  EXPECT_EQ(iter[4].ind, 8);
}

TYPED_TEST(RegionIndexTest, Iteration) {
  typename Region<TypeParam>::RegionIndices region{
      TypeParam{0}, TypeParam{2},  TypeParam{4}, TypeParam{6},
      TypeParam{8}, TypeParam{10}, TypeParam{12}};
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
  typename Region<TypeParam>::RegionIndices region{
      TypeParam{0}, TypeParam{2},  TypeParam{4}, TypeParam{6},
      TypeParam{8}, TypeParam{10}, TypeParam{12}};
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
  ~FieldIndexTest() override = default;
};

using FieldIndexTypes = ::testing::Types<Ind2D, Ind3D>;
TYPED_TEST_SUITE(FieldIndexTest, FieldIndexTypes);

TYPED_TEST(FieldIndexTest, Constructor) {
  TypeParam index(1);
  EXPECT_EQ(index.ind, 1);
}

TYPED_TEST(FieldIndexTest, CopyConstructor) {
  TypeParam index(2);
  TypeParam index2(index);
  EXPECT_EQ(index2.ind, 2);
}

TYPED_TEST(FieldIndexTest, Assignment) {
  TypeParam index;
  index = TypeParam{3};
  EXPECT_EQ(index.ind, 3);
}

TYPED_TEST(FieldIndexTest, CopyAssignment) {
  TypeParam index(4);
  TypeParam index2;
  index2 = index;
  EXPECT_EQ(index2.ind, 4);
}

TYPED_TEST(FieldIndexTest, PreIncrement) {
  TypeParam index(5);
  TypeParam index2;
  index2 = ++index;
  EXPECT_EQ(index2.ind, 6);
  EXPECT_EQ(index.ind, 6);
}

TYPED_TEST(FieldIndexTest, PostIncrement) {
  TypeParam index(7);
  TypeParam index2;
  index2 = index++;
  EXPECT_EQ(index2.ind, 7);
  EXPECT_EQ(index.ind, 8);
}

TYPED_TEST(FieldIndexTest, PreDecrement) {
  TypeParam index(9);
  TypeParam index2;
  index2 = --index;
  EXPECT_EQ(index2.ind, 8);
  EXPECT_EQ(index.ind, 8);
}

TYPED_TEST(FieldIndexTest, PostDecrement) {
  TypeParam index(10);
  TypeParam index2;
  index2 = index--;
  EXPECT_EQ(index2.ind, 10);
  EXPECT_EQ(index.ind, 9);
}

TYPED_TEST(FieldIndexTest, Equality) {
  TypeParam index(11);
  TypeParam index2(11);
  EXPECT_EQ(index, index2);
}

TYPED_TEST(FieldIndexTest, EqualityInt) {
  TypeParam index(11);
  EXPECT_EQ(index.ind, 11);
}

TYPED_TEST(FieldIndexTest, Inequality) {
  TypeParam index(12);
  TypeParam index2(13);
  EXPECT_NE(index, index2);
}

TYPED_TEST(FieldIndexTest, InequalityInt) {
  TypeParam index(12);
  EXPECT_NE(index.ind, 13);
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
  EXPECT_EQ(index.ind, 22);
  EXPECT_EQ(index2.ind, 3);
  EXPECT_EQ(index3.ind, 25);
}

TYPED_TEST(FieldIndexTest, AdditionInt) {
  TypeParam index(22);
  TypeParam index2, index3;
  index2 = index + 1;
  index3 = 1 + index;
  EXPECT_EQ(index2.ind, 23);
  EXPECT_EQ(index3.ind, 23);
}

TYPED_TEST(FieldIndexTest, Subtraction) {
  TypeParam index(23);
  TypeParam index2(3);
  TypeParam index3;
  index3 = index - index2;
  EXPECT_EQ(index.ind, 23);
  EXPECT_EQ(index2.ind, 3);
  EXPECT_EQ(index3.ind, 20);
}

TYPED_TEST(FieldIndexTest, SubtractionInt) {
  TypeParam index(24);
  TypeParam index2;
  index2 = index - 1;
  EXPECT_EQ(index2.ind, 23);
}

TYPED_TEST(FieldIndexTest, InPlaceAddition) {
  TypeParam index(25);
  index += 1;
  EXPECT_EQ(index.ind, 26);
}

TYPED_TEST(FieldIndexTest, InPlaceSubtraction) {
  TypeParam index(27);
  index -= 1;
  EXPECT_EQ(index.ind, 26);
}

TYPED_TEST(FieldIndexTest, Modulus) {
  TypeParam index(28);
  TypeParam index2;
  index2 = index % 27;
  EXPECT_EQ(index.ind, 28);
  EXPECT_EQ(index2.ind, 1);
}

/// Test fixture to make sure the global mesh is our fake one
class IndexOffsetTest : public ::testing::Test {
protected:
  IndexOffsetTest() {
    WithQuietOutput quiet{output_info};
    delete mesh;
    mesh = new FakeMesh(nx, ny, nz);
    mesh->createDefaultRegions();
  }

  ~IndexOffsetTest() override {
    delete mesh;
    mesh = nullptr;
  }

public:
  static const int nx;
  static const int ny;
  static const int nz;
};

const int IndexOffsetTest::nx = 3;
const int IndexOffsetTest::ny = 5;
const int IndexOffsetTest::nz = 7;

TEST_F(IndexOffsetTest, X) {
  const auto &region = mesh->getRegion3D("RGN_ALL");

  auto index = region.cbegin();

  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < ny; ++j) {
      for (int k = 0; k < nz; ++k) {
        EXPECT_EQ(index->x(), i);
        ++index;
      }
    }
  }
}

TEST_F(IndexOffsetTest, Y) {
  const auto &region = mesh->getRegion3D("RGN_ALL");

  auto index = region.cbegin();

  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < ny; ++j) {
      for (int k = 0; k < nz; ++k) {
        EXPECT_EQ(index->y(), j);
        ++index;
      }
    }
  }
}

TEST_F(IndexOffsetTest, Z) {
  const auto &region = mesh->getRegion3D("RGN_ALL");

  auto index = region.cbegin();

  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < ny; ++j) {
      for (int k = 0; k < nz; ++k) {
        EXPECT_EQ(index->z(), k);
        ++index;
      }
    }
  }
}

TEST_F(IndexOffsetTest, XPlusOne) {
  const auto &region = mesh->getRegion3D("RGN_ALL");

  auto index = region.cbegin();

  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < ny; ++j) {
      for (int k = 0; k < nz; ++k) {
        EXPECT_EQ(index->x(), i);
        EXPECT_EQ(index->y(), j);
        EXPECT_EQ(index->z(), k);

        if (i >= (nx - 1)) {
          // skip this point
        } else {
          EXPECT_EQ(index->xp().x(), i + 1);
          EXPECT_EQ(index->xp().y(), j);
          EXPECT_EQ(index->xp().z(), k);
        }
        ++index;
      }
    }
  }
}

TEST_F(IndexOffsetTest, YPlusOne) {
  const auto &region = mesh->getRegion3D("RGN_ALL");

  auto index = region.cbegin();

  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < ny; ++j) {
      for (int k = 0; k < nz; ++k) {
        EXPECT_EQ(index->x(), i);
        EXPECT_EQ(index->y(), j);
        EXPECT_EQ(index->z(), k);

        if (j >= (ny - 1)) {
#if CHECK > 3
          EXPECT_THROW(index->yp(), BoutException);
#endif
        } else {
          EXPECT_EQ(index->yp().x(), i);
          EXPECT_EQ(index->yp().y(), j + 1);
          EXPECT_EQ(index->yp().z(), k);
        }
        ++index;
      }
    }
  }
}

TEST_F(IndexOffsetTest, ZPlusOne) {
  const auto &region = mesh->getRegion3D("RGN_ALL");

  auto index = region.cbegin();

  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < ny; ++j) {
      for (int k = 0; k < nz; ++k) {
        EXPECT_EQ(index->x(), i);
        EXPECT_EQ(index->y(), j);
        EXPECT_EQ(index->z(), k);

        EXPECT_EQ(index->zp().x(), i);
        EXPECT_EQ(index->zp().y(), j);
        EXPECT_EQ(index->zp().z(), (k + 1) % nz);
        ++index;
      }
    }
  }
}

TEST_F(IndexOffsetTest, XPlusOneGeneric) {
  const auto &region = mesh->getRegion3D("RGN_ALL");

  auto index = region.cbegin();

  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < ny; ++j) {
      for (int k = 0; k < nz; ++k) {
        EXPECT_EQ(index->x(), i);
        EXPECT_EQ(index->y(), j);
        EXPECT_EQ(index->z(), k);

        if (i >= (nx - 1)) {
          // skip this point
        } else {
          EXPECT_EQ((index->plus<1, DIRECTION::X>().x()), i + 1);
          EXPECT_EQ((index->plus<1, DIRECTION::X>().y()), j);
          EXPECT_EQ((index->plus<1, DIRECTION::X>().z()), k);
        }
        ++index;
      }
    }
  }
}

TEST_F(IndexOffsetTest, YPlusOneGeneric) {
  const auto &region = mesh->getRegion3D("RGN_ALL");

  auto index = region.cbegin();

  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < ny; ++j) {
      for (int k = 0; k < nz; ++k) {
        EXPECT_EQ(index->x(), i);
        EXPECT_EQ(index->y(), j);
        EXPECT_EQ(index->z(), k);

        if (j >= (ny - 1)) {
#if CHECK > 3
          EXPECT_THROW(index->yp(), BoutException);
#endif
        } else {
          EXPECT_EQ((index->plus<1, DIRECTION::Y>().x()), i);
          EXPECT_EQ((index->plus<1, DIRECTION::Y>().y()), j + 1);
          EXPECT_EQ((index->plus<1, DIRECTION::Y>().z()), k);
        }
        ++index;
      }
    }
  }
}

TEST_F(IndexOffsetTest, ZPlusOneGeneric) {
  const auto &region = mesh->getRegion3D("RGN_ALL");

  auto index = region.cbegin();

  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < ny; ++j) {
      for (int k = 0; k < nz; ++k) {
        EXPECT_EQ(index->x(), i);
        EXPECT_EQ(index->y(), j);
        EXPECT_EQ(index->z(), k);

        EXPECT_EQ((index->plus<1, DIRECTION::Z>().x()), i);
        EXPECT_EQ((index->plus<1, DIRECTION::Z>().y()), j);
        EXPECT_EQ((index->plus<1, DIRECTION::Z>().z()), (k + 1) % nz);
        ++index;
      }
    }
  }
}

TEST_F(IndexOffsetTest, XMinusOne) {
  const auto &region = mesh->getRegion3D("RGN_ALL");

  auto index = region.cbegin();

  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < ny; ++j) {
      for (int k = 0; k < nz; ++k) {
        EXPECT_EQ(index->x(), i);
        EXPECT_EQ(index->y(), j);
        EXPECT_EQ(index->z(), k);

        if (i < 1) {
          // skip this point
        } else {
          EXPECT_EQ(index->xm().x(), i - 1);
          EXPECT_EQ(index->xm().y(), j);
          EXPECT_EQ(index->xm().z(), k);
        }
        ++index;
      }
    }
  }
}

TEST_F(IndexOffsetTest, YMinusOne) {
  const auto &region = mesh->getRegion3D("RGN_ALL");

  auto index = region.cbegin();

  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < ny; ++j) {
      for (int k = 0; k < nz; ++k) {
        EXPECT_EQ(index->x(), i);
        EXPECT_EQ(index->y(), j);
        EXPECT_EQ(index->z(), k);

        if (j < 1) {
#if CHECK > 3
          EXPECT_THROW(index->ym(), BoutException);
#endif
        } else {
          EXPECT_EQ(index->ym().x(), i);
          EXPECT_EQ(index->ym().y(), j - 1);
          EXPECT_EQ(index->ym().z(), k);
        }
        ++index;
      }
    }
  }
}

TEST_F(IndexOffsetTest, ZMinusOne) {
  const auto &region = mesh->getRegion3D("RGN_ALL");

  auto index = region.cbegin();

  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < ny; ++j) {
      for (int k = 0; k < nz; ++k) {
        EXPECT_EQ(index->x(), i);
        EXPECT_EQ(index->y(), j);
        EXPECT_EQ(index->z(), k);

        EXPECT_EQ(index->zm().x(), i);
        EXPECT_EQ(index->zm().y(), j);
        EXPECT_EQ(index->zm().z(), (k - 1 + nz) % nz);
        ++index;
      }
    }
  }
}

TEST_F(IndexOffsetTest, XMinusOneGeneric) {
  const auto &region = mesh->getRegion3D("RGN_ALL");

  auto index = region.cbegin();

  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < ny; ++j) {
      for (int k = 0; k < nz; ++k) {
        EXPECT_EQ(index->x(), i);
        EXPECT_EQ(index->y(), j);
        EXPECT_EQ(index->z(), k);

        if (i < 1) {
          // skip this point
        } else {
          EXPECT_EQ((index->minus<1, DIRECTION::X>().x()), i - 1);
          EXPECT_EQ((index->minus<1, DIRECTION::X>().y()), j);
          EXPECT_EQ((index->minus<1, DIRECTION::X>().z()), k);
        }
        ++index;
      }
    }
  }
}

TEST_F(IndexOffsetTest, YMinusOneGeneric) {
  const auto &region = mesh->getRegion3D("RGN_ALL");

  auto index = region.cbegin();

  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < ny; ++j) {
      for (int k = 0; k < nz; ++k) {
        EXPECT_EQ(index->x(), i);
        EXPECT_EQ(index->y(), j);
        EXPECT_EQ(index->z(), k);

        if (j < 1) {
#if CHECK > 3
          EXPECT_THROW(index->ym(), BoutException);
#endif
        } else {
          EXPECT_EQ((index->minus<1, DIRECTION::Y>().x()), i);
          EXPECT_EQ((index->minus<1, DIRECTION::Y>().y()), j - 1);
          EXPECT_EQ((index->minus<1, DIRECTION::Y>().z()), k);
        }
        ++index;
      }
    }
  }
}

TEST_F(IndexOffsetTest, ZMinusOneGeneric) {
  const auto &region = mesh->getRegion3D("RGN_ALL");

  auto index = region.cbegin();

  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < ny; ++j) {
      for (int k = 0; k < nz; ++k) {
        EXPECT_EQ(index->x(), i);
        EXPECT_EQ(index->y(), j);
        EXPECT_EQ(index->z(), k);

        EXPECT_EQ((index->minus<1, DIRECTION::Z>().x()), i);
        EXPECT_EQ((index->minus<1, DIRECTION::Z>().y()), j);
        EXPECT_EQ((index->minus<1, DIRECTION::Z>().z()), (k - 1 + nz) % nz);
        ++index;
      }
    }
  }
}

TEST_F(IndexOffsetTest, XPlusTwo) {
  const auto &region = mesh->getRegion3D("RGN_ALL");

  auto index = region.cbegin();

  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < ny; ++j) {
      for (int k = 0; k < nz; ++k) {
        EXPECT_EQ(index->x(), i);
        EXPECT_EQ(index->y(), j);
        EXPECT_EQ(index->z(), k);

        if (i >= (nx - 2)) {
          // skip this point
        } else {
          EXPECT_EQ(index->xpp().x(), i + 2);
          EXPECT_EQ(index->xpp().y(), j);
          EXPECT_EQ(index->xpp().z(), k);
        }
        ++index;
      }
    }
  }
}

TEST_F(IndexOffsetTest, YPlusTwo) {
  const auto &region = mesh->getRegion3D("RGN_ALL");

  auto index = region.cbegin();

  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < ny; ++j) {
      for (int k = 0; k < nz; ++k) {
        EXPECT_EQ(index->x(), i);
        EXPECT_EQ(index->y(), j);
        EXPECT_EQ(index->z(), k);

        if (j >= (ny - 2)) {
#if CHECK > 3
          EXPECT_THROW(index->ypp(), BoutException);
#endif
        } else {
          EXPECT_EQ(index->ypp().x(), i);
          EXPECT_EQ(index->ypp().y(), j + 2);
          EXPECT_EQ(index->ypp().z(), k);
        }
        ++index;
      }
    }
  }
}

TEST_F(IndexOffsetTest, ZPlusTwo) {
  const auto &region = mesh->getRegion3D("RGN_ALL");

  auto index = region.cbegin();

  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < ny; ++j) {
      for (int k = 0; k < nz; ++k) {
        EXPECT_EQ(index->x(), i);
        EXPECT_EQ(index->y(), j);
        EXPECT_EQ(index->z(), k);

        EXPECT_EQ(index->zpp().x(), i);
        EXPECT_EQ(index->zpp().y(), j);
        EXPECT_EQ(index->zpp().z(), (k + 2 + nz) % nz);
        ++index;
      }
    }
  }
}

TEST_F(IndexOffsetTest, XMinusTwo) {
  const auto &region = mesh->getRegion3D("RGN_ALL");

  auto index = region.cbegin();

  for (int j = 0; j < ny; ++j) {
    for (int k = 0; k < nz; ++k) {
      ++index;
      ++index;
    }
  }
  for (int i = 2; i < nx; ++i) {
    for (int j = 0; j < ny; ++j) {
      for (int k = 0; k < nz; ++k) {
        EXPECT_EQ(index->x(), i);
        EXPECT_EQ(index->y(), j);
        EXPECT_EQ(index->z(), k);

        EXPECT_EQ(index->xmm().x(), i - 2);
        EXPECT_EQ(index->xmm().y(), j);
        EXPECT_EQ(index->xmm().z(), k);
        ++index;
      }
    }
  }
}

TEST_F(IndexOffsetTest, YMinusTwo) {
  const auto &region = mesh->getRegion3D("RGN_ALL");

  auto index = region.cbegin();

  for (int i = 0; i < nx; ++i) {
    for (int k = 0; k < nz; ++k) {
      ++index;
      ++index;
    }
    for (int j = 2; j < ny; ++j) {
      for (int k = 0; k < nz; ++k) {
        EXPECT_EQ(index->x(), i);
        EXPECT_EQ(index->y(), j);
        EXPECT_EQ(index->z(), k);

        EXPECT_EQ(index->ymm().x(), i);
        EXPECT_EQ(index->ymm().y(), j - 2);
        EXPECT_EQ(index->ymm().z(), k);
        ++index;
      }
    }
  }
}

TEST_F(IndexOffsetTest, ZMinusTwo) {
  const auto &region = mesh->getRegion3D("RGN_ALL");

  auto index = region.cbegin();

  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < ny; ++j) {
      for (int k = 0; k < nz; ++k) {
        EXPECT_EQ(index->x(), i);
        EXPECT_EQ(index->y(), j);
        EXPECT_EQ(index->z(), k);

        EXPECT_EQ(index->zmm().x(), i);
        EXPECT_EQ(index->zmm().y(), j);
        EXPECT_EQ(index->zmm().z(), (k - 2 + nz) % nz);
        ++index;
      }
    }
  }
}

TEST_F(IndexOffsetTest, Offset111) {
  const auto &region = mesh->getRegion3D("RGN_ALL");

  auto index = region.cbegin();

  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < ny; ++j) {
      for (int k = 0; k < nz; ++k) {
        EXPECT_EQ(index->x(), i);
        EXPECT_EQ(index->y(), j);
        EXPECT_EQ(index->z(), k);

        if (i >= (nx - 1) or j >= (ny - 1)) {
          // skip this point
        } else {
          EXPECT_EQ(index->offset(1, 1, 1).x(), i + 1);
          EXPECT_EQ(index->offset(1, 1, 1).y(), j + 1);
          EXPECT_EQ(index->offset(1, 1, 1).z(), (k + 1 + nz) % nz);
        }
        ++index;
      }
    }
  }
}

TEST_F(IndexOffsetTest, Offsetm1m1m1) {
  const auto &region = mesh->getRegion3D("RGN_ALL");

  auto index = region.cbegin();

  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < ny; ++j) {
      for (int k = 0; k < nz; ++k) {
        EXPECT_EQ(index->x(), i);
        EXPECT_EQ(index->y(), j);
        EXPECT_EQ(index->z(), k);

        if (i < 1 or j < 1) {
          // skip this point
        } else {
          EXPECT_EQ(index->offset(-1, -1, -1).x(), i - 1);
          EXPECT_EQ(index->offset(-1, -1, -1).y(), j - 1);
          EXPECT_EQ(index->offset(-1, -1, -1).z(), (k - 1 + nz) % nz);
        }
        ++index;
      }
    }
  }
}

#if CHECK > 2
TEST_F(IndexOffsetTest, ZNegativeOffsetInd3D) {
  const auto &region = mesh->getRegion3D("RGN_ALL");

  auto index = region.cbegin();

  EXPECT_THROW(index->zp(-1), BoutException);
  EXPECT_THROW(index->zm(-1), BoutException);
}
#endif

TEST_F(IndexOffsetTest, ZOffsetZeroInd3D) {
  const auto &region = mesh->getRegion3D("RGN_ALL");

  auto index = region.cbegin();

  EXPECT_EQ(index->zp(0), *index);
  EXPECT_EQ(index->zm(0), *index);
  EXPECT_EQ(index->offset(0, 0, 0), *index);
}

TEST_F(IndexOffsetTest, XInd2D) {
  const auto &region = mesh->getRegion2D("RGN_ALL");

  auto index = region.cbegin();

  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < ny; ++j) {
      EXPECT_EQ(index->x(), i);
      ++index;
    }
  }
}

TEST_F(IndexOffsetTest, YInd2D) {
  const auto &region = mesh->getRegion2D("RGN_ALL");

  auto index = region.cbegin();

  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < ny; ++j) {
      EXPECT_EQ(index->y(), j);
      ++index;
    }
  }
}

TEST_F(IndexOffsetTest, ZInd2D) {
  const auto &region = mesh->getRegion2D("RGN_ALL");

  auto index = region.cbegin();

  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < ny; ++j) {
      EXPECT_EQ(index->z(), 0);
      ++index;
    }
  }
}

TEST_F(IndexOffsetTest, XPlusOneInd2D) {
  const auto &region = mesh->getRegion2D("RGN_ALL");

  auto index = region.cbegin();

  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < ny; ++j) {
      EXPECT_EQ(index->x(), i);
      EXPECT_EQ(index->y(), j);
      EXPECT_EQ(index->z(), 0);

      if (i >= (nx - 1)) {
        // skip this point
      } else {
        EXPECT_EQ(index->xp().x(), i + 1);
        EXPECT_EQ(index->xp().y(), j);
        EXPECT_EQ(index->xp().z(), 0);
      }
      ++index;
    }
  }
}

TEST_F(IndexOffsetTest, YPlusOneInd2D) {
  const auto &region = mesh->getRegion2D("RGN_ALL");

  auto index = region.cbegin();

  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < ny; ++j) {
      EXPECT_EQ(index->x(), i);
      EXPECT_EQ(index->y(), j);
      EXPECT_EQ(index->z(), 0);

      if (j >= (ny - 1)) {
#if CHECK > 3
        EXPECT_THROW(index->yp(), BoutException);
#endif
      } else {
        EXPECT_EQ(index->yp().x(), i);
        EXPECT_EQ(index->yp().y(), j + 1);
        EXPECT_EQ(index->yp().z(), 0);
      }
      ++index;
    }
  }
}

TEST_F(IndexOffsetTest, ZPlusOneInd2D) {
  const auto &region = mesh->getRegion2D("RGN_ALL");

  auto index = region.cbegin();

  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < ny; ++j) {
      EXPECT_EQ(index->x(), i);
      EXPECT_EQ(index->y(), j);
      EXPECT_EQ(index->z(), 0);

      EXPECT_EQ(index->zp().x(), i);
      EXPECT_EQ(index->zp().y(), j);
      EXPECT_EQ(index->zp().z(), 0);
      ++index;
    }
  }
}

TEST_F(IndexOffsetTest, XPlusOneInd2DGeneric) {
  const auto &region = mesh->getRegion2D("RGN_ALL");

  auto index = region.cbegin();

  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < ny; ++j) {
      EXPECT_EQ(index->x(), i);
      EXPECT_EQ(index->y(), j);
      EXPECT_EQ(index->z(), 0);

      if (i >= (nx - 1)) {
        // skip this point
      } else {
        EXPECT_EQ((index->plus<1, DIRECTION::X>().x()), i + 1);
        EXPECT_EQ((index->plus<1, DIRECTION::X>().y()), j);
        EXPECT_EQ((index->plus<1, DIRECTION::X>().z()), 0);
      }
      ++index;
    }
  }
}

TEST_F(IndexOffsetTest, YPlusOneInd2DGeneric) {
  const auto &region = mesh->getRegion2D("RGN_ALL");

  auto index = region.cbegin();

  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < ny; ++j) {
      EXPECT_EQ(index->x(), i);
      EXPECT_EQ(index->y(), j);
      EXPECT_EQ(index->z(), 0);

      if (j >= (ny - 1)) {
#if CHECK > 3
        EXPECT_THROW(index->yp(), BoutException);
#endif
      } else {
        EXPECT_EQ((index->plus<1, DIRECTION::Y>().x()), i);
        EXPECT_EQ((index->plus<1, DIRECTION::Y>().y()), j + 1);
        EXPECT_EQ((index->plus<1, DIRECTION::Y>().z()), 0);
      }
      ++index;
    }
  }
}

TEST_F(IndexOffsetTest, ZPlusOneInd2DGeneric) {
  const auto &region = mesh->getRegion2D("RGN_ALL");

  auto index = region.cbegin();

  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < ny; ++j) {
      EXPECT_EQ(index->x(), i);
      EXPECT_EQ(index->y(), j);
      EXPECT_EQ(index->z(), 0);

      EXPECT_EQ((index->plus<1, DIRECTION::Z>().x()), i);
      EXPECT_EQ((index->plus<1, DIRECTION::Z>().y()), j);
      EXPECT_EQ((index->plus<1, DIRECTION::Z>().z()), 0);
      ++index;
    }
  }
}

TEST_F(IndexOffsetTest, XMinusOneInd2D) {
  const auto &region = mesh->getRegion2D("RGN_ALL");

  auto index = region.cbegin();

  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < ny; ++j) {
      EXPECT_EQ(index->x(), i);
      EXPECT_EQ(index->y(), j);
      EXPECT_EQ(index->z(), 0);

      if (i < 1) {
        // skip this point
      } else {
        EXPECT_EQ(index->xm().x(), i - 1);
        EXPECT_EQ(index->xm().y(), j);
        EXPECT_EQ(index->xm().z(), 0);
      }
      ++index;
    }
  }
}

TEST_F(IndexOffsetTest, YMinusOneInd2D) {
  const auto &region = mesh->getRegion2D("RGN_ALL");

  auto index = region.cbegin();

  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < ny; ++j) {
      EXPECT_EQ(index->x(), i);
      EXPECT_EQ(index->y(), j);
      EXPECT_EQ(index->z(), 0);

      if (j < 1) {
#if CHECK > 3
        EXPECT_THROW(index->ym(), BoutException);
#endif
      } else {
        EXPECT_EQ(index->ym().x(), i);
        EXPECT_EQ(index->ym().y(), j - 1);
        EXPECT_EQ(index->ym().z(), 0);
      }
      ++index;
    }
  }
}

TEST_F(IndexOffsetTest, ZMinusOneInd2D) {
  const auto &region = mesh->getRegion2D("RGN_ALL");

  auto index = region.cbegin();

  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < ny; ++j) {
      EXPECT_EQ(index->x(), i);
      EXPECT_EQ(index->y(), j);
      EXPECT_EQ(index->z(), 0);

      EXPECT_EQ(index->zm().x(), i);
      EXPECT_EQ(index->zm().y(), j);
      EXPECT_EQ(index->zm().z(), 0);
      ++index;
    }
  }
}

TEST_F(IndexOffsetTest, XMinusOneInd2DGeneric) {
  const auto &region = mesh->getRegion2D("RGN_ALL");

  auto index = region.cbegin();

  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < ny; ++j) {
      EXPECT_EQ(index->x(), i);
      EXPECT_EQ(index->y(), j);
      EXPECT_EQ(index->z(), 0);

      if (i < 1) {
        // skip this point
      } else {
        EXPECT_EQ((index->minus<1, DIRECTION::X>().x()), i - 1);
        EXPECT_EQ((index->minus<1, DIRECTION::X>().y()), j);
        EXPECT_EQ((index->minus<1, DIRECTION::X>().z()), 0);
      }
      ++index;
    }
  }
}

TEST_F(IndexOffsetTest, YMinusOneInd2DGeneric) {
  const auto &region = mesh->getRegion2D("RGN_ALL");

  auto index = region.cbegin();

  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < ny; ++j) {
      EXPECT_EQ(index->x(), i);
      EXPECT_EQ(index->y(), j);
      EXPECT_EQ(index->z(), 0);

      if (j < 1) {
#if CHECK > 3
        EXPECT_THROW(index->ym(), BoutException);
#endif
      } else {
        EXPECT_EQ((index->minus<1, DIRECTION::Y>().x()), i);
        EXPECT_EQ((index->minus<1, DIRECTION::Y>().y()), j - 1);
        EXPECT_EQ((index->minus<1, DIRECTION::Y>().z()), 0);
      }
      ++index;
    }
  }
}

TEST_F(IndexOffsetTest, ZMinusOneInd2DGeneric) {
  const auto &region = mesh->getRegion2D("RGN_ALL");

  auto index = region.cbegin();

  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < ny; ++j) {
      EXPECT_EQ(index->x(), i);
      EXPECT_EQ(index->y(), j);
      EXPECT_EQ(index->z(), 0);

      EXPECT_EQ((index->minus<1, DIRECTION::Z>().x()), i);
      EXPECT_EQ((index->minus<1, DIRECTION::Z>().y()), j);
      EXPECT_EQ((index->minus<1, DIRECTION::Z>().z()), 0);
      ++index;
    }
  }
}

TEST_F(IndexOffsetTest, XPlusTwoInd2D) {
  const auto &region = mesh->getRegion2D("RGN_ALL");

  auto index = region.cbegin();

  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < ny; ++j) {
      EXPECT_EQ(index->x(), i);
      EXPECT_EQ(index->y(), j);
      EXPECT_EQ(index->z(), 0);

      if (i >= (nx - 2)) {
        // skip this point
      } else {
        EXPECT_EQ(index->xpp().x(), i + 2);
        EXPECT_EQ(index->xpp().y(), j);
        EXPECT_EQ(index->xpp().z(), 0);
      }
      ++index;
    }
  }
}

TEST_F(IndexOffsetTest, YPlusTwoInd2D) {
  const auto &region = mesh->getRegion2D("RGN_ALL");

  auto index = region.cbegin();

  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < ny; ++j) {
      EXPECT_EQ(index->x(), i);
      EXPECT_EQ(index->y(), j);
      EXPECT_EQ(index->z(), 0);

      if (j >= (ny - 2)) {
#if CHECK > 3
        EXPECT_THROW(index->ypp(), BoutException);
#endif
      } else {
        EXPECT_EQ(index->ypp().x(), i);
        EXPECT_EQ(index->ypp().y(), j + 2);
        EXPECT_EQ(index->ypp().z(), 0);
      }
      ++index;
    }
  }
}

TEST_F(IndexOffsetTest, ZPlusTwoInd2D) {
  const auto &region = mesh->getRegion2D("RGN_ALL");

  auto index = region.cbegin();

  EXPECT_EQ(index->zpp(), *index);
}

#if CHECK > 2
TEST_F(IndexOffsetTest, ZNegativeOffsetInd2D) {
  const auto &region = mesh->getRegion2D("RGN_ALL");

  auto index = region.cbegin();

  EXPECT_THROW(index->zp(-1), BoutException);
  EXPECT_THROW(index->zm(-1), BoutException);
}
#endif

TEST_F(IndexOffsetTest, ZOffsetZeroInd2D) {
  const auto &region = mesh->getRegion2D("RGN_ALL");

  auto index = region.cbegin();

  EXPECT_EQ(index->zp(0), *index);
  EXPECT_EQ(index->zm(0), *index);
  EXPECT_EQ(index->offset(0, 0, 0), *index);
}

TEST_F(IndexOffsetTest, XMinusTwoInd2D) {
  const auto &region = mesh->getRegion2D("RGN_ALL");

  auto index = region.cbegin();

  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < ny; ++j) {
      EXPECT_EQ(index->x(), i);
      EXPECT_EQ(index->y(), j);
      EXPECT_EQ(index->z(), 0);

      if (i < 2) {
        // skip this point
      } else {
        EXPECT_EQ(index->xmm().x(), i - 2);
        EXPECT_EQ(index->xmm().y(), j);
        EXPECT_EQ(index->xmm().z(), 0);
      }
      ++index;
    }
  }
}

TEST_F(IndexOffsetTest, YMinusTwoInd2D) {
  const auto &region = mesh->getRegion2D("RGN_ALL");

  auto index = region.cbegin();

  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < ny; ++j) {
      EXPECT_EQ(index->x(), i);
      EXPECT_EQ(index->y(), j);
      EXPECT_EQ(index->z(), 0);

      if (j < 2) {
#if CHECK > 3
        EXPECT_THROW(index->ymm(), BoutException);
#endif
      } else {
        EXPECT_EQ(index->ymm().x(), i);
        EXPECT_EQ(index->ymm().y(), j - 2);
        EXPECT_EQ(index->ymm().z(), 0);
      }
      ++index;
    }
  }
}

TEST_F(IndexOffsetTest, ZMinusTwoInd2D) {
  const auto &region = mesh->getRegion2D("RGN_ALL");

  auto index = region.cbegin();
  EXPECT_EQ(index->zmm(), *index);
}

TEST_F(IndexOffsetTest, Offset111Ind2D) {
  const auto &region = mesh->getRegion2D("RGN_ALL");

  auto index = region.cbegin();

  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < ny; ++j) {
      EXPECT_EQ(index->x(), i);
      EXPECT_EQ(index->y(), j);
      EXPECT_EQ(index->z(), 0);

      if (i >= (nx - 1) or j >= (ny - 1)) {
        // skip this point
      } else {
        EXPECT_EQ(index->offset(1, 1, 1).x(), i + 1);
        EXPECT_EQ(index->offset(1, 1, 1).y(), j + 1);
        EXPECT_EQ(index->offset(1, 1, 1).z(), 0);
      }
      ++index;
    }
  }
}

TEST_F(IndexOffsetTest, Offsetm1m1m1Ind2D) {
  const auto &region = mesh->getRegion2D("RGN_ALL");

  auto index = region.cbegin();

  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < ny; ++j) {
      EXPECT_EQ(index->x(), i);
      EXPECT_EQ(index->y(), j);
      EXPECT_EQ(index->z(), 0);

      if (i < 1 or j < 1) {
        // skip this point
      } else {
        EXPECT_EQ(index->offset(-1, -1, -1).x(), i - 1);
        EXPECT_EQ(index->offset(-1, -1, -1).y(), j - 1);
        EXPECT_EQ(index->offset(-1, -1, -1).z(), 0);
      }
      ++index;
    }
  }
}
