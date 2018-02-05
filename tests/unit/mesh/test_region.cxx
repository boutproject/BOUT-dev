#include "gtest/gtest.h"

#include "output.hxx"
#include "bout/constants.hxx"
#include "bout/mesh.hxx"
#include "boutexception.hxx"
#include "bout/region.hxx"
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

TEST_F(RegionTest, maxBlockSize) {
  EXPECT_TRUE(MAXREGIONBLOCKSIZE>0);
}

TEST_F(RegionTest, indicesFromRange) {
  Region<ind3D> region(0,mesh->LocalNx-1,0,mesh->LocalNy-1,
		       0,mesh->LocalNz-1,mesh->LocalNy,mesh->LocalNz);

  int nmesh = RegionTest::nx*RegionTest::ny*RegionTest::nz;
  auto regionIndices = region.getIndices();

  //Check index vector is expected length
  EXPECT_EQ(regionIndices.size(), nmesh);
  
  for(int i=0; i < nmesh; i++){
    EXPECT_EQ(regionIndices[i].ind, i);
  }
}

TEST_F(RegionTest, numberOfBlocks) {
  Region<ind3D> region(0,mesh->LocalNx-1,0,mesh->LocalNy-1,
		       0,mesh->LocalNz-1,mesh->LocalNy,mesh->LocalNz);

  auto blocks = region.getBlocks();
  int nmesh = RegionTest::nx*RegionTest::ny*RegionTest::nz;
  int nblocks = blocks.size();

  //Calculate expected number of blocks, assuming region is entirely
  //contiguous, as is the case here.
  int nexpected = nmesh/MAXREGIONBLOCKSIZE; //Integer division required
  //Above is only correct when MAXREGIONBLOCKSIZE is a factor of nmesh,
  //so check if there needs to be an extra block for remainder
  if((1.0*nmesh)/MAXREGIONBLOCKSIZE > 0){
    nexpected++;
  }
  
  EXPECT_EQ(nexpected, nblocks);
}

TEST_F(RegionTest, contiguousBlockSize) {
  Region<ind3D> region(0,mesh->LocalNx-1,0,mesh->LocalNy-1,
		       0,mesh->LocalNz-1,mesh->LocalNy,mesh->LocalNz);

  const int nmesh = RegionTest::nx*RegionTest::ny*RegionTest::nz;
  
  auto blocks = region.getBlocks();
  const int nblocks = blocks.size();
  auto firstBlock = blocks[0];
  auto lastBlock = blocks[nblocks-1];
  
  //Calculate expected size of blocks, assuming region is entirely
  //contiguous, as is the case here.
  const int expectedFirstBlockSize = nmesh >= MAXREGIONBLOCKSIZE ? MAXREGIONBLOCKSIZE : nmesh;
  EXPECT_EQ(1 + firstBlock.second.ind-firstBlock.first.ind,expectedFirstBlockSize);
  const int expectedLastBlockSize = nmesh % MAXREGIONBLOCKSIZE;
  if (expectedLastBlockSize != 0){
    EXPECT_EQ(1 + lastBlock.second.ind-lastBlock.first.ind, expectedLastBlockSize);
  }else{
    //If no remainder then expect same block size as in first block
    EXPECT_EQ(1 + lastBlock.second.ind - lastBlock.first.ind, expectedFirstBlockSize);
  }
}

TEST_F(RegionTest, defaultRegions) {
  const int nmesh = RegionTest::nx*RegionTest::ny*RegionTest::nz;
  EXPECT_EQ(mesh->getRegion("RGN_ALL").getIndices().size(), nmesh);

  const int nnobndry = RegionTest::nz*(mesh->yend-mesh->ystart+1)*(mesh->xend-mesh->xstart+1);
  EXPECT_EQ(mesh->getRegion("RGN_NOBNDRY").getIndices().size(), nnobndry);

  const int nnox = RegionTest::nz*RegionTest::ny*(mesh->xend-mesh->xstart+1);
  EXPECT_EQ(mesh->getRegion("RGN_NOX").getIndices().size(), nnox);

  const int nnoy = RegionTest::nz*RegionTest::nx*(mesh->yend-mesh->ystart+1);
  EXPECT_EQ(mesh->getRegion("RGN_NOY").getIndices().size(), nnoy);
  
}
