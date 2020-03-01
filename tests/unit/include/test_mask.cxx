#include "gtest/gtest.h"

#include "mask.hxx"
#include "boutexception.hxx"
#include "test_extras.hxx"

/// Global mesh
namespace bout{
namespace globals{
extern Mesh *mesh;
} // namespace globals
} // namespace bout

using MaskTest = FakeMeshFixture;

TEST_F(MaskTest, Indexing) {
  constexpr int nx{2};
  constexpr int ny{2};
  constexpr int nz{2};

  BoutMask mask{nx, ny, nz};

  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < ny; ++j) {
      for (int k = 0; k < nz; ++k) {
        EXPECT_FALSE(mask(i, j, k));
      }
    }
  }

  mask(0, 1, 0) = true;
  EXPECT_TRUE(mask(0, 1, 0));
}

TEST_F(MaskTest, ConstIndexing) {
  constexpr int nx{2};
  constexpr int ny{2};
  constexpr int nz{2};

  const BoutMask mask{nx, ny, nz};

  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < ny; ++j) {
      for (int k = 0; k < nz; ++k) {
        EXPECT_FALSE(mask(i, j, k));
      }
    }
  }
}

#if CHECK >= 2
TEST_F(MaskTest, BoundsChecking) {
  constexpr int nx{2};
  constexpr int ny{2};
  constexpr int nz{2};

  BoutMask mask{nx, ny, nz};

  EXPECT_THROW(mask(nx + 1, 0, 0), BoutException);
  EXPECT_THROW(mask(0, ny + 1, 0), BoutException);
  EXPECT_THROW(mask(0, 0, nz + 1), BoutException);
  EXPECT_THROW(mask(-1, 0, 0), BoutException);
  EXPECT_THROW(mask(0, -1, 0), BoutException);
  EXPECT_THROW(mask(0, 0, -1), BoutException);
}

TEST_F(MaskTest, ConstBoundsChecking) {
  constexpr int nx{2};
  constexpr int ny{2};
  constexpr int nz{2};

  const BoutMask mask{nx, ny, nz};

  EXPECT_THROW(mask(nx + 1, 0, 0), BoutException);
  EXPECT_THROW(mask(0, ny + 1, 0), BoutException);
  EXPECT_THROW(mask(0, 0, nz + 1), BoutException);
  EXPECT_THROW(mask(-1, 0, 0), BoutException);
  EXPECT_THROW(mask(0, -1, 0), BoutException);
  EXPECT_THROW(mask(0, 0, -1), BoutException);
}
#endif

TEST_F(MaskTest, Assignment) {
  constexpr int nx{2};
  constexpr int ny{2};
  constexpr int nz{2};

  BoutMask mask{nx, ny, nz};

  mask = true;

  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < ny; ++j) {
      for (int k = 0; k < nz; ++k) {
        EXPECT_TRUE(mask(i, j, k));
      }
    }
  }
}

TEST_F(MaskTest, CreateOnGlobalMesh) {
  BoutMask mask{};

  for (int i = 0; i < MaskTest::nx; ++i) {
    for (int j = 0; j < MaskTest::ny; ++j) {
      for (int k = 0; k < MaskTest::nz; ++k) {
        EXPECT_FALSE(mask(i, j, k));
      }
    }
  }
}

TEST_F(MaskTest, CreateOnMesh) {
  constexpr int nx{2};
  constexpr int ny{3};
  constexpr int nz{4};
  
  FakeMesh mesh{nx, ny, nz};

  BoutMask mask{mesh};

  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < ny; ++j) {
      for (int k = 0; k < nz; ++k) {
        EXPECT_FALSE(mask(i, j, k));
      }
    }
  }
}

TEST_F(MaskTest, CreateOnMeshWithValueTrue) {
  constexpr int nx{2};
  constexpr int ny{3};
  constexpr int nz{4};
  
  FakeMesh mesh{nx, ny, nz};

  BoutMask mask{mesh, true};

  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < ny; ++j) {
      for (int k = 0; k < nz; ++k) {
        EXPECT_TRUE(mask(i, j, k));
      }
    }
  }
}

TEST_F(MaskTest, CreateOnMeshWithValueFalse) {
  constexpr int nx{2};
  constexpr int ny{3};
  constexpr int nz{4};
  
  FakeMesh mesh{nx, ny, nz};

  BoutMask mask{mesh, false};

  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < ny; ++j) {
      for (int k = 0; k < nz; ++k) {
        EXPECT_FALSE(mask(i, j, k));
      }
    }
  }
}

TEST_F(MaskTest, CreateFromIndicesWithValueTrue) {
  constexpr int nx{2};
  constexpr int ny{2};
  constexpr int nz{2};

  BoutMask mask{nx, ny, nz, true};

  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < ny; ++j) {
      for (int k = 0; k < nz; ++k) {
        EXPECT_TRUE(mask(i, j, k));
      }
    }
  }

  mask(0, 1, 0) = false;
  EXPECT_FALSE(mask(0, 1, 0));
}

TEST_F(MaskTest, CreateFromIndicesWithValueFalse) {
  constexpr int nx{2};
  constexpr int ny{2};
  constexpr int nz{2};

  BoutMask mask{nx, ny, nz, false};

  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < ny; ++j) {
      for (int k = 0; k < nz; ++k) {
        EXPECT_FALSE(mask(i, j, k));
      }
    }
  }

  mask(0, 1, 0) = true;
  EXPECT_TRUE(mask(0, 1, 0));
}
