#include "gtest/gtest.h"

#include "bout/indexoffset.hxx"
#include "field3d.hxx"
#include "test_extras.hxx"

#include <iostream>

/// Global mesh
extern Mesh *mesh;

/// Test fixture to make sure the global mesh is our fake one
class IndexOffsetTest : public ::testing::Test {
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

const int IndexOffsetTest::nx = 3;
const int IndexOffsetTest::ny = 5;
const int IndexOffsetTest::nz = 7;


TEST_F(IndexOffsetTest, X) {
  auto region = mesh->getRegion("RGN_ALL");

  IndexOffset<Ind3D> offset(*mesh);

  auto index = std::begin(region);

  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < ny; ++j) {
      for (int k = 0; k < nz; ++k) {
        EXPECT_EQ(offset.x(*index), i);
        ++index;
      }
    }
  }
}

TEST_F(IndexOffsetTest, Y) {
  auto region = mesh->getRegion("RGN_ALL");

  IndexOffset<Ind3D> offset(*mesh);

  auto index = std::begin(region);

  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < ny; ++j) {
      for (int k = 0; k < nz; ++k) {
        EXPECT_EQ(offset.y(*index), j);
        ++index;
      }
    }
  }
}

TEST_F(IndexOffsetTest, Z) {
  auto region = mesh->getRegion("RGN_ALL");

  IndexOffset<Ind3D> offset(*mesh);

  auto index = std::begin(region);

  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < ny; ++j) {
      for (int k = 0; k < nz; ++k) {
        EXPECT_EQ(offset.z(*index), k);
        ++index;
      }
    }
  }
}

TEST_F(IndexOffsetTest, XPlusOne) {
  auto region = mesh->getRegion("RGN_ALL");

  IndexOffset<Ind3D> offset(*mesh);

  auto index = std::begin(region);

  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < ny; ++j) {
      for (int k = 0; k < nz; ++k) {
        EXPECT_EQ(offset.x(*index), i);
        EXPECT_EQ(offset.y(*index), j);
        EXPECT_EQ(offset.z(*index), k);

        EXPECT_EQ(offset.x(offset.xp(*index)), i + 1);
        EXPECT_EQ(offset.y(offset.xp(*index)), j);
        EXPECT_EQ(offset.z(offset.xp(*index)), k);
        ++index;
      }
    }
  }
}

TEST_F(IndexOffsetTest, YPlusOne) {
  auto region = mesh->getRegion("RGN_ALL");

  IndexOffset<Ind3D> offset(*mesh);

  auto index = std::begin(region);

  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < ny - 1; ++j) {
      for (int k = 0; k < nz; ++k) {
        EXPECT_EQ(offset.x(*index), i);
        EXPECT_EQ(offset.y(*index), j);
        EXPECT_EQ(offset.z(*index), k);

        EXPECT_EQ(offset.x(offset.yp(*index)), i);
        EXPECT_EQ(offset.y(offset.yp(*index)), j + 1);
        EXPECT_EQ(offset.z(offset.yp(*index)), k);
        ++index;
      }
    }
    for (int k = 0; k < nz; ++k) {
      ++index;
    }
  }
}

TEST_F(IndexOffsetTest, ZPlusOne) {
  auto region = mesh->getRegion("RGN_ALL");

  IndexOffset<Ind3D> offset(*mesh);

  auto index = std::begin(region);

  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < ny; ++j) {
      for (int k = 0; k < nz - 1; ++k) {
        EXPECT_EQ(offset.x(*index), i);
        EXPECT_EQ(offset.y(*index), j);
        EXPECT_EQ(offset.z(*index), k);

        EXPECT_EQ(offset.x(offset.zp(*index)), i);
        EXPECT_EQ(offset.y(offset.zp(*index)), j);
        EXPECT_EQ(offset.z(offset.zp(*index)), k + 1);
        ++index;
      }
      ++index;
    }
  }
}

TEST_F(IndexOffsetTest, XMinusOne) {
  auto region = mesh->getRegion("RGN_ALL");

  IndexOffset<Ind3D> offset(*mesh);

  auto index = std::begin(region);

  for (int j = 0; j < ny; ++j) {
    for (int k = 0; k < nz; ++k) {
      ++index;
    }
  }
  for (int i = 1; i < nx; ++i) {
    for (int j = 0; j < ny; ++j) {
      for (int k = 0; k < nz; ++k) {
        EXPECT_EQ(offset.x(*index), i);
        EXPECT_EQ(offset.y(*index), j);
        EXPECT_EQ(offset.z(*index), k);

        EXPECT_EQ(offset.x(offset.xm(*index)), i - 1);
        EXPECT_EQ(offset.y(offset.xm(*index)), j);
        EXPECT_EQ(offset.z(offset.xm(*index)), k);
        ++index;
      }
    }
  }
}

TEST_F(IndexOffsetTest, YMinusOne) {
  auto region = mesh->getRegion("RGN_ALL");

  IndexOffset<Ind3D> offset(*mesh);

  auto index = std::begin(region);

  for (int i = 0; i < nx; ++i) {
    for (int k = 0; k < nz; ++k) {
      ++index;
    }
    for (int j = 1; j < ny; ++j) {
      for (int k = 0; k < nz; ++k) {
        EXPECT_EQ(offset.x(*index), i);
        EXPECT_EQ(offset.y(*index), j);
        EXPECT_EQ(offset.z(*index), k);

        EXPECT_EQ(offset.x(offset.ym(*index)), i);
        EXPECT_EQ(offset.y(offset.ym(*index)), j - 1);
        EXPECT_EQ(offset.z(offset.ym(*index)), k);
        ++index;
      }
    }
  }
}

TEST_F(IndexOffsetTest, ZMinusOne) {
  auto region = mesh->getRegion("RGN_ALL");

  IndexOffset<Ind3D> offset(*mesh);

  auto index = std::begin(region);

  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < ny; ++j) {
      ++index;
      for (int k = 1; k < nz; ++k) {
        EXPECT_EQ(offset.x(*index), i);
        EXPECT_EQ(offset.y(*index), j);
        EXPECT_EQ(offset.z(*index), k);

        EXPECT_EQ(offset.x(offset.zm(*index)), i);
        EXPECT_EQ(offset.y(offset.zm(*index)), j);
        EXPECT_EQ(offset.z(offset.zm(*index)), k - 1);
        ++index;
      }
    }
  }
}
