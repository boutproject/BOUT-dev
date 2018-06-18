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
      for (int k = 0; k < nz; ++k) {
        EXPECT_EQ(offset.x(*index), i);
        EXPECT_EQ(offset.y(*index), j);
        EXPECT_EQ(offset.z(*index), k);

        EXPECT_EQ(offset.x(offset.zp(*index)), i);
        EXPECT_EQ(offset.y(offset.zp(*index)), j);
        EXPECT_EQ(offset.z(offset.zp(*index)), (k + 1) % nz);
        ++index;
      }
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
      for (int k = 0; k < nz; ++k) {
        EXPECT_EQ(offset.x(*index), i);
        EXPECT_EQ(offset.y(*index), j);
        EXPECT_EQ(offset.z(*index), k);

        EXPECT_EQ(offset.x(offset.zm(*index)), i);
        EXPECT_EQ(offset.y(offset.zm(*index)), j);
        EXPECT_EQ(offset.z(offset.zm(*index)), (k - 1 + nz) % nz);
        ++index;
      }
    }
  }
}

TEST_F(IndexOffsetTest, XPlusTwo) {
  auto region = mesh->getRegion("RGN_ALL");

  IndexOffset<Ind3D> offset(*mesh);

  auto index = std::begin(region);

  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < ny; ++j) {
      for (int k = 0; k < nz; ++k) {
        EXPECT_EQ(offset.x(*index), i);
        EXPECT_EQ(offset.y(*index), j);
        EXPECT_EQ(offset.z(*index), k);

        EXPECT_EQ(offset.x(offset.xpp(*index)), i + 2);
        EXPECT_EQ(offset.y(offset.xpp(*index)), j);
        EXPECT_EQ(offset.z(offset.xpp(*index)), k);
        ++index;
      }
    }
  }
}

TEST_F(IndexOffsetTest, YPlusTwo) {
  auto region = mesh->getRegion("RGN_ALL");

  IndexOffset<Ind3D> offset(*mesh);

  auto index = std::begin(region);

  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < ny - 2; ++j) {
      for (int k = 0; k < nz; ++k) {
        EXPECT_EQ(offset.x(*index), i);
        EXPECT_EQ(offset.y(*index), j);
        EXPECT_EQ(offset.z(*index), k);

        EXPECT_EQ(offset.x(offset.ypp(*index)), i);
        EXPECT_EQ(offset.y(offset.ypp(*index)), j + 2);
        EXPECT_EQ(offset.z(offset.ypp(*index)), k);
        ++index;
      }
    }
    for (int k = 0; k < nz; ++k) {
      ++index;
      ++index;
    }
  }
}

TEST_F(IndexOffsetTest, ZPlusTwo) {
  auto region = mesh->getRegion("RGN_ALL");

  IndexOffset<Ind3D> offset(*mesh);

  auto index = std::begin(region);

  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < ny; ++j) {
      for (int k = 0; k < nz; ++k) {
        EXPECT_EQ(offset.x(*index), i);
        EXPECT_EQ(offset.y(*index), j);
        EXPECT_EQ(offset.z(*index), k);

        EXPECT_EQ(offset.x(offset.zpp(*index)), i);
        EXPECT_EQ(offset.y(offset.zpp(*index)), j);
        EXPECT_EQ(offset.z(offset.zpp(*index)), (k + 2 + nz) % nz);
        ++index;
      }
    }
  }
}

TEST_F(IndexOffsetTest, XMinusTwo) {
  auto region = mesh->getRegion("RGN_ALL");

  IndexOffset<Ind3D> offset(*mesh);

  auto index = std::begin(region);

  for (int j = 0; j < ny; ++j) {
    for (int k = 0; k < nz; ++k) {
      ++index;
      ++index;
    }
  }
  for (int i = 2; i < nx; ++i) {
    for (int j = 0; j < ny; ++j) {
      for (int k = 0; k < nz; ++k) {
        EXPECT_EQ(offset.x(*index), i);
        EXPECT_EQ(offset.y(*index), j);
        EXPECT_EQ(offset.z(*index), k);

        EXPECT_EQ(offset.x(offset.xmm(*index)), i - 2);
        EXPECT_EQ(offset.y(offset.xmm(*index)), j);
        EXPECT_EQ(offset.z(offset.xmm(*index)), k);
        ++index;
      }
    }
  }
}

TEST_F(IndexOffsetTest, YMinusTwo) {
  auto region = mesh->getRegion("RGN_ALL");

  IndexOffset<Ind3D> offset(*mesh);

  auto index = std::begin(region);

  for (int i = 0; i < nx; ++i) {
    for (int k = 0; k < nz; ++k) {
      ++index;
      ++index;
    }
    for (int j = 2; j < ny; ++j) {
      for (int k = 0; k < nz; ++k) {
        EXPECT_EQ(offset.x(*index), i);
        EXPECT_EQ(offset.y(*index), j);
        EXPECT_EQ(offset.z(*index), k);

        EXPECT_EQ(offset.x(offset.ymm(*index)), i);
        EXPECT_EQ(offset.y(offset.ymm(*index)), j - 2);
        EXPECT_EQ(offset.z(offset.ymm(*index)), k);
        ++index;
      }
    }
  }
}

TEST_F(IndexOffsetTest, ZMinusTwo) {
  auto region = mesh->getRegion("RGN_ALL");

  IndexOffset<Ind3D> offset(*mesh);

  auto index = std::begin(region);

  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < ny; ++j) {
      for (int k = 0; k < nz; ++k) {
        EXPECT_EQ(offset.x(*index), i);
        EXPECT_EQ(offset.y(*index), j);
        EXPECT_EQ(offset.z(*index), k);

        EXPECT_EQ(offset.x(offset.zmm(*index)), i);
        EXPECT_EQ(offset.y(offset.zmm(*index)), j);
        EXPECT_EQ(offset.z(offset.zmm(*index)), (k - 2 + nz) % nz);
        ++index;
      }
    }
  }
}

TEST_F(IndexOffsetTest, XInd2D) {
  auto region = mesh->getRegion2D("RGN_ALL");

  IndexOffset<Ind2D> offset(*mesh);

  auto index = std::begin(region);

  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < ny; ++j) {
      EXPECT_EQ(offset.x(*index), i);
      ++index;
    }
  }
}

TEST_F(IndexOffsetTest, YInd2D) {
  auto region = mesh->getRegion2D("RGN_ALL");

  IndexOffset<Ind2D> offset(*mesh);

  auto index = std::begin(region);

  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < ny; ++j) {
      EXPECT_EQ(offset.y(*index), j);
      ++index;
    }
  }
}

TEST_F(IndexOffsetTest, ZInd2D) {
  auto region = mesh->getRegion2D("RGN_ALL");

  IndexOffset<Ind2D> offset(*mesh);

  auto index = std::begin(region);

  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < ny; ++j) {
      EXPECT_EQ(offset.z(*index), 0);
      ++index;
    }
  }
}

TEST_F(IndexOffsetTest, XPlusOneInd2D) {
  auto region = mesh->getRegion2D("RGN_ALL");

  IndexOffset<Ind2D> offset(*mesh);

  auto index = std::begin(region);

  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < ny; ++j) {
      EXPECT_EQ(offset.x(*index), i);
      EXPECT_EQ(offset.y(*index), j);
      EXPECT_EQ(offset.z(*index), 0);

      EXPECT_EQ(offset.x(offset.xp(*index)), i + 1);
      EXPECT_EQ(offset.y(offset.xp(*index)), j);
      EXPECT_EQ(offset.z(offset.xp(*index)), 0);
      ++index;
    }
  }
}

TEST_F(IndexOffsetTest, YPlusOneInd2D) {
  auto region = mesh->getRegion2D("RGN_ALL");

  IndexOffset<Ind2D> offset(*mesh);

  auto index = std::begin(region);

  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < ny - 1; ++j) {
      EXPECT_EQ(offset.x(*index), i);
      EXPECT_EQ(offset.y(*index), j);
      EXPECT_EQ(offset.z(*index), 0);

      EXPECT_EQ(offset.x(offset.yp(*index)), i);
      EXPECT_EQ(offset.y(offset.yp(*index)), j + 1);
      EXPECT_EQ(offset.z(offset.yp(*index)), 0);
      ++index;
    }
    ++index;
  }
}

TEST_F(IndexOffsetTest, ZPlusOneInd2D) {
  auto region = mesh->getRegion2D("RGN_ALL");

  IndexOffset<Ind2D> offset(*mesh);

  auto index = std::begin(region);

  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < ny; ++j) {
      EXPECT_EQ(offset.x(*index), i);
      EXPECT_EQ(offset.y(*index), j);
      EXPECT_EQ(offset.z(*index), 0);

      EXPECT_EQ(offset.x(offset.zp(*index)), i);
      EXPECT_EQ(offset.y(offset.zp(*index)), j);
      EXPECT_EQ(offset.z(offset.zp(*index)), 0);
      ++index;
    }
  }
}

TEST_F(IndexOffsetTest, XMinusOneInd2D) {
  auto region = mesh->getRegion2D("RGN_ALL");

  IndexOffset<Ind2D> offset(*mesh);

  auto index = std::begin(region);

  for (int j = 0; j < ny; ++j) {
    ++index;
  }
  for (int i = 1; i < nx; ++i) {
    for (int j = 0; j < ny; ++j) {
      EXPECT_EQ(offset.x(*index), i);
      EXPECT_EQ(offset.y(*index), j);
      EXPECT_EQ(offset.z(*index), 0);

      EXPECT_EQ(offset.x(offset.xm(*index)), i - 1);
      EXPECT_EQ(offset.y(offset.xm(*index)), j);
      EXPECT_EQ(offset.z(offset.xm(*index)), 0);
      ++index;
    }
  }
}

TEST_F(IndexOffsetTest, YMinusOneInd2D) {
  auto region = mesh->getRegion2D("RGN_ALL");

  IndexOffset<Ind2D> offset(*mesh);

  auto index = std::begin(region);

  for (int i = 0; i < nx; ++i) {
    ++index;
    for (int j = 1; j < ny; ++j) {
      EXPECT_EQ(offset.x(*index), i);
      EXPECT_EQ(offset.y(*index), j);
      EXPECT_EQ(offset.z(*index), 0);

      EXPECT_EQ(offset.x(offset.ym(*index)), i);
      EXPECT_EQ(offset.y(offset.ym(*index)), j - 1);
      EXPECT_EQ(offset.z(offset.ym(*index)), 0);
      ++index;
    }
  }
}

TEST_F(IndexOffsetTest, ZMinusOneInd2D) {
  auto region = mesh->getRegion2D("RGN_ALL");

  IndexOffset<Ind2D> offset(*mesh);

  auto index = std::begin(region);

  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < ny; ++j) {
      EXPECT_EQ(offset.x(*index), i);
      EXPECT_EQ(offset.y(*index), j);
      EXPECT_EQ(offset.z(*index), 0);

      EXPECT_EQ(offset.x(offset.zm(*index)), i);
      EXPECT_EQ(offset.y(offset.zm(*index)), j);
      EXPECT_EQ(offset.z(offset.zm(*index)), 0);
      ++index;
    }
  }
}

TEST_F(IndexOffsetTest, XPlusTwoInd2D) {
  auto region = mesh->getRegion2D("RGN_ALL");

  IndexOffset<Ind2D> offset(*mesh);

  auto index = std::begin(region);

  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < ny; ++j) {
      EXPECT_EQ(offset.x(*index), i);
      EXPECT_EQ(offset.y(*index), j);
      EXPECT_EQ(offset.z(*index), 0);

      EXPECT_EQ(offset.x(offset.xpp(*index)), i + 2);
      EXPECT_EQ(offset.y(offset.xpp(*index)), j);
      EXPECT_EQ(offset.z(offset.xpp(*index)), 0);
      ++index;
    }
  }
}

TEST_F(IndexOffsetTest, YPlusTwoInd2D) {
  auto region = mesh->getRegion2D("RGN_ALL");

  IndexOffset<Ind2D> offset(*mesh);

  auto index = std::begin(region);

  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < ny - 2; ++j) {
      EXPECT_EQ(offset.x(*index), i);
      EXPECT_EQ(offset.y(*index), j);
      EXPECT_EQ(offset.z(*index), 0);

      EXPECT_EQ(offset.x(offset.ypp(*index)), i);
      EXPECT_EQ(offset.y(offset.ypp(*index)), j + 2);
      EXPECT_EQ(offset.z(offset.ypp(*index)), 0);
      ++index;
    }
    ++index;
    ++index;
  }
}

TEST_F(IndexOffsetTest, ZPlusTwoInd2D) {
  auto region = mesh->getRegion2D("RGN_ALL");

  IndexOffset<Ind2D> offset(*mesh);

  auto index = std::begin(region);

  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < ny; ++j) {
      EXPECT_EQ(offset.x(*index), i);
      EXPECT_EQ(offset.y(*index), j);
      EXPECT_EQ(offset.z(*index), 0);

      EXPECT_EQ(offset.x(offset.zpp(*index)), i);
      EXPECT_EQ(offset.y(offset.zpp(*index)), j);
      EXPECT_EQ(offset.z(offset.zpp(*index)), 0);
      ++index;
    }
  }
}

TEST_F(IndexOffsetTest, XMinusTwoInd2D) {
  auto region = mesh->getRegion2D("RGN_ALL");

  IndexOffset<Ind2D> offset(*mesh);

  auto index = std::begin(region);

  for (int j = 0; j < ny; ++j) {
    ++index;
    ++index;
  }
  for (int i = 2; i < nx; ++i) {
    for (int j = 0; j < ny; ++j) {
      EXPECT_EQ(offset.x(*index), i);
      EXPECT_EQ(offset.y(*index), j);
      EXPECT_EQ(offset.z(*index), 0);

      EXPECT_EQ(offset.x(offset.xmm(*index)), i - 2);
      EXPECT_EQ(offset.y(offset.xmm(*index)), j);
      EXPECT_EQ(offset.z(offset.xmm(*index)), 0);
      ++index;
    }
  }
}

TEST_F(IndexOffsetTest, YMinusTwoInd2D) {
  auto region = mesh->getRegion2D("RGN_ALL");

  IndexOffset<Ind2D> offset(*mesh);

  auto index = std::begin(region);

  for (int i = 0; i < nx; ++i) {
    ++index;
    ++index;
    for (int j = 2; j < ny; ++j) {
      EXPECT_EQ(offset.x(*index), i);
      EXPECT_EQ(offset.y(*index), j);
      EXPECT_EQ(offset.z(*index), 0);

      EXPECT_EQ(offset.x(offset.ymm(*index)), i);
      EXPECT_EQ(offset.y(offset.ymm(*index)), j - 2);
      EXPECT_EQ(offset.z(offset.ymm(*index)), 0);
      ++index;
    }
  }
}

TEST_F(IndexOffsetTest, ZMinusTwoInd2D) {
  auto region = mesh->getRegion2D("RGN_ALL");

  IndexOffset<Ind2D> offset(*mesh);

  auto index = std::begin(region);

  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < ny; ++j) {
      EXPECT_EQ(offset.x(*index), i);
      EXPECT_EQ(offset.y(*index), j);
      EXPECT_EQ(offset.z(*index), 0);

      EXPECT_EQ(offset.x(offset.zmm(*index)), i);
      EXPECT_EQ(offset.y(offset.zmm(*index)), j);
      EXPECT_EQ(offset.z(offset.zmm(*index)), 0);
      ++index;
    }
  }
}
