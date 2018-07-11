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

        if (i >= (nx - 1)) {
#if CHECK > 3
          EXPECT_THROW(offset.xp(*index), BoutException);
#endif
        } else {
          EXPECT_EQ(offset.x(offset.xp(*index)), i + 1);
          EXPECT_EQ(offset.y(offset.xp(*index)), j);
          EXPECT_EQ(offset.z(offset.xp(*index)), k);
        }
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
    for (int j = 0; j < ny; ++j) {
      for (int k = 0; k < nz; ++k) {
        EXPECT_EQ(offset.x(*index), i);
        EXPECT_EQ(offset.y(*index), j);
        EXPECT_EQ(offset.z(*index), k);

        if (j >= (ny - 1)) {
#if CHECK > 3
          EXPECT_THROW(offset.yp(*index), BoutException);
#endif
        } else {
          EXPECT_EQ(offset.x(offset.yp(*index)), i);
          EXPECT_EQ(offset.y(offset.yp(*index)), j + 1);
          EXPECT_EQ(offset.z(offset.yp(*index)), k);
        }
        ++index;
      }
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

  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < ny; ++j) {
      for (int k = 0; k < nz; ++k) {
        EXPECT_EQ(offset.x(*index), i);
        EXPECT_EQ(offset.y(*index), j);
        EXPECT_EQ(offset.z(*index), k);

        if (i < 1) {
#if CHECK > 3
          EXPECT_THROW(offset.xm(*index), BoutException);
#endif
        } else {
          EXPECT_EQ(offset.x(offset.xm(*index)), i - 1);
          EXPECT_EQ(offset.y(offset.xm(*index)), j);
          EXPECT_EQ(offset.z(offset.xm(*index)), k);
        }
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
    for (int j = 0; j < ny; ++j) {
      for (int k = 0; k < nz; ++k) {
        EXPECT_EQ(offset.x(*index), i);
        EXPECT_EQ(offset.y(*index), j);
        EXPECT_EQ(offset.z(*index), k);

        if (j < 1) {
#if CHECK > 3
          EXPECT_THROW(offset.ym(*index), BoutException);
#endif
        } else {
          EXPECT_EQ(offset.x(offset.ym(*index)), i);
          EXPECT_EQ(offset.y(offset.ym(*index)), j - 1);
          EXPECT_EQ(offset.z(offset.ym(*index)), k);
        }
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

        if (i >= (nx - 2)) {
#if CHECK > 3
          EXPECT_THROW(offset.xpp(*index), BoutException);
#endif
        } else {
          EXPECT_EQ(offset.x(offset.xpp(*index)), i + 2);
          EXPECT_EQ(offset.y(offset.xpp(*index)), j);
          EXPECT_EQ(offset.z(offset.xpp(*index)), k);
        }
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
    for (int j = 0; j < ny; ++j) {
      for (int k = 0; k < nz; ++k) {
        EXPECT_EQ(offset.x(*index), i);
        EXPECT_EQ(offset.y(*index), j);
        EXPECT_EQ(offset.z(*index), k);

        if (j >= (ny - 2)) {
#if CHECK > 3
          EXPECT_THROW(offset.ypp(*index), BoutException);
#endif
        } else {
          EXPECT_EQ(offset.x(offset.ypp(*index)), i);
          EXPECT_EQ(offset.y(offset.ypp(*index)), j + 2);
          EXPECT_EQ(offset.z(offset.ypp(*index)), k);
        }
        ++index;
      }
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

TEST_F(IndexOffsetTest, Offset111) {
  auto region = mesh->getRegion("RGN_ALL");

  IndexOffset<Ind3D> offset(*mesh);

  auto index = std::begin(region);

  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < ny; ++j) {
      for (int k = 0; k < nz; ++k) {
        EXPECT_EQ(offset.x(*index), i);
        EXPECT_EQ(offset.y(*index), j);
        EXPECT_EQ(offset.z(*index), k);

        if (i >= (nx - 1) or j >= (ny - 1)) {
#if CHECK > 3
          EXPECT_THROW(offset.offset(*index, 1, 1, 1), BoutException);
#endif
        } else {
          EXPECT_EQ(offset.x(offset.offset(*index, 1, 1, 1)), i + 1);
          EXPECT_EQ(offset.y(offset.offset(*index, 1, 1, 1)), j + 1);
          EXPECT_EQ(offset.z(offset.offset(*index, 1, 1, 1)), (k + 1 + nz) % nz);
        }
        ++index;
      }
    }
  }
}

TEST_F(IndexOffsetTest, Offsetm1m1m1) {
  auto region = mesh->getRegion("RGN_ALL");

  IndexOffset<Ind3D> offset(*mesh);

  auto index = std::begin(region);

  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < ny; ++j) {
      for (int k = 0; k < nz; ++k) {
        EXPECT_EQ(offset.x(*index), i);
        EXPECT_EQ(offset.y(*index), j);
        EXPECT_EQ(offset.z(*index), k);

        if (i < 1 or j < 1) {
#if CHECK > 3
          EXPECT_THROW(offset.offset(*index, -1, -1, -1), BoutException);
#endif
        } else {
          EXPECT_EQ(offset.x(offset.offset(*index, -1, -1, -1)), i - 1);
          EXPECT_EQ(offset.y(offset.offset(*index, -1, -1, -1)), j - 1);
          EXPECT_EQ(offset.z(offset.offset(*index, -1, -1, -1)), (k - 1 + nz) % nz);
        }
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

      if (i >= (nx - 1)) {
#if CHECK > 3
        EXPECT_THROW(offset.xp(*index), BoutException);
#endif
      } else {
        EXPECT_EQ(offset.x(offset.xp(*index)), i + 1);
        EXPECT_EQ(offset.y(offset.xp(*index)), j);
        EXPECT_EQ(offset.z(offset.xp(*index)), 0);
      }
      ++index;
    }
  }
}

TEST_F(IndexOffsetTest, YPlusOneInd2D) {
  auto region = mesh->getRegion2D("RGN_ALL");

  IndexOffset<Ind2D> offset(*mesh);

  auto index = std::begin(region);

  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < ny; ++j) {
      EXPECT_EQ(offset.x(*index), i);
      EXPECT_EQ(offset.y(*index), j);
      EXPECT_EQ(offset.z(*index), 0);

      if (j >= (ny - 1)) {
#if CHECK > 3
        EXPECT_THROW(offset.yp(*index), BoutException);
#endif
      } else {
        EXPECT_EQ(offset.x(offset.yp(*index)), i);
        EXPECT_EQ(offset.y(offset.yp(*index)), j + 1);
        EXPECT_EQ(offset.z(offset.yp(*index)), 0);
      }
      ++index;
    }
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

  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < ny; ++j) {
      EXPECT_EQ(offset.x(*index), i);
      EXPECT_EQ(offset.y(*index), j);
      EXPECT_EQ(offset.z(*index), 0);

      if (i < 1) {
#if CHECK > 3
        EXPECT_THROW(offset.xm(*index), BoutException);
#endif
      } else {
        EXPECT_EQ(offset.x(offset.xm(*index)), i - 1);
        EXPECT_EQ(offset.y(offset.xm(*index)), j);
        EXPECT_EQ(offset.z(offset.xm(*index)), 0);
      }
      ++index;
    }
  }
}

TEST_F(IndexOffsetTest, YMinusOneInd2D) {
  auto region = mesh->getRegion2D("RGN_ALL");

  IndexOffset<Ind2D> offset(*mesh);

  auto index = std::begin(region);

  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < ny; ++j) {
      EXPECT_EQ(offset.x(*index), i);
      EXPECT_EQ(offset.y(*index), j);
      EXPECT_EQ(offset.z(*index), 0);

      if (j < 1) {
#if CHECK > 3
        EXPECT_THROW(offset.ym(*index), BoutException);
#endif
      } else {
        EXPECT_EQ(offset.x(offset.ym(*index)), i);
        EXPECT_EQ(offset.y(offset.ym(*index)), j - 1);
        EXPECT_EQ(offset.z(offset.ym(*index)), 0);
      }
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

      if (i >= (nx - 2)) {
#if CHECK > 3
        EXPECT_THROW(offset.xpp(*index), BoutException);
#endif
      } else {
        EXPECT_EQ(offset.x(offset.xpp(*index)), i + 2);
        EXPECT_EQ(offset.y(offset.xpp(*index)), j);
        EXPECT_EQ(offset.z(offset.xpp(*index)), 0);
      }
      ++index;
    }
  }
}

TEST_F(IndexOffsetTest, YPlusTwoInd2D) {
  auto region = mesh->getRegion2D("RGN_ALL");

  IndexOffset<Ind2D> offset(*mesh);

  auto index = std::begin(region);

  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < ny; ++j) {
      EXPECT_EQ(offset.x(*index), i);
      EXPECT_EQ(offset.y(*index), j);
      EXPECT_EQ(offset.z(*index), 0);

      if (j >= (ny - 2)) {
#if CHECK > 3
        EXPECT_THROW(offset.ypp(*index), BoutException);
#endif
      } else {
        EXPECT_EQ(offset.x(offset.ypp(*index)), i);
        EXPECT_EQ(offset.y(offset.ypp(*index)), j + 2);
        EXPECT_EQ(offset.z(offset.ypp(*index)), 0);
      }
      ++index;
    }
  }
}

#if CHECK > 2
TEST_F(IndexOffsetTest, ZPlusTwoInd2D) {
  auto region = mesh->getRegion2D("RGN_ALL");

  IndexOffset<Ind2D> offset(*mesh);

  auto index = std::begin(region);

  EXPECT_THROW(offset.zpp(*index), BoutException);
}
#endif

TEST_F(IndexOffsetTest, XMinusTwoInd2D) {
  auto region = mesh->getRegion2D("RGN_ALL");

  IndexOffset<Ind2D> offset(*mesh);

  auto index = std::begin(region);

  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < ny; ++j) {
      EXPECT_EQ(offset.x(*index), i);
      EXPECT_EQ(offset.y(*index), j);
      EXPECT_EQ(offset.z(*index), 0);

      if (i < 2) {
#if CHECK > 3
        EXPECT_THROW(offset.xmm(*index), BoutException);
#endif
      } else {
        EXPECT_EQ(offset.x(offset.xmm(*index)), i - 2);
        EXPECT_EQ(offset.y(offset.xmm(*index)), j);
        EXPECT_EQ(offset.z(offset.xmm(*index)), 0);
      }
      ++index;
    }
  }
}

TEST_F(IndexOffsetTest, YMinusTwoInd2D) {
  auto region = mesh->getRegion2D("RGN_ALL");

  IndexOffset<Ind2D> offset(*mesh);

  auto index = std::begin(region);

  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < ny; ++j) {
      EXPECT_EQ(offset.x(*index), i);
      EXPECT_EQ(offset.y(*index), j);
      EXPECT_EQ(offset.z(*index), 0);

      if (j < 2) {
#if CHECK > 3
        EXPECT_THROW(offset.ymm(*index), BoutException);
#endif
      } else {
        EXPECT_EQ(offset.x(offset.ymm(*index)), i);
        EXPECT_EQ(offset.y(offset.ymm(*index)), j - 2);
        EXPECT_EQ(offset.z(offset.ymm(*index)), 0);
      }
      ++index;
    }
  }
}

#if CHECK > 2
TEST_F(IndexOffsetTest, ZMinusTwoInd2D) {
  auto region = mesh->getRegion2D("RGN_ALL");

  IndexOffset<Ind2D> offset(*mesh);

  auto index = std::begin(region);

  EXPECT_THROW(offset.zmm(*index), BoutException);
}
#endif

TEST_F(IndexOffsetTest, Offset111Ind2D) {
  auto region = mesh->getRegion2D("RGN_ALL");

  IndexOffset<Ind2D> offset(*mesh);

  auto index = std::begin(region);

  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < ny; ++j) {
      EXPECT_EQ(offset.x(*index), i);
      EXPECT_EQ(offset.y(*index), j);
      EXPECT_EQ(offset.z(*index), 0);

      if (i >= (nx - 1) or j >= (ny - 1)) {
#if CHECK > 3
        EXPECT_THROW(offset.offset(*index, 1, 1, 1), BoutException);
#endif
      } else {
        EXPECT_EQ(offset.x(offset.offset(*index, 1, 1, 1)), i + 1);
        EXPECT_EQ(offset.y(offset.offset(*index, 1, 1, 1)), j + 1);
        EXPECT_EQ(offset.z(offset.offset(*index, 1, 1, 1)), 0);
      }
      ++index;
    }
  }
}

TEST_F(IndexOffsetTest, Offsetm1m1m1Ind2D) {
  auto region = mesh->getRegion2D("RGN_ALL");

  IndexOffset<Ind2D> offset(*mesh);

  auto index = std::begin(region);

  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < ny; ++j) {
      EXPECT_EQ(offset.x(*index), i);
      EXPECT_EQ(offset.y(*index), j);
      EXPECT_EQ(offset.z(*index), 0);

      if (i < 1 or j < 1) {
#if CHECK > 3
        EXPECT_THROW(offset.offset(*index, -1, -1, -1), BoutException);
#endif
      } else {
        EXPECT_EQ(offset.x(offset.offset(*index, -1, -1, -1)), i - 1);
        EXPECT_EQ(offset.y(offset.offset(*index, -1, -1, -1)), j - 1);
        EXPECT_EQ(offset.z(offset.offset(*index, -1, -1, -1)), 0);
      }
      ++index;
    }
  }
}

TEST_F(IndexOffsetTest, XIndPerp) {
  auto region = mesh->getRegionPerp("RGN_ALL");

  IndexOffset<IndPerp> offset(*mesh);

  auto index = std::begin(region);

  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < nz; ++j) {
      EXPECT_EQ(offset.x(*index), i);
      ++index;
    }
  }
}

TEST_F(IndexOffsetTest, YIndPerp) {
  auto region = mesh->getRegionPerp("RGN_ALL");

  IndexOffset<IndPerp> offset(*mesh);

  auto index = std::begin(region);

  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < nz; ++j) {
      EXPECT_EQ(offset.y(*index), 0);
      ++index;
    }
  }
}

TEST_F(IndexOffsetTest, ZIndPerp) {
  auto region = mesh->getRegionPerp("RGN_ALL");

  IndexOffset<IndPerp> offset(*mesh);

  auto index = std::begin(region);

  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < nz; ++j) {
      EXPECT_EQ(offset.z(*index), j);
      ++index;
    }
  }
}

TEST_F(IndexOffsetTest, XPlusOneIndPerp) {
  auto region = mesh->getRegionPerp("RGN_ALL");

  IndexOffset<IndPerp> offset(*mesh);

  auto index = std::begin(region);

  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < nz; ++j) {
      EXPECT_EQ(offset.x(*index), i);
      EXPECT_EQ(offset.y(*index), 0);
      EXPECT_EQ(offset.z(*index), j);

      if (i >= (nx - 1)) {
#if CHECK > 3
        EXPECT_THROW(offset.xp(*index), BoutException);
#endif
      } else {
        EXPECT_EQ(offset.x(offset.xp(*index)), i + 1);
        EXPECT_EQ(offset.y(offset.xp(*index)), 0);
        EXPECT_EQ(offset.z(offset.xp(*index)), j);
      }
      ++index;
    }
  }
}

#if CHECK > 2
TEST_F(IndexOffsetTest, YPlusOneIndPerp) {
  auto region = mesh->getRegionPerp("RGN_ALL");

  IndexOffset<IndPerp> offset(*mesh);

  auto index = std::begin(region);

  EXPECT_THROW(offset.yp(*index), BoutException);
}
#endif

TEST_F(IndexOffsetTest, ZPlusOneIndPerp) {
  auto region = mesh->getRegionPerp("RGN_ALL");

  IndexOffset<IndPerp> offset(*mesh);

  auto index = std::begin(region);

  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < nz; ++j) {
      EXPECT_EQ(offset.x(*index), i);
      EXPECT_EQ(offset.y(*index), 0);
      EXPECT_EQ(offset.z(*index), j);

      EXPECT_EQ(offset.x(offset.zp(*index)), i);
      EXPECT_EQ(offset.y(offset.zp(*index)), 0);
      EXPECT_EQ(offset.z(offset.zp(*index)), (j + 1) % nz);

      ++index;
    }
  }
}

TEST_F(IndexOffsetTest, XMinusOneIndPerp) {
  auto region = mesh->getRegionPerp("RGN_ALL");

  IndexOffset<IndPerp> offset(*mesh);

  auto index = std::begin(region);

  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < nz; ++j) {
      EXPECT_EQ(offset.x(*index), i);
      EXPECT_EQ(offset.y(*index), 0);
      EXPECT_EQ(offset.z(*index), j);

      if (i < 1) {
#if CHECK > 3
        EXPECT_THROW(offset.xm(*index), BoutException);
#endif
      } else {
        EXPECT_EQ(offset.x(offset.xm(*index)), i - 1);
        EXPECT_EQ(offset.y(offset.xm(*index)), 0);
        EXPECT_EQ(offset.z(offset.xm(*index)), j);
      }
      ++index;
    }
  }
}

#if CHECK > 2
TEST_F(IndexOffsetTest, YMinusOneIndPerp) {
  auto region = mesh->getRegionPerp("RGN_ALL");

  IndexOffset<IndPerp> offset(*mesh);

  auto index = std::begin(region);

  EXPECT_THROW(offset.ym(*index), BoutException);
}
#endif

TEST_F(IndexOffsetTest, ZMinusOneIndPerp) {
  auto region = mesh->getRegionPerp("RGN_ALL");

  IndexOffset<IndPerp> offset(*mesh);

  auto index = std::begin(region);

  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < nz; ++j) {
      EXPECT_EQ(offset.x(*index), i);
      EXPECT_EQ(offset.y(*index), 0);
      EXPECT_EQ(offset.z(*index), j);

      EXPECT_EQ(offset.x(offset.zm(*index)), i);
      EXPECT_EQ(offset.y(offset.zm(*index)), 0);
      EXPECT_EQ(offset.z(offset.zm(*index)), (j - 1 + nz) % nz);
      ++index;
    }
  }
}

TEST_F(IndexOffsetTest, XPlusTwoIndPerp) {
  auto region = mesh->getRegionPerp("RGN_ALL");

  IndexOffset<IndPerp> offset(*mesh);

  auto index = std::begin(region);

  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < nz; ++j) {
      EXPECT_EQ(offset.x(*index), i);
      EXPECT_EQ(offset.y(*index), 0);
      EXPECT_EQ(offset.z(*index), j);

      if (i >= (nx - 2)) {
#if CHECK > 3
        EXPECT_THROW(offset.xpp(*index), BoutException);
#endif
      } else {
        EXPECT_EQ(offset.x(offset.xpp(*index)), i + 2);
        EXPECT_EQ(offset.y(offset.xpp(*index)), 0);
        EXPECT_EQ(offset.z(offset.xpp(*index)), j);
      }
      ++index;
    }
  }
}

#if CHECK > 2
TEST_F(IndexOffsetTest, YPlusTwoIndPerp) {
  auto region = mesh->getRegionPerp("RGN_ALL");

  IndexOffset<IndPerp> offset(*mesh);

  auto index = std::begin(region);

  EXPECT_THROW(offset.ypp(*index), BoutException);
}
#endif

TEST_F(IndexOffsetTest, ZPlusTwoIndPerp) {
  auto region = mesh->getRegionPerp("RGN_ALL");

  IndexOffset<IndPerp> offset(*mesh);

  auto index = std::begin(region);

  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < nz; ++j) {
      EXPECT_EQ(offset.x(*index), i);
      EXPECT_EQ(offset.y(*index), 0);
      EXPECT_EQ(offset.z(*index), j);

      EXPECT_EQ(offset.x(offset.zpp(*index)), i);
      EXPECT_EQ(offset.y(offset.zpp(*index)), 0);
      EXPECT_EQ(offset.z(offset.zpp(*index)), (j + 2 + nz) % nz);
      ++index;
    }
  }
}

TEST_F(IndexOffsetTest, XMinusTwoIndPerp) {
  auto region = mesh->getRegionPerp("RGN_ALL");

  IndexOffset<IndPerp> offset(*mesh);

  auto index = std::begin(region);

  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < nz; ++j) {
      EXPECT_EQ(offset.x(*index), i);
      EXPECT_EQ(offset.y(*index), 0);
      EXPECT_EQ(offset.z(*index), j);

      if (i < 2) {
#if CHECK > 3
        EXPECT_THROW(offset.xmm(*index), BoutException);
#endif
      } else {
        EXPECT_EQ(offset.x(offset.xmm(*index)), i - 2);
        EXPECT_EQ(offset.y(offset.xmm(*index)), 0);
        EXPECT_EQ(offset.z(offset.xmm(*index)), j);
      }
      ++index;
    }
  }
}

#if CHECK > 2
TEST_F(IndexOffsetTest, YMinusTwoIndPerp) {
  auto region = mesh->getRegionPerp("RGN_ALL");

  IndexOffset<IndPerp> offset(*mesh);

  auto index = std::begin(region);

  EXPECT_THROW(offset.ymm(*index), BoutException);
}
#endif

TEST_F(IndexOffsetTest, ZMinusTwoIndPerp) {
  auto region = mesh->getRegionPerp("RGN_ALL");

  IndexOffset<IndPerp> offset(*mesh);

  auto index = std::begin(region);

  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < nz; ++j) {

      EXPECT_EQ(offset.x(*index), i);
      EXPECT_EQ(offset.y(*index), 0);
      EXPECT_EQ(offset.z(*index), j);

      EXPECT_EQ(offset.x(offset.zmm(*index)), i);
      EXPECT_EQ(offset.y(offset.zmm(*index)), 0);
      EXPECT_EQ(offset.z(offset.zmm(*index)), (j - 2 + nz) % nz);

      ++index;
    }
  }
}

#if CHECK > 2
TEST_F(IndexOffsetTest, Offset111IndPerp) {
  auto region = mesh->getRegionPerp("RGN_ALL");

  IndexOffset<IndPerp> offset(*mesh);

  auto index = std::begin(region);

  EXPECT_THROW(offset.offset(*index, 1, 1, 1), BoutException);
}
#endif

#if CHECK > 2
TEST_F(IndexOffsetTest, Offsetm1m1m1IndPerp) {
  auto region = mesh->getRegionPerp("RGN_ALL");

  IndexOffset<IndPerp> offset(*mesh);

  auto index = std::begin(region);

  EXPECT_THROW(offset.offset(*index, -1, -1, -1), BoutException);
}
#endif
