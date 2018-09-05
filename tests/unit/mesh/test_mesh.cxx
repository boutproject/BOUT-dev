#include "gtest/gtest.h"

#include "bout/mesh.hxx"
#include "bout/region.hxx"
#include "boutexception.hxx"
#include "output.hxx"

#include "test_extras.hxx"

/// Test fixture to make sure the global mesh is our fake one
class MeshTest : public ::testing::Test {
protected:
  static void SetUpTestCase() { output_info.disable(); }

  static void TearDownTestCase() { output_info.enable(); }

public:
  MeshTest() : localmesh(nx, ny, nz) {}
  static const int nx = 3;
  static const int ny = 5;
  static const int nz = 7;
  FakeMesh localmesh;
};

TEST_F(MeshTest, CreateDefaultRegions) {
  EXPECT_NO_THROW(localmesh.createDefaultRegions());
  EXPECT_THROW(localmesh.createDefaultRegions(), BoutException);
}

TEST_F(MeshTest, GetRegionFromMesh) {
  localmesh.createDefaultRegions();
  EXPECT_NO_THROW(localmesh.getRegion("RGN_ALL"));
  EXPECT_NO_THROW(localmesh.getRegion("RGN_NOBNDRY"));
  EXPECT_NO_THROW(localmesh.getRegion("RGN_NOX"));
  EXPECT_NO_THROW(localmesh.getRegion("RGN_NOY"));
  EXPECT_THROW(localmesh.getRegion("SOME_MADE_UP_REGION_NAME"), BoutException);
}

TEST_F(MeshTest, GetRegion3DFromMesh) {
  localmesh.createDefaultRegions();
  EXPECT_NO_THROW(localmesh.getRegion3D("RGN_ALL"));
  EXPECT_NO_THROW(localmesh.getRegion3D("RGN_NOBNDRY"));
  EXPECT_THROW(localmesh.getRegion3D("SOME_MADE_UP_REGION_NAME"), BoutException);
}

TEST_F(MeshTest, GetRegion2DFromMesh) {
  localmesh.createDefaultRegions();
  EXPECT_NO_THROW(localmesh.getRegion2D("RGN_ALL"));
  EXPECT_NO_THROW(localmesh.getRegion2D("RGN_NOBNDRY"));
  EXPECT_THROW(localmesh.getRegion2D("SOME_MADE_UP_REGION_NAME"), BoutException);
}

TEST_F(MeshTest, GetRegionPerpFromMesh) {
  localmesh.createDefaultRegions();
  EXPECT_NO_THROW(localmesh.getRegionPerp("RGN_ALL"));
  EXPECT_NO_THROW(localmesh.getRegionPerp("RGN_NOBNDRY"));
  EXPECT_THROW(localmesh.getRegionPerp("SOME_MADE_UP_REGION_NAME"), BoutException);
}

TEST_F(MeshTest, AddRegionToMesh) {
  Region<Ind3D> junk(0, 0, 0, 0, 0, 0, 1, 1);
  EXPECT_NO_THROW(localmesh.addRegion("RGN_JUNK", junk));
  EXPECT_THROW(localmesh.addRegion("RGN_JUNK", junk), BoutException);
  Region<Ind2D> junk2D(0, 0, 0, 0, 0, 0, 1, 1);
  EXPECT_NO_THROW(localmesh.addRegion("RGN_JUNK", junk2D));
  Region<IndPerp> junkPerp(0, 0, 0, 0, 0, 0, 1, 1);
  EXPECT_NO_THROW(localmesh.addRegion("RGN_JUNK", junkPerp));
}

TEST_F(MeshTest, AddRegion3DToMesh) {
  Region<Ind3D> junk(0, 0, 0, 0, 0, 0, 1, 1);
  EXPECT_NO_THROW(localmesh.addRegion3D("RGN_JUNK_3D", junk));
  EXPECT_THROW(localmesh.addRegion3D("RGN_JUNK_3D", junk), BoutException);
}

TEST_F(MeshTest, AddRegion2DToMesh) {
  Region<Ind2D> junk(0, 0, 0, 0, 0, 0, 1, 1);
  EXPECT_NO_THROW(localmesh.addRegion2D("RGN_JUNK_2D", junk));
  EXPECT_THROW(localmesh.addRegion2D("RGN_JUNK_2D", junk), BoutException);
}

TEST_F(MeshTest, AddRegionPerpToMesh) {
  Region<IndPerp> junk(0, 0, 0, 0, 0, 0, 1, 1);
  EXPECT_NO_THROW(localmesh.addRegionPerp("RGN_JUNK_Perp", junk));
  EXPECT_THROW(localmesh.addRegionPerp("RGN_JUNK_Perp", junk), BoutException);
}

TEST_F(MeshTest, Ind2DTo3D) {
  Ind2D index2d_0(0);
  Ind2D index2d_7(7);
  Ind2D index2d_14(14);

  EXPECT_EQ(localmesh.ind2Dto3D(index2d_0, 0), Ind3D(0));
  EXPECT_EQ(localmesh.ind2Dto3D(index2d_0, 1), Ind3D(1));

  EXPECT_EQ(localmesh.ind2Dto3D(index2d_7, 0), Ind3D(49));
  EXPECT_EQ(localmesh.ind2Dto3D(index2d_7, 1), Ind3D(50));

  EXPECT_EQ(localmesh.ind2Dto3D(index2d_14, 0), Ind3D(98));
  EXPECT_EQ(localmesh.ind2Dto3D(index2d_14, 1), Ind3D(99));
}

TEST_F(MeshTest, Ind3DTo2D) {
  Ind3D index3d_0(0);
  Ind3D index3d_49(49);
  Ind3D index3d_98(98);

  EXPECT_EQ(localmesh.ind3Dto2D(index3d_0), Ind2D(0));
  EXPECT_EQ(localmesh.ind3Dto2D(index3d_49), Ind2D(7));
  EXPECT_EQ(localmesh.ind3Dto2D(index3d_98), Ind2D(14));
}

TEST_F(MeshTest, MapInd3DTo2D) {
  Ind3D index3d_0(0);
  Ind3D index3d_49(49);
  Ind3D index3d_98(98);

  EXPECT_EQ(localmesh.ind3Dto2D(index3d_0).ind, 0);
  EXPECT_EQ(localmesh.ind3Dto2D(index3d_49).ind, 7);
  EXPECT_EQ(localmesh.ind3Dto2D(index3d_98).ind, 14);
}

TEST_F(MeshTest, IndPerpTo3D) {
  std::vector<int> globalInds = {0, 1, 49, 50, 98, 99};

  for (const auto &i : globalInds) {
    const auto tmp3D = Ind3D(i, ny, nz);
    const auto tmpPerp = IndPerp(tmp3D.x() * nz + tmp3D.z(), 1, nz);

    EXPECT_EQ(tmp3D.x(), tmpPerp.x());
    EXPECT_EQ(tmp3D.z(), tmpPerp.z());

    const auto converted = localmesh.indPerpto3D(tmpPerp, tmp3D.y());

    EXPECT_EQ(tmp3D.ind, converted.ind);
  }
}

TEST_F(MeshTest, Ind3DToPerp) {
  std::vector<int> perpInds{0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10,
                            11, 12, 13, 14, 15, 16, 17, 18, 19, 20};
  std::vector<int> jyVals{0, 1, 2, 3, 4};

  for (const auto &jy : jyVals) {
    for (const auto &i : perpInds) {
      const auto tmpPerp = IndPerp(i, 1, nz);
      const auto tmp3D = Ind3D(tmpPerp.z() + nz * (jy + ny * tmpPerp.x()), ny, nz);

      EXPECT_EQ(tmp3D.x(), tmpPerp.x());
      EXPECT_EQ(tmp3D.y(), jy);
      EXPECT_EQ(tmp3D.z(), tmpPerp.z());

      const auto converted = localmesh.ind3DtoPerp(tmp3D);

      EXPECT_EQ(tmpPerp.ind, converted.ind);
    }
  }
}

extern Mesh::deriv_func fDDX, fDDY, fDDZ;       ///< Differencing methods for each dimension
extern Mesh::deriv_func fD2DX2, fD2DY2, fD2DZ2; ///< second differential operators
extern Mesh::upwind_func fVDDX, fVDDY, fVDDZ;   ///< Upwind functions in the three directions
extern Mesh::flux_func fFDDX, fFDDY, fFDDZ;     ///< Default flux functions
extern Mesh::deriv_func sfDDX, sfDDY, sfDDZ;
extern Mesh::deriv_func sfD2DX2, sfD2DY2, sfD2DZ2;
extern Mesh::flux_func sfVDDX, sfVDDY, sfVDDZ;
extern Mesh::flux_func sfFDDX, sfFDDY, sfFDDZ;

#define RESET                                   \
  fDDX=nullptr; fDDY=nullptr; fDDZ=nullptr;     \
  fD2DX2=nullptr; fD2DY2=nullptr; fD2DZ2=nullptr;       \
  fVDDX=nullptr; fVDDY=nullptr; fVDDZ=nullptr;          \
  fFDDX=nullptr; fFDDY=nullptr; fFDDZ=nullptr;          \
  sfDDX=nullptr; sfDDY=nullptr; sfDDZ=nullptr;          \
  sfD2DX2=nullptr; sfD2DY2=nullptr; sfD2DZ2=nullptr;    \
  sfVDDX=nullptr; sfVDDY=nullptr; sfVDDZ=nullptr;       \
  sfFDDX=nullptr; sfFDDY=nullptr; sfFDDZ=nullptr;


BoutReal DDX_C2(stencil &f);
BoutReal DDX_C4(stencil &f);
BoutReal DDX_CWENO2(stencil &f);
BoutReal DDX_S2(stencil &f);
BoutReal VDDX_C2(BoutReal vc, stencil &f);
BoutReal VDDX_C4(BoutReal vc, stencil &f);
BoutReal VDDX_U1(BoutReal vc, stencil &f);
BoutReal VDDX_U2(BoutReal vc, stencil &f);
BoutReal VDDX_U3(BoutReal vc, stencil &f);
BoutReal VDDX_WENO3(BoutReal vc, stencil &f);
BoutReal DDX_CWENO3(stencil &f);
BoutReal FDDX_U1(stencil &v, stencil &f);
BoutReal FDDX_C2(stencil &v, stencil &f);
BoutReal FDDX_C4(stencil &v, stencil &f);
BoutReal DDX_C2_stag(stencil &f);
BoutReal DDX_C4_stag(stencil &f);
BoutReal VDDX_U1_stag(stencil &v, stencil &f);
BoutReal VDDX_U2_stag(stencil &v, stencil &f);
BoutReal VDDX_C2_stag(stencil &v, stencil &f);
BoutReal VDDX_C4_stag(stencil &v, stencil &f);
BoutReal FDDX_U1_stag(stencil &v, stencil &f);

TEST_F(MeshTest, SetDerivativesDefault) {
  RESET
  Options opt;
  localmesh.initDerivs(&opt);
  EXPECT_EQ(fDDX,&DDX_C2);
  EXPECT_EQ(sfDDX,&DDX_C2_stag);
}

TEST_F(MeshTest, SetDerivativesDiff) {
  RESET
  Options opt;
  opt.getSection("diff")->set("first","C4","test");
  localmesh.initDerivs(&opt);
  EXPECT_EQ(fDDY,&DDX_C4);
  EXPECT_EQ(sfDDY,&DDX_C4_stag);
}

TEST_F(MeshTest, SetDerivativesDiffStag) {
  RESET
  Options opt;
  opt.getSection("diff")->set("first","C2","test");
  opt.getSection("diff")->set("firstStag","C4","test");
  localmesh.initDerivs(&opt);
  EXPECT_EQ(fDDZ,&DDX_C2);
  EXPECT_EQ(sfDDZ,&DDX_C4_stag);
}

TEST_F(MeshTest, SetDerivativesDdxBeforeDiff) {
  RESET
  Options opt;
  opt.getSection("diff")->set("firstStag","C4","test");
  opt.getSection("ddx")->set("firstStag","C2","test");
  localmesh.initDerivs(&opt);
  EXPECT_EQ(sfDDZ,&DDX_C4_stag);
  EXPECT_EQ(sfDDX,&DDX_C2_stag);
}

TEST_F(MeshTest, SetDerivativesDiffStagBeforeDDXNone) {
  RESET
  Options opt;
  opt.getSection("diff")->set("firstStag","C4","test");
  opt.getSection("ddx")->set("first","C2","test");
  localmesh.initDerivs(&opt);
  EXPECT_EQ(sfDDX,&DDX_C4_stag);
}

TEST_F(MeshTest, SetDerivativesInvalid) {
  RESET
  Options opt;
  // An invalid but unused option is fine
  opt.getSection("diff")->set("firstStag","XXX","test");
  opt.getSection("ddx")->set("firstStag","C2","test");
  opt.getSection("ddy")->set("firstStag","C2","test");
  opt.getSection("ddz")->set("firstStag","C2","test");
  localmesh.initDerivs(&opt);
  EXPECT_EQ(sfDDX,&DDX_C2_stag);
}

TEST_F(MeshTest, SetDerivativesInvalid2) {
  RESET
  Options opt;
  opt.getSection("diff")->set("firstStag","XXX","test");
  EXPECT_THROW(localmesh.initDerivs(&opt),BoutException);
}

TEST_F(MeshTest, SetDerivativesInvalid3) {
  RESET
  Options opt;
  // Invalid for this option - expect error
  opt.getSection("diff")->set("SecondStag","C4","test");
  EXPECT_THROW(localmesh.initDerivs(&opt),BoutException);
}


#undef RESET
