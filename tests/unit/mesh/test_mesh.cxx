#include "gtest/gtest.h"
#include "gmock/gmock.h"

#include "bout/mesh.hxx"
#include "bout/region.hxx"
#include "boutexception.hxx"
#include "output.hxx"

#include "test_extras.hxx"

/// Test fixture to make sure the global mesh is our fake one
class MeshTest : public ::testing::Test {
public:
  MeshTest() : localmesh(nx, ny, nz) {}
  static const int nx = 3;
  static const int ny = 5;
  static const int nz = 7;
  FakeMesh localmesh;

  WithQuietOutput quiet_warn{output_warn};
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

TEST_F(MeshTest, HasRegion3D) {
  localmesh.createDefaultRegions();
  EXPECT_TRUE(localmesh.hasRegion3D("RGN_ALL"));
  EXPECT_TRUE(localmesh.hasRegion3D("RGN_NOBNDRY"));
  EXPECT_FALSE(localmesh.hasRegion3D("SOME_MADE_UP_REGION_NAME"));
}

TEST_F(MeshTest, HasRegion2D) {
  localmesh.createDefaultRegions();
  EXPECT_TRUE(localmesh.hasRegion2D("RGN_ALL"));
  EXPECT_TRUE(localmesh.hasRegion2D("RGN_NOBNDRY"));
  EXPECT_FALSE(localmesh.hasRegion2D("SOME_MADE_UP_REGION_NAME"));
}

TEST_F(MeshTest, HasRegionPMesh) {
  localmesh.createDefaultRegions();
  EXPECT_TRUE(localmesh.hasRegionPerp("RGN_ALL"));
  EXPECT_TRUE(localmesh.hasRegionPerp("RGN_NOBNDRY"));
  EXPECT_FALSE(localmesh.hasRegionPerp("SOME_MADE_UP_REGION_NAME"));
}

TEST_F(MeshTest, GetRegionTemplatedFromMesh) {
  using namespace ::testing;
  localmesh.createDefaultRegions();

  const auto& region3d = localmesh.getRegion3D("RGN_ALL");
  const auto& regionT_3d = localmesh.getRegion<Field3D>("RGN_ALL");
  EXPECT_THAT(regionT_3d, ElementsAreArray(region3d));

  const auto& region2d = localmesh.getRegion2D("RGN_ALL");
  const auto& regionT_2d = localmesh.getRegion<Field2D>("RGN_ALL");
  EXPECT_THAT(regionT_2d, ElementsAreArray(region2d));

  const auto& regionPerp = localmesh.getRegionPerp("RGN_ALL");
  const auto& regionT_Perp = localmesh.getRegion<FieldPerp>("RGN_ALL");
  EXPECT_THAT(regionT_Perp, ElementsAreArray(regionPerp));
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

TEST_F(MeshTest, GetStringNoSource) {
  std::string string_value;
  EXPECT_NE(localmesh.get(string_value, "no_source"), 0);
  EXPECT_EQ(string_value, "");
}

TEST_F(MeshTest, GetStringNoSourceWithDefault) {
  std::string string_value;
  const std::string default_value = "some default";
  EXPECT_NE(localmesh.get(string_value, "no_source", default_value), 0);
  EXPECT_EQ(string_value, default_value);
}

TEST_F(MeshTest, GetIntNoSource) {
  int int_value;
  EXPECT_NE(localmesh.get(int_value, "no_source"), 0);
  EXPECT_EQ(int_value, 0);
}

TEST_F(MeshTest, GetIntNoSourceWithDefault) {
  int int_value;
  constexpr int default_value = 42;
  EXPECT_NE(localmesh.get(int_value, "no_source", default_value), 0);
  EXPECT_EQ(int_value, default_value);
}

TEST_F(MeshTest, GetBoutRealNoSource) {
  BoutReal boutreal_value;
  EXPECT_NE(localmesh.get(boutreal_value, "no_source"), 0);
  EXPECT_EQ(boutreal_value, 0.0);
}

TEST_F(MeshTest, GetBoutRealNoSourceWithDefault) {
  BoutReal boutreal_value;
  constexpr BoutReal default_value = 3.14;
  EXPECT_NE(localmesh.get(boutreal_value, "no_source", default_value), 0);
  EXPECT_DOUBLE_EQ(boutreal_value, default_value);
}

TEST_F(MeshTest, GetField2DNoSource) {
  WithQuietOutput warn{output_warn};

  localmesh.createDefaultRegions();
  localmesh.setCoordinates(nullptr);

  Field2D field2d_value{&localmesh};
  EXPECT_NE(localmesh.get(field2d_value, "no_source"), 0);
  EXPECT_TRUE(IsFieldEqual(field2d_value, 0.0));
}

TEST_F(MeshTest, GetField2DNoSourceWithDefault) {
  WithQuietOutput warn{output_warn};

  localmesh.createDefaultRegions();
  localmesh.setCoordinates(nullptr);

  Field2D field2d_value{&localmesh};
  constexpr BoutReal default_value = 4.2;
  EXPECT_NE(localmesh.get(field2d_value, "no_source", default_value), 0);
  EXPECT_TRUE(IsFieldEqual(field2d_value, default_value));
}

TEST_F(MeshTest, GetField3DNoSource) {
  WithQuietOutput warn{output_warn};

  localmesh.createDefaultRegions();
  localmesh.setCoordinates(nullptr);

  Field3D field3d_value{&localmesh};
  EXPECT_NE(localmesh.get(field3d_value, "no_source"), 0);
  EXPECT_TRUE(IsFieldEqual(field3d_value, 0.0));
}

TEST_F(MeshTest, GetField3DNoSourceWithDefault) {
  WithQuietOutput warn{output_warn};

  localmesh.createDefaultRegions();
  localmesh.setCoordinates(nullptr);

  Field3D field3d_value{&localmesh};
  constexpr BoutReal default_value = 4.2;
  EXPECT_NE(localmesh.get(field3d_value, "no_source", default_value), 0);
  EXPECT_TRUE(IsFieldEqual(field3d_value, default_value));
}

TEST_F(MeshTest, MsgLen) {
  localmesh.createDefaultRegions();
  localmesh.setCoordinates(nullptr);

  Field3D f3D_1(0., &localmesh);
  Field3D f3D_2(0., &localmesh);
  Field2D f2D_1(0., &localmesh);
  Field2D f2D_2(0., &localmesh);

  std::vector<FieldData*> var_list {&f3D_1, &f2D_1, &f3D_2, &f2D_2};

  const int len = localmesh.msg_len(var_list, 0, nx, 0, ny);

  EXPECT_EQ(len, 2 * (nx * ny * nz) + 2 * (nx * ny));
}
