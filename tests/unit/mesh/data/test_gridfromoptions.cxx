#include "gtest/gtest.h"

#include "options.hxx"
#include "output.hxx"
#include "test_extras.hxx"
#include "bout/constants.hxx"
#include "bout/griddata.hxx"
#include "bout/mesh.hxx"

#include <numeric>
#include <string>
#include <vector>

/// Global mesh
namespace bout {
namespace globals {
extern Mesh* mesh;
} // namespace globals
} // namespace bout

// The unit tests use the global mesh
using namespace bout::globals;

class GridFromOptionsTest : public FakeMeshFixture {
public:
  GridFromOptionsTest() : FakeMeshFixture(), options(), griddata(&options) {
    // We need a parallel transform as FieldFactory::create3D wants to
    // un-field-align the result
    mesh->setParallelTransform(
        bout::utils::make_unique<ParallelTransformIdentity>(*mesh));

    output_info.disable();
    output_warn.disable();
    options["f"] = expected_string;
  }

  ~GridFromOptionsTest() {
    Options::cleanup();
    output_info.enable();
    output_warn.enable();
  }

  Options options;
  GridFromOptions griddata;
  std::string expected_string{"x + y + z + 3"};
};

TEST_F(GridFromOptionsTest, HasVar) {
  EXPECT_TRUE(griddata.hasVar("f"));
  EXPECT_FALSE(griddata.hasVar("non-existent"));
}

TEST_F(GridFromOptionsTest, GetString) {
  std::string result{"wrong"};

  EXPECT_TRUE(griddata.get(mesh, result, "f"));
  EXPECT_EQ(result, expected_string);
}

TEST_F(GridFromOptionsTest, GetStringNone) {
  std::string result{"wrong"};

  EXPECT_FALSE(griddata.get(mesh, result, "non-existent"));
  EXPECT_EQ(result, std::string{});
}

TEST_F(GridFromOptionsTest, GetInt) {
  int result{-1};
  int expected{3};

  EXPECT_TRUE(griddata.get(mesh, result, "f"));
  EXPECT_EQ(result, expected);
}

TEST_F(GridFromOptionsTest, GetIntNone) {
  int result{-1};
  int expected{0};

  EXPECT_FALSE(griddata.get(mesh, result, "non-existent"));
  EXPECT_EQ(result, expected);
}

TEST_F(GridFromOptionsTest, GetBoutReal) {
  BoutReal result{-1.};
  BoutReal expected{3.};

  EXPECT_TRUE(griddata.get(mesh, result, "f"));
  EXPECT_EQ(result, expected);
}

TEST_F(GridFromOptionsTest, GetBoutRealNone) {
  BoutReal result{-1.};
  BoutReal expected{0};

  EXPECT_FALSE(griddata.get(mesh, result, "non-existent"));
  EXPECT_EQ(result, expected);
}

TEST_F(GridFromOptionsTest, GetField2D) {
  Field2D result{mesh};
  Field2D expected{makeField<Field2D>(
      [](Field2D::ind_type& index) {
        return index.x() + (TWOPI * index.y()) + (TWOPI * index.z() / nz) + 3;
      },
      mesh)};

  EXPECT_TRUE(griddata.get(mesh, result, "f"));
  EXPECT_TRUE(IsFieldEqual(result, expected));
}

TEST_F(GridFromOptionsTest, GetField2DNone) {
  Field2D result{mesh};
  BoutReal expected{0.};

  EXPECT_FALSE(griddata.get(mesh, result, "non-existent"));
  EXPECT_TRUE(IsFieldEqual(result, expected));
}

TEST_F(GridFromOptionsTest, GetField2DNoneWithDefault) {
  Field2D result{mesh};
  BoutReal default_value{-32};

  EXPECT_FALSE(griddata.get(mesh, result, "non-existent", default_value));
  EXPECT_TRUE(IsFieldEqual(result, default_value));
}

TEST_F(GridFromOptionsTest, GetField3D) {
  Field3D result{mesh};
  Field3D expected{makeField<Field3D>(
      [](Field3D::ind_type& index) {
        return index.x() + (TWOPI * index.y()) + (TWOPI * index.z() / nz) + 3;
      },
      mesh)};

  EXPECT_TRUE(griddata.get(mesh, result, "f"));
  EXPECT_TRUE(IsFieldEqual(result, expected));

  EXPECT_FALSE(griddata.get(mesh, result, "non-existent"));
  EXPECT_TRUE(IsFieldEqual(result, 0.));

  BoutReal default_value{-64};
  EXPECT_FALSE(griddata.get(mesh, result, "non-existent", default_value));
  EXPECT_TRUE(IsFieldEqual(result, default_value));
}

TEST_F(GridFromOptionsTest, GetField3DNone) {
  Field3D result{mesh};
  BoutReal expected{0.};

  EXPECT_FALSE(griddata.get(mesh, result, "non-existent"));
  EXPECT_TRUE(IsFieldEqual(result, expected));
}

TEST_F(GridFromOptionsTest, GetField3DNoneWithDefault) {
  Field3D result{mesh};
  BoutReal default_value{-64};

  EXPECT_FALSE(griddata.get(mesh, result, "non-existent", default_value));
  EXPECT_TRUE(IsFieldEqual(result, default_value));
}

TEST_F(GridFromOptionsTest, GetVectorInt) {
  std::vector<int> result{};
  std::vector<int> expected{3, 3, 3};

  EXPECT_TRUE(griddata.get(mesh, result, "f", 3));
  EXPECT_EQ(result, expected);
}

TEST_F(GridFromOptionsTest, GetVectorIntNone) {
  std::vector<int> result{-1, -1, -1};
  std::vector<int> expected{};

  EXPECT_FALSE(griddata.get(mesh, result, "non-existent", 3));
  EXPECT_EQ(result, expected);
}

TEST_F(GridFromOptionsTest, GetVectorBoutRealX) {
  std::vector<BoutReal> result{};
  std::vector<BoutReal> expected{3., 4., 5.};

  EXPECT_TRUE(griddata.get(mesh, result, "f", nx));
  EXPECT_EQ(result, expected);
}

TEST_F(GridFromOptionsTest, GetVectorBoutRealXOffset) {
  std::vector<BoutReal> result{};
  std::vector<BoutReal> expected{4., 5., 6.};

  EXPECT_TRUE(griddata.get(mesh, result, "f", nx, 1, GridDataSource::Direction::X));
  EXPECT_EQ(result, expected);
}

TEST_F(GridFromOptionsTest, GetVectorBoutRealXMeshOffset) {
  std::vector<BoutReal> result{};
  std::vector<BoutReal> expected{2., 3., 4.};

  mesh->OffsetX = 1;
  mesh->OffsetY = 100;
  mesh->OffsetZ = 100;

  EXPECT_TRUE(griddata.get(mesh, result, "f", nx, 0, GridDataSource::Direction::X));
  EXPECT_EQ(result, expected);
}

TEST_F(GridFromOptionsTest, GetVectorBoutRealXNone) {
  std::vector<BoutReal> result{};
  std::vector<BoutReal> default_expected{};
  EXPECT_FALSE(griddata.get(mesh, result, "non-existent", nx));
  EXPECT_EQ(result, default_expected);
}

TEST_F(GridFromOptionsTest, GetVectorBoutRealY) {
  std::vector<BoutReal> result{};
  std::vector<BoutReal> expected{3., 3. + TWOPI, 3. + (2. * TWOPI), 3. + (3. * TWOPI),
                                 3. + (4. * TWOPI)};

  EXPECT_TRUE(griddata.get(mesh, result, "f", ny, 0, GridDataSource::Direction::Y));
  EXPECT_EQ(result, expected);
}

TEST_F(GridFromOptionsTest, GetVectorBoutRealYOffset) {
  std::vector<BoutReal> result{};
  std::vector<BoutReal> expected{3. + TWOPI, 3. + (2. * TWOPI), 3. + (3. * TWOPI),
                                 3. + (4. * TWOPI), 3. + (5. * TWOPI)};

  EXPECT_TRUE(griddata.get(mesh, result, "f", ny, 1, GridDataSource::Direction::Y));
  EXPECT_EQ(result, expected);
}

TEST_F(GridFromOptionsTest, GetVectorBoutRealYMeshOffset) {
  std::vector<BoutReal> result{};
  std::vector<BoutReal> expected{3. - TWOPI, 3., 3. + TWOPI, 3. + (2. * TWOPI),
                                 3. + (3. * TWOPI)};

  mesh->OffsetX = 100;
  mesh->OffsetY = 1;
  mesh->OffsetZ = 100;

  EXPECT_TRUE(griddata.get(mesh, result, "f", ny, 0, GridDataSource::Direction::Y));
  EXPECT_EQ(result, expected);
}

TEST_F(GridFromOptionsTest, GetVectorBoutRealYNone) {
  std::vector<BoutReal> result{};
  std::vector<BoutReal> default_expected{};
  EXPECT_FALSE(griddata.get(mesh, result, "non-existent", ny));
  EXPECT_EQ(result, default_expected);
}

TEST_F(GridFromOptionsTest, GetVectorBoutRealZ) {
  std::vector<BoutReal> result{};
  std::vector<BoutReal> expected{3.,
                                 3. + (1. * TWOPI / nz),
                                 3. + (2. * TWOPI / nz),
                                 3. + (3. * TWOPI / nz),
                                 3. + (4. * TWOPI / nz),
                                 3. + (5. * TWOPI / nz),
                                 3. + (6. * TWOPI / nz)};

  EXPECT_TRUE(griddata.get(mesh, result, "f", nz, 0, GridDataSource::Direction::Z));
  EXPECT_EQ(result, expected);
}

TEST_F(GridFromOptionsTest, GetVectorBoutRealZOffset) {
  std::vector<BoutReal> result{};
  std::vector<BoutReal> expected{3. + (1. * TWOPI / nz), 3. + (2. * TWOPI / nz),
                                 3. + (3. * TWOPI / nz), 3. + (4. * TWOPI / nz),
                                 3. + (5. * TWOPI / nz), 3. + (6. * TWOPI / nz),
                                 3. + (7. * TWOPI / nz)};

  EXPECT_TRUE(griddata.get(mesh, result, "f", nz, 1, GridDataSource::Direction::Z));
  EXPECT_EQ(result, expected);
}

TEST_F(GridFromOptionsTest, GetVectorBoutRealZMeshOffset) {
  std::vector<BoutReal> result{};
  std::vector<BoutReal> expected{3. + (-1. * TWOPI / nz), 3.,
                                 3. + (1. * TWOPI / nz),  3. + (2. * TWOPI / nz),
                                 3. + (3. * TWOPI / nz),  3. + (4. * TWOPI / nz),
                                 3. + (5. * TWOPI / nz)};

  mesh->OffsetX = 100;
  mesh->OffsetY = 100;
  mesh->OffsetZ = 1;

  EXPECT_TRUE(griddata.get(mesh, result, "f", nz, 0, GridDataSource::Direction::Z));
  EXPECT_EQ(result, expected);
}

TEST_F(GridFromOptionsTest, GetVectorBoutRealZNone) {
  std::vector<BoutReal> result{};
  std::vector<BoutReal> default_expected{};
  EXPECT_FALSE(griddata.get(mesh, result, "non-existent", nz));
  EXPECT_EQ(result, default_expected);
}
