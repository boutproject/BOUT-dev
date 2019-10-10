#include "gtest/gtest.h"

#include "boutexception.hxx"
#include "field.hxx"
#include "initialprofiles.hxx"
#include "output.hxx"
#include "test_extras.hxx"
#include "bout/constants.hxx"
#include "bout/mesh.hxx"

/// Global mesh
namespace bout {
namespace globals {
extern Mesh* mesh;
} // namespace globals
} // namespace bout

// The unit tests use the global mesh
using namespace bout::globals;

class InitialProfileTest : public FakeMeshFixture {
public:
  InitialProfileTest() : FakeMeshFixture() {
    // We need Coordinates so a parallel transform is available as
    // FieldFactory::create3D wants to un-field-align the result
    static_cast<FakeMesh*>(mesh)->setCoordinates(test_coords);

    mesh->getCoordinates()->setParallelTransform(
        bout::utils::make_unique<ParallelTransformIdentity>(*mesh,
          std::vector<BoutReal>()));
  }

  virtual ~InitialProfileTest() { Options::cleanup(); }
  WithQuietOutput quiet{output_info};
};

TEST_F(InitialProfileTest, Field2D) {
  Field2D result{mesh};

  Options::root()["f"]["function"] = "pi";
  initial_profile("f", result);

  EXPECT_TRUE(IsFieldEqual(result, PI));
}

TEST_F(InitialProfileTest, Field3D) {
  Field3D result{mesh};

  Options::root()["g"]["function"] = "2*pi";
  initial_profile("g", result);

  EXPECT_TRUE(IsFieldEqual(result, TWOPI));
}

TEST_F(InitialProfileTest, Vector2DCovariant) {
  Vector2D result{mesh};

  Options::root()["f_x"]["function"] = "1";
  Options::root()["f_y"]["function"] = "2";
  Options::root()["f_z"]["function"] = "3";

  initial_profile("f", result);

  EXPECT_TRUE(IsFieldEqual(result.x, 1));
  EXPECT_TRUE(IsFieldEqual(result.y, 2));
  EXPECT_TRUE(IsFieldEqual(result.z, 3));
}

TEST_F(InitialProfileTest, Vector3DCovariant) {
  Vector3D result{mesh};

  Options::root()["g_x"]["function"] = "4";
  Options::root()["g_y"]["function"] = "5";
  Options::root()["g_z"]["function"] = "6";

  initial_profile("g", result);

  EXPECT_TRUE(IsFieldEqual(result.x, 4));
  EXPECT_TRUE(IsFieldEqual(result.y, 5));
  EXPECT_TRUE(IsFieldEqual(result.z, 6));
}

TEST_F(InitialProfileTest, Vector2DContravariant) {
  Vector2D result{mesh};
  result.covariant = false;

  Options::root()["fx"]["function"] = "7";
  Options::root()["fy"]["function"] = "8";
  Options::root()["fz"]["function"] = "9";

  initial_profile("f", result);

  EXPECT_TRUE(IsFieldEqual(result.x, 7));
  EXPECT_TRUE(IsFieldEqual(result.y, 8));
  EXPECT_TRUE(IsFieldEqual(result.z, 9));
}

TEST_F(InitialProfileTest, Vector3DContravariant) {
  Vector3D result{mesh};
  result.covariant = false;

  Options::root()["gx"]["function"] = "10";
  Options::root()["gy"]["function"] = "11";
  Options::root()["gz"]["function"] = "12";

  initial_profile("g", result);

  EXPECT_TRUE(IsFieldEqual(result.x, 10));
  EXPECT_TRUE(IsFieldEqual(result.y, 11));
  EXPECT_TRUE(IsFieldEqual(result.z, 12));
}

TEST_F(InitialProfileTest, Field2DWithScale) {
  Field2D result{mesh};

  Options::root()["f"]["function"] = "pi";
  Options::root()["f"]["scale"] = 2;
  initial_profile("f", result);

  EXPECT_TRUE(IsFieldEqual(result, TWOPI));
}

TEST_F(InitialProfileTest, Field3DWithScale) {
  Field3D result{mesh};

  Options::root()["g"]["function"] = "2*pi";
  Options::root()["g"]["scale"] = 3;
  initial_profile("g", result);

  EXPECT_TRUE(IsFieldEqual(result, 3 * TWOPI));
}

TEST_F(InitialProfileTest, Vector2DWithScale) {
  Vector2D result{mesh};

  Options::root()["f_x"]["function"] = "1";
  Options::root()["f_x"]["scale"] = 4;
  Options::root()["f_y"]["function"] = "2";
  Options::root()["f_y"]["scale"] = 5;
  Options::root()["f_z"]["function"] = "3";
  Options::root()["f_z"]["scale"] = 6;

  initial_profile("f", result);

  EXPECT_TRUE(IsFieldEqual(result.x, 4));
  EXPECT_TRUE(IsFieldEqual(result.y, 10));
  EXPECT_TRUE(IsFieldEqual(result.z, 18));
}

TEST_F(InitialProfileTest, Vector3DWithScale) {
  Vector3D result{mesh};

  Options::root()["g_x"]["function"] = "4";
  Options::root()["g_x"]["scale"] = 7;
  Options::root()["g_y"]["function"] = "5";
  Options::root()["g_y"]["scale"] = 8;
  Options::root()["g_z"]["function"] = "6";
  Options::root()["g_z"]["scale"] = 9;

  initial_profile("g", result);

  EXPECT_TRUE(IsFieldEqual(result.x, 28));
  EXPECT_TRUE(IsFieldEqual(result.y, 40));
  EXPECT_TRUE(IsFieldEqual(result.z, 54));
}
