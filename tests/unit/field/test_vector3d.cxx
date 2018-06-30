#include "gtest/gtest.h"

#include "bout/constants.hxx"
#include "bout/mesh.hxx"
#include "boutexception.hxx"
#include "vector3d.hxx"
#include "test_extras.hxx"
#include "unused.hxx"

/// Global mesh
extern Mesh *mesh;

/// Test fixture to make sure the global mesh is our fake one
class Vector3DTest : public ::testing::Test {
protected:
  static void SetUpTestCase() {
    // Delete any existing mesh
    if (mesh != nullptr) {
      // Delete boundary regions
      for (auto &r : mesh->getBoundaries()) {
        delete r;
      }

      delete mesh;
      mesh = nullptr;
    }
    mesh = new FakeMesh(nx, ny, nz);
    mesh->createDefaultRegions();

    mesh->addBoundary(new BoundaryRegionXIn("core", 1, ny - 2));
    mesh->addBoundary(new BoundaryRegionXOut("sol", 1, ny - 2));
    mesh->addBoundary(new BoundaryRegionYUp("upper_target", 1, nx - 2));
    mesh->addBoundary(new BoundaryRegionYDown("lower_target", 1, nx - 2));
  }

  static void TearDownTestCase() {
    if (mesh != nullptr) {
      // Delete boundary regions
      for (auto &r : mesh->getBoundaries()) {
        delete r;
      }
    }
    delete mesh;
    mesh = nullptr;
  }

public:
  static const int nx;
  static const int ny;
  static const int nz;
};

const int Vector3DTest::nx = 5;
const int Vector3DTest::ny = 5;
const int Vector3DTest::nz = 3;

TEST_F(Vector3DTest, ApplyBoundaryString) {
  Vector3D v;
  v = 0.0;
  v.applyBoundary("dirichlet(1.0)");

  // boundary cell in x
  EXPECT_DOUBLE_EQ(v.x(0,2,0), 2.0);
  EXPECT_DOUBLE_EQ(v.y(4,2,1), 2.0);
  
  // boundary cell in y
  EXPECT_DOUBLE_EQ(v.x(2,0,2), 2.0);
  EXPECT_DOUBLE_EQ(v.z(2,4,0), 2.0);

  // Middle cell not changed
  EXPECT_DOUBLE_EQ(v.x(2,2,1), 0.0);
}

TEST_F(Vector3DTest, IsReal) {
  Vector3D vector;

  EXPECT_TRUE(vector.isReal());
}

TEST_F(Vector3DTest, Is3D) {
  Vector3D vector;

  EXPECT_TRUE(vector.is3D());
}

TEST_F(Vector3DTest, ByteSize) {
  Vector3D vector;

  EXPECT_EQ(vector.byteSize(), 3 * sizeof(BoutReal));
}

TEST_F(Vector3DTest, BoutRealSize) {
  Vector3D vector;

  EXPECT_EQ(vector.BoutRealSize(), 3);
}

TEST_F(Vector3DTest, AssignFromBoutReal) {
  Vector3D vector;

  vector = 0.0;

  EXPECT_TRUE(IsField3DEqualBoutReal(vector.x, 0.0));
  EXPECT_TRUE(IsField3DEqualBoutReal(vector.y, 0.0));
  EXPECT_TRUE(IsField3DEqualBoutReal(vector.z, 0.0));
}

TEST_F(Vector3DTest, AssignFromVector3D) {
  Vector3D vector1, vector2;
  vector1.x = 1.0;
  vector1.y = 2.0;
  vector1.z = 3.0;

  vector2 = vector1;

  EXPECT_TRUE(IsField3DEqualBoutReal(vector2.x, 1.0));
  EXPECT_TRUE(IsField3DEqualBoutReal(vector2.y, 2.0));
  EXPECT_TRUE(IsField3DEqualBoutReal(vector2.z, 3.0));
}

TEST_F(Vector3DTest, CreateFromVector3D) {
  Vector3D vector1;
  vector1.x = 4.0;
  vector1.y = 5.0;
  vector1.z = 6.0;

  Vector3D vector2{vector1};

  EXPECT_TRUE(IsField3DEqualBoutReal(vector2.x, 4.0));
  EXPECT_TRUE(IsField3DEqualBoutReal(vector2.y, 5.0));
  EXPECT_TRUE(IsField3DEqualBoutReal(vector2.z, 6.0));
}

TEST_F(Vector3DTest, UnaryMinus) {
  Vector3D vector1;
  vector1.x = 7.0;
  vector1.y = 8.0;
  vector1.z = 9.0;

  Vector3D vector2;

  vector2 = -vector1;

  EXPECT_TRUE(IsField3DEqualBoutReal(vector2.x, -7.0));
  EXPECT_TRUE(IsField3DEqualBoutReal(vector2.y, -8.0));
  EXPECT_TRUE(IsField3DEqualBoutReal(vector2.z, -9.0));
}

TEST_F(Vector3DTest, AddEqualsVector3D) {
  Vector3D vector1;
  vector1.x = 10.0;
  vector1.y = 11.0;
  vector1.z = 12.0;

  Vector3D vector2;
  vector2.x = 1.0;
  vector2.y = 2.0;
  vector2.z = 3.0;

  vector2 += vector1;

  EXPECT_TRUE(IsField3DEqualBoutReal(vector2.x, 11.0));
  EXPECT_TRUE(IsField3DEqualBoutReal(vector2.y, 13.0));
  EXPECT_TRUE(IsField3DEqualBoutReal(vector2.z, 15.0));
}

TEST_F(Vector3DTest, AddVector3DVector2D) {
  Vector3D vector1;
  vector1.x = 13.0;
  vector1.y = 14.0;
  vector1.z = 15.0;

  Vector2D vector2;
  vector2.x = 4.0;
  vector2.y = 5.0;
  vector2.z = 6.0;

  Vector3D result = vector1 + vector2;

  EXPECT_TRUE(IsField3DEqualBoutReal(result.x, 17.0));
  EXPECT_TRUE(IsField3DEqualBoutReal(result.y, 19.0));
  EXPECT_TRUE(IsField3DEqualBoutReal(result.z, 21.0));
}

TEST_F(Vector3DTest, AddVector3DVector3D) {
  Vector3D vector1;
  vector1.x = 3.0;
  vector1.y = 4.0;
  vector1.z = 5.0;

  Vector3D vector2;
  vector2.x = 4.0;
  vector2.y = 5.0;
  vector2.z = 6.0;

  Vector3D result = vector1 + vector2;

  EXPECT_TRUE(IsField3DEqualBoutReal(result.x, 7.0));
  EXPECT_TRUE(IsField3DEqualBoutReal(result.y, 9.0));
  EXPECT_TRUE(IsField3DEqualBoutReal(result.z, 11.0));
}

TEST_F(Vector3DTest, MinusEqualsVector3D) {
  Vector3D vector1;
  vector1.x = 100.0;
  vector1.y = 101.0;
  vector1.z = 102.0;

  Vector3D vector2;
  vector2.x = 3.0;
  vector2.y = 2.0;
  vector2.z = 1.0;

  vector2 -= vector1;

  EXPECT_TRUE(IsField3DEqualBoutReal(vector2.x, -97.0));
  EXPECT_TRUE(IsField3DEqualBoutReal(vector2.y, -99.0));
  EXPECT_TRUE(IsField3DEqualBoutReal(vector2.z, -101.0));
}

TEST_F(Vector3DTest, MinusVector3DVector2D) {
  Vector3D vector1;
  vector1.x = 0.0;
  vector1.y = 1.0;
  vector1.z = 2.0;

  Vector2D vector2;
  vector2.x = 3.0;
  vector2.y = 2.0;
  vector2.z = 1.0;

  Vector3D result = vector1 - vector2;

  EXPECT_TRUE(IsField3DEqualBoutReal(result.x, -3.0));
  EXPECT_TRUE(IsField3DEqualBoutReal(result.y, -1.0));
  EXPECT_TRUE(IsField3DEqualBoutReal(result.z, 1.0));
}

TEST_F(Vector3DTest, MinusVector3DVector3D) {
  Vector3D vector1;
  vector1.x = 10.0;
  vector1.y = 11.0;
  vector1.z = 12.0;

  Vector3D vector2;
  vector2.x = 3.0;
  vector2.y = 2.0;
  vector2.z = 1.0;

  Vector3D result = vector1 - vector2;

  EXPECT_TRUE(IsField3DEqualBoutReal(result.x, 7.0));
  EXPECT_TRUE(IsField3DEqualBoutReal(result.y, 9.0));
  EXPECT_TRUE(IsField3DEqualBoutReal(result.z, 11.0));
}

TEST_F(Vector3DTest, MultiplyEqualsBoutReal) {
  Vector3D vector;
  vector.x = 4.0;
  vector.y = 5.0;
  vector.z = 6.0;

  BoutReal real {4.0};

  vector *= real;

  EXPECT_TRUE(IsField3DEqualBoutReal(vector.x, 16.0));
  EXPECT_TRUE(IsField3DEqualBoutReal(vector.y, 20.0));
  EXPECT_TRUE(IsField3DEqualBoutReal(vector.z, 24.0));
}

TEST_F(Vector3DTest, MultiplyEqualsField2D) {
  Vector3D vector;
  vector.x = 4.0;
  vector.y = 5.0;
  vector.z = 6.0;

  Field2D field{40.0};

  vector *= field;

  EXPECT_TRUE(IsField3DEqualBoutReal(vector.x, 160.0));
  EXPECT_TRUE(IsField3DEqualBoutReal(vector.y, 200.0));
  EXPECT_TRUE(IsField3DEqualBoutReal(vector.z, 240.0));
}

TEST_F(Vector3DTest, MultiplyVector3DBoutReal) {
  Vector3D vector;
  vector.x = 1.0;
  vector.y = 2.0;
  vector.z = 3.0;

  BoutReal real {2.0};

  Vector3D result = vector * real;

  EXPECT_TRUE(IsField3DEqualBoutReal(result.x, 2.0));
  EXPECT_TRUE(IsField3DEqualBoutReal(result.y, 4.0));
  EXPECT_TRUE(IsField3DEqualBoutReal(result.z, 6.0));
}

TEST_F(Vector3DTest, MultiplyVector3DField2D) {
  Vector3D vector;
  vector.x = 1.0;
  vector.y = 2.0;
  vector.z = 3.0;

  Field2D field{3.0};

  Vector3D result = vector * field;

  EXPECT_TRUE(IsField3DEqualBoutReal(result.x, 3.0));
  EXPECT_TRUE(IsField3DEqualBoutReal(result.y, 6.0));
  EXPECT_TRUE(IsField3DEqualBoutReal(result.z, 9.0));
}

TEST_F(Vector3DTest, MultiplyVector3DField3D) {
  Vector3D vector;
  vector.x = 1.0;
  vector.y = 2.0;
  vector.z = 3.0;

  Field3D field{4.0};

  Vector3D result = vector * field;

  EXPECT_TRUE(IsField3DEqualBoutReal(result.x, 4.0));
  EXPECT_TRUE(IsField3DEqualBoutReal(result.y, 8.0));
  EXPECT_TRUE(IsField3DEqualBoutReal(result.z, 12.0));
}

TEST_F(Vector3DTest, DivideEqualsBoutReal) {
  Vector3D vector;
  vector.x = 1.0;
  vector.y = 2.0;
  vector.z = 3.0;

  BoutReal real{10.0};

  vector /= real;

  EXPECT_TRUE(IsField3DEqualBoutReal(vector.x, 0.1));
  EXPECT_TRUE(IsField3DEqualBoutReal(vector.y, 0.2));
  EXPECT_TRUE(IsField3DEqualBoutReal(vector.z, 0.3));
}

TEST_F(Vector3DTest, DivideEqualsConstField2D) {
  Vector3D vector;
  vector.x = 1.0;
  vector.y = 2.0;
  vector.z = 3.0;

  Field2D field{5.0};

  vector /= field;

  EXPECT_TRUE(IsField3DEqualBoutReal(vector.x, 0.2));
  EXPECT_TRUE(IsField3DEqualBoutReal(vector.y, 0.4));
  EXPECT_TRUE(IsField3DEqualBoutReal(vector.z, 0.6));
}

TEST_F(Vector3DTest, DivideVector3DBoutReal) {
  Vector3D vector;
  vector.x = 1.0;
  vector.y = 2.0;
  vector.z = 3.0;

  BoutReal real {2.0};

  Vector3D result = vector / real;

  EXPECT_TRUE(IsField3DEqualBoutReal(result.x, 0.5));
  EXPECT_TRUE(IsField3DEqualBoutReal(result.y, 1.0));
  EXPECT_TRUE(IsField3DEqualBoutReal(result.z, 1.5));
}

TEST_F(Vector3DTest, DivideVector3DField2D) {
  Vector3D vector;
  vector.x = 1.0;
  vector.y = 2.0;
  vector.z = 3.0;

  Field2D field{4.0};

  Vector3D result = vector / field;

  EXPECT_TRUE(IsField3DEqualBoutReal(result.x, 0.25));
  EXPECT_TRUE(IsField3DEqualBoutReal(result.y, 0.5));
  EXPECT_TRUE(IsField3DEqualBoutReal(result.z, 0.75));
}

TEST_F(Vector3DTest, DivideVector3DField3D) {
  Vector3D vector;
  vector.x = 2.0;
  vector.y = 4.0;
  vector.z = 6.0;

  Field3D field{2.0};

  Vector3D result = vector / field;

  EXPECT_TRUE(IsField3DEqualBoutReal(result.x, 1.0));
  EXPECT_TRUE(IsField3DEqualBoutReal(result.y, 2.0));
  EXPECT_TRUE(IsField3DEqualBoutReal(result.z, 3.0));
}
