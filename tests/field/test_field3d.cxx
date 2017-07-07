#include "gtest/gtest.h"

#include "field3d.hxx"
#include "bout/mesh.hxx"
#include "test_extras.hxx"

// Global mesh
extern Mesh *mesh;

class Field3DTest : public ::testing::Test {
public:
  Field3DTest() {
    // Make sure the global mesh is our fake one
    if (mesh == nullptr) {
      delete mesh;
      mesh = nullptr;
    }
    mesh = new FakeMesh(nx, ny, nz);
  }
  ~Field3DTest() {
    delete mesh;
    mesh = nullptr;
  }

  int nx = 3;
  int ny = 5;
  int nz = 7;
};

TEST_F(Field3DTest, Allocate) {
  Field3D field;

  field.allocate();

  EXPECT_TRUE(field.isAllocated());
}

TEST_F(Field3DTest, GetGridSizes) {
  Field3D field;

  field.allocate();

  EXPECT_EQ(field.getNx(), nx);
  EXPECT_EQ(field.getNy(), ny);
  EXPECT_EQ(field.getNz(), nz);
}

TEST_F(Field3DTest, CreateFromBoutReal) {
  Field3D field(1.0);

  EXPECT_TRUE(IsField3DEqualBoutReal(field, 1.0));
}

TEST_F(Field3DTest, AssignFromBoutReal) {
  Field3D field;

  field = 2.0;

  EXPECT_TRUE(IsField3DEqualBoutReal(field, 2.0));
}

TEST_F(Field3DTest, Add) {
  Field3D a, b, c;

  a = 1.0;
  b = 2.0;
  c = a + b;

  EXPECT_TRUE(IsField3DEqualBoutReal(c, 3.0));
}
