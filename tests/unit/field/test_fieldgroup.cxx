#include "gtest/gtest.h"
#include "test_extras.hxx"

#include "bout/fieldgroup.hxx"
#include "bout/mesh.hxx"
#include "field2d.hxx"
#include "field3d.hxx"
#include "unused.hxx"

/// Global mesh
extern Mesh *mesh;

/// Test fixture to make sure the global mesh is our fake one
class FieldGroupTest : public ::testing::Test {
protected:
  static void SetUpTestCase() {
    // Delete any existing mesh
    if (mesh != nullptr) {
      delete mesh;
      mesh = nullptr;
    }
    mesh = new FakeMesh(nx, ny, nz);
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

const int FieldGroupTest::nx = 3;
const int FieldGroupTest::ny = 5;
const int FieldGroupTest::nz = 7;

TEST_F(FieldGroupTest, CreateWithField2D) {
  Field2D a;
  FieldGroup group(a);

  EXPECT_EQ(group.size(), 1);
  EXPECT_EQ(group.size_field3d(), 0);
}

TEST_F(FieldGroupTest, CreateWithField3D) {
  Field3D a;
  FieldGroup group(a);

  EXPECT_EQ(group.size(), 1);
  EXPECT_EQ(group.size_field3d(), 1);
}

TEST_F(FieldGroupTest, CreateWithField2DField3D) {
  Field2D a;
  Field3D b;
  FieldGroup group(a, b);

  EXPECT_EQ(group.size(), 2);
  EXPECT_EQ(group.size_field3d(), 1);
}

TEST_F(FieldGroupTest, AddField2D) {
  Field2D a;
  FieldGroup group;

  group.add(a);

  EXPECT_EQ(group.size(), 1);
  EXPECT_EQ(group.size_field3d(), 0);
}

TEST_F(FieldGroupTest, AddField3D) {
  Field3D a;
  FieldGroup group;

  group.add(a);

  EXPECT_EQ(group.size(), 1);
  EXPECT_EQ(group.size_field3d(), 1);
}

TEST_F(FieldGroupTest, AddFieldGroup) {
  Field3D a;
  FieldGroup group, other_group;

  group.add(a);

  other_group.add(group);

  EXPECT_EQ(other_group.size(), 1);
  EXPECT_EQ(other_group.size_field3d(), 1);
}

TEST_F(FieldGroupTest, AddEqualsFieldGroup) {
  Field3D a;
  FieldGroup group, other_group;

  group.add(a);

  other_group += group;

  EXPECT_EQ(other_group.size(), 1);
  EXPECT_EQ(other_group.size_field3d(), 1);
}

TEST_F(FieldGroupTest, AddOperatorFieldGroup) {
  Field3D a;
  FieldGroup group, other_group, result_group;

  group.add(a);

  result_group = other_group + group;

  EXPECT_EQ(result_group.size(), 1);
  EXPECT_EQ(result_group.size_field3d(), 1);
}

TEST_F(FieldGroupTest, AddField2DField3D) {
  Field2D a;
  Field3D b;
  FieldGroup group;

  group.add(a, b);

  EXPECT_EQ(group.size(), 2);
  EXPECT_EQ(group.size_field3d(), 1);
}

TEST_F(FieldGroupTest, Empty) {
  FieldGroup group;

  EXPECT_TRUE(group.empty());
}

TEST_F(FieldGroupTest, NotEmpty) {
  Field2D a;
  FieldGroup group(a);

  EXPECT_FALSE(group.empty());
}

TEST_F(FieldGroupTest, NotEmptyField3D) {
  Field3D a;
  FieldGroup group(a);

  EXPECT_FALSE(group.empty());
}

TEST_F(FieldGroupTest, Clear) {
  Field2D a;
  Field3D b;
  FieldGroup group(a, b);

  group.clear();
  EXPECT_TRUE(group.empty());
}

TEST_F(FieldGroupTest, MakeUnique) {
  Field3D a;
  FieldGroup group;

  group.add(a);
  group.add(a);

  EXPECT_EQ(group.size(), 2);

  group.makeUnique();

  EXPECT_EQ(group.size(), 1);
  EXPECT_EQ(group.size_field3d(), 1);
}

