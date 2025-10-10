#include "gtest/gtest.h"

#include "bout/field2d.hxx"
#include "bout/field3d.hxx"
#include "bout/fieldgroup.hxx"
#include "bout/vector2d.hxx"
#include "bout/vector3d.hxx"

#include <type_traits>

TEST(FieldGroupTest, CreateWithField2D) {
  Field2D a;
  FieldGroup group(a);

  EXPECT_EQ(group.size(), 1);
  EXPECT_EQ(group.size_field3d(), 0);
}

TEST(FieldGroupTest, CreateWithField3D) {
  Field3D a;
  FieldGroup group(a);

  EXPECT_EQ(group.size(), 1);
  EXPECT_EQ(group.size_field3d(), 1);
}

TEST(FieldGroupTest, CreateWithField2DField3D) {
  Field2D a;
  Field3D b;
  FieldGroup group(a, b);

  EXPECT_EQ(group.size(), 2);
  EXPECT_EQ(group.size_field3d(), 1);
}

TEST(FieldGroupTest, CreateWithVector2D) {
  Vector2D a;
  FieldGroup group(a);

  EXPECT_EQ(group.size(), 3);
  EXPECT_EQ(group.size_field3d(), 0);
}

TEST(FieldGroupTest, CreateWithVector3D) {
  Vector3D a;
  FieldGroup group(a);

  EXPECT_EQ(group.size(), 3);
  EXPECT_EQ(group.size_field3d(), 3);
}

TEST(FieldGroupTest, CreateWithVector2DVector3D) {
  Vector2D a;
  Vector3D b;
  FieldGroup group(a, b);

  EXPECT_EQ(group.size(), 6);
  EXPECT_EQ(group.size_field3d(), 3);
}

TEST(FieldGroupTest, CreateWithField2DField3DVector2DVector3D) {
  Field2D a;
  Field3D b;
  Vector2D c;
  Vector3D d;
  FieldGroup group(a, b, c, d);

  EXPECT_EQ(group.size(), 8);
  EXPECT_EQ(group.size_field3d(), 4);
}

TEST(FieldGroupTest, AddField2D) {
  Field2D a;
  FieldGroup group;

  group.add(a);

  EXPECT_EQ(group.size(), 1);
  EXPECT_EQ(group.size_field3d(), 0);
}

TEST(FieldGroupTest, AddField3D) {
  Field3D a;
  FieldGroup group;

  group.add(a);

  EXPECT_EQ(group.size(), 1);
  EXPECT_EQ(group.size_field3d(), 1);
}

TEST(FieldGroupTest, AddVector2D) {
  Vector2D a;
  FieldGroup group;

  group.add(a);

  EXPECT_EQ(group.size(), 3);
  EXPECT_EQ(group.size_field3d(), 0);
}

TEST(FieldGroupTest, AddVector3D) {
  Vector3D a;
  FieldGroup group;

  group.add(a);

  EXPECT_EQ(group.size(), 3);
  EXPECT_EQ(group.size_field3d(), 3);
}

TEST(FieldGroupTest, AddFieldGroup) {
  Field3D a;
  FieldGroup group, other_group;

  group.add(a);

  other_group.add(group);

  EXPECT_EQ(other_group.size(), 1);
  EXPECT_EQ(other_group.size_field3d(), 1);
}

TEST(FieldGroupTest, AddEqualsFieldGroup) {
  Field3D a;
  FieldGroup group, other_group;

  group.add(a);

  other_group += group;

  EXPECT_EQ(other_group.size(), 1);
  EXPECT_EQ(other_group.size_field3d(), 1);
}

TEST(FieldGroupTest, AddOperatorFieldGroup) {
  Field3D a;
  FieldGroup group, other_group, result_group;

  group.add(a);

  result_group = other_group + group;

  EXPECT_EQ(result_group.size(), 1);
  EXPECT_EQ(result_group.size_field3d(), 1);
}

TEST(FieldGroupTest, AddField2DField3D) {
  Field2D a;
  Field3D b;
  FieldGroup group;

  group.add(a, b);

  EXPECT_EQ(group.size(), 2);
  EXPECT_EQ(group.size_field3d(), 1);
}

TEST(FieldGroupTest, AddVector2DVector3D) {
  Vector2D a;
  Vector3D b;
  FieldGroup group;

  group.add(a, b);

  EXPECT_EQ(group.size(), 6);
  EXPECT_EQ(group.size_field3d(), 3);
}

TEST(FieldGroupTest, AddField2DField3DVector2DVector3D) {
  Field2D a;
  Field3D b;
  Vector2D c;
  Vector3D d;
  FieldGroup group;

  group.add(a, b, c, d);

  EXPECT_EQ(group.size(), 8);
  EXPECT_EQ(group.size_field3d(), 4);
}

TEST(FieldGroupTest, Empty) {
  FieldGroup group;

  EXPECT_TRUE(group.empty());
}

TEST(FieldGroupTest, NotEmpty) {
  Field2D a;
  FieldGroup group(a);

  EXPECT_FALSE(group.empty());
}

TEST(FieldGroupTest, NotEmptyField3D) {
  Field3D a;
  FieldGroup group(a);

  EXPECT_FALSE(group.empty());
}

TEST(FieldGroupTest, Clear) {
  Field2D a;
  Field3D b;
  FieldGroup group(a, b);

  group.clear();
  EXPECT_TRUE(group.empty());
}

TEST(FieldGroupTest, MakeUnique) {
  Field3D a;
  FieldGroup group;

  group.add(a);
  group.add(a);

  EXPECT_EQ(group.size(), 2);

  group.makeUnique();

  EXPECT_EQ(group.size(), 1);
  EXPECT_EQ(group.size_field3d(), 1);
}

TEST(FieldGroupTest, Get) {
  Field3D a;
  Field2D b;
  FieldGroup group;

  group.add(a);
  group.add(b);

  auto vec = group.get();

  EXPECT_EQ(vec.size(), 2);
}

TEST(FieldGroupTest, Field3d) {
  Field3D a;
  Field2D b;
  FieldGroup group;

  group.add(a);
  group.add(b);

  auto vec = group.field3d();

  EXPECT_EQ(vec.size(), 1);
}

TEST(FieldGroupTest, Iterator) {
  Field3D a;
  Field2D b;
  FieldGroup group;

  group.add(a);
  group.add(b);

  int count = 0;

  for (auto iter = group.begin(); iter < group.end(); ++iter) {
    ++count;
  }

  EXPECT_EQ(count, 2);
}

TEST(FieldGroupTest, ConstIterator) {
  Field3D a;
  Field2D b;
  const FieldGroup group(a, b);

  int count = 0;

  for (auto iter = group.begin(); iter < group.end(); ++iter) {
    ++count;
  }

  EXPECT_EQ(count, 2);
}

TEST(FieldGroupTest, NotConstructableFromInt) {
  EXPECT_FALSE((std::is_constructible_v<FieldGroup, int>));
}

TEST(FieldGroupTest, NotConstructableFromBoutReal) {
  EXPECT_FALSE((std::is_constructible_v<FieldGroup, BoutReal>));
}
