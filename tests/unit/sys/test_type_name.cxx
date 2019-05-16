#include "gtest/gtest.h"
#include "bout/sys/type_name.hxx"
#include "field2d.hxx"
#include "field3d.hxx"

#include <string>

using bout::utils::typeName;

TEST(TypeNameTest, BoolName) {
  EXPECT_EQ(typeName<bool>(), "bool");
}

TEST(TypeNameTest, IntName) {
  EXPECT_EQ(typeName<int>(), "int");
}

TEST(TypeNameTest, StringName) {
  EXPECT_EQ(typeName<std::string>(), "string");
}

TEST(TypeNameTest, BoutRealName) {
  EXPECT_EQ(typeName<BoutReal>(), "BoutReal");
}

TEST(TypeNameTest, Field2DName) {
  EXPECT_EQ(typeName<Field2D>(), "Field2D");
}

TEST(TypeNameTest, Field3DName) {
  EXPECT_EQ(typeName<Field3D>(), "Field3D");
}
