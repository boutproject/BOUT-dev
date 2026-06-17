#include "bout/field2d.hxx"
#include "bout/field3d.hxx"
#include <bout/mesh.hxx>
#include "gtest/gtest.h"

class NVectorFixture : public ::testing::Test {
public:
  Mesh mesh;
  Field2D field1;

  NVectorFixture() : field1{0} {}
}

TEST_F(NVectorFixture, Constructor) {
}
