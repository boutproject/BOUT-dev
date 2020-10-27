#include "gtest/gtest.h"

#include "bout/constants.hxx"
#include "bout/mesh.hxx"
#include "boutexception.hxx"
#include "field.hxx"
#include "output.hxx"
#include "test_extras.hxx"

/// Global mesh
namespace bout{
namespace globals{
extern Mesh *mesh;
} // namespace globals
} // namespace bout

// The unit tests use the global mesh
using namespace bout::globals;

// Reuse the "standard" fixture for FakeMesh
using FieldTest = FakeMeshFixture;


TEST_F(FieldTest, GetGlobalMesh) {
  Field field;

  auto localmesh = field.getMesh();

  EXPECT_EQ(localmesh, mesh);
}

TEST_F(FieldTest, GetLocalMesh) {
  FakeMesh myMesh{nx + 1, ny + 2, nz + 3};
  myMesh.setCoordinates(nullptr);
  Field field(&myMesh, CELL_CENTRE, {YDirectionType::Standard, ZDirectionType::Standard});

  auto localmesh = field.getMesh();

  EXPECT_EQ(localmesh, &myMesh);
}

TEST_F(FieldTest, GetGridSizes) {
  Field field;

  EXPECT_EQ(field.getNx(), nx);
  EXPECT_EQ(field.getNy(), ny);
  EXPECT_EQ(field.getNz(), nz);
}

TEST_F(FieldTest, AreFieldsCompatibleTrue) {
  // Create a field with non-default members
  Field field{mesh_staggered, CELL_XLOW, {YDirectionType::Aligned, ZDirectionType::Average}};

  Field field2{field};

  EXPECT_TRUE(areFieldsCompatible(field, field2));
  EXPECT_EQ(field.getMesh(), field2.getMesh());
  EXPECT_EQ(field.getLocation(), field2.getLocation());
  EXPECT_EQ(field.getDirectionY(), field2.getDirectionY());
  EXPECT_EQ(field.getDirectionZ(), field2.getDirectionZ());
}

TEST_F(FieldTest, AreFieldsCompatibleFalseMesh) {
  // Create a field with default members
  Field field;

  FakeMesh myMesh{nx + 1, ny + 2, nz + 3};
  myMesh.setCoordinates(nullptr);

  // Create a field with all members set explicitly, and a non-default mesh
  Field field2{&myMesh, CELL_CENTRE, {YDirectionType::Standard, ZDirectionType::Standard}};

  EXPECT_FALSE(areFieldsCompatible(field, field2));
  EXPECT_NE(field.getMesh(), field2.getMesh());
  EXPECT_EQ(field.getLocation(), field2.getLocation());
  EXPECT_EQ(field.getDirectionY(), field2.getDirectionY());
  EXPECT_EQ(field.getDirectionZ(), field2.getDirectionZ());
}

TEST_F(FieldTest, AreFieldsCompatibleFalseLocation) {
  // Create a field with default members
  Field field{mesh_staggered, CELL_CENTRE, {YDirectionType::Standard, ZDirectionType::Standard}};

  // Create a field with all members set explicitly, and a non-default location
  Field field2{mesh_staggered, CELL_XLOW, {YDirectionType::Standard, ZDirectionType::Standard}};

  EXPECT_FALSE(areFieldsCompatible(field, field2));
  EXPECT_EQ(field.getMesh(), field2.getMesh());
  EXPECT_NE(field.getLocation(), field2.getLocation());
  EXPECT_EQ(field.getDirectionY(), field2.getDirectionY());
  EXPECT_EQ(field.getDirectionZ(), field2.getDirectionZ());
}

TEST_F(FieldTest, AreFieldsCompatibleFalseDirectionY) {
  // Create a field with default members
  Field field;

  // Create a field with all members set explicitly, and a non-default mesh
  Field field2{mesh, CELL_CENTRE, {YDirectionType::Aligned, ZDirectionType::Standard}};

  EXPECT_FALSE(areFieldsCompatible(field, field2));
  EXPECT_EQ(field.getMesh(), field2.getMesh());
  EXPECT_EQ(field.getLocation(), field2.getLocation());
  EXPECT_NE(field.getDirectionY(), field2.getDirectionY());
  EXPECT_EQ(field.getDirectionZ(), field2.getDirectionZ());
}

TEST_F(FieldTest, AreFieldsCompatibleTrueZAverage) {
  // Create a field with default members
  Field field;

  // Create a field with all members set explicitly, and a non-default mesh
  Field field2{mesh, CELL_CENTRE, {YDirectionType::Standard, ZDirectionType::Average}};

  EXPECT_TRUE(areFieldsCompatible(field, field2));
  EXPECT_EQ(field.getMesh(), field2.getMesh());
  EXPECT_EQ(field.getLocation(), field2.getLocation());
  EXPECT_EQ(field.getDirectionY(), field2.getDirectionY());
  EXPECT_NE(field.getDirectionZ(), field2.getDirectionZ());
}

TEST_F(FieldTest, AreFieldsCompatibleTrueYAlignedZAverage) {
  // Create a field with y aligned
  Field field{mesh, CELL_CENTRE, {YDirectionType::Aligned, ZDirectionType::Standard}};

  // Create a field with all members set explicitly, and a non-default mesh
  Field field2{mesh, CELL_CENTRE, {YDirectionType::Standard, ZDirectionType::Average}};

  EXPECT_TRUE(areFieldsCompatible(field, field2));
  EXPECT_EQ(field.getMesh(), field2.getMesh());
  EXPECT_EQ(field.getLocation(), field2.getLocation());
  EXPECT_NE(field.getDirectionY(), field2.getDirectionY());
  EXPECT_NE(field.getDirectionZ(), field2.getDirectionZ());
}

TEST_F(FieldTest, AreFieldsCompatibleFalseYAlignedZAverage2) {
  // Create a field with default members
  Field field;

  // Create a field with all members set explicitly, and a non-default mesh
  // Note it doesn't make sense for a field to be y-aligned and z-average,
  // because for a z-average field there is no difference between y-standard
  // and y-aligned.
  Field field2{mesh, CELL_CENTRE, {YDirectionType::Aligned, ZDirectionType::Average}};

  EXPECT_FALSE(areFieldsCompatible(field, field2));
  EXPECT_EQ(field.getMesh(), field2.getMesh());
  EXPECT_EQ(field.getLocation(), field2.getLocation());
  EXPECT_NE(field.getDirectionY(), field2.getDirectionY());
  EXPECT_NE(field.getDirectionZ(), field2.getDirectionZ());
}

TEST_F(FieldTest, filledFromAuto) {

  Field3D f;
  Field3D result = filledFrom(f, [](auto i) {
                                   return i.x() * i.y();
                                 });

  // Note: Serial so compiles with OpenMP
  BOUT_FOR_SERIAL(i, result.getRegion("RGN_ALL")) {
    ASSERT_DOUBLE_EQ( result[i], i.x() * i.y());
  }
}

TEST_F(FieldTest, filledFromInd3D) {
  Field3D f;
  Field3D result = filledFrom(f, [](Ind3D& i) {
                                   return i.x() + i.z() - 2*i.y();
                                 });

  BOUT_FOR_SERIAL(i, result.getRegion("RGN_ALL")) {
    ASSERT_DOUBLE_EQ( result[i], i.x() + i.z() - 2*i.y());
  }
}

TEST_F(FieldTest, filledFromConstInd3D) {
  Field3D f;
  Field3D result = filledFrom(f, [](const Ind3D& i) {
                                   return i.x() * i.y();
                                 });

  BOUT_FOR_SERIAL(i, result.getRegion("RGN_ALL")) {
    ASSERT_DOUBLE_EQ( result[i], i.x() * i.y());
  }
}

TEST_F(FieldTest, IsUniformTrue3D) {
  Field3D uniform{1.0};
  EXPECT_TRUE(isUniform(uniform));
}

TEST_F(FieldTest, IsUniformFalse3D) {
  Field3D nonuniform{2.0};
  nonuniform(1, 2, 3) = 3.0;
  EXPECT_FALSE(isUniform(nonuniform));
}

TEST_F(FieldTest, IsUniformTrue2D) {
  Field2D uniform{4.0};
  EXPECT_TRUE(isUniform(uniform));
}

TEST_F(FieldTest, IsUniformFalse2D) {
  Field2D nonuniform{5.0};
  nonuniform(2, 1) = 6.0;
  EXPECT_FALSE(isUniform(nonuniform));
}

TEST_F(FieldTest, GetUniformTrue3D) {
  Field3D uniform{10.0};
  EXPECT_EQ(getUniform(uniform), 10.0);
}

TEST_F(FieldTest, GetUniformFalse3D) {
  Field3D nonuniform{20.0};
  nonuniform(1, 2, 3) = 30.0;
#if CHECK > 1
  EXPECT_THROW(getUniform(nonuniform), BoutException);
#endif
}

TEST_F(FieldTest, GetUniformTrue2D) {
  Field2D uniform{40.0};
  EXPECT_EQ(getUniform(uniform), 40.0);
}

TEST_F(FieldTest, GetUniformFalse2D) {
  Field2D nonuniform{50.0};
  nonuniform(2, 1) = 60.0;
#if CHECK > 1
  EXPECT_THROW(getUniform(nonuniform), BoutException);
#endif
}
