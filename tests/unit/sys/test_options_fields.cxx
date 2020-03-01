// Test the Options class methods which involve Field2D/3D
// These need a Mesh fixture, unlike the tests with scalars

#include "gtest/gtest.h"

#include "bout/mesh.hxx"
#include "field3d.hxx"
#include "test_extras.hxx"
#include "unused.hxx"

/// Global mesh
namespace bout {
namespace globals {
extern Mesh* mesh;
}
} // namespace bout

// Reuse the "standard" fixture for FakeMesh
class OptionsFieldTest : public FakeMeshFixture {
  WithQuietOutput quiet_warn{output_warn};
};

TEST_F(OptionsFieldTest, StoreField3D) {
  Field3D field = 1.0;
  Options options;

  EXPECT_FALSE(options.isValue());
  
  options = field;

  EXPECT_TRUE(options.isValue());
}

TEST_F(OptionsFieldTest, StoreField2D) {
  Field2D field = 1.0;
  Options options;

  EXPECT_FALSE(options.isValue());
  
  options = field;

  EXPECT_TRUE(options.isValue());
}

TEST_F(OptionsFieldTest, RetrieveField3D) {
  Field3D field = 1.0;
  field(0,1,1) = 2.0;
  
  Options options;
  options = field;

  Field3D other = options;

  EXPECT_DOUBLE_EQ(other(0,1,0), 1.0);
  EXPECT_DOUBLE_EQ(other(0,1,1), 2.0);
}

TEST_F(OptionsFieldTest, RetrieveField3DfromBoutReal) {
  Options options;
  options = 1.2;

  Field3D other = options;

  EXPECT_DOUBLE_EQ(other(0,1,0), 1.2);
  EXPECT_DOUBLE_EQ(other(0,0,1), 1.2);
}

TEST_F(OptionsFieldTest, RetrieveBoutRealfromField3D) {
  Options options;
  Field3D field = 1.2;
  options = field;

  EXPECT_THROW(BoutReal UNUSED(value) = options, BoutException);
}

TEST_F(OptionsFieldTest, RetrieveField2DfromField3D) {
  Options options;
  Field3D field = 1.2;
  options = field;

  EXPECT_THROW(Field2D value = options.as<Field2D>(), BoutException);
}

TEST_F(OptionsFieldTest, RetrieveField3DfromField2D) {
  Options options;
  Field2D field = 1.2;
  options = field;

  Field3D value = options.as<Field3D>();
  EXPECT_DOUBLE_EQ(value(0,1,0), 1.2);
}

TEST_F(OptionsFieldTest, RetrieveStringfromField3D) {
  Options options;
  Field3D field = 1.2;
  options = field;

  WithQuietOutput quiet{output_info};
  EXPECT_EQ(options.as<std::string>(), "<Field3D>");
}

TEST_F(OptionsFieldTest, RetrieveStringfromField2D) {
  Options options;
  Field2D field = 1.2;
  options = field;

  WithQuietOutput quiet{output_info};
  EXPECT_EQ(options.as<std::string>(), "<Field2D>");
}

TEST_F(OptionsFieldTest, RetrieveField3DfromString) {
  Options options;
  options = "1 + 2";

  WithQuietOutput quiet{output_info};

  Field3D other = options.as<Field3D>();

  EXPECT_DOUBLE_EQ(other(0,1,0), 3.0);
  EXPECT_DOUBLE_EQ(other(0,0,1), 3.0);
}

TEST_F(OptionsFieldTest, RetrieveField2DfromString) {
  Options options;
  options = "1 + 2";

  Field2D other = options.as<Field2D>();

  EXPECT_DOUBLE_EQ(other(0,1,0), 3.0);
  EXPECT_DOUBLE_EQ(other(0,0,1), 3.0);
}

TEST_F(OptionsFieldTest, RetrieveField2DfromBadString) {
  Options options;
  options = "1 + ";

  EXPECT_THROW(Field2D other = options.as<Field2D>(), ParseException);
}

