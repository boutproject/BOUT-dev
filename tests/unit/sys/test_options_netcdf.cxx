// Test reading and writing to NetCDF

#ifdef NCDF4

#include "gtest/gtest.h"

#include "bout/mesh.hxx"
#include "field3d.hxx"
#include "test_extras.hxx"
#include "options_netcdf.hxx"

/// Global mesh
extern Mesh *mesh;

// Reuse the "standard" fixture for FakeMesh
using OptionsNetCDFTest = FakeMeshFixture;



TEST_F(OptionsNetCDFTest, ReadWriteInt) {
  // Temporary file
  OptionsNetCDF file(std::tmpnam(nullptr));
  
  Options options;
  options["test"] = 42;

  // Write the file
  file.write(options);

  // Read again
  Options data = file.read();

  EXPECT_EQ(data["test"], 42);
}

TEST_F(OptionsNetCDFTest, ReadWriteString) {
  std::string filename = std::tmpnam(nullptr);

  Options options;
  options["test"] = "hello";

  // Write file
  OptionsNetCDF(filename).write(options);

  // Read file
  Options data = OptionsNetCDF(filename).read();

  EXPECT_EQ(data["test"], std::string("hello"));
}

TEST_F(OptionsNetCDFTest, ReadWriteField2D) {
  std::string filename = std::tmpnam(nullptr);

  Options options;
  options["test"] = Field2D(1.0);

  // Write file
  OptionsNetCDF(filename).write(options);

  // Read file
  Options data = OptionsNetCDF(filename).read();

  Field2D value = data["test"].as<Field2D>(mesh);
  
  EXPECT_DOUBLE_EQ(value(0,1), 1.0);
  EXPECT_DOUBLE_EQ(value(1,0), 1.0);
}

TEST_F(OptionsNetCDFTest, ReadWriteField3D) {
  std::string filename = std::tmpnam(nullptr);
  
  Options options;
  options["test"] = Field3D(2.4);
  
  // Write file
  OptionsNetCDF(filename).write(options);

  // Read file
  Options data = OptionsNetCDF(filename).read();

  Field3D value = data["test"].as<Field3D>(mesh);
  
  EXPECT_DOUBLE_EQ(value(0,1,0), 2.4);
  EXPECT_DOUBLE_EQ(value(1,0,1), 2.4);
  EXPECT_DOUBLE_EQ(value(1,1,1), 2.4);
}

TEST_F(OptionsNetCDFTest, AttributeInt) {
  std::string filename = std::tmpnam(nullptr);
  
  Options options;
  options["test"] = 3;
  options["test"].attributes["thing"] = 4;

  // Write file
  OptionsNetCDF(filename).write(options);

  // Read file
  Options data = OptionsNetCDF(filename).read();
  EXPECT_EQ(options["test"].attributes["thing"].as<int>(), 4);
}

TEST_F(OptionsNetCDFTest, AttributeBoutReal) {
  std::string filename = std::tmpnam(nullptr);
  
  Options options;
  options["test"] = 3;
  options["test"].attributes["thing"] = 3.14;

  // Write file
  OptionsNetCDF(filename).write(options);

  // Read file
  Options data = OptionsNetCDF(filename).read();
  EXPECT_DOUBLE_EQ(options["test"].attributes["thing"].as<BoutReal>(), 3.14);
}

TEST_F(OptionsNetCDFTest, AttributeString) {
  std::string filename = std::tmpnam(nullptr);
  
  Options options;
  options["test"] = 3;
  options["test"].attributes["thing"] = "hello";

  // Write file
  OptionsNetCDF(filename).write(options);

  // Read file
  Options data = OptionsNetCDF(filename).read();
  EXPECT_EQ(options["test"].attributes["thing"].as<std::string>(), "hello");
}


#endif // NCDF4
