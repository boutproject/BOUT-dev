// Test reading and writing to NetCDF

#ifdef NCDF4

#include "gtest/gtest.h"

#include "bout/mesh.hxx"
#include "field3d.hxx"
#include "test_extras.hxx"
#include "options_netcdf.hxx"

/// Global mesh
namespace bout {
namespace globals {
extern Mesh* mesh;
}
} // namespace bout

// Reuse the "standard" fixture for FakeMesh
using OptionsNetCDFTest = FakeMeshFixture;

TEST_F(OptionsNetCDFTest, ReadWriteInt) {
  // Temporary file
  OptionsNetCDF file(std::tmpnam(nullptr));

  {
    Options options;
    options["test"] = 42;
    
    // Write the file
    file.write(options);
  }

  // Read again
  Options data = file.read();

  EXPECT_EQ(data["test"], 42);
}

TEST_F(OptionsNetCDFTest, ReadWriteString) {
  std::string filename = std::tmpnam(nullptr);

  {
    Options options;
    options["test"] = "hello";

    // Write file
    OptionsNetCDF(filename).write(options);
  }
  
  // Read file
  Options data = OptionsNetCDF(filename).read();

  EXPECT_EQ(data["test"], std::string("hello"));
}

TEST_F(OptionsNetCDFTest, ReadWriteField2D) {
  std::string filename = std::tmpnam(nullptr);

  {
    Options options;
    options["test"] = Field2D(1.0);
    
    // Write file
    OptionsNetCDF(filename).write(options);
  }
  
  // Read file
  Options data = OptionsNetCDF(filename).read();

  Field2D value = data["test"].as<Field2D>(bout::globals::mesh);
  
  EXPECT_DOUBLE_EQ(value(0,1), 1.0);
  EXPECT_DOUBLE_EQ(value(1,0), 1.0);
}

TEST_F(OptionsNetCDFTest, ReadWriteField3D) {
  std::string filename = std::tmpnam(nullptr);

  {
    Options options;
    options["test"] = Field3D(2.4);
    
    // Write file
    OptionsNetCDF(filename).write(options);
  }
  
  // Read file
  Options data = OptionsNetCDF(filename).read();

  Field3D value = data["test"].as<Field3D>(bout::globals::mesh);
  
  EXPECT_DOUBLE_EQ(value(0,1,0), 2.4);
  EXPECT_DOUBLE_EQ(value(1,0,1), 2.4);
  EXPECT_DOUBLE_EQ(value(1,1,1), 2.4);
}

TEST_F(OptionsNetCDFTest, Groups) {
  std::string filename = std::tmpnam(nullptr);

  {
    Options options;
    options["test"]["key"] = 42;
    
    // Write file
    OptionsNetCDF(filename).write(options);
  }
  
  // Read file
  Options data = OptionsNetCDF(filename).read();
  EXPECT_EQ(data["test"]["key"], 42);
}

TEST_F(OptionsNetCDFTest, AttributeInt) {
  std::string filename = std::tmpnam(nullptr);

  {
    Options options;
    options["test"] = 3;
    options["test"].attributes["thing"] = 4;

    // Write file
    OptionsNetCDF(filename).write(options);
  }
  
  // Read file
  Options data = OptionsNetCDF(filename).read();
  EXPECT_EQ(data["test"].attributes["thing"].as<int>(), 4);
}

TEST_F(OptionsNetCDFTest, AttributeBoutReal) {
  std::string filename = std::tmpnam(nullptr);

  {
    Options options;
    options["test"] = 3;
    options["test"].attributes["thing"] = 3.14;
    
    // Write file
    OptionsNetCDF(filename).write(options);
  }
  
  // Read file
  Options data = OptionsNetCDF(filename).read();
  EXPECT_DOUBLE_EQ(data["test"].attributes["thing"].as<BoutReal>(), 3.14);
}

TEST_F(OptionsNetCDFTest, AttributeString) {
  std::string filename = std::tmpnam(nullptr);

  {
    Options options;
    options["test"] = 3;
    options["test"].attributes["thing"] = "hello";
    
    // Write file
    OptionsNetCDF(filename).write(options);
  }
  
  // Read file
  Options data = OptionsNetCDF(filename).read();
  EXPECT_EQ(data["test"].attributes["thing"].as<std::string>(), "hello");
}

TEST_F(OptionsNetCDFTest, Field2DWriteCellCentre) {
  std::string filename = std::tmpnam(nullptr);

  {
    Options options;
    options["f2d"] = Field2D(2.0);
    
    // Write file
    OptionsNetCDF(filename).write(options);
  }
  
  // Read file
  Options data = OptionsNetCDF(filename).read();

  EXPECT_EQ(data["f2d"].attributes["cell_location"].as<std::string>(), toString(CELL_CENTRE));
}

TEST_F(OptionsNetCDFTest, Field2DWriteCellYLow) {
  std::string filename = std::tmpnam(nullptr);

  {
    Options options;
    options["f2d"] = Field2D(2.0, mesh_staggered).setLocation(CELL_YLOW);
    
    // Write file
    OptionsNetCDF(filename).write(options);
  }
  
  // Read file
  Options data = OptionsNetCDF(filename).read();

  EXPECT_EQ(data["f2d"].attributes["cell_location"].as<std::string>(), toString(CELL_YLOW));
}

TEST_F(OptionsNetCDFTest, Field3DWriteCellCentre) {
  std::string filename = std::tmpnam(nullptr);

  {
    Options options;
    options["f3d"] = Field3D(2.0);
    
    // Write file
    OptionsNetCDF(filename).write(options);
  }
  
  // Read file
  Options data = OptionsNetCDF(filename).read();

  EXPECT_EQ(data["f3d"].attributes["cell_location"].as<std::string>(), toString(CELL_CENTRE));
}

TEST_F(OptionsNetCDFTest, Field3DWriteCellYLow) {
  std::string filename = std::tmpnam(nullptr);

  {
    Options options;
    options["f3d"] = Field3D(2.0, mesh_staggered).setLocation(CELL_YLOW);
    
    // Write file
    OptionsNetCDF(filename).write(options);
  }
  
  // Read file
  Options data = OptionsNetCDF(filename).read();

  EXPECT_EQ(data["f3d"].attributes["cell_location"].as<std::string>(), toString(CELL_YLOW));
}


#endif // NCDF4
