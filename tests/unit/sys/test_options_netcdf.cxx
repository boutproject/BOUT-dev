// Test reading and writing to NetCDF

#include "bout/build_defines.hxx"

#if BOUT_HAS_NETCDF && !BOUT_HAS_LEGACY_NETCDF

#include "gtest/gtest.h"

#include "test_extras.hxx"
#include "test_tmpfiles.hxx"
#include "bout/field3d.hxx"
#include "bout/mesh.hxx"
#include "bout/options_io.hxx"

using bout::OptionsIO;

#include <string>

#include "fake_mesh_fixture.hxx"

// Reuse the "standard" fixture for FakeMesh
class OptionsNetCDFTest : public FakeMeshFixture {
public:
  // A temporary filename
  bout::testing::TempFile filename;
  WithQuietOutput quiet{output_info};
};

TEST_F(OptionsNetCDFTest, ReadWriteInt) {
  {
    Options options;
    options["test"] = 42;

    // Write the file
    OptionsIO::create(filename)->write(options);
  }

  // Read again
  Options data = OptionsIO::create(filename)->read();

  EXPECT_EQ(data["test"], 42);
}

TEST_F(OptionsNetCDFTest, ReadWriteString) {
  {
    Options options;
    options["test"] = std::string{"hello"};

    // Write file
    OptionsIO::create(filename)->write(options);
  }

  // Read file
  Options data = OptionsIO::create(filename)->read();

  EXPECT_EQ(data["test"], std::string("hello"));
}

TEST_F(OptionsNetCDFTest, ReadWriteField2D) {
  {
    Options options;
    options["test"] = Field2D(1.0);

    // Write file
    OptionsIO::create(filename)->write(options);
  }

  // Read file
  Options data = OptionsIO::create(filename)->read();

  Field2D value = data["test"].as<Field2D>(bout::globals::mesh);

  EXPECT_DOUBLE_EQ(value(0, 1), 1.0);
  EXPECT_DOUBLE_EQ(value(1, 0), 1.0);
}

TEST_F(OptionsNetCDFTest, ReadWriteField3D) {
  {
    Options options;
    options["test"] = Field3D(2.4);

    // Write file
    OptionsIO::create(filename)->write(options);
  }

  // Read file
  Options data = OptionsIO::create(filename)->read();

  Field3D value = data["test"].as<Field3D>(bout::globals::mesh);

  EXPECT_DOUBLE_EQ(value(0, 1, 0), 2.4);
  EXPECT_DOUBLE_EQ(value(1, 0, 1), 2.4);
  EXPECT_DOUBLE_EQ(value(1, 1, 1), 2.4);
}

TEST_F(OptionsNetCDFTest, Groups) {
  {
    Options options;
    options["test"]["key"] = 42;

    // Write file
    OptionsIO::create(filename)->write(options);
  }

  // Read file
  Options data = OptionsIO::create(filename)->read();
  EXPECT_EQ(data["test"]["key"], 42);
}

TEST_F(OptionsNetCDFTest, AttributeInt) {
  {
    Options options;
    options["test"] = 3;
    options["test"].attributes["thing"] = 4;

    // Write file
    OptionsIO::create(filename)->write(options);
  }

  // Read file
  Options data = OptionsIO::create(filename)->read();
  EXPECT_EQ(data["test"].attributes["thing"].as<int>(), 4);
}

TEST_F(OptionsNetCDFTest, AttributeBoutReal) {
  {
    Options options;
    options["test"] = 3;
    options["test"].attributes["thing"] = 3.14;

    // Write file
    OptionsIO::create(filename)->write(options);
  }

  // Read file
  Options data = OptionsIO::create(filename)->read();
  EXPECT_DOUBLE_EQ(data["test"].attributes["thing"].as<BoutReal>(), 3.14);
}

TEST_F(OptionsNetCDFTest, AttributeString) {
  {
    Options options;
    options["test"] = 3;
    options["test"].attributes["thing"] = "hello";

    // Write file
    OptionsIO::create(filename)->write(options);
  }

  // Read file
  Options data = OptionsIO::create(filename)->read();
  EXPECT_EQ(data["test"].attributes["thing"].as<std::string>(), "hello");
}

TEST_F(OptionsNetCDFTest, Field2DWriteCellCentre) {
  {
    Options options;
    options["f2d"] = Field2D(2.0);

    // Write file
    OptionsIO::create(filename)->write(options);
  }

  // Read file
  Options data = OptionsIO::create(filename)->read();

  EXPECT_EQ(data["f2d"].attributes["cell_location"].as<std::string>(),
            toString(CELL_CENTRE));
}

TEST_F(OptionsNetCDFTest, Field2DWriteCellYLow) {
  {
    Options options;
    options["f2d"] = Field2D(2.0, mesh_staggered).setLocation(CELL_YLOW);

    // Write file
    OptionsIO::create(filename)->write(options);
  }

  // Read file
  Options data = OptionsIO::create(filename)->read();

  EXPECT_EQ(data["f2d"].attributes["cell_location"].as<std::string>(),
            toString(CELL_YLOW));
}

TEST_F(OptionsNetCDFTest, Field3DWriteCellCentre) {
  {
    Options options;
    options["f3d"] = Field3D(2.0);

    // Write file
    OptionsIO::create(filename)->write(options);
  }

  // Read file
  Options data = OptionsIO::create(filename)->read();

  EXPECT_EQ(data["f3d"].attributes["cell_location"].as<std::string>(),
            toString(CELL_CENTRE));
}

TEST_F(OptionsNetCDFTest, Field3DWriteCellYLow) {
  {
    Options options;
    options["f3d"] = Field3D(2.0, mesh_staggered).setLocation(CELL_YLOW);

    // Write file
    OptionsIO::create(filename)->write(options);
  }

  // Read file
  Options data = OptionsIO::create(filename)->read();

  EXPECT_EQ(data["f3d"].attributes["cell_location"].as<std::string>(),
            toString(CELL_YLOW));
}

TEST_F(OptionsNetCDFTest, FieldPerpWriteCellCentre) {
  {
    Options options;
    FieldPerp fperp(3.0);
    fperp.setIndex(2);
    options["fperp"] = fperp;

    // Ensure MPI is initialised, otherwise we end up creating threads while
    // the file is open, and the lock is not removed on closing.
    fperp.getMesh()->getXcomm();

    // Write file
    OptionsIO::create(filename)->write(options);
  }

  // Read file
  Options data = OptionsIO::create(filename)->read();

  EXPECT_EQ(data["fperp"].attributes["cell_location"].as<std::string>(),
            toString(CELL_CENTRE));
  EXPECT_EQ(data["fperp"].attributes["yindex_global"].as<int>(), 2);
}

TEST_F(OptionsNetCDFTest, VerifyTimesteps) {
  {
    Options options;
    options["thing1"] = 1.0;
    options["thing1"].attributes["time_dimension"] = "t";

    OptionsIO::create(filename)->write(options);
  }

  EXPECT_NO_THROW(OptionsIO::create(filename)->verifyTimesteps());

  {
    Options options;
    options["thing1"] = 2.0;
    options["thing1"].attributes["time_dimension"] = "t";

    options["thing2"] = 3.0;
    options["thing2"].attributes["time_dimension"] = "t";

    OptionsIO::create({{"type", "netcdf"}, {"file", filename.string()}, {"append", true}})
        ->write(options);
  }

  EXPECT_THROW(OptionsIO::create(filename)->verifyTimesteps(), BoutException);
}

TEST_F(OptionsNetCDFTest, WriteTimeDimension) {
  {
    Options options;
    options["thing1"].assignRepeat(1.0);       // default time dim
    options["thing2"].assignRepeat(2.0, "t2"); // non-default

    // Only write non-default time dim
    OptionsIO::create(filename)->write(options, "t2");
  }

  Options data = OptionsIO::create(filename)->read();

  EXPECT_FALSE(data.isSet("thing1"));
  EXPECT_TRUE(data.isSet("thing2"));
}

TEST_F(OptionsNetCDFTest, WriteMultipleTimeDimensions) {
  {
    Options options;
    options["thing1_t1"].assignRepeat(1.0); // default time dim
    options["thing2_t1"].assignRepeat(1.0); // default time dim

    options["thing3_t2"].assignRepeat(2.0, "t2"); // non-default
    options["thing4_t2"].assignRepeat(2.0, "t2"); // non-default

    // Write the non-default time dim twice
    OptionsIO::create(filename)->write(options, "t2");
    OptionsIO::create(filename)->write(options, "t2");
    OptionsIO::create(filename)->write(options, "t");
  }

  EXPECT_NO_THROW(OptionsIO::create(filename)->verifyTimesteps());
}

#endif // BOUT_HAS_NETCDF
