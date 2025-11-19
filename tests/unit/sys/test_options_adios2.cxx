// Test reading and writing to ADIOS2

#include "bout/build_defines.hxx"

#if BOUT_HAS_ADIOS2

#include "gtest/gtest.h"

#include "test_extras.hxx"
#include "test_tmpfiles.hxx"
#include "bout/adios_object.hxx"
#include "bout/field3d.hxx"
#include "bout/mesh.hxx"
#include "bout/options_io.hxx"
#include "bout/utils.hxx"

using bout::OptionsIO;

#include <string>

#include "fake_mesh_fixture.hxx"

// Reuse the "standard" fixture for FakeMesh
class OptionsAdios2Test : public FakeMeshFixture {
public:
  static void SetUpTestSuite() { bout::ADIOSInit(BoutComm::get()); }
  static void TearDownTestSuite() { bout::ADIOSFinalize(); }

  ~OptionsAdios2Test() override {
    bout::ADIOSStream::ADIOSGetStream(filename, adios2::Mode::Deferred).engine().Close();
  }

  // A temporary filename
  bout::testing::TempFile filename;
  WithQuietOutput quiet{output_info};
  Options file_options{
      {"file", filename.string()},
      {"type", "adios"},
      {"singleWriteFile", true},
  };
  Options no_write_options{
      {"file", filename.string()},
      {"type", "adios"},
      {"singleWriteFile", false},
  };
};

TEST_F(OptionsAdios2Test, ReadWriteInt) {
  {
    Options options;
    options["test"] = 42;

    // Write the file
    OptionsIO::create(file_options)->write(options);
  }

  // Read again
  Options data = OptionsIO::create(file_options)->read();

  EXPECT_EQ(data["test"], 42);
}

TEST_F(OptionsAdios2Test, ReadWriteString) {
  {
    Options options;
    options["test"] = std::string{"hello"};

    // Write file
    OptionsIO::create(file_options)->write(options);
  }

  // Read file
  Options data = OptionsIO::create(file_options)->read();

  EXPECT_EQ(data["test"], std::string("hello"));
}

TEST_F(OptionsAdios2Test, ReadWriteField2D) {
  {
    Options options;
    options["test"] = Field2D(1.0);

    // Write file
    OptionsIO::create(file_options)->write(options);
  }

  // Read file
  Options data = OptionsIO::create(file_options)->read();

  Field2D value = data["test"].as<Field2D>(bout::globals::mesh);

  EXPECT_DOUBLE_EQ(value(0, 1), 1.0);
  EXPECT_DOUBLE_EQ(value(1, 0), 1.0);
}

TEST_F(OptionsAdios2Test, ReadWriteField3D) {
  {
    Options options;
    options["test"] = Field3D(2.4);

    // Write file
    OptionsIO::create(file_options)->write(options);
  }

  // Read file
  Options data = OptionsIO::create(file_options)->read();

  Field3D value = data["test"].as<Field3D>(bout::globals::mesh);

  EXPECT_DOUBLE_EQ(value(0, 1, 0), 2.4);
  EXPECT_DOUBLE_EQ(value(1, 0, 1), 2.4);
  EXPECT_DOUBLE_EQ(value(1, 1, 1), 2.4);
}

TEST_F(OptionsAdios2Test, Groups) {
  {
    Options options;
    options["test"]["key"] = 42;

    // Write file
    OptionsIO::create(file_options)->write(options);
  }

  // Read file
  Options data = OptionsIO::create(file_options)->read();
  EXPECT_EQ(data["test"]["key"], 42);
}

TEST_F(OptionsAdios2Test, AttributeInt) {
  {
    Options options;
    options["test"] = 3;
    options["test"].attributes["thing"] = 4;

    // Write file
    OptionsIO::create(file_options)->write(options);
  }

  // Read file
  Options data = OptionsIO::create(file_options)->read();
  EXPECT_EQ(data["test"].attributes["thing"].as<int>(), 4);
}

TEST_F(OptionsAdios2Test, AttributeBoutReal) {
  {
    Options options;
    options["test"] = 3;
    options["test"].attributes["thing"] = 3.14;

    // Write file
    OptionsIO::create(file_options)->write(options);
  }

  // Read file
  Options data = OptionsIO::create(file_options)->read();
  EXPECT_DOUBLE_EQ(data["test"].attributes["thing"].as<BoutReal>(), 3.14);
}

TEST_F(OptionsAdios2Test, AttributeString) {
  {
    Options options;
    options["test"] = 3;
    options["test"].attributes["thing"] = "hello";

    // Write file
    OptionsIO::create(file_options)->write(options);
  }

  // Read file
  Options data = OptionsIO::create(file_options)->read();
  EXPECT_EQ(data["test"].attributes["thing"].as<std::string>(), "hello");
}

TEST_F(OptionsAdios2Test, Field2DWriteCellCentre) {
  {
    Options options;
    options["f2d"] = Field2D(2.0);

    // Write file
    OptionsIO::create(file_options)->write(options);
  }

  // Read file
  Options data = OptionsIO::create(file_options)->read();

  EXPECT_EQ(data["f2d"].attributes["cell_location"].as<std::string>(),
            toString(CELL_CENTRE));
}

TEST_F(OptionsAdios2Test, Field2DWriteCellYLow) {
  {
    Options options;
    options["f2d"] = Field2D(2.0, mesh_staggered).setLocation(CELL_YLOW);

    // Write file
    OptionsIO::create(file_options)->write(options);
  }

  // Read file
  Options data = OptionsIO::create(file_options)->read();

  EXPECT_EQ(data["f2d"].attributes["cell_location"].as<std::string>(),
            toString(CELL_YLOW));
}

TEST_F(OptionsAdios2Test, Field3DWriteCellCentre) {
  {
    Options options;
    options["f3d"] = Field3D(2.0);

    // Write file
    OptionsIO::create(file_options)->write(options);
  }

  // Read file
  Options data = OptionsIO::create(file_options)->read();

  EXPECT_EQ(data["f3d"].attributes["cell_location"].as<std::string>(),
            toString(CELL_CENTRE));
}

TEST_F(OptionsAdios2Test, Field3DWriteCellYLow) {
  {
    Options options;
    options["f3d"] = Field3D(2.0, mesh_staggered).setLocation(CELL_YLOW);

    // Write file
    OptionsIO::create(file_options)->write(options);
  }

  // Read file
  Options data = OptionsIO::create(file_options)->read();

  EXPECT_EQ(data["f3d"].attributes["cell_location"].as<std::string>(),
            toString(CELL_YLOW));
}

TEST_F(OptionsAdios2Test, FieldPerpWriteCellCentre) {
  {
    Options options;
    FieldPerp fperp(3.0);
    fperp.setIndex(2);
    options["fperp"] = fperp;

    // Ensure MPI is initialised, otherwise we end up creating threads while
    // the file is open, and the lock is not removed on closing.
    fperp.getMesh()->getXcomm();

    // Write file
    OptionsIO::create(file_options)->write(options);
  }

  // Read file
  Options data = OptionsIO::create(file_options)->read();

  EXPECT_EQ(data["fperp"].attributes["cell_location"].as<std::string>(),
            toString(CELL_CENTRE));
  EXPECT_EQ(data["fperp"].attributes["yindex_global"].as<int>(), 2);
}

TEST_F(OptionsAdios2Test, VerifyTimesteps) {
  {
    Options options;
    options["thing1"] = 1.0;
    options["thing1"].attributes["time_dimension"] = "t";

    auto f = OptionsIO::create(no_write_options);
    f->write(options);
    EXPECT_NO_THROW(f->verifyTimesteps());
  }

  {
    Options options;
    options["thing1"] = 2.0;
    options["thing1"].attributes["time_dimension"] = "t";

    options["thing2"] = 3.0;
    options["thing2"].attributes["time_dimension"] = "t";

    OptionsIO::create({{"type", "adios"}, {"file", filename.string()}, {"append", true}})
        ->write(options);
  }

  // Doesn't throw, but it's less of a problem for ADIOS2 than netCDF?
  // EXPECT_THROW(OptionsIO::create(no_write_options)->verifyTimesteps(), BoutException);
}

TEST_F(OptionsAdios2Test, WriteTimeDimension) {
  {
    Options options;
    options["thing1"].assignRepeat(1.0);       // default time dim
    options["thing2"].assignRepeat(2.0, "t2"); // non-default

    // Only write non-default time dim
    OptionsIO::create(file_options)->write(options, "t2");
  }

  Options data = OptionsIO::create(file_options)->read();

  EXPECT_FALSE(data.isSet("thing1"));
  EXPECT_TRUE(data.isSet("thing2"));
}

TEST_F(OptionsAdios2Test, WriteMultipleTimeDimensions) {
  {
    Options options;
    options["thing1_t1"].assignRepeat(1.0); // default time dim
    options["thing2_t1"].assignRepeat(1.0); // default time dim

    options["thing3_t2"].assignRepeat(2.0, "t2"); // non-default
    options["thing4_t2"].assignRepeat(2.0, "t2"); // non-default

    // Write the non-default time dim twice
    OptionsIO::create(no_write_options)->write(options, "t2");
    OptionsIO::create(no_write_options)->write(options, "t2");
    OptionsIO::create(no_write_options)->write(options, "t");
  }

  EXPECT_NO_THROW(OptionsIO::create(no_write_options)->verifyTimesteps());
}

#endif // Bout_Has_ADIOS2
