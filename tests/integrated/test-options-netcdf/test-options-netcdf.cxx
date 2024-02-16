
#include "bout/bout.hxx"

#include "bout/options_io.hxx"
#include "bout/optionsreader.hxx"

using bout::OptionsIO;

int main(int argc, char** argv) {
  BoutInitialise(argc, argv);

  // Read values from a NetCDF file
  auto file = OptionsIO::create("test.nc");

  auto values = file->read();

  values.printUnused();

  // Write to an INI text file
  OptionsReader* reader = OptionsReader::getInstance();
  reader->write(&values, "test-out.ini");

  // Write to a NetCDF file
  OptionsIO::create("test-out.nc")->write(values);

  ///////////////////////////

  // Write the BOUT.inp settings to NetCDF file
  OptionsIO::create("settings.nc")->write(Options::root());

  // Read back in
  auto settings = OptionsIO::create("settings.nc")->read();

  // Write to INI file
  reader->write(&settings, "settings.ini");

  ///////////////////////////
  // Write fields

  Options fields;
  fields["f2d"] = Field2D(1.0);
  fields["f3d"] = Field3D(2.0);
  fields["fperp"] = FieldPerp(3.0);
  OptionsIO::create("fields.nc")->write(fields);

  ///////////////////////////
  // Read fields

  Options fields_in = OptionsIO::create("fields.nc")->read();

  auto f2d = fields_in["f2d"].as<Field2D>(bout::globals::mesh);
  auto f3d = fields_in["f3d"].as<Field3D>(bout::globals::mesh);
  auto fperp = fields_in["fperp"].as<FieldPerp>(bout::globals::mesh);

  Options fields2;
  fields2["f2d"] = f2d;
  fields2["f3d"] = f3d;
  fields2["fperp"] = fperp;

  // Write out again
  OptionsIO::create("fields2.nc")->write(fields2);

  ///////////////////////////
  // Time dependent values

  Options data;
  data["scalar"] = 1.0;
  data["scalar"].attributes["time_dimension"] = "t";

  data["field"] = Field3D(2.0);
  data["field"].attributes["time_dimension"] = "t";

  OptionsIO::create("time.nc")->write(data);

  // Update time-dependent values
  data["scalar"] = 2.0;
  data["field"] = Field3D(3.0);

  // Append data to file
  OptionsIO::create({{"file", "time.nc"}, {"append", true}})->write(data);

  BoutFinalise();
};
