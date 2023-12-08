
#include "bout/bout.hxx"

#include "bout/options_adios.hxx"
#include "bout/options_io.hxx"
#include "bout/options_netcdf.hxx"
#include "bout/optionsreader.hxx"

using bout::OptionsADIOS;

int main(int argc, char** argv) {
  BoutInitialise(argc, argv);

  // Read values from a NetCDF file
  bout::OptionsNetCDF file("test.nc");

  auto values = file.read();

  values.printUnused();

  // Write to an INI text file
  OptionsReader* reader = OptionsReader::getInstance();
  reader->write(&values, "test-out.ini");

  // Write to ADIOS file
  OptionsADIOS("test-out.bp", bout::OptionsIO::FileMode::replace, true).write(values);

  ///////////////////////////

  // Write the BOUT.inp settings to NetCDF file
  OptionsADIOS("settings.bp", bout::OptionsIO::FileMode::replace, true)
      .write(Options::root());

  // Read back in
  auto settings = OptionsADIOS("settings.bp").read();

  // Write to INI file
  reader->write(&settings, "settings.ini");

  ///////////////////////////
  // Write fields

  Options fields;
  fields["f2d"] = Field2D(1.0);
  fields["f3d"] = Field3D(2.0);
  fields["fperp"] = FieldPerp(3.0);
  auto f = OptionsADIOS("fields.bp");
  /*
     write() for adios only buffers data but does not guarantee writing to disk
     unless singleWriteFile is set to true 
     */
  f.write(fields);
  // indicate completion of step, required to get data on disk
  f.verifyTimesteps();

  ///////////////////////////
  // Read fields

  Options fields_in = OptionsADIOS("fields.bp").read();

  auto f2d = fields_in["f2d"].as<Field2D>(bout::globals::mesh);
  auto f3d = fields_in["f3d"].as<Field3D>(bout::globals::mesh);
  auto fperp = fields_in["fperp"].as<FieldPerp>(bout::globals::mesh);

  Options fields2;
  fields2["f2d"] = f2d;
  fields2["f3d"] = f3d;
  fields2["fperp"] = fperp;

  // Write out again
  auto f2 = bout::OptionsIOFactory("fields2.bp", bout::OptionsIO::FileMode::replace,
                                   bout::OptionsIO::Library::ADIOS, true);
  f2->write(fields2);

  ///////////////////////////
  // Time dependent values

  Options data;
  data["scalar"] = 1.0;
  data["scalar"].attributes["time_dimension"] = "t";

  data["field"] = Field3D(2.0);
  data["field"].attributes["time_dimension"] = "t";

  OptionsADIOS("time.bp", bout::OptionsIO::FileMode::replace, true).write(data);

  // Update time-dependent values
  data["scalar"] = 2.0;
  data["field"] = Field3D(3.0);

  // Append data to file
  OptionsADIOS("time.bp", bout::OptionsIO::FileMode::append, true).write(data);

  BoutFinalise();
};