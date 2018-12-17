
#include "bout.hxx"

#include "options_netcdf.hxx"
#include "optionsreader.hxx"

int main(int argc, char** argv) {
  BoutInitialise(argc, argv);

  // Read values from a NetCDF file
  OptionsNetCDF file("test.nc");

  auto values = file.read();

  values.printUnused();

  // Write to an INI text file
  OptionsReader *reader = OptionsReader::getInstance();
  reader->write(&values, "test-out.ini");

  // Write to a NetCDF file
  OptionsNetCDF("test-out.nc").write(values);

  ///////////////////////////

  Options::root()["f2d"] = Field2D(1.0);
  Options::root()["f3d"] = Field3D(2.0);

  Options::root()["time_test"] = 1.0;
  Options::root()["time_test"].attributes["time_dimension"] = std::string("t");
  
  // Write the BOUT.inp settings to NetCDF file
  OptionsNetCDF("settings.nc").write(Options::root());

  // Read back in
  auto settings = OptionsNetCDF("settings.nc").read();

  // Write to INI file
  reader->write(&settings, "settings.ini");
  
  BoutFinalise();
};
