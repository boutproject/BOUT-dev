
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
  reader->write(&values, "test.settings");
  
  BoutFinalise();
};
