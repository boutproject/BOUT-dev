
#include "bout/bout.hxx"

#include "bout/options_io.hxx"
#include "bout/optionsreader.hxx"

using bout::OptionsIO;

int main(int argc, char** argv) {
  BoutInitialise(argc, argv);

  // Read values from a NetCDF file
  auto file = bout::OptionsIO::create("test.nc");

  auto values = file->read();

  values.printUnused();

  // Write to an INI text file
  OptionsReader* reader = OptionsReader::getInstance();
  reader->write(&values, "test-out.ini");

  // Write to ADIOS file, by setting file type "adios"
  OptionsIO::create({{"file", "test-out.bp"},
                     {"type", "adios"},
                     {"append", false},
                     {"singleWriteFile", true}})
      ->write(values);

  ///////////////////////////

  // Write the BOUT.inp settings to ADIOS file
  OptionsIO::create({{"file", "settings.bp"},
                     {"type", "adios"},
                     {"append", false},
                     {"singleWriteFile", true}})
      ->write(Options::root());

  // Read back in
  auto settings = OptionsIO::create({{"file", "settings.bp"}, {"type", "adios"}})->read();

  // Write to INI file
  reader->write(&settings, "settings.ini");

  ///////////////////////////
  // Write fields

  Options fields;
  fields["f2d"] = Field2D(1.0);
  fields["f3d"] = Field3D(2.0);
  fields["fperp"] = FieldPerp(3.0);
  auto f = OptionsIO::create({{"file", "fields.bp"}, {"type", "adios"}});
  /*
     write() for adios only buffers data but does not guarantee writing to disk
     unless singleWriteFile is set to true 
     */
  f->write(fields);
  // indicate completion of step, required to get data on disk
  f->verifyTimesteps();

  ///////////////////////////
  // Read fields

  Options fields_in =
      OptionsIO::create({{"file", "fields.bp"}, {"type", "adios"}})->read();

  auto f2d = fields_in["f2d"].as<Field2D>(bout::globals::mesh);
  auto f3d = fields_in["f3d"].as<Field3D>(bout::globals::mesh);
  auto fperp = fields_in["fperp"].as<FieldPerp>(bout::globals::mesh);

  Options fields2;
  fields2["f2d"] = f2d;
  fields2["f3d"] = f3d;
  fields2["fperp"] = fperp;

  // Write out again
  auto f2 = OptionsIO::create({{"file", "fields2.bp"},
                               {"type", "adios"},
                               {"append", false},
                               {"singleWriteFile", true}});
  f2->write(fields2);

  ///////////////////////////
  // Time dependent values

  Options data;
  data["scalar"] = 1.0;
  data["scalar"].attributes["time_dimension"] = "t";

  data["field"] = Field3D(2.0);
  data["field"].attributes["time_dimension"] = "t";

  OptionsIO::create({{"file", "time.bp"},
                     {"type", "adios"},
                     {"append", false},
                     {"singleWriteFile", true}})
      ->write(data);

  // Update time-dependent values
  data["scalar"] = 2.0;
  data["field"] = Field3D(3.0);

  // Append data to file
  OptionsIO::create({{"file", "time.bp"},
                     {"type", "adios"},
                     {"append", true},
                     {"singleWriteFile", true}})
      ->write(data);

  BoutFinalise();
};
