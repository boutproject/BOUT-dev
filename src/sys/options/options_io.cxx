#include "bout/options_io.hxx"
#include "bout/globals.hxx"
#include "bout/mesh.hxx"

#include "options_adios.hxx"
#include "options_netcdf.hxx"

namespace bout {
std::unique_ptr<OptionsIO> OptionsIO::create(const std::string& file) {
  return OptionsIOFactory::getInstance().createFile(file);
}

std::unique_ptr<OptionsIO> OptionsIO::create(Options& config) {
  auto& factory = OptionsIOFactory::getInstance();
  return factory.create(factory.getType(&config), config);
}

OptionsIOFactory::ReturnType OptionsIOFactory::createRestart(Options* optionsptr) const {
  Options& options = optionsptr ? *optionsptr : Options::root()["restart_files"];

  // Set defaults
  options["path"].overrideDefault(
      Options::root()["datadir"].withDefault<std::string>("data"));
  options["prefix"].overrideDefault("BOUT.restart");
  options["append"].overrideDefault(false);
  options["singleWriteFile"].overrideDefault(true);
  return create(getType(&options), options);
}

OptionsIOFactory::ReturnType OptionsIOFactory::createOutput(Options* optionsptr) const {
  Options& options = optionsptr ? *optionsptr : Options::root()["output"];

  // Set defaults
  options["path"].overrideDefault(
      Options::root()["datadir"].withDefault<std::string>("data"));
  options["prefix"].overrideDefault("BOUT.dmp");
  options["append"].overrideDefault(Options::root()["append"]
                                        .doc("Add output data to existing (dump) files?")
                                        .withDefault<bool>(false));
  return create(getType(&options), options);
}

OptionsIOFactory::ReturnType OptionsIOFactory::createFile(const std::string& file) const {
  Options options{{"file", file}};
  return create(getDefaultType(), options);
}

void writeDefaultOutputFile(Options& data) {
  // Add mesh information
  bout::globals::mesh->outputVars(data);
  // Write to the default output file
  OptionsIOFactory::getInstance().createOutput()->write(data);
}

} // namespace bout
