#include "bout/options_io.hxx"
#include "bout/bout.hxx"
#include "bout/globals.hxx"
#include "bout/mesh.hxx"
#include "bout/options.hxx"
#include "bout/version.hxx"

#include "options_adios.hxx"
#include "options_netcdf.hxx"

#include <memory>
#include <string>

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
  options["replace"].overrideDefault(true); // Restart files are overwritten each output
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
  options["replace"].overrideDefault(
      Options::root()["replace"]
          .doc("Replace output (dump) files if they exist and not appending?")
          .withDefault<bool>(false));
  return create(getType(&options), options);
}

OptionsIOFactory::ReturnType OptionsIOFactory::createFile(const std::string& file) const {
  Options options{{"file", file}};
  return create(getDefaultType(), options);
}

void writeDefaultOutputFile(Options& data) {
  // Add BOUT++ version and flags
  bout::experimental::addBuildFlagsToOptions(data);
  // Add mesh information
  bout::globals::mesh->outputVars(data);
  // Write to the default output file
  OptionsIOFactory::getInstance().createOutput()->write(data);
}

void OptionsIO::write(const std::string& prefix, Options& data, Mesh* mesh) {
  Options file_options = {{"prefix", prefix}};
  data["BOUT_VERSION"].force(bout::version::as_double);
  if (mesh != nullptr) {
    mesh->outputVars(data);
  }
  OptionsIOFactory::getInstance().createOutput(&file_options)->write(data);
}
} // namespace bout
