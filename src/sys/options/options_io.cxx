#include "bout/build_config.hxx"

#include "bout/options_io.hxx"

#include "bout/bout.hxx"
#include "bout/globals.hxx"
#include "bout/mesh.hxx"
#include "bout/sys/timer.hxx"

#include <exception>
#include <iostream>
#include <vector>

#include "bout/options_adios.hxx"
#include "bout/options_netcdf.hxx"

namespace bout {

OptionsIO::OptionsIO() {}

OptionsIO::OptionsIO(std::string filename, FileMode mode, bool singleWriteFile)
    : filename(std::move(filename)), file_mode(mode), singleWriteFile(singleWriteFile) {}

OptionsIO::~OptionsIO() = default;
OptionsIO::OptionsIO(OptionsIO&&) noexcept = default;
OptionsIO& OptionsIO::operator=(OptionsIO&&) noexcept = default;

OptionsIO::Library getIOLibrary(Options& options) {
  if (options["iolibrary"].isSet()) {
    // Solver-specific IO library
    std::string iolib = options["iolibrary"];
    std::transform(iolib.begin(), iolib.end(), iolib.begin(), ::tolower);
    if (iolib == "adios")
      return OptionsIO::Library::ADIOS;
    else if (iolib == "netcdf")
      return OptionsIO::Library::NetCDF;
    else
      return OptionsIO::Library::Invalid;
  } else {
    return OptionsIO::defaultIOLibrary;
  }
}

std::shared_ptr<OptionsIO> OptionsIOFactory(std::string filename,
                                            OptionsIO::FileMode mode,
                                            const OptionsIO::Library library,
                                            const bool singleWriteFile) {
  if (library == OptionsIO::Library::ADIOS) {
    return std::make_shared<OptionsADIOS>(OptionsADIOS(filename, mode, singleWriteFile));
  } else if (library == OptionsIO::Library::NetCDF) {
    return std::make_shared<OptionsNetCDF>(
        OptionsNetCDF(filename, mode, singleWriteFile));
  } else {
    return nullptr;
  }
}

std::string getRestartDirectoryName(Options& options) {
  if (options["restartdir"].isSet()) {
    // Solver-specific restart directory
    return options["restartdir"].withDefault<std::string>("data");
  }
  // Use the root data directory
  return options["datadir"].withDefault<std::string>("data");
}

std::string getRestartFilename(Options& options, const OptionsIO::Library library) {
  return getRestartFilename(options, BoutComm::rank(), library);
}

std::string getRestartFilename(Options& options, int rank,
                               const OptionsIO::Library library) {
  if (library == OptionsIO::Library::ADIOS)
    return fmt::format("{}/BOUT.restart.bp", bout::getRestartDirectoryName(options));
  else if (library == OptionsIO::Library::NetCDF)
    return fmt::format("{}/BOUT.restart.{}.nc", bout::getRestartDirectoryName(options),
                       rank);
  else
    return fmt::format("{}/BOUT.restart.{}.data", bout::getRestartDirectoryName(options),
                       rank);
}

std::string getOutputFilename(Options& options, const OptionsIO::Library library) {
  return getOutputFilename(options, BoutComm::rank(), library);
}

std::string getOutputFilename(Options& options, int rank,
                              const OptionsIO::Library library) {
  if (library == OptionsIO::Library::ADIOS)
    return fmt::format("{}/BOUT.dmp.bp",
                       options["datadir"].withDefault<std::string>("data"));
  else if (library == OptionsIO::Library::NetCDF)
    return fmt::format("{}/BOUT.dmp.{}.nc",
                       options["datadir"].withDefault<std::string>("data"), rank);
  else
    return fmt::format("{}/BOUT.dmp.{}.data",
                       options["datadir"].withDefault<std::string>("data"), rank);
}

void writeDefaultOutputFile(const OptionsIO::Library library) {
  writeDefaultOutputFile(Options::root(), library);
}

void writeDefaultOutputFile(Options& options, const OptionsIO::Library library) {
  bout::experimental::addBuildFlagsToOptions(options);
  bout::globals::mesh->outputVars(options);
  auto mode = options["append"]
                      .doc("Add output data to existing (dump) files?")
                      .withDefault(false)
                  ? bout::OptionsIO::FileMode::append
                  : bout::OptionsIO::FileMode::replace;
  auto io = OptionsIOFactory(getOutputFilename(options, library), mode, library);
  io->write(options);
}

} // namespace bout
