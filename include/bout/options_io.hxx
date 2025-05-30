/// Parent class for IO to binary files and streams
///
///
/// Usage:
///
/// 1. Dump files, containing time history:
///
///       auto dump = OptionsIOFactory::getInstance().createOutput();
///       dump->write(data);
///
///    where data is an Options tree. By default dump files are configured
///    with the root `output` section, or an Option tree can be passed to
///    `createOutput`.
///
/// 2. Restart files:
///
///       auto restart = OptionsIOFactory::getInstance().createOutput();
///       restart->write(data);
///
///    where data is an Options tree. By default restart files are configured
///    with the root `restart_files` section, or an Option tree can be passed to
///    `createRestart`.
///
/// 3. Ad-hoc single files
///   Note: The caller should consider how multiple processors interact with the file.
///
///       auto file = OptionsIOFactory::getInstance().createFile("some_file.nc");
///   or
///       auto file = OptionsIO::create("some_file.nc");
///
///

#pragma once

#ifndef OPTIONS_IO_H
#define OPTIONS_IO_H

#include "bout/build_defines.hxx"
#include "bout/generic_factory.hxx"
#include "bout/options.hxx"

#include <memory>
#include <string>

namespace bout {

class OptionsIO {
public:
  /// No default constructor, as settings are required
  OptionsIO() = delete;

  /// Constructor specifies the kind of file, and options to control
  /// the name of file, mode of operation etc.
  OptionsIO(Options&) {}

  virtual ~OptionsIO() = default;

  OptionsIO(const OptionsIO&) = delete;
  OptionsIO(OptionsIO&&) noexcept = default;
  OptionsIO& operator=(const OptionsIO&) = delete;
  OptionsIO& operator=(OptionsIO&&) noexcept = default;

  /// Read options from file
  virtual Options read(bool lazy = true) = 0;

  /// Write options to file
  void write(const Options& options) { write(options, "t"); }
  virtual void write(const Options& options, const std::string& time_dim) = 0;

  /// NetCDF: Check that all variables with the same time dimension have the
  /// same size in that dimension. Throws BoutException if there are
  /// any differences, otherwise is silent.
  /// ADIOS: Indicate completion of an output step.
  virtual void verifyTimesteps() const = 0;

  /// Create an OptionsIO for I/O to the given file.
  /// This uses the default file type and default options.
  static std::unique_ptr<OptionsIO> create(const std::string& file);

  /// Create an OptionsIO for I/O to the given file.
  /// The file will be configured using the given `config` options:
  ///  - "type" : string   The file type e.g. "netcdf" or "adios"
  ///  - "file" : string   Name of the file
  ///  - "append" : bool   Append to existing data (Default is false)
  static std::unique_ptr<OptionsIO> create(Options& config);

  /// Create an OptionsIO for I/O to the given file.
  /// The file will be configured using the given `config` options:
  ///  - "type" : string   The file type e.g. "netcdf" or "adios"
  ///  - "file" : string   Name of the file
  ///  - "append" : bool   Append to existing data (Default is false)
  ///
  /// Example:
  ///
  ///     auto file = OptionsIO::create({
  ///         {"file", "some_file.nc"},
  ///         {"type", "netcdf"},
  ///         {"append", false}
  ///     });
  static std::unique_ptr<OptionsIO> create(Options::InitializerList config_list) {
    Options config(config_list); // Construct an Options to pass by reference
    return create(config);
  }
};

class OptionsIOFactory : public Factory<OptionsIO, OptionsIOFactory, Options&> {
public:
  static constexpr auto type_name = "OptionsIO";
  static constexpr auto section_name = "io";
  static constexpr auto option_name = "type";
  static constexpr auto default_type =
#if BOUT_HAS_NETCDF
      "netcdf";
#elif BOUT_HAS_ADIOS2
      "adios";
#else
      "invalid";
#endif

  /// Create a restart file, configured with options (if given),
  /// or root "restart_files" section.
  ///
  /// Options:
  ///  - "type"    The type of file e.g "netcdf" or "adios"
  ///  - "file"    Name of the file. Default is <path>/<prefix>.[type-dependent]
  ///  - "path"    Path to restart files. Default is root "datadir" option,
  ///              that defaults to "data"
  ///  - "prefix"  Default is "BOUT.restart"
  ReturnType createRestart(Options* optionsptr = nullptr) const;

  /// Create an output file for writing time history.
  /// Configure with options (if given), or root "output" section.
  ///
  /// Options:
  ///  - "type"    The type of file e.g "netcdf" or "adios"
  ///  - "file"    Name of the file. Default is <path>/<prefix>.[type]
  ///  - "path"    Path to output files. Default is root "datadir" option,
  ///              that defaults to "data"
  ///  - "prefix"  Default is "BOUT.dmp"
  ///  - "append"  Append to existing file? Default is root "append" option,
  ///              that defaults to false.
  ReturnType createOutput(Options* optionsptr = nullptr) const;

  /// Create a single file (e.g. mesh file) of the default type
  ReturnType createFile(const std::string& file) const;
};

/// Simpler name for Factory registration helper class
///
/// Usage:
///
///     #include <bout/options_io.hxx>
///     namespace {
///     RegisterOptionsIO<MyOptionsIO> registeroptionsiomine("myoptionsio");
///     }
template <typename DerivedType>
using RegisterOptionsIO = OptionsIOFactory::RegisterInFactory<DerivedType>;

/// Simpler name for indicating that an OptionsIO implementation
/// is unavailable.
///
/// Usage:
///
///     namespace {
///     RegisterUnavailableOptionsIO
///         unavailablemyoptionsio("myoptiosio", "BOUT++ was not configured with MyOptionsIO");
///     }
using RegisterUnavailableOptionsIO = OptionsIOFactory::RegisterUnavailableInFactory;

/// Convenient wrapper function around OptionsIOFactory::createOutput
/// Opens a dump file configured with the `output` root section,
/// and writes the given `data` to the file.
void writeDefaultOutputFile(Options& data);

} // namespace bout

#endif //  OPTIONS_IO_H
