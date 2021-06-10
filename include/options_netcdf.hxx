
#pragma once

#ifndef __OPTIONS_NETCDF_H__
#define __OPTIONS_NETCDF_H__

#include "bout/build_config.hxx"

#if !BOUT_HAS_NETCDF || BOUT_HAS_LEGACY_NETCDF

#include <string>

#include "boutexception.hxx"
#include "options.hxx"

namespace bout {

class OptionsNetCDF {
public:
  enum class FileMode {
                       replace, ///< Overwrite file when writing
                       append   ///< Append to file when writing
  };
  
  OptionsNetCDF(const std::string &filename, FileMode mode = FileMode::replace) {}
  OptionsNetCDF(const OptionsNetCDF&) = default;
  OptionsNetCDF(OptionsNetCDF&&) = default;
  OptionsNetCDF& operator=(const OptionsNetCDF&) = default;
  OptionsNetCDF& operator=(OptionsNetCDF&&) = default;

  /// Read options from file
  Options read() {
    throw BoutException("OptionsNetCDF not available\n");
  }

  /// Write options to file
  void write(const Options &options) {
    throw BoutException("OptionsNetCDF not available\n");
  }
};

}

#else

#include <string>

#include "options.hxx"

namespace bout {

class OptionsNetCDF {
public:
  enum class FileMode {
                       replace, ///< Overwrite file when writing
                       append   ///< Append to file when writing
  };

  OptionsNetCDF() {}
  explicit OptionsNetCDF(std::string filename, FileMode mode = FileMode::replace)
      : filename(std::move(filename)), file_mode(mode) {}
  OptionsNetCDF(const OptionsNetCDF&) = default;
  OptionsNetCDF(OptionsNetCDF&&) = default;
  OptionsNetCDF& operator=(const OptionsNetCDF&) = default;
  OptionsNetCDF& operator=(OptionsNetCDF&&) = default;

  /// Read options from file
  Options read();

  /// Write options to file
  void write(const Options& options) { write(options, "t"); }
  void write(const Options& options, const std::string& time_dim);

  /// Check that all variables with the same time dimension have the
  /// same size in that dimension. Throws BoutException if there are
  /// any differences, otherwise is silent
  void verifyTimesteps() const;
private:
  std::string filename;
  FileMode file_mode{FileMode::replace};
};

} // namespace bout

#endif

namespace bout {
/// Name of the directory for restart files
std::string getRestartDirectoryName(Options& options);
/// Name of the restart file on this rank
std::string getRestartFilename(Options& options);
/// Name of the restart file on \p rank
std::string getRestartFilename(Options& options, int rank);
/// Name of the main output file on this rank
std::string getOutputFilename(Options& options);
/// Name of the main output file on \p rank
std::string getOutputFilename(Options& options, int rank);
/// Write `Options::root()` to the main output file, overwriting any
/// existing files
void writeDefaultOutputFile();
/// Write \p options to the main output file, overwriting any existing
/// files
void writeDefaultOutputFile(Options& options);
} // namespace bout

#endif //  __OPTIONS_NETCDF_H__
