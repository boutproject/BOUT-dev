
#pragma once

#ifndef OPTIONS_NETCDF_H
#define OPTIONS_NETCDF_H

#include "bout/build_defines.hxx"

#include "bout/options.hxx"
#include "bout/options_io.hxx"

#if !BOUT_HAS_NETCDF || BOUT_HAS_LEGACY_NETCDF

namespace {
RegisterUnavailableOptionsIO
    registerunavailableoptionsnetcdf("netcdf", "BOUT++ was not configured with NetCDF");
}

#else

#include <memory>
#include <netcdf>
#include <string>

namespace bout {

class OptionsNetCDF : public OptionsIO {
public:
  // Constructors need to be defined in implementation due to forward
  // declaration of NcFile
  OptionsNetCDF() = delete;

  /// Create an OptionsNetCDF
  ///
  /// Options:
  ///  - "file"   The name of the file
  ///    If not set then "path" and "prefix" options must be set,
  ///    and file is set to {path}/{prefix}.{rank}.nc
  ///  - "append"  File mode, default is false
  OptionsNetCDF(Options& options);

  ~OptionsNetCDF() {}

  OptionsNetCDF(const OptionsNetCDF&) = delete;
  OptionsNetCDF(OptionsNetCDF&&) noexcept = default;
  OptionsNetCDF& operator=(const OptionsNetCDF&) = delete;
  OptionsNetCDF& operator=(OptionsNetCDF&&) noexcept = default;

  /// Read options from file
  Options read(bool lazy = true) override;

  /// Write options to file
  void write(const Options& options) { write(options, "t"); }
  void write(const Options& options, const std::string& time_dim);

  /// Check that all variables with the same time dimension have the
  /// same size in that dimension. Throws BoutException if there are
  /// any differences, otherwise is silent
  void verifyTimesteps() const;

private:
  enum class FileMode {
    replace, ///< Overwrite file when writing
    append   ///< Append to file when writing
  };

  /// Pointer to netCDF file so we don't introduce direct dependence
  std::unique_ptr<netCDF::NcFile> data_file = nullptr;

  /// Name of the file on disk
  std::string filename;
  /// How to open the file for writing
  FileMode file_mode{FileMode::replace};
};

namespace {
RegisterOptionsIO<OptionsNetCDF> registeroptionsnetcdf("netcdf");
}

} // namespace bout

#endif

#endif //  OPTIONS_NETCDF_H
