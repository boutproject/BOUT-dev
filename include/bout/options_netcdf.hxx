
#pragma once

#ifndef __OPTIONS_NETCDF_H__
#define __OPTIONS_NETCDF_H__

#include "bout/build_config.hxx"

#include "bout/options.hxx"
#include "bout/options_io.hxx"

#if !BOUT_HAS_NETCDF || BOUT_HAS_LEGACY_NETCDF

#include <string>

#include "bout/boutexception.hxx"
#include "bout/options.hxx"

namespace bout {

class OptionsNetCDF : public OptionsIO {
public:
  OptionsNetCDF(const std::string& filename,
                bout::OptionsIO::FileMode mode = bout::OptionsIO::FileMode::replace,
                bool singleWriteFile = false) {}
  OptionsNetCDF(const OptionsNetCDF&) = default;
  OptionsNetCDF(OptionsNetCDF&&) = default;
  OptionsNetCDF& operator=(const OptionsNetCDF&) = default;
  OptionsNetCDF& operator=(OptionsNetCDF&&) = default;

  /// Read options from file
  Options read() { throw BoutException("OptionsNetCDF not available\n"); }

  /// Write options to file
  void write(const Options& options, const std::string& time_dim) {
    throw BoutException("OptionsNetCDF not available\n");
  }

  void verifyTimesteps() const { throw BoutException("OptionsADIOS not available\n"); }
};

} // namespace bout

#else

#include <memory>
#include <string>

/// Forward declare netCDF file type so we don't need to depend
/// directly on netCDF
namespace netCDF {
class NcFile;
}

namespace bout {

class OptionsNetCDF : public OptionsIO {
public:
  // Constructors need to be defined in implementation due to forward
  // declaration of NcFile
  OptionsNetCDF();
  OptionsNetCDF(std::string filename, FileMode mode = FileMode::replace,
                bool singleWriteFile = false);
  ~OptionsNetCDF();
  OptionsNetCDF(const OptionsNetCDF&) = delete;
  OptionsNetCDF(OptionsNetCDF&&) noexcept;
  OptionsNetCDF& operator=(const OptionsNetCDF&) = delete;
  OptionsNetCDF& operator=(OptionsNetCDF&&) noexcept;

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
  /// Pointer to netCDF file so we don't introduce direct dependence
  std::unique_ptr<netCDF::NcFile> data_file;
};

} // namespace bout

#endif

#endif //  __OPTIONS_NETCDF_H__
