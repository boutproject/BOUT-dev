
#pragma once

#ifndef __OPTIONS_NETCDF_H__
#define __OPTIONS_NETCDF_H__

#include "bout/build_config.hxx"

#if !BOUT_HAS_NETCDF || BOUT_HAS_LEGACY_NETCDF

#include <string>

#include "boutexception.hxx"
#include "options.hxx"

namespace bout {
namespace experimental {

class OptionsNetCDF {
public:
  enum class FileMode {
                       replace, ///< Overwrite file when writing
                       append   ///< Append to file when writing
  };
  
  OptionsNetCDF(const std::string &filename, FileMode mode = FileMode::replace) {}

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
}

#else

#include <string>

#include "options.hxx"

namespace bout {
namespace experimental {

class OptionsNetCDF {
public:
  enum class FileMode {
                       replace, ///< Overwrite file when writing
                       append   ///< Append to file when writing
  };

  OptionsNetCDF() {}
  explicit OptionsNetCDF(std::string filename, FileMode mode = FileMode::replace)
      : filename(std::move(filename)), file_mode(mode) {}

  /// Read options from file
  Options read();

  /// Write options to file
  void write(const Options &options);

  /// Check that all variables with the same time dimension have the
  /// same size in that dimension. Throws BoutException if there are
  /// any differences, otherwise is silent
  void verifyTimesteps() const;
private:
  std::string filename;
  FileMode file_mode;
};

}
}

#endif

#endif //  __OPTIONS_NETCDF_H__
