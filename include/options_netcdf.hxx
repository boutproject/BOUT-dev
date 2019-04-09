
#pragma once

#ifndef __OPTIONS_NETCDF_H__
#define __OPTIONS_NETCDF_H__

#ifndef NCDF4

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
  
  OptionsNetCDF(const std::string &filename, FileMode mode = FileMode::replace) : filename(filename), file_mode(mode) {}

  /// Read options from file
  Options read();

  /// Write options to file
  void write(const Options &options);

private:
  std::string filename;
  FileMode file_mode;

  /// NetCDF doesn't seem to keep track of the current index
  /// for each variable. This map keeps track of the current
  /// index being written for each time dimension
  std::map<int, size_t> time_index;
};

}
}

#endif

#endif //  __OPTIONS_NETCDF_H__
