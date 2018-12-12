
#pragma once

#ifndef __OPTIONS_NETCDF_H__
#define __OPTIONS_NETCDF_H__

#ifndef NCDF4

#include <string>

#include "boutexception.hxx"
#include "options.hxx"

class OptionsNetCDF {
public:
  OptionsNetCDF(const std::string &filename) {}

  /// Read options from file
  Options read() {
    throw BoutException("OptionNetCDF not available\n");
  }

  /// Write options to file
  void write(const Options &options) {
    throw BoutException("OptionNetCDF not available\n");
  }
};

#else

#include <string>

#include "options.hxx"

class OptionsNetCDF {
public:
  OptionsNetCDF(const std::string &filename) : filename(filename) {}

  /// Read options from file
  Options read();

  /// Write options to file
  void write(const Options &options);

private:
  std::string filename;
};

#endif

#endif //  __OPTIONS_NETCDF_H__
