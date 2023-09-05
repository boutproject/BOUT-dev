
#pragma once

#ifndef __OPTIONS_ADIOS_H__
#define __OPTIONS_ADIOS_H__

#include "bout/build_config.hxx"

#if !BOUT_HAS_ADIOS

#include <string>

#include "bout/boutexception.hxx"
#include "bout/options.hxx"

namespace bout {

class OptionsADIOS {
public:
  enum class FileMode {
    replace, ///< Overwrite file when writing
    append   ///< Append to file when writing
  };

  OptionsADIOS(const std::string& filename, FileMode mode = FileMode::replace) {}
  OptionsADIOS(const OptionsADIOS&) = default;
  OptionsADIOS(OptionsADIOS&&) = default;
  OptionsADIOS& operator=(const OptionsADIOS&) = default;
  OptionsADIOS& operator=(OptionsADIOS&&) = default;

  /// Read options from file
  Options read() { throw BoutException("OptionsADIOS not available\n"); }

  /// Write options to file
  void write(const Options& options) {
    throw BoutException("OptionsADIOS not available\n");
  }
};

} // namespace bout

#else

#include <memory>
#include <string>

#include "bout/adios_object.hxx"
#include "bout/options.hxx"

namespace bout {

/// Forward declare ADIOS file type so we don't need to depend
/// directly on ADIOS
struct ADIOSStream;

class OptionsADIOS {
public:
  enum class FileMode {
    replace, ///< Overwrite file when writing
    append   ///< Append to file when writing
  };

  // Constructors need to be defined in implementation due to forward
  // declaration of NcFile
  OptionsADIOS();
  explicit OptionsADIOS(std::string filename, FileMode mode = FileMode::replace);
  ~OptionsADIOS();
  OptionsADIOS(const OptionsADIOS&) = delete;
  OptionsADIOS(OptionsADIOS&&) noexcept;
  OptionsADIOS& operator=(const OptionsADIOS&) = delete;
  OptionsADIOS& operator=(OptionsADIOS&&) noexcept;

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
  /// Name of the file on disk
  std::string filename;
  /// How to open the file for writing
  FileMode file_mode{FileMode::replace};
  /// Pointer to ADIOS file so we don't introduce direct dependence
  std::unique_ptr<ADIOSStream> data_file;
};

} // namespace bout

#endif // BOUT_HAS_ADIOS

#endif //  __OPTIONS_ADIOS_H__
