
#pragma once

#ifndef OPTIONS_ADIOS_H
#define OPTIONS_ADIOS_H

#include "bout/build_defines.hxx"
#include "bout/options.hxx"
#include "bout/options_io.hxx"

#if !BOUT_HAS_ADIOS2

namespace {
bout::RegisterUnavailableOptionsIO
    registerunavailableoptionsadios("adios", "BOUT++ was not configured with ADIOS2");
}

#else

#include <memory>
#include <string>

namespace bout {

/// Forward declare ADIOS file type so we don't need to depend
/// directly on ADIOS
struct ADIOSStream;

class OptionsADIOS : public OptionsIO {
public:
  // Constructors need to be defined in implementation due to forward
  // declaration of ADIOSStream
  OptionsADIOS() = delete;

  /// Create an OptionsADIOS
  ///
  /// Options:
  ///  - "file"   The name of the file
  ///    If not set then "path" and "prefix" must be set,
  ///    and file is set to {path}/{prefix}.bp
  ///  - "append"
  ///  - "singleWriteFile"
  OptionsADIOS(Options& options);

  OptionsADIOS(const OptionsADIOS&) = delete;
  OptionsADIOS(OptionsADIOS&&) noexcept = default;
  ~OptionsADIOS() = default;

  OptionsADIOS& operator=(const OptionsADIOS&) = delete;
  OptionsADIOS& operator=(OptionsADIOS&&) noexcept = default;

  /// Read options from file
  Options read() override;

  /// Write options to file
  void write(const Options& options, const std::string& time_dim) override;

  /// Check that all variables with the same time dimension have the
  /// same size in that dimension. Throws BoutException if there are
  /// any differences, otherwise is silent
  void verifyTimesteps() const override;

private:
  enum class FileMode {
    replace, ///< Overwrite file when writing
    append   ///< Append to file when writing
  };

  /// Name of the file on disk
  std::string filename;
  /// How to open the file for writing
  FileMode file_mode{FileMode::replace};
  bool singleWriteFile = false;
};

namespace {
RegisterOptionsIO<OptionsADIOS> registeroptionsadios("adios");
}

} // namespace bout

#endif // BOUT_HAS_ADIOS2
#endif // OPTIONS_ADIOS_H
