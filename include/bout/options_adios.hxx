
#pragma once

#ifndef __OPTIONS_ADIOS_H__
#define __OPTIONS_ADIOS_H__

#include "bout/build_config.hxx"
#include "bout/options.hxx"
#include "bout/options_io.hxx"

#if !BOUT_HAS_ADIOS

#include <string>

#include "bout/boutexception.hxx"

namespace bout {

class OptionsADIOS : public OptionsIO {
public:
  OptionsADIOS() {}
  explicit OptionsADIOS(
      const std::string& filename [[maybe_unused]],
      [[maybe_unused]] bout::OptionsIO::FileMode mode = bout::OptionsIO::FileMode::replace,
      [[maybe_unused]] bool singleWriteFile = false) {}
  OptionsADIOS(const OptionsADIOS&) = delete;
  OptionsADIOS(OptionsADIOS&&) noexcept = default;
  ~OptionsADIOS() = default;

  OptionsADIOS& operator=(const OptionsADIOS&) = delete;
  OptionsADIOS& operator=(OptionsADIOS&&) noexcept = default;

  /// Read options from file
  Options read() override { throw BoutException("OptionsADIOS not available\n"); }

  /// Write options to file
  void write([[maybe_unused]] const Options& options, [[maybe_unused]] const std::string& time_dim) override {
    throw BoutException("OptionsADIOS not available\n");
  }

  void verifyTimesteps() const override { throw BoutException("OptionsADIOS not available\n"); }
};

} // namespace bout

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
  OptionsADIOS() {}
  OptionsADIOS(std::string filename,
               bout::OptionsIO::FileMode mode = bout::OptionsIO::FileMode::replace,
               bool singleWriteFile = false) : OptionsIO(filename, mode, singleWriteFile) {}
  OptionsADIOS(const OptionsADIOS&) = delete;
  OptionsADIOS(OptionsADIOS&&) noexcept = default;
  ~OptionsADIOS() = default;

  OptionsADIOS& operator=(const OptionsADIOS&) = delete;
  OptionsADIOS& operator=(OptionsADIOS&&) noexcept = default;

  /// Read options from file
  Options read() override;

  /// Write options to file
  void write(const Options& options) override { write(options, "t"); }
  void write(const Options& options, const std::string& time_dim) override;

  /// Check that all variables with the same time dimension have the
  /// same size in that dimension. Throws BoutException if there are
  /// any differences, otherwise is silent
  void verifyTimesteps() const override;

private:
};

} // namespace bout

#endif // BOUT_HAS_ADIOS

#endif //  __OPTIONS_ADIOS_H__
