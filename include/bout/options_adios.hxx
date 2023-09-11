
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
  explicit OptionsADIOS(const std::string& filename, bout::OptionsIO::FileMode mode = bout::OptionsIO::FileMode::replace) {}
  OptionsADIOS(const OptionsADIOS&) = default;
  OptionsADIOS(OptionsADIOS&&) = default;
  OptionsADIOS& operator=(const OptionsADIOS&) = default;
  OptionsADIOS& operator=(OptionsADIOS&&) = default;

  /// Read options from file
  Options read() { throw BoutException("OptionsADIOS not available\n"); }

  /// Write options to file
  void write(const Options& options, const std::string& time_dim) {
    throw BoutException("OptionsADIOS not available\n");
  }

  void verifyTimesteps() const{
    throw BoutException("OptionsADIOS not available\n");
  }


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
  OptionsADIOS();
  OptionsADIOS(std::string filename, bout::OptionsIO::FileMode mode =
                                                  bout::OptionsIO::FileMode::replace);
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
  /// Pointer to ADIOS stream so we don't introduce direct dependence
  std::unique_ptr<ADIOSStream> stream;
};

} // namespace bout

#endif // BOUT_HAS_ADIOS

#endif //  __OPTIONS_ADIOS_H__
