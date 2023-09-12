/* Parent class for IO/ADIOS option classes */

#pragma once

#ifndef __OPTIONS_IO_H__
#define __OPTIONS_IO_H__

#include "bout/build_config.hxx"

#include <memory>
#include <string>

#include "bout/options.hxx"

namespace bout {

class OptionsIO {
public:
  enum class FileMode {
    replace, ///< Overwrite file when writing
    append   ///< Append to file when writing
  };

  enum class Library { ADIOS, NetCDF, Invalid };

  static const Library defaultIOLibrary =
#if BOUT_HAS_ADIOS
      Library::ADIOS;
#elif BOUT_HAS_NETCDF
      Library::NetCDF;
#else
      Library::Invalid;
#endif

  OptionsIO();
  OptionsIO(std::string filename, FileMode mode = FileMode::replace);
  ~OptionsIO();
  OptionsIO(const OptionsIO&) = delete;
  OptionsIO(OptionsIO&&) noexcept;
  OptionsIO& operator=(const OptionsIO&) = delete;
  OptionsIO& operator=(OptionsIO&&) noexcept;

  /// Read options from file
  virtual Options read() = 0;

  /// Write options to file
  void write(const Options& options) { write(options, "t"); }
  virtual void write(const Options& options, const std::string& time_dim) = 0;

  /// Check that all variables with the same time dimension have the
  /// same size in that dimension. Throws BoutException if there are
  /// any differences, otherwise is silent
  virtual void verifyTimesteps() const = 0;

protected:
  /// Name of the file on disk
  std::string filename;
  /// How to open the file for writing
  FileMode file_mode{FileMode::replace};
  Library library = Library::Invalid;
};

std::shared_ptr<OptionsIO>
OptionsIOFactory(std::string filename,
                 OptionsIO::FileMode mode = OptionsIO::FileMode::replace,
                 const OptionsIO::Library library = OptionsIO::defaultIOLibrary);

OptionsIO::Library getIOLibrary(Options& options);

/// Name of the directory for restart files
std::string getRestartDirectoryName(Options& options);
/// Name of the restart file on this rank
std::string getRestartFilename(Options& options, const OptionsIO::Library library);
/// Name of the restart file on \p rank
std::string getRestartFilename(Options& options, int rank,
                               const OptionsIO::Library library);
/// Name of the main output file on this rank
std::string getOutputFilename(Options& options, const OptionsIO::Library library);
/// Name of the main output file on \p rank
std::string getOutputFilename(Options& options, int rank,
                              const OptionsIO::Library library);
/// Write `Options::root()` to the main output file, overwriting any
/// existing files
void writeDefaultOutputFile(
    const OptionsIO::Library library = OptionsIO::defaultIOLibrary);
/// Write \p options to the main output file, overwriting any existing
/// files
void writeDefaultOutputFile(
    Options& options, const OptionsIO::Library library = OptionsIO::defaultIOLibrary);

} // namespace bout

#endif //  __OPTIONS_IO_H__
