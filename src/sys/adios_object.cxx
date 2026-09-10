#include "bout/build_defines.hxx"

#if BOUT_HAS_ADIOS2

#include "bout/adios_object.hxx"
#include "bout/boutexception.hxx"

#include <adios2.h>

#include <memory>
#include <string>
#include <unordered_map>

namespace bout {

static ADIOSPtr adios = nullptr;
static std::unordered_map<std::string, std::unique_ptr<ADIOSStream>> adiosStreams;

namespace {
std::string GetIOName(const std::string& fname, adios2::Mode mode) {
  switch (mode) {
  case adios2::Mode::Read:
    return "read_" + fname;
  case adios2::Mode::ReadRandomAccess:
    return "read_random_access_" + fname;
  case adios2::Mode::Append:
    return "append_" + fname;
  case adios2::Mode::Write:
    return "write_" + fname;
  default:
    return fname;
  }
}

adios2::IO GetIO(const std::string& fname, adios2::Mode mode,
                 const std::string& engineType) {
  auto io = GetADIOSPtr()->DeclareIO(GetIOName(fname, mode));
  io.SetEngine(engineType);
  return io;
}
} // namespace

void ADIOSInit(MPI_Comm comm) { adios = std::make_shared<adios2::ADIOS>(comm); }

void ADIOSInit(const std::string configFile, MPI_Comm comm) {
  adios = std::make_shared<adios2::ADIOS>(configFile, comm);
}

void ADIOSFinalize() {
  if (adios == nullptr) {
    throw BoutException(
        "ADIOS needs to be initialized first before calling ADIOSFinalize()");
  }
  adiosStreams.clear();
  adios.reset();
}

ADIOSPtr GetADIOSPtr() {
  if (adios == nullptr) {
    throw BoutException(
        "ADIOS needs to be initialized first before calling GetADIOSPtr()");
  }
  return adios;
}

ADIOSStream::~ADIOSStream() { finish(); }

ADIOSStream::ADIOSStream(const std::string& fname, adios2::Mode mode,
                         const std::string& engineType)
    : io(GetIO(fname, mode, engineType)), fname(fname), io_name(GetIOName(fname, mode)),
      file_mode(mode) {}

ADIOSStream& ADIOSStream::ADIOSGetStream(const std::string& fname, adios2::Mode mode,
                                         const std::string& engineType) {
  const auto key = GetIOName(fname, mode);
  auto it = adiosStreams.find(key);
  if (it == adiosStreams.end() || it->second->isRetired()) {
    it = adiosStreams
             .insert_or_assign(key, std::unique_ptr<ADIOSStream>(
                                        new ADIOSStream(fname, mode, engineType)))
             .first;
  }
  return *it->second;
}

} // namespace bout
#endif //BOUT_HAS_ADIOS2
