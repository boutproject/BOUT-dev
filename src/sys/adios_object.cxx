#include "bout/build_defines.hxx"

#if BOUT_HAS_ADIOS2

#include "bout/adios_object.hxx"
#include "bout/boutexception.hxx"

#include <adios2.h>

#include <string>
#include <unordered_map>

namespace bout {

static ADIOSPtr adios = nullptr;
static std::unordered_map<std::string, ADIOSStream> adiosStreams;

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

ADIOSStream::~ADIOSStream() {
  if (engine_) {
    if (isInStep) {
      engine_.EndStep();
      isInStep = false;
    }
    engine_.Close();
  }
}

ADIOSStream::ADIOSStream(const std::string& fname, adios2::Mode mode,
                         const std::string& engineType)
    : io(GetIO(fname, mode, engineType)), fname(fname), file_mode(mode) {}

ADIOSStream& ADIOSStream::ADIOSGetStream(const std::string& fname, adios2::Mode mode,
                                         const std::string& engineType) {
  const auto key = GetIOName(fname, mode);
  auto it = adiosStreams.find(key);
  if (it == adiosStreams.end()) {
    it = adiosStreams.emplace(key, ADIOSStream(fname, mode, engineType)).first;
  }
  return it->second;
}

void ADIOSSetParameters(const std::string& input, char delimKeyValue, char delimItem,
                        adios2::IO& io) {
  auto lf_Trim = [](std::string& input) {
    input.erase(0, input.find_first_not_of(" \n\r\t")); // prefixing spaces
    input.erase(input.find_last_not_of(" \n\r\t") + 1); // suffixing spaces
  };

  std::istringstream inputSS(input);
  std::string parameter;
  while (std::getline(inputSS, parameter, delimItem)) {
    const size_t position = parameter.find(delimKeyValue);
    if (position == std::string::npos) {
      throw BoutException("ADIOSSetParameters(): wrong format for IO parameter "
                          + parameter + ", format must be key" + delimKeyValue
                          + "value for each entry");
    }

    std::string key = parameter.substr(0, position);
    lf_Trim(key);
    std::string value = parameter.substr(position + 1);
    lf_Trim(value);
    if (value.length() == 0) {
      throw BoutException("ADIOS2SetParameters: empty value in IO parameter " + parameter
                          + ", format must be key" + delimKeyValue + "value");
    }
    io.SetParameter(key, value);
  }
}

} // namespace bout
#endif //BOUT_HAS_ADIOS2
