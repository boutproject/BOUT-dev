#include "bout/build_defines.hxx"

#if BOUT_HAS_ADIOS2

#include "bout/adios_object.hxx"
#include "bout/boutexception.hxx"

#include <adios2.h>

#include <unordered_map>

namespace bout {

static ADIOSPtr adios = nullptr;
static std::unordered_map<std::string, ADIOSStream> adiosStreams;

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

IOPtr GetIOPtr(const std::string IOName) {
  auto adios = GetADIOSPtr();
  IOPtr io = nullptr;
  try {
    io = std::make_shared<adios2::IO>(adios->AtIO(IOName));
  } catch (std::invalid_argument& e) {
  }
  return io;
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

ADIOSStream& ADIOSStream::ADIOSGetStream(const std::string& fname, adios2::Mode mode) {
  auto it = adiosStreams.find(fname);
  if (it == adiosStreams.end()) {
    it = adiosStreams.emplace(fname, ADIOSStream(fname, mode)).first;
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
