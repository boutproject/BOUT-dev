#include "bout/build_defines.hxx"

#if BOUT_HAS_ADIOS2

#include "bout/adios_object.hxx"
#include "bout/boutexception.hxx"

#include <exception>
#include <iostream>
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
  if (engine) {
    if (isInStep) {
      engine.EndStep();
      isInStep = false;
    }
    engine.Close();
  }
}

ADIOSStream& ADIOSStream::ADIOSGetStream(const std::string& fname) {
  auto it = adiosStreams.find(fname);
  if (it == adiosStreams.end()) {
    it = adiosStreams.emplace(fname, ADIOSStream(fname)).first;
  }
  return it->second;
}

void ADIOSSetParameters(const std::string& input, const char delimKeyValue,
                        const char delimItem, adios2::IO& io) {
  auto lf_Trim = [](std::string& input) {
    input.erase(0, input.find_first_not_of(" \n\r\t")); // prefixing spaces
    input.erase(input.find_last_not_of(" \n\r\t") + 1); // suffixing spaces
  };

  std::istringstream inputSS(input);
  std::string parameter;
  while (std::getline(inputSS, parameter, delimItem)) {
    const size_t position = parameter.find(delimKeyValue);
    if (position == parameter.npos) {
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
