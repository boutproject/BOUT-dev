#include "bout/build_config.hxx"

#if BOUT_HAS_ADIOS

#include "bout/adios_object.hxx"

#include <exception>
#include <iostream>
#include <unordered_map>

namespace bout {

static ADIOSPtr adios = nullptr;
ADIOSStruct adios_restart;
ADIOSStruct adios_dump;
//ADIOS2Param* adios2params;

void ADIOSInit(MPI_Comm comm) { adios = std::make_shared<adios2::ADIOS>(comm); }

void ADIOSInit(const std::string configFile, MPI_Comm comm) {
  adios = std::make_shared<adios2::ADIOS>(configFile, comm);
}

void ADIOSFinalize() {
  if (adios == nullptr) {
    throw std::runtime_error(
        "ADIOS needs to be initialized first before calling ADIOSFinalize()");
  }
  if (adios_dump.adiosStep > 0 && adios_dump.engine) {
    adios_dump.engine.Close();
  }
  adios.reset();
}

ADIOSPtr GetADIOSPtr() {
  if (adios == nullptr) {
    throw std::runtime_error(
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

ADIOSStruct& ADIOSGetStruct(const std::string& fname) {
  if (fname.find(".restart") != std::string::npos)
    return adios_restart;
  //if (fname.find(".dmp") != std::string::npos)
  //  return adios_dump;
  return adios_dump;
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
      throw std::invalid_argument("ADIOSSetParameters(): wrong format for IO parameter "
                                  + parameter + ", format must be key" + delimKeyValue
                                  + "value for each entry");
    }

    std::string key = parameter.substr(0, position);
    lf_Trim(key);
    std::string value = parameter.substr(position + 1);
    lf_Trim(value);
    if (value.length() == 0) {
      throw std::invalid_argument("ADIOS2SetParameters: empty value in IO parameter "
                                  + parameter + ", format must be key" + delimKeyValue
                                  + "value");
    }
    io.SetParameter(key, value);
  }
}

} // namespace bout
#endif //BOUT_HAS_ADIOS
