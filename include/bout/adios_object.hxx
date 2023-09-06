#ifndef ADIOS_OBJECT_HXX
#define ADIOS_OBJECT_HXX

#include "bout/build_config.hxx"

#if BOUT_HAS_ADIOS

#include <adios2.h>
#include <memory>
#include <mpi.h>

namespace bout {

void ADIOSInit(MPI_Comm comm);
void ADIOSInit(const std::string configFile, MPI_Comm comm);
void ADIOSFinalize();

using ADIOSPtr = std::shared_ptr<adios2::ADIOS>;
using EnginePtr = std::shared_ptr<adios2::Engine>;
using IOPtr = std::shared_ptr<adios2::IO>;

ADIOSPtr GetADIOSPtr();
IOPtr GetIOPtr(const std::string IOName);

class ADIOSStream {
public:
  adios2::IO io;
  adios2::Engine engine;
  adios2::Variable<double> vTime;
  adios2::Variable<int> vStep;
  int adiosStep = 0;
};

extern ADIOSStream adios_restart;
extern ADIOSStream adios_dump;
//extern ADIOS2Param *adios2params;

/** return one of the extern variable based on the target file name */
ADIOSStream& ADIOSGetStream(const std::string& fname);

/** Set user parameters for an IO group */
void ADIOSSetParameters(const std::string& input, const char delimKeyValue,
                        const char delimItem, adios2::IO& io);

} // namespace bout

#endif //BOUT_HAS_ADIOS
#endif //ADIOS_OBJECT_HXX
