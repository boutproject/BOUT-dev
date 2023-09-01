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

struct ADIOSStruct {
  adios2::IO io;
  adios2::Engine engine;
  adios2::Variable<int> vCellDim, vPhysDim, vNTotalElem, vNTotalNode;
  adios2::Variable<double> vPhysTime;
  adios2::Variable<int> vStep;
  adios2::Variable<int64_t> vElementConnectivity, vElementRange;
  adios2::Variable<double> vXco, vYco, vZco;
  std::vector<adios2::Variable<double>> vFlowVars;
  std::vector<adios2::Variable<int64_t>> vBoundariesConnectivity;
  std::vector<adios2::Variable<int64_t>> vBoundariesRange;
  int adiosStep = 0;
};

extern ADIOSStruct adios_restart;
extern ADIOSStruct adios_dump;
//extern ADIOS2Param *adios2params;

/** return one of the extern variable based on the target file name */
ADIOSStruct& ADIOSGetStruct(const std::string& fname);

/** Set user parameters for an IO group */
void ADIOSSetParameters(const std::string& input, const char delimKeyValue,
                        const char delimItem, adios2::IO& io);

} // namespace bout

#endif //BOUT_HAS_ADIOS
#endif //ADIOS_OBJECT_HXX
