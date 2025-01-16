/*!************************************************************************
 * Provides access to the ADIOS library, handling initialisation and
 * finalisation.
 *
 * Usage
 * -----
 *
 * #include <bout/adios_object.hxx>
 *
 **************************************************************************/

#ifndef ADIOS_OBJECT_HXX
#define ADIOS_OBJECT_HXX

#include "bout/build_defines.hxx"

#if BOUT_HAS_ADIOS2

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
  bool isInStep = false; // true if BeginStep was called and EndStep was not yet called

  /** create or return the ADIOSStream based on the target file name */
  static ADIOSStream& ADIOSGetStream(const std::string& fname);

  ~ADIOSStream();

  template <class T>
  adios2::Variable<T> GetValueVariable(const std::string& varname) {
    auto v = io.InquireVariable<T>(varname);
    if (!v) {
      v = io.DefineVariable<T>(varname);
    }
    return v;
  }

  template <class T>
  adios2::Variable<T> GetArrayVariable(const std::string& varname, adios2::Dims& shape) {
    adios2::Variable<T> v = io.InquireVariable<T>(varname);
    if (!v) {
      adios2::Dims start(shape.size());
      v = io.DefineVariable<T>(varname, shape, start, shape);
    } else {
      v.SetShape(shape);
    }
    return v;
  }

private:
  ADIOSStream(const std::string fname) : fname(fname){};
  std::string fname;
};

/** Set user parameters for an IO group */
void ADIOSSetParameters(const std::string& input, const char delimKeyValue,
                        const char delimItem, adios2::IO& io);

} // namespace bout

#endif //BOUT_HAS_ADIOS2
#endif //ADIOS_OBJECT_HXX
