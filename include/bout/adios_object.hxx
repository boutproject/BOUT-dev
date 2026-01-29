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

#include "bout/boutexception.hxx"

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
  adios2::Variable<double> vTime;
  adios2::Variable<int> vStep;
  int adiosStep = 0;

  /** create or return the ADIOSStream based on the target file name */
  static ADIOSStream& ADIOSGetStream(const std::string& fname, adios2::Mode mode);

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
  adios2::Variable<T>
  GetArrayVariable(const std::string& varname, const adios2::Dims& shape,
                   const std::vector<std::string>& dimNames, int rank) {
    adios2::Variable<T> v = io.InquireVariable<T>(varname);
    if (!v) {
      adios2::Dims start(shape.size());
      v = io.DefineVariable<T>(varname, shape, start, shape);
      if (!rank && dimNames.size()) {
        io.DefineAttribute<std::string>("__xarray_dimensions__", dimNames.data(),
                                        dimNames.size(), varname, "/", true);
      }
    } else {
      v.SetShape(shape);
    }
    return v;
  }

  auto engine() -> adios2::Engine& {
    if (not engine_) {
      engine_ = io.Open(fname, file_mode);
      if (not engine_) {
        throw BoutException("Could not open ADIOS file '{:s}' for writing", fname);
      }
    }
    return engine_;
  }

  void beginStep() {
    if (not isInStep) {
      engine().BeginStep();
      isInStep = true;
      adiosStep = static_cast<int>(engine().CurrentStep());
    }
  }

  void endStep() {
    if (isInStep) {
      engine().EndStep();
      isInStep = false;
    }
  }

  void finish() {
    if (engine_) {
      engine().EndStep();
      engine().Close();
    }
  }

private:
  ADIOSStream(const std::string& fname, adios2::Mode mode)
      : fname(fname), file_mode(mode) {

    ADIOSPtr adiosp = GetADIOSPtr();
    std::string ioname = "write_" + fname;
    try {
      io = adiosp->AtIO(ioname);
    } catch (const std::invalid_argument& e) {
      io = adiosp->DeclareIO(ioname);
      io.SetEngine("BP5");
    }
  };

  std::string fname;
  adios2::Mode file_mode;
  adios2::Engine engine_;

  /// true if BeginStep was called and EndStep was not yet called
  bool isInStep = false;
};

/** Set user parameters for an IO group */
void ADIOSSetParameters(const std::string& input, char delimKeyValue, char delimItem,
                        adios2::IO& io);

} // namespace bout

#endif //BOUT_HAS_ADIOS2
#endif //ADIOS_OBJECT_HXX
