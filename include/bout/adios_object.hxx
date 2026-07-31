/*!\file*******************************************************************
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

#include "bout/assert.hxx"
#include "bout/boutcomm.hxx"
#include "bout/boutexception.hxx"
#include "bout/utils.hxx"

#include <adios2.h>
#include <memory>
#include <mpi.h>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

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
  static ADIOSStream& ADIOSGetStream(const std::string& fname, adios2::Mode mode,
                                     const std::string& engineType = "BP5");

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
  void Get(const std::string& varname, T& value, adios2::Mode mode = adios2::Mode::Sync) {
    auto variable = io.InquireVariable<T>(varname);
    ASSERT1(variable);
    ASSERT1(variable.ShapeID() == adios2::ShapeID::GlobalValue);
    engine().Get(variable, &value, mode);
  }

  template <class T>
  void Get(const std::string& varname, Array<T>& value,
           adios2::Mode mode = adios2::Mode::Sync) {
    GetArrayLike(varname, value, mode);
  }

  template <class T>
  void Get(const std::string& varname, Matrix<T>& value,
           adios2::Mode mode = adios2::Mode::Sync) {
    GetArrayLike(varname, value, mode);
  }

  template <class T>
  void Get(const std::string& varname, Tensor<T>& value,
           adios2::Mode mode = adios2::Mode::Sync) {
    GetArrayLike(varname, value, mode);
  }

  template <class T>
  void Put(const std::string& varname, T value) {
    if (BoutComm::rank() != 0) {
      return;
    }
    engine().Put(GetValueVariable<T>(varname), value);
  }

  template <class T>
  void Put(const std::string& varname, const Array<T>& value) {
    PutArrayLike(varname, value);
  }

  template <class T>
  void Put(const std::string& varname, const Matrix<T>& value) {
    PutArrayLike(varname, value);
  }

  template <class T>
  void Put(const std::string& varname, const Tensor<T>& value) {
    PutArrayLike(varname, value);
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
        throw BoutException("Could not open ADIOS file '{:s}'", fname);
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
      if (isInStep) {
        engine().EndStep();
        isInStep = false;
      }
      engine().Close();
      engine_ = adios2::Engine();
    }
  }

private:
  ADIOSStream(const std::string& fname, adios2::Mode mode, const std::string& engineType);

  template <class Container>
  void GetArrayLike(const std::string& varname, Container& value, adios2::Mode mode) {
    using T = typename Container::data_type;
    auto variable = io.InquireVariable<T>(varname);
    ASSERT1(variable);
    ASSERT1(variable.ShapeID() == adios2::ShapeID::GlobalArray);

    const auto shape = variable.Shape();
    auto dims_attr = io.InquireAttribute<std::string>(varname + "/__xarray_dimensions__");
    std::size_t offset = 0;

    if (dims_attr) {
      const auto dim_names = dims_attr.Data();
      if (!dim_names.empty() && dim_names[0] == "rank") {
        ASSERT1(!shape.empty());
        ASSERT1(static_cast<std::size_t>(BoutComm::rank()) < shape[0]);

        adios2::Dims start{static_cast<std::size_t>(BoutComm::rank())};
        adios2::Dims count{1};
        for (std::size_t i = 1; i < shape.size(); i++) {
          start.push_back(0);
          count.push_back(shape[i]);
        }

        variable.SetSelection(adios2::Box<adios2::Dims>{start, count});
        variable.SetMemorySelection(adios2::Box<adios2::Dims>{start, count});
        offset = 1;
      }
    }

    constexpr auto ndims = std::tuple_size_v<decltype(value.shape())>;
    ASSERT1(shape.size() == ndims + offset);

    resizeForShape(value, shape, offset, std::make_index_sequence<ndims>{});
    engine().Get(variable, value.begin(), mode);
  }

  template <class Container>
  void PutArrayLike(const std::string& varname, const Container& value) {
    using T = typename Container::data_type;
    auto var = GetArrayVariable<T>(varname, makeShape(BoutComm::size(), value),
                                   makeDimNames(value), BoutComm::rank());
    var.SetSelection(adios2::Box<adios2::Dims>{makeStart(value), makeShape(1, value)});
    engine().Put<T>(var, value.begin());
  }

  template <class Container>
  adios2::Dims makeShape(std::size_t first, const Container& value) const {
    return std::apply(
        [first](auto... sizes) {
          return adios2::Dims{first, static_cast<std::size_t>(sizes)...};
        },
        value.shape());
  }

  template <class Container>
  adios2::Dims makeStart(const Container& value) const {
    constexpr auto ndims = std::tuple_size_v<decltype(value.shape())>;
    adios2::Dims start(ndims + 1, 0);
    start[0] = static_cast<std::size_t>(BoutComm::rank());
    return start;
  }

  template <class Container>
  std::vector<std::string> makeDimNames(const Container& value) const {
    constexpr auto ndims = std::tuple_size_v<decltype(value.shape())>;
    std::vector<std::string> dim_names{"rank"};
    dim_names.reserve(ndims + 1);
    for (std::size_t i = 0; i < ndims; i++) {
      dim_names.push_back("dim_" + std::to_string(i));
    }
    return dim_names;
  }

  template <class Container, std::size_t... I>
  void resizeForShape(Container& value, const adios2::Dims& shape, std::size_t offset,
                      std::index_sequence<I...> /*indices*/) {
    value.reallocate(static_cast<typename Container::size_type>(shape[offset + I])...);
  }

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
