#include "bout/build_defines.hxx"

#if BOUT_HAS_ADIOS2

#include "options_adios.hxx"

#include "bout/adios_object.hxx"
#include "bout/array.hxx"
#include "bout/assert.hxx"
#include "bout/bout_types.hxx"
#include "bout/boutexception.hxx"
#include "bout/field2d.hxx"
#include "bout/globals.hxx"
#include "bout/mesh.hxx"
#include "bout/options.hxx"
#include "bout/options_io.hxx"
#include "bout/output.hxx"
#include "bout/sys/timer.hxx"
#include "bout/sys/variant.hxx"
#include "bout/traits.hxx"
#include "bout/utils.hxx"

#include <adios2.h> // IWYU pragma: keep
#include <fmt/format.h>
#include <fmt/ranges.h>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <exception>
#include <iterator>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {
auto to_int_dims(const adios2::Dims& dims) {
  std::vector<int> int_dims;
  int_dims.reserve(dims.size());
  std::transform(dims.begin(), dims.end(), std::back_inserter(int_dims),
                 [](auto dim) { return static_cast<int>(dim); });
  return int_dims;
}

// Helper class to construct Adios hyperslices
struct Selection {
  // Offset of this processor's data into the global array
  adios2::Dims start;
  // The size of the mapped region
  adios2::Dims count;
  // Where the actual data starts in data pointer (to exclude ghost cells)
  adios2::Dims mem_start;
  // The actual size of data pointer in memory (including ghost cells)
  adios2::Dims mem_count;
  // Global shape, including boundaries but not guard cells
  adios2::Dims shape;
  // Shape of the local variable to read into
  std::vector<int> dims;

  // Distributed Field/Array/Matrix/Tensor
  bool should_set_selection{false};

  Selection(const std::vector<std::string>& dim_names, const std::vector<int>& dim_sizes,
            const Mesh& mesh) {
    const auto ndims = dim_names.size();
    const bool dim0_is_rank = ndims > 0 ? (dim_names[0] == "rank") : false;
    const bool dim0_is_x = ndims > 0 ? (dim_names[0] == "x") : false;
    const bool dim1_is_y = ndims > 1 ? (dim_names[1] == "y") : false;
    const bool dim1_is_z = ndims > 1 ? (dim_names[1] == "z") : false;
    const bool dim2_is_z = ndims > 2 ? (dim_names[2] == "z") : false;

    should_set_selection = dim0_is_rank or (ndims == 2 and dim0_is_x and dim1_is_y)
                           or (ndims == 2 and dim0_is_x and dim1_is_z)
                           or (ndims == 3 and dim0_is_x and dim1_is_y and dim2_is_z);

    if (dim0_is_rank) {
      const auto ndim_sizes = dim_sizes.size();
      ASSERT3(ndim_sizes > 1);

      // This is a distributed array, so the local variable is going
      // to be shape (dim_sizes[1]...) (that is, drop the rank)
      dims.push_back(dim_sizes[1]);
      // but we tell adios to read our rank's bit with the full ndims
      start = {static_cast<std::size_t>(BoutComm::rank()), 0};
      count = {std::size_t{1}, static_cast<std::size_t>(dim_sizes[0])};

      if (ndim_sizes > 2) {
        dims.push_back(dim_sizes[2]);
        start.push_back(0);
        count.push_back(dim_sizes[1]);
      }
      if (ndim_sizes > 3) {
        dims.push_back(dim_sizes[3]);
        start.push_back(0);
        count.push_back(dim_sizes[2]);
      }
      mem_count = count;
      mem_start = start;
      return;
    }

    if (ndims > 0) {
      dims.push_back(dim0_is_x ? mesh.LocalNx : dim_sizes[0]);
    }
    if (ndims > 1) {
      if (dim1_is_y) {
        dims.push_back(mesh.LocalNy);
      } else if (dim1_is_z) {
        dims.push_back(mesh.LocalNz);
      } else {
        dims.push_back(dim_sizes[1]);
      }
    }
    if (ndims > 2) {
      dims.push_back(dim2_is_z ? mesh.LocalNz : dim_sizes[2]);
    }

    shape.push_back(static_cast<std::size_t>(mesh.GlobalNx));
    start.push_back(static_cast<std::size_t>(mesh.MapGlobalX));
    count.push_back(static_cast<std::size_t>(mesh.MapCountX));
    mem_start.push_back(static_cast<std::size_t>(mesh.MapLocalX));
    mem_count.push_back(static_cast<std::size_t>(mesh.LocalNx));

    if (dim1_is_y) {
      shape.push_back(static_cast<std::size_t>(mesh.GlobalNy));
      start.push_back(static_cast<std::size_t>(mesh.MapGlobalY));
      count.push_back(static_cast<std::size_t>(mesh.MapCountY));
      mem_start.push_back(static_cast<std::size_t>(mesh.MapLocalY));
      mem_count.push_back(static_cast<std::size_t>(mesh.LocalNy));
    } else if (dim1_is_z) {
      shape.push_back(static_cast<std::size_t>(mesh.GlobalNz));
      start.push_back(static_cast<std::size_t>(mesh.MapGlobalZ));
      count.push_back(static_cast<std::size_t>(mesh.MapCountZ));
      mem_start.push_back(static_cast<std::size_t>(mesh.MapLocalZ));
      mem_count.push_back(static_cast<std::size_t>(mesh.LocalNz));
    }

    if (dim2_is_z) {
      shape.push_back(static_cast<std::size_t>(mesh.GlobalNz));
      start.push_back(static_cast<std::size_t>(mesh.MapGlobalZ));
      count.push_back(static_cast<std::size_t>(mesh.MapCountZ));
      mem_start.push_back(static_cast<std::size_t>(mesh.MapLocalZ));
      mem_count.push_back(static_cast<std::size_t>(mesh.LocalNz));
    }
  }

  auto selection() const { return adios2::Box<adios2::Dims>{start, count}; }
  auto memorySelection() const { return adios2::Box<adios2::Dims>{mem_start, mem_count}; }
};

template <template <class> class T, class U>
auto read_variable(adios2::IO& io, adios2::Engine& reader, const std::string& name,
                   const Selection& selection, T<U>& value) {
  auto variable = io.InquireVariable<U>(name);

  if (selection.should_set_selection) {
    variable.SetSelection(selection.selection());
    variable.SetMemorySelection(selection.memorySelection());
  }

  reader.Get<U>(variable, value.begin(), adios2::Mode::Sync);
  return Options(value);
}

template <class T>
Options readVariable(adios2::Engine& reader, adios2::IO& io, const std::string& name,
                     const std::string& type) {
  auto variable = io.InquireVariable<T>(name);

  using bout::globals::mesh;

  if (variable.ShapeID() == adios2::ShapeID::GlobalValue) {
    T value;
    reader.Get<T>(variable, &value, adios2::Mode::Sync);
    return Options(value);
  }

  if (variable.ShapeID() == adios2::ShapeID::LocalArray) {
    throw BoutException(
        "ADIOS reader did not implement reading local arrays like `{}` '{}' in file '{}'",
        type, name, reader.Name());
  }

  if (type != "double" && type != "float" && type != "int32_t") {
    throw BoutException("ADIOS reader did not implement reading arrays that are not "
                        "`double`/`float`/`int32_t` type. "
                        "Found `{}` '{}' in file '{}'",
                        type, name, reader.Name());
  }

  if (type == "double" && sizeof(BoutReal) != sizeof(double)) {
    throw BoutException(
        "ADIOS does not allow for implicit type conversions. BoutReal type is "
        "float but found `{}` '{}' in file '{}'",
        type, name, reader.Name());
  }

  if (type == "float" && sizeof(BoutReal) != sizeof(float)) {
    throw BoutException(
        "ADIOS reader does not allow for implicit type conversions. BoutReal type is "
        "double but found `{}` '{}' in file '{}'",
        type, name, reader.Name());
  }

  const auto dims = to_int_dims(variable.Shape());
  auto dims_attr =
      io.InquireAttribute<std::string>(fmt::format("{}/__xarray_dimensions__", name))
          .Data();

  const auto selection = Selection(dims_attr, dims, *mesh);
  const auto ndims = selection.dims.size();

  switch (ndims) {
  case 1: {
    if (type == "float" or type == "double") {
      Array<BoutReal> value(dims[0]);
      return read_variable(io, reader, name, selection, value);
    }
    if (type == "int32_t") {
      Array<int> value(dims[0]);
      return read_variable(io, reader, name, selection, value);
    }
    break;
  }
  case 2: {
    if (type == "float" or type == "double") {
      Matrix<BoutReal> value(selection.dims[0], selection.dims[1]);
      return read_variable(io, reader, name, selection, value);
    }
    if (type == "int32_t") {
      Matrix<int> value(selection.dims[0], selection.dims[1]);
      return read_variable(io, reader, name, selection, value);
    }
    break;
  }
  case 3: {
    if (type == "float" or type == "double") {
      Tensor<BoutReal> value(selection.dims[0], selection.dims[1], selection.dims[2]);
      return read_variable(io, reader, name, selection, value);
    }
    if (type == "int32_t") {
      Tensor<int> value(selection.dims[0], selection.dims[1], selection.dims[2]);
      return read_variable(io, reader, name, selection, value);
    }
    break;
  }
  // Throw below
  default:
    break;
  }
  auto dims_str = fmt::format("[{}]", fmt::join(dims, ", "));
  throw BoutException(
      "ADIOS reader failed to read '{}' (shape: {}, type: '{}') in file '{}'", name,
      dims_str, type, reader.Name());
}

Options readVariable(adios2::Engine& reader, adios2::IO& io, const std::string& name,
                     const std::string& type) {
  if (type == adios2::GetType<std::string>()) {
    return readVariable<std::string>(reader, io, name, type);
  }
  if (type == adios2::GetType<char>()) {
    return readVariable<char>(reader, io, name, type);
  }
  if (type == adios2::GetType<int8_t>()) {
    return readVariable<int8_t>(reader, io, name, type);
  }
  if (type == adios2::GetType<int16_t>()) {
    return readVariable<int16_t>(reader, io, name, type);
  }
  if (type == adios2::GetType<int32_t>()) {
    return readVariable<int32_t>(reader, io, name, type);
  }
  if (type == adios2::GetType<int64_t>()) {
    return readVariable<int64_t>(reader, io, name, type);
  }
  if (type == adios2::GetType<uint8_t>()) {
    return readVariable<uint8_t>(reader, io, name, type);
  }
  if (type == adios2::GetType<uint16_t>()) {
    return readVariable<uint16_t>(reader, io, name, type);
  }
  if (type == adios2::GetType<uint32_t>()) {
    return readVariable<uint32_t>(reader, io, name, type);
  }
  if (type == adios2::GetType<uint64_t>()) {
    return readVariable<uint64_t>(reader, io, name, type);
  }
  if (type == adios2::GetType<float>()) {
    return readVariable<float>(reader, io, name, type);
  }
  if (type == adios2::GetType<double>()) {
    return readVariable<double>(reader, io, name, type);
  }
  if (type == adios2::GetType<long double>()) {
    return readVariable<long double>(reader, io, name, type);
  }

  output_warn.write("ADIOS readVariable can't read type '{}' (variable '{}')", type,
                    name);
  return Options{};
}

bool readAttribute(adios2::IO& io, const std::string& name, const std::string& type,
                   Options& result) {
  // Attribute is the part of 'name' after the last '/' separator
  std::string attrname = name;
  auto pos = name.find_last_of('/');
  if (pos != std::string::npos) {
    attrname = name.substr(pos + 1);
  }

  if (type == adios2::GetType<int>()) {
    result.attributes[attrname] = *io.InquireAttribute<int>(name).Data().data();
    return true;
  }
  if (type == adios2::GetType<BoutReal>()) {
    result.attributes[attrname] = *io.InquireAttribute<BoutReal>(name).Data().data();
    return true;
  }
  if (type == adios2::GetType<std::string>()) {
    result.attributes[attrname] = *io.InquireAttribute<std::string>(name).Data().data();
    return true;
  }

  output_warn.write("ADIOS readAttribute can't read type '{}' (variable '{}')", type,
                    name);
  return false;
}
} // namespace

namespace bout {
OptionsADIOS::OptionsADIOS(Options& options) : OptionsIO(options) {
  if (options["file"].doc("File name. Defaults to <path>/<prefix>.pb").isSet()) {
    filename = options["file"].as<std::string>();
  } else {
    // Both path and prefix must be set
    filename = fmt::format("{}/{}.bp", options["path"].as<std::string>(),
                           options["prefix"].as<std::string>());
  }

  file_mode = (options["append"].doc("Append to existing file?").withDefault<bool>(false))
                  ? adios2::Mode::Append
                  : adios2::Mode::Write;

  singleWriteFile = options["singleWriteFile"].withDefault<bool>(false);
}

Options OptionsADIOS::read([[maybe_unused]] bool lazy) {
  const Timer timer("io");

  // Open file
  ADIOSPtr const adiosp = GetADIOSPtr();
  adios2::IO io;
  const std::string ioname = "read_" + filename;
  try {
    io = adiosp->AtIO(ioname);
  } catch (const std::invalid_argument& e) {
    io = adiosp->DeclareIO(ioname);
  }

  adios2::Engine reader = io.Open(filename, adios2::Mode::ReadRandomAccess);
  if (!reader) {
    throw BoutException("Could not open ADIOS file '{:s}' for reading", filename);
  }

  Options result;

  // Iterate over all variables
  for (const auto& varpair : io.AvailableVariables()) {
    const auto& var_name = varpair.first; // Name of the variable

    auto it = varpair.second.find("Type");
    const std::string& var_type = it->second;

    Options* varptr = &result;
    for (const auto& piece : strsplit(var_name, '/')) {
      varptr = &(*varptr)[piece]; // Navigate to subsection if needed
    }
    Options& var = *varptr;

    // Note: Copying the value rather than simple assignment is used
    // because the Options assignment operator overwrites full_name.
    var = 0; // Setting is_section to false
    var.value = readVariable(reader, io, var_name, var_type).value;
    var.attributes["source"] = filename;

    // Get variable attributes
    for (const auto& attpair : io.AvailableAttributes(var_name, "/", true)) {
      const auto& att_name = attpair.first; // Attribute name
      const auto& att = attpair.second;     // attribute params

      auto it = att.find("Type");
      const std::string& att_type = it->second;
      readAttribute(io, att_name, att_type, var);
    }
  }

  reader.Close();

  return result;
}

void OptionsADIOS::verifyTimesteps() const {
  // This doesn't _verify_ the timesteps, but does flush to disk.
  // Maybe this is fine, because ADIOS2 doesn't require every variable
  // to be in sync?
  if (singleWriteFile) {
    return;
  }

  ADIOSStream& stream = ADIOSStream::ADIOSGetStream(filename, file_mode);

  stream.endStep();
}
} // namespace bout

const std::vector<std::string> DIMS_X = {"x"};
const std::vector<std::string> DIMS_XY = {"x", "y"};
const std::vector<std::string> DIMS_XZ = {"x", "z"};
const std::vector<std::string> DIMS_XYZ = {"x", "y", "z"};

namespace {
using bout::utils::tuple_index_sequence;

template <class Tuple, std::size_t... I>
auto make_shape_impl(std::size_t first, const Tuple& t,
                     std::index_sequence<I...> /* index */) {
  return adios2::Dims{first, static_cast<std::size_t>(std::get<I>(t))...};
}
// Return an `adios2::Dims` with value ``{first, value.shape()[0]...}``
template <class T>
auto make_shape(std::size_t first, const T& value) {
  const auto shape = value.shape();
  return make_shape_impl(first, shape, tuple_index_sequence<decltype(shape)>{});
}

template <class Tuple, std::size_t... I>
auto make_start_impl(const Tuple& t, std::index_sequence<I...> /* index */) {
  // Hey look, a legitimate use of the comma operator to get a bunch
  // of zeros the length of the index_sequence!
  return adios2::Dims{static_cast<std::size_t>(BoutComm::rank()),
                      (std::get<I>(t), std::size_t{0})...};
}
// Return an `adios2::Dims` with value ``{rank, 0...}``, with as many zeros as the dimension of ``T``
template <class T>
auto make_start(const T& value) {
  const auto shape = value.shape();
  return make_start_impl(shape, tuple_index_sequence<decltype(shape)>{});
}

template <std::size_t... I>
auto make_dims_impl(std::index_sequence<I...> /*index*/) {
  using namespace std::string_literals;
  return std::vector{"rank"s, (fmt::format("dim_{}", I))...};
}
// Return vector of dimension names: `{"rank", "dim_0", ...}`
template <class T>
auto make_dims(const T& value) {
  return make_dims_impl(tuple_index_sequence<decltype(value.shape())>{});
}

/// Visit a variant type, and put the data into a NcVar
struct ADIOSPutVarVisitor {
  ADIOSPutVarVisitor(const std::string& name, bout::ADIOSStream& stream)
      : varname(name), stream(stream) {}
  template <typename T>
  void operator()(const T& value) {
    adios2::Variable<T> var = stream.GetValueVariable<T>(varname);
    stream.engine().Put(var, value);
  }

  void operator()(const Array<int>& value) { amt_put_helper(value); }
  void operator()(const Array<BoutReal>& value) { amt_put_helper(value); }
  void operator()(const Matrix<int>& value) { amt_put_helper(value); }
  void operator()(const Matrix<BoutReal>& value) { amt_put_helper(value); }
  void operator()(const Tensor<int>& value) { amt_put_helper(value); }
  void operator()(const Tensor<BoutReal>& value) { amt_put_helper(value); }
  void operator()(bool value) {
    // Scalars are only written from processor 0
    if (BoutComm::rank() != 0) {
      return;
    }
    stream.engine().Put(stream.GetValueVariable<int>(varname), static_cast<int>(value));
  }
  void operator()(int value) {
    // Scalars are only written from processor 0
    if (BoutComm::rank() != 0) {
      return;
    }
    stream.engine().Put(stream.GetValueVariable<int>(varname), value);
  }
  void operator()(BoutReal value) {
    // Scalars are only written from processor 0
    if (BoutComm::rank() != 0) {
      return;
    }
    stream.engine().Put(stream.GetValueVariable<BoutReal>(varname), value);
  }

  void operator()(const std::string& value) {
    // Scalars are only written from processor 0
    if (BoutComm::rank() != 0) {
      return;
    }
    stream.engine().Put<std::string>(stream.GetValueVariable<std::string>(varname), value,
                                     adios2::Mode::Sync);
  }

  void operator()(const Field2D& value) {
    // Empty dim_sizes is fine because we provide full list of dim_names
    auto selection = Selection(DIMS_XY, {}, *value.getMesh());
    auto var = stream.GetArrayVariable<BoutReal>(varname, selection.shape, DIMS_XY,
                                                 BoutComm::rank());
    var.SetSelection(selection.selection());
    var.SetMemorySelection(selection.memorySelection());
    stream.engine().Put(var, &value(0, 0));
  }

  void operator()(const Field3D& value) {
    // Empty dim_sizes is fine because we provide full list of dim_names
    auto selection = Selection(DIMS_XYZ, {}, *value.getMesh());
    auto var = stream.GetArrayVariable<BoutReal>(varname, selection.shape, DIMS_XYZ,
                                                 BoutComm::rank());
    var.SetSelection(selection.selection());
    var.SetMemorySelection(selection.memorySelection());
    stream.engine().Put(var, &value(0, 0, 0));
  }

  void operator()(const FieldPerp& value) {
    // Empty dim_sizes is fine because we provide full list of dim_names
    auto selection = Selection(DIMS_XZ, {}, *value.getMesh());
    auto var = stream.GetArrayVariable<BoutReal>(varname, selection.shape, DIMS_XZ,
                                                 BoutComm::rank());
    var.SetSelection(selection.selection());
    var.SetMemorySelection(selection.memorySelection());
    stream.engine().Put<BoutReal>(var, &value(0, 0));
  }

private:
  // Ok to keep const/refs here as visitor is only a temporary
  const std::string& varname; // NOLINT(*-avoid-const-or-ref-data-members)
  bout::ADIOSStream& stream;  // NOLINT(*-avoid-const-or-ref-data-members)

  // helper for `Array`, `Matrix`, `Tensor`
  template <template <class> class T, class U>
  void amt_put_helper(const T<U>& value) {
    const auto shape = make_shape(static_cast<std::size_t>(BoutComm::size()), value);
    auto var =
        stream.GetArrayVariable<U>(varname, shape, make_dims(value), BoutComm::rank());
    var.SetSelection({make_start(value), make_shape(1, value)});
    stream.engine().Put<U>(var, value.begin());
  }
};

/// Visit a variant type, and put the data into a NcVar
struct ADIOSPutAttVisitor {
  ADIOSPutAttVisitor(const std::string& varname, const std::string& attrname,
                     bout::ADIOSStream& stream)
      : varname(varname), attrname(attrname), stream(stream) {}
  template <typename T>
  void operator()(const T& value) {
    stream.io.DefineAttribute<T>(attrname, value, varname, "/", false);
  }

  void operator()(bool value) {
    stream.io.DefineAttribute<int>(attrname, static_cast<int>(value), varname, "/",
                                   false);
  }

private:
  // Ok to keep const/refs here as visitor is only a temporary
  const std::string& varname;  // NOLINT(*-avoid-const-or-ref-data-members)
  const std::string& attrname; // NOLINT(*-avoid-const-or-ref-data-members)
  bout::ADIOSStream& stream;   // NOLINT(*-avoid-const-or-ref-data-members)
};

void writeGroup(const Options& options, bout::ADIOSStream& stream,
                const std::string& groupname, const std::string& time_dimension) {

  for (const auto& childpair : options.getChildren()) {
    const auto& name = childpair.first;
    const auto& child = childpair.second;

    if (child.isSection()) {
      TRACE("Writing group '{:s}'", name);
      writeGroup(child, stream, name, time_dimension);
      continue;
    }

    if (child.isValue()) {
      try {
        auto time_it = child.attributes.find("time_dimension");
        if (time_it == child.attributes.end()) {
          if (stream.adiosStep > 0) {
            // we should only write the non-varying values in the first step
            continue;
          }
        } else {
          // Has a time dimension

          const auto& time_name = bout::utils::get<std::string>(time_it->second);

          // Only write time-varying values that match current time
          // dimension being written
          if (time_name != time_dimension) {
            continue;
          }
        }

        // Write the variable
        // Note: ADIOS2 uses '/' to as a group separator; BOUT++ uses ':'
        const std::string varname =
            groupname.empty() ? name : fmt::format("{}/{}", groupname, name);
        bout::utils::visit(ADIOSPutVarVisitor(varname, stream), child.value);

        // Write attributes
        if (BoutComm::rank() == 0) {
          for (const auto& attribute : child.attributes) {
            const std::string& att_name = attribute.first;
            const auto& att = attribute.second;

            bout::utils::visit(ADIOSPutAttVisitor(varname, att_name, stream), att);
          }
        }

      } catch (const std::exception& e) {
        throw BoutException("Error while writing value '{:s}' : {:s}", name, e.what());
      }
    }
  }
}
} // namespace

namespace bout {
/// Write options to file
void OptionsADIOS::write(const Options& options, const std::string& time_dim) {
  const Timer timer("io");

  // ADIOSStream is just a BOUT++ object, it does not create anything inside ADIOS
  ADIOSStream& stream = ADIOSStream::ADIOSGetStream(filename, file_mode);

  // Multiple write() calls allowed in a single adios step to output multiple
  // Options objects in the same step. verifyTimesteps() will indicate the
  // completion of the step (and adios will publish the step).
  stream.beginStep();

  writeGroup(options, stream, "", time_dim);

  // In singleWriteFile mode, we complete the step and close the file
  if (singleWriteFile) {
    stream.finish();
  }
}

} // namespace bout

#endif // BOUT_HAS_ADIOS2
