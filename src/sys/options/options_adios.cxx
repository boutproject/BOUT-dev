#include "bout/build_defines.hxx"

#if BOUT_HAS_ADIOS2

#include "options_adios.hxx"

#include "bout/adios_object.hxx"
#include "bout/array.hxx"
#include "bout/assert.hxx"
#include "bout/bout_types.hxx"
#include "bout/boutexception.hxx"
#include "bout/field2d.hxx"
#include "bout/options.hxx"
#include "bout/options_io.hxx"
#include "bout/output.hxx"
#include "bout/sys/timer.hxx"
#include "bout/sys/variant.hxx"
#include "bout/utils.hxx"

#include <adios2.h> // IWYU pragma: keep
#include <fmt/format.h>
#include <fmt/ranges.h>

#include <cstddef>
#include <cstdint>
#include <exception>
#include <stdexcept>
#include <string>

namespace {
std::size_t getLocalNDims(adios2::IO& io, const std::string& name,
                          const adios2::Dims& shape) {
  auto dims_attr = io.InquireAttribute<std::string>(name + "/__xarray_dimensions__");
  if (dims_attr) {
    const auto dim_names = dims_attr.Data();
    if (!dim_names.empty() && dim_names[0] == "rank") {
      ASSERT1(!shape.empty());
      return shape.size() - 1;
    }
  }
  return shape.size();
}

template <class T>
Options readArrayVariable(bout::ADIOSStream& stream, const std::string& name, T& value) {
  try {
    stream.Get(name, value);
  } catch (const std::exception& e) {
    output_warn.write("ADIOS exception while reading '{}': {}\n", name, e.what());
    return Options{};
  }
  return Options(value);
}

template <class T>
Options readVariable(bout::ADIOSStream& stream, const std::string& name,
                     const std::string& type) {
  auto variable = stream.io.InquireVariable<T>(name);

  if (variable.ShapeID() == adios2::ShapeID::GlobalValue) {
    T value;
    stream.Get(name, value);
    return Options(value);
  }

  if (variable.ShapeID() == adios2::ShapeID::LocalArray) {
    throw BoutException("ADIOS reader did not implement reading local arrays like `{}` "
                        "'{}' in file '{}'\n",
                        type, name, stream.engine().Name());
  }

  if (type != "double" && type != "float" && type != "int32_t" && type != "int64_t") {
    throw BoutException("ADIOS reader did not implement reading arrays that are not "
                        "`double`/`float`/`int32_t`/`int64_t` type. "
                        "Found `{}` '{}' in file '{}'\n",
                        type, name, stream.engine().Name());
  }

  if (type == "double" && sizeof(BoutReal) != sizeof(double)) {
    throw BoutException(
        "ADIOS does not allow for implicit type conversions. BoutReal type is "
        "float but found `{}` '{}' in file '{}'\n",
        type, name, stream.engine().Name());
  }

  if (type == "float" && sizeof(BoutReal) != sizeof(float)) {
    throw BoutException(
        "ADIOS reader does not allow for implicit type conversions. BoutReal type is "
        "double but found `{}` '{}' in file '{}'",
        type, name, stream.engine().Name());
  }

  const auto shape = variable.Shape();
  const auto ndims = getLocalNDims(stream.io, name, shape);

  switch (ndims) {
  case 1: {
    if (type == "float" or type == "double") {
      Array<BoutReal> value;
      return readArrayVariable(stream, name, value);
    }
    if (type == "int32_t" or type == "int64_t") {
      Array<int> value;
      return readArrayVariable(stream, name, value);
    }
    break;
  }
  case 2: {
    if (type == "float" or type == "double") {
      Matrix<BoutReal> value;
      return readArrayVariable(stream, name, value);
    }
    if (type == "int32_t" or type == "int64_t") {
      Matrix<int> value;
      return readArrayVariable(stream, name, value);
    }
    break;
  }
  case 3: {
    if (type == "float" or type == "double") {
      Tensor<BoutReal> value;
      return readArrayVariable(stream, name, value);
    }
    if (type == "int32_t" or type == "int64_t") {
      Tensor<int> value;
      return readArrayVariable(stream, name, value);
    }
    break;
  }
  // Throw below
  default:
    break;
  }
  auto dims_str = fmt::format("[{}]", fmt::join(shape, ", "));
  throw BoutException(
      "ADIOS reader failed to read '{}' (shape: {}, type: '{}') in file '{}'\n", name,
      dims_str, type, stream.engine().Name());
}

Options readVariable(bout::ADIOSStream& stream, const std::string& name,
                     const std::string& type) {
  if (type == adios2::GetType<std::string>()) {
    return readVariable<std::string>(stream, name, type);
  }
  if (type == adios2::GetType<char>()) {
    return readVariable<char>(stream, name, type);
  }
  if (type == adios2::GetType<int8_t>()) {
    return readVariable<int8_t>(stream, name, type);
  }
  if (type == adios2::GetType<int16_t>()) {
    return readVariable<int16_t>(stream, name, type);
  }
  if (type == adios2::GetType<int32_t>()) {
    return readVariable<int32_t>(stream, name, type);
  }
  if (type == adios2::GetType<int64_t>()) {
    return readVariable<int64_t>(stream, name, type);
  }
  if (type == adios2::GetType<uint8_t>()) {
    return readVariable<uint8_t>(stream, name, type);
  }
  if (type == adios2::GetType<uint16_t>()) {
    return readVariable<uint16_t>(stream, name, type);
  }
  if (type == adios2::GetType<uint32_t>()) {
    return readVariable<uint32_t>(stream, name, type);
  }
  if (type == adios2::GetType<uint64_t>()) {
    return readVariable<uint64_t>(stream, name, type);
  }
  if (type == adios2::GetType<float>()) {
    return readVariable<float>(stream, name, type);
  }
  if (type == adios2::GetType<double>()) {
    return readVariable<double>(stream, name, type);
  }
  if (type == adios2::GetType<long double>()) {
    return readVariable<long double>(stream, name, type);
  }

  output_warn.write("ADIOS readVariable can't read type '{}' (variable '{}')\n", type,
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

  output_warn.write("ADIOS readAttribute can't read type '{}' (variable '{}')\n", type,
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
  ADIOSStream& stream =
      ADIOSStream::ADIOSGetStream(filename, adios2::Mode::ReadRandomAccess);
  adios2::Engine& reader = stream.engine();
  if (!reader) {
    throw BoutException("Could not open ADIOS file '{:s}' for reading\n", filename);
  }

  Options result;

  // Iterate over all variables
  for (const auto& varpair : stream.io.AvailableVariables()) {
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
    var.value = readVariable(stream, var_name, var_type).value;
    var.attributes["source"] = filename;

    // Get variable attributes
    for (const auto& attpair : stream.io.AvailableAttributes(var_name, "/", true)) {
      const auto& att_name = attpair.first; // Attribute name
      const auto& att = attpair.second;     // attribute params

      auto it = att.find("Type");
      const std::string& att_type = it->second;
      readAttribute(stream.io, att_name, att_type, var);
    }
  }

  stream.finish();

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

namespace {
/// Visit a variant type, and put the data into a NcVar
struct ADIOSPutVarVisitor {
  ADIOSPutVarVisitor(const std::string& name, bout::ADIOSStream& stream)
      : varname(name), stream(stream) {}
  template <typename T>
  void operator()(const T& value) {
    stream.Put(varname, value);
  }

  void operator()(bool value) { stream.Put(varname, static_cast<int>(value)); }

private:
  // Ok to keep const/refs here as visitor is only a temporary
  const std::string& varname; // NOLINT(*-avoid-const-or-ref-data-members)
  bout::ADIOSStream& stream;  // NOLINT(*-avoid-const-or-ref-data-members)
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
