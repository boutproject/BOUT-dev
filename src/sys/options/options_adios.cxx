#include "bout/build_config.hxx"

#if BOUT_HAS_ADIOS

#include "bout/adios_object.hxx"
#include "bout/options_adios.hxx"

#include "bout/bout.hxx"
#include "bout/globals.hxx"
#include "bout/mesh.hxx"
#include "bout/sys/timer.hxx"

#include "adios2.h"
#include <exception>
#include <iostream>
#include <vector>

namespace bout {

OptionsADIOS::OptionsADIOS() : stream(nullptr) {}

OptionsADIOS::OptionsADIOS(std::string filename, FileMode mode)
    : OptionsIO(filename, mode), stream(nullptr) {}

OptionsADIOS::~OptionsADIOS() = default;
OptionsADIOS::OptionsADIOS(OptionsADIOS&&) noexcept = default;
OptionsADIOS& OptionsADIOS::operator=(OptionsADIOS&&) noexcept = default;

/// Name of the attribute used to track individual variable's time indices
constexpr auto current_time_index_name = "current_time_index";

template <class T>
bool readVariable(adios2::Engine& reader, adios2::IO& io, const std::string& name,
                  Options& result) {
  std::vector<T> data;
  adios2::Variable<T> variable = io.InquireVariable<T>(name);

  if (variable.ShapeID() == adios2::ShapeID::GlobalValue) {
    T value;
    reader.Get<T>(variable, &value, adios2::Mode::Sync);
    result[name] = value;
    return true;
  }

  if (variable.ShapeID() == adios2::ShapeID::LocalArray) {
    if (!BoutComm::rank()) {
      std::cout << "    LocalArray not supported" << std::endl;
    }
    return false;
  }

  auto dims = variable.Shape();
  auto ndims = dims.size();

  switch (ndims) {
  case 1: {
    Array<T> value(static_cast<int>(dims[0]));
    reader.Get<T>(variable, value.begin(), adios2::Mode::Sync);
    //result[name] = value;
    break;
  }
  case 2: {
    Matrix<T> value(static_cast<int>(dims[0]), static_cast<int>(dims[1]));
    reader.Get<T>(variable, value.begin(), adios2::Mode::Sync);
    //[name] = value;
    break;
  }
  case 3: {
    Tensor<T> value(static_cast<int>(dims[0]), static_cast<int>(dims[1]),
                    static_cast<int>(dims[2]));
    reader.Get<T>(variable, value.begin(), adios2::Mode::Sync);
    //result[name] = value;
    break;
  }
  }

  if (!BoutComm::rank())
    std::cout << "    array has " << variable.Steps() << " steps" << std::endl;

  /* Need to read the data here */
  result[name] = data.data();
  return true;
}

bool readVariable(adios2::Engine& reader, adios2::IO& io, const std::string& name,
                  const std::string& type, Options& result) {
  bool ret;
#define declare_template_instantiation(T)            \
  if (type == adios2::GetType<T>()) {                \
    ret = readVariable<T>(reader, io, name, result); \
  }
  ADIOS2_FOREACH_ATTRIBUTE_PRIMITIVE_STDTYPE_1ARG(declare_template_instantiation)
#undef declare_template_instantiation
  return ret;
}

bool readAttribute(adios2::IO& io, const std::string& name, const std::string& type,
                   Options& result) {
#define declare_template_instantiation(T)                  \
  if (type == adios2::GetType<T>()) {                      \
    adios2::Attribute<T> a = io.InquireAttribute<T>(name); \
    result[name] = a.Data().data();                        \
  }
  ADIOS2_FOREACH_ATTRIBUTE_PRIMITIVE_STDTYPE_1ARG(declare_template_instantiation)
#undef declare_template_instantiation
  return true;
}

Options OptionsADIOS::read() {
  Timer timer("io");

  // Open file
  ADIOSPtr adiosp = GetADIOSPtr();
  adios2::IO io;
  try {
    io = adiosp->AtIO(filename);
  } catch (const std::invalid_argument& e) {
    std::cerr << e.what() << '\n';
    io = adiosp->DeclareIO(filename);
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
    readVariable(reader, io, var_name, var_type, result);

    result[var_name].attributes["source"] = filename;

    // Get variable attributes
    for (const auto& attpair : io.AvailableAttributes(var_name, "/", true)) {
      const auto& att_name = attpair.first; // Attribute name
      const auto& att = attpair.second;     // attribute params
      auto it = att.find("Type");
      const std::string& att_type = it->second;
      readAttribute(io, att_name, att_type, result);
    }
  }

  reader.Close();

  return result;
}

/// Helper struct for returning errors from verifyTimesteps(NcGroup)
struct TimeDimensionError {
  std::string variable_name;
  std::string time_name;
  std::size_t expected_size;
  std::size_t current_size;
};

std::vector<TimeDimensionError> verifyTimesteps(adios2::Engine& reader) {
  reader.CurrentStep(); // just to avoid warning on unused variable

  // Variables with mismatched dimension sizes. Note that this might
  // be a little odd: if the first variable we come across has the
  // "wrong" dimension size, we will actually list all the others as
  // being wrong!
  std::vector<TimeDimensionError> errors;

  return errors;
}

void OptionsADIOS::verifyTimesteps() const {

  // Open file
  ADIOSPtr adiosp = GetADIOSPtr();
  adios2::IO io;
  try {
    io = adiosp->AtIO(filename);
  } catch (const std::invalid_argument& e) {
    std::cerr << e.what() << '\n';
    io = adiosp->DeclareIO(filename);
  }

  adios2::Engine reader = io.Open(filename, adios2::Mode::ReadRandomAccess);
  if (!reader) {
    throw BoutException("Could not open ADIOS file '{:s}' for reading", filename);
  }

  auto errors = bout::verifyTimesteps(reader);

  if (errors.empty()) {
    // No errors
    return;
  }

  std::string error_string;
  for (const auto& error : errors) {
    error_string += fmt::format(
        "  variable: {}; dimension: {}; expected size: {}; actual size: {}\n",
        error.variable_name, error.time_name, error.expected_size, error.current_size);
  }
  throw BoutException("ERROR: When checking timesteps in file '{}', some ({}) variables "
                      "did not have the expected size(s):\n{}",
                      filename, errors.size(), error_string);
}

/// Visit a variant type, and put the data into a NcVar
struct ADIOSPutVarVisitor {
  ADIOSPutVarVisitor(const std::string& name, ADIOSStream& stream)
      : varname(name), stream(stream) {}
  template <typename T>
  void operator()(const T& value) {
    adios2::Variable<T> var = stream.io.DefineVariable<T>(varname);
    stream.engine.Put(var, value);
  }

private:
  const std::string& varname;
  ADIOSStream& stream;
};

//template <>
//void ADIOSPutVarVisitor::operator()<bool>(const bool& value) {
//  int int_val = value ? 1 : 0;
//}

/// Visit a variant type, and put the data into a NcVar
struct ADIOSPutAttVisitor {
  ADIOSPutAttVisitor(const std::string& varname, const std::string& attrname,
                     ADIOSStream& stream)
      : varname(varname), attrname(attrname), stream(stream) {}
  template <typename T>
  void operator()(const T& value) {
    stream.io.DefineAttribute(attrname, value, varname, "/", false);
  }

private:
  const std::string& varname;
  const std::string& attrname;
  ADIOSStream& stream;
};

void writeGroup(const Options& options, ADIOSStream& stream, const std::string& groupname,
                const std::string& time_dimension) {

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
        if (time_it != child.attributes.end()) {
          // Has a time dimension

          const auto& time_name = bout::utils::get<std::string>(time_it->second);

          // Only write time-varying values that match current time
          // dimension being written
          if (time_name != time_dimension) {
            continue;
          }
        }

        // Write the variable
        std::string varname = groupname.empty() ? name : groupname + "/" + name;
        bout::utils::visit(ADIOSPutVarVisitor(varname, stream), child.value);

        // Write attributes
        for (const auto& attribute : child.attributes) {
          const std::string& att_name = attribute.first;
          const auto& att = attribute.second;

          bout::utils::visit(ADIOSPutAttVisitor(varname, att_name, stream), att);
        }

      } catch (const std::exception& e) {
        throw BoutException("Error while writing value '{:s}' : {:s}", name, e.what());
      }
    }
  }
}

/// Write options to file
void OptionsADIOS::write(const Options& options, const std::string& time_dim) {
  Timer timer("io");

  adios2::Mode mode =
      (file_mode == FileMode::replace) ? adios2::Mode::Write : adios2::Mode::Append;

  if (not stream) {
    stream = std::make_unique<ADIOSStream>();
  }

  // Open file
  ADIOSPtr adiosp = GetADIOSPtr();
  try {
    stream->io = adiosp->AtIO(filename);
  } catch (const std::invalid_argument& e) {
    std::cerr << e.what() << '\n';
    stream->io = adiosp->DeclareIO(filename);
  }

  stream->engine = stream->io.Open(filename, mode);
  if (!stream->engine) {
    throw BoutException("Could not open ADIOS file '{:s}' for writing", filename);
  }

  stream->engine.BeginStep();
  writeGroup(options, *stream, "", time_dim);

  stream->engine.EndStep();
}

} // namespace bout

#endif // BOUT_HAS_ADIOS
