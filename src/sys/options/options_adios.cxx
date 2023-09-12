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

OptionsADIOS::OptionsADIOS() {}

OptionsADIOS::OptionsADIOS(std::string filename, FileMode mode)
    : OptionsIO(filename, mode) {}

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
  bool ret = false;
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
  std::string ioname = "read_" + filename;
  try {
    io = adiosp->AtIO(ioname);
  } catch (const std::invalid_argument& e) {
    io = adiosp->DeclareIO(ioname);
  }

  std::cout << "OptionsADIOS::read: open " << filename << std::endl;
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
  ADIOSStream& stream = ADIOSStream::ADIOSGetStream(filename);
  std::cout << "OptionsADIOS::write: END adios step = " << stream.engine.CurrentStep()
            << std::endl;
  stream.engine.EndStep();
  stream.isInStep = false;
  return;
}

/// Visit a variant type, and put the data into a NcVar
struct ADIOSPutVarVisitor {
  ADIOSPutVarVisitor(const std::string& name, ADIOSStream& stream)
      : varname(name), stream(stream) {}
  template <typename T>
  void operator()(const T& value) {
    adios2::Variable<T> var = stream.GetValueVariable<T>(varname);
    stream.engine.Put(var, value);
  }

private:
  const std::string& varname;
  ADIOSStream& stream;
};

template <>
void ADIOSPutVarVisitor::operator()<bool>(const bool& value) {
  // Scalars are only written from processor 0
  if (BoutComm::rank() != 0) {
    return;
  }
  adios2::Variable<int> var = stream.GetValueVariable<int>(varname);
  stream.engine.Put<int>(var, static_cast<int>(value));
}

template <>
void ADIOSPutVarVisitor::operator()<int>(const int& value) {
  // Scalars are only written from processor 0
  if (BoutComm::rank() != 0) {
    return;
  }
  adios2::Variable<int> var = stream.GetValueVariable<int>(varname);
  stream.engine.Put<int>(var, value);
}

template <>
void ADIOSPutVarVisitor::operator()<BoutReal>(const BoutReal& value) {
  // Scalars are only written from processor 0
  if (BoutComm::rank() != 0) {
    return;
  }
  adios2::Variable<BoutReal> var = stream.GetValueVariable<BoutReal>(varname);
  stream.engine.Put<BoutReal>(var, value);
}

template <>
void ADIOSPutVarVisitor::operator()<std::string>(const std::string& value) {
  // Scalars are only written from processor 0
  if (BoutComm::rank() != 0) {
    return;
  }
  adios2::Variable<std::string> var = stream.GetValueVariable<std::string>(varname);
  stream.engine.Put<std::string>(var, value, adios2::Mode::Sync);
}

template <>
void ADIOSPutVarVisitor::operator()<Field2D>(const Field2D& value) {
  // Get the mesh that describes how local data relates to global arrays
  auto mesh = value.getMesh();

  // The global size of this array includes boundary cells but not communication guard cells.
  // In general this array will be sparse because it may have gaps.
  adios2::Dims shape = {static_cast<size_t>(mesh->GlobalNx),
                        static_cast<size_t>(mesh->GlobalNy)};

  // Offset of this processor's data into the global array
  adios2::Dims start = {static_cast<size_t>(mesh->MapGlobalX),
                        static_cast<size_t>(mesh->MapGlobalY)};

  // The size of the mapped region
  adios2::Dims count = {static_cast<size_t>(mesh->MapCountX),
                        static_cast<size_t>(mesh->MapCountY)};
  adios2::Variable<BoutReal> var = stream.GetArrayVariable<BoutReal>(varname, shape);
  var.SetSelection({start, count});

  // Get a Span of an internal engine buffer that we can fill with data
  auto span = stream.engine.Put<BoutReal>(var);

  // Iterate over the span and the Field2D array
  // Note: The Field2D includes communication guard cells that are not written
  auto it = span.begin();
  for (int x = mesh->MapLocalX; x < (mesh->MapLocalX + mesh->MapCountX); ++x) {
    for (int y = mesh->MapLocalY; y < (mesh->MapLocalY + mesh->MapCountY); ++y) {
      *it = value(x, y);
      ++it;
    }
  }
  ASSERT1(it == span.end());
}

template <>
void ADIOSPutVarVisitor::operator()<Field3D>(const Field3D& value) {
  // Get the mesh that describes how local data relates to global arrays
  auto mesh = value.getMesh();

  // The global size of this array includes boundary cells but not communication guard cells.
  // In general this array will be sparse because it may have gaps.
  adios2::Dims shape = {static_cast<size_t>(mesh->GlobalNx),
                        static_cast<size_t>(mesh->GlobalNy),
                        static_cast<size_t>(mesh->GlobalNz)};

  // Offset of this processor's data into the global array
  adios2::Dims start = {static_cast<size_t>(mesh->MapGlobalX),
                        static_cast<size_t>(mesh->MapGlobalY),
                        static_cast<size_t>(mesh->MapGlobalZ)};

  // The size of the mapped region
  adios2::Dims count = {static_cast<size_t>(mesh->MapCountX),
                        static_cast<size_t>(mesh->MapCountY),
                        static_cast<size_t>(mesh->MapCountZ)};

  adios2::Variable<BoutReal> var = stream.GetArrayVariable<BoutReal>(varname, shape);
  var.SetSelection({start, count});
  // Get a Span of an internal engine buffer.
  auto span = stream.engine.Put<BoutReal>(var);

  // Iterate over the span and the Field3D array
  // Note: The Field3D includes communication guard cells that are not written
  auto it = span.begin();
  for (int x = mesh->MapLocalX; x < (mesh->MapLocalX + mesh->MapCountX); ++x) {
    for (int y = mesh->MapLocalY; y < (mesh->MapLocalY + mesh->MapCountY); ++y) {
      for (int z = mesh->MapLocalZ; z < (mesh->MapLocalZ + mesh->MapCountZ); ++z) {
        *it = value(x, y, z);
        ++it;
      }
    }
  }
  ASSERT1(it == span.end());
}

template <>
void ADIOSPutVarVisitor::operator()<FieldPerp>(const FieldPerp& value) {
  // Get the mesh that describes how local data relates to global arrays
  auto mesh = value.getMesh();

  // The global size of this array includes boundary cells but not communication guard cells.
  // In general this array will be sparse because it may have gaps.
  adios2::Dims shape = {static_cast<size_t>(mesh->GlobalNx),
                        static_cast<size_t>(mesh->GlobalNz)};

  // Offset of this processor's data into the global array
  adios2::Dims start = {static_cast<size_t>(mesh->MapGlobalX),
                        static_cast<size_t>(mesh->MapGlobalZ)};

  // The size of the mapped region
  adios2::Dims count = {static_cast<size_t>(mesh->MapCountX),
                        static_cast<size_t>(mesh->MapCountZ)};

  adios2::Variable<BoutReal> var = stream.GetArrayVariable<BoutReal>(varname, shape);
  var.SetSelection({start, count});

  // Get a Span of an internal engine buffer.
  auto span = stream.engine.Put<BoutReal>(var);

  // Iterate over the span and the Field3D array
  // Note: The Field3D includes communication guard cells that are not written
  auto it = span.begin();
  for (int x = mesh->MapLocalX; x < (mesh->MapLocalX + mesh->MapCountX); ++x) {
    for (int z = mesh->MapLocalZ; z < (mesh->MapLocalZ + mesh->MapCountZ); ++z) {
      *it = value(x, z);
      ++it;
    }
  }
  ASSERT1(it == span.end());
}

template <>
void ADIOSPutVarVisitor::operator()<Array<BoutReal>>(const Array<BoutReal>& value) {
  // Pointer to data. Assumed to be contiguous array
  adios2::Dims shape = {(size_t)BoutComm::size(), (size_t)value.size()};
  adios2::Dims start = {(size_t)BoutComm::rank(), 0};
  adios2::Dims count = {1, shape[1]};
  adios2::Variable<BoutReal> var = stream.GetArrayVariable<BoutReal>(varname, shape);
  var.SetSelection({start, count});
  stream.engine.Put<BoutReal>(var, value.begin());
}

template <>
void ADIOSPutVarVisitor::operator()<Matrix<BoutReal>>(const Matrix<BoutReal>& value) {
  // Pointer to data. Assumed to be contiguous array
  auto s = value.shape();
  adios2::Dims shape = {(size_t)BoutComm::size(), (size_t)std::get<0>(s),
                        (size_t)std::get<1>(s)};
  adios2::Dims start = {(size_t)BoutComm::rank(), 0, 0};
  adios2::Dims count = {1, shape[1], shape[2]};
  adios2::Variable<BoutReal> var = stream.GetArrayVariable<BoutReal>(varname, shape);
  var.SetSelection({start, count});
  stream.engine.Put<BoutReal>(var, value.begin());
}

template <>
void ADIOSPutVarVisitor::operator()<Tensor<BoutReal>>(const Tensor<BoutReal>& value) {
  // Pointer to data. Assumed to be contiguous array
  auto s = value.shape();
  adios2::Dims shape = {(size_t)BoutComm::size(), (size_t)std::get<0>(s),
                        (size_t)std::get<1>(s), (size_t)std::get<2>(s)};
  adios2::Dims start = {(size_t)BoutComm::rank(), 0, 0, 0};
  adios2::Dims count = {1, shape[1], shape[2], shape[3]};
  adios2::Variable<BoutReal> var = stream.GetArrayVariable<BoutReal>(varname, shape);
  var.SetSelection({start, count});
  stream.engine.Put<BoutReal>(var, value.begin());
}

/// Visit a variant type, and put the data into a NcVar
struct ADIOSPutAttVisitor {
  ADIOSPutAttVisitor(const std::string& varname, const std::string& attrname,
                     ADIOSStream& stream)
      : varname(varname), attrname(attrname), stream(stream) {}
  template <typename T>
  void operator()(const T& value) {
    stream.io.DefineAttribute<T>(attrname, value, varname, "/", false);
  }

private:
  const std::string& varname;
  const std::string& attrname;
  ADIOSStream& stream;
};

template <>
void ADIOSPutAttVisitor::operator()<bool>(const bool& value) {
  stream.io.DefineAttribute<int>(attrname, (int)value, varname, "/", false);
}

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
        std::string varname = groupname.empty() ? name : groupname + "/" + name;
        bout::utils::visit(ADIOSPutVarVisitor(varname, stream), child.value);

        // Write attributes
        if (!BoutComm::rank()) {
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

/// Write options to file
void OptionsADIOS::write(const Options& options, const std::string& time_dim) {
  Timer timer("io");

  // ADIOSStream is just a BOUT++ object, it does not create anything inside ADIOS
  ADIOSStream& stream = ADIOSStream::ADIOSGetStream(filename);

  // Need to have an adios2::IO object first, which can only be created once.
  if (!stream.io) {
    ADIOSPtr adiosp = GetADIOSPtr();
    std::string ioname = "write_" + filename;
    try {
      stream.io = adiosp->AtIO(ioname);
    } catch (const std::invalid_argument& e) {
      stream.io = adiosp->DeclareIO(ioname);
      stream.io.SetEngine("BP4");
    }
  }

  if (file_mode == FileMode::append) {
    // Open file once and keep it open, close in stream desctructor
    if (!stream.engine) {
      std::cout << "OptionsADIOS::write: open for append " << filename
                << " timedim = " << time_dim << std::endl;
      stream.engine = stream.io.Open(filename, adios2::Mode::Append);
      if (!stream.engine) {
        throw BoutException("Could not open ADIOS file '{:s}' for writing", filename);
      }
    } else {
      std::cout << "OptionsADIOS::write: revisiting append  " << filename
                << " timedim = " << time_dim << std::endl;
    }
  } else {
    if (!stream.engine) {
      std::cout << "OptionsADIOS::write: create " << filename << " timedim = " << time_dim
                << std::endl;
      stream.engine = stream.io.Open(filename, adios2::Mode::Write);
    } else {
      std::cout << "OptionsADIOS::write: revisiting write  " << filename
                << " timedim = " << time_dim << std::endl;
    }
  }

  std::cout << "          BetweenStepPairs = " << stream.isInStep << std::endl;

  if (!stream.isInStep) {
    stream.engine.BeginStep();
    stream.isInStep = true;
    stream.adiosStep = stream.engine.CurrentStep();
    std::cout << "OptionsADIOS::write: BEGIN adios step = " << stream.adiosStep
              << "  BetweenStepPairs = " << stream.isInStep << std::endl;
  } else {
    std::cout << "OptionsADIOS::write: continue writing adios step = " << stream.adiosStep
              << std::endl;
  }
  writeGroup(options, stream, "", time_dim);
}

} // namespace bout

#endif // BOUT_HAS_ADIOS
