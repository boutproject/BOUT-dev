#include "bout/build_config.hxx"

#if BOUT_HAS_ADIOS

#include "options_adios.hxx"
#include "bout/adios_object.hxx"

#include "bout/bout.hxx"
#include "bout/globals.hxx"
#include "bout/mesh.hxx"
#include "bout/sys/timer.hxx"

#include "adios2.h"
#include <exception>
#include <iostream>
#include <vector>

namespace bout {
/// Name of the attribute used to track individual variable's time indices
constexpr auto current_time_index_name = "current_time_index";

OptionsADIOS::OptionsADIOS(Options& options) : OptionsIO(options) {
  if (options["file"].doc("File name. Defaults to <path>/<prefix>.pb").isSet()) {
    filename = options["file"].as<std::string>();
  } else {
    // Both path and prefix must be set
    filename = fmt::format("{}/{}.bp", options["path"].as<std::string>(),
                           options["prefix"].as<std::string>());
  }

  file_mode = (options["append"].doc("Append to existing file?").withDefault<bool>(false))
                  ? FileMode::append
                  : FileMode::replace;

  singleWriteFile = options["singleWriteFile"].withDefault<bool>(false);
}

template <class T>
bool readVariable(adios2::Engine& reader, adios2::IO& io, const std::string& name,
                  const std::string& type, Options& result) {
  std::vector<T> data;
  adios2::Variable<T> variable = io.InquireVariable<T>(name);

  if (variable.ShapeID() == adios2::ShapeID::GlobalValue) {
    T value;
    reader.Get<T>(variable, &value, adios2::Mode::Sync);
    result[name] = value;
    return true;
  }

  if (variable.ShapeID() == adios2::ShapeID::LocalArray) {
    throw std::invalid_argument(
        "ADIOS reader did not implement reading local arrays like " + type + " " + name
        + " in file " + reader.Name());
  }

  if (type != "double" && type != "float") {
    throw std::invalid_argument(
        "ADIOS reader did not implement reading arrays that are not double/float type. "
        "Found "
        + type + " " + name + " in file " + reader.Name());
  }

  if (type == "double" && sizeof(BoutReal) != sizeof(double)) {
    throw std::invalid_argument(
        "ADIOS does not allow for implicit type conversions. BoutReal type is "
        "float but found "
        + type + " " + name + " in file " + reader.Name());
  }

  if (type == "float" && sizeof(BoutReal) != sizeof(float)) {
    throw std::invalid_argument(
        "ADIOS reader does not allow for implicit type conversions. BoutReal type is "
        "double but found "
        + type + " " + name + " in file " + reader.Name());
  }

  auto dims = variable.Shape();
  auto ndims = dims.size();
  adios2::Variable<BoutReal> variableD = io.InquireVariable<BoutReal>(name);

  switch (ndims) {
  case 1: {
    Array<BoutReal> value(static_cast<int>(dims[0]));
    BoutReal* data = value.begin();
    reader.Get<BoutReal>(variableD, data, adios2::Mode::Sync);
    result[name] = value;
    break;
  }
  case 2: {
    Matrix<BoutReal> value(static_cast<int>(dims[0]), static_cast<int>(dims[1]));
    BoutReal* data = value.begin();
    reader.Get<BoutReal>(variableD, data, adios2::Mode::Sync);
    result[name] = value;
    break;
  }
  case 3: {
    Tensor<BoutReal> value(static_cast<int>(dims[0]), static_cast<int>(dims[1]),
                           static_cast<int>(dims[2]));
    BoutReal* data = value.begin();
    reader.Get<BoutReal>(variableD, data, adios2::Mode::Sync);
    result[name] = value;
    break;
  }
  }

  return true;
}

bool readVariable(adios2::Engine& reader, adios2::IO& io, const std::string& name,
                  const std::string& type, Options& result) {
  bool ret = false;
#define declare_template_instantiation(T)                  \
  if (type == adios2::GetType<T>()) {                      \
    ret = readVariable<T>(reader, io, name, type, result); \
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

void OptionsADIOS::verifyTimesteps() const {
  ADIOSStream& stream = ADIOSStream::ADIOSGetStream(filename);
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

  // Where the actual data starts in data pointer (to exclude ghost cells)
  adios2::Dims memStart = {static_cast<size_t>(mesh->MapLocalX),
                           static_cast<size_t>(mesh->MapLocalY)};

  // The actual size of data pointer in memory (including ghost cells)
  adios2::Dims memCount = {static_cast<size_t>(value.getNx()),
                           static_cast<size_t>(value.getNy())};

  adios2::Variable<BoutReal> var = stream.GetArrayVariable<BoutReal>(varname, shape);
  /* std::cout << "PutVar Field2D rank " << BoutComm::rank() << " var = " << varname
            << " shape = " << shape[0] << "x" << shape[1] << " count = " << count[0]
            << "x" << count[1] << " Nx*Ny = " << value.getNx() << "x" << value.getNy()
            << " memStart = " << memStart[0] << "x" << memStart[1]
            << " memCount = " << memCount[0] << "x" << memCount[1] << std::endl;*/
  var.SetSelection({start, count});
  var.SetMemorySelection({memStart, memCount});
  stream.engine.Put<BoutReal>(var, &value(0, 0));
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

  // Where the actual data starts in data pointer (to exclude ghost cells)
  adios2::Dims memStart = {static_cast<size_t>(mesh->MapLocalX),
                           static_cast<size_t>(mesh->MapLocalY),
                           static_cast<size_t>(mesh->MapLocalZ)};

  // The actual size of data pointer in memory (including ghost cells)
  adios2::Dims memCount = {static_cast<size_t>(value.getNx()),
                           static_cast<size_t>(value.getNy()),
                           static_cast<size_t>(value.getNz())};

  adios2::Variable<BoutReal> var = stream.GetArrayVariable<BoutReal>(varname, shape);
  /*std::cout << "PutVar Field3D rank " << BoutComm::rank() << " var = " << varname
            << " shape = " << shape[0] << "x" << shape[1] << "x" << shape[2]
            << " count = " << count[0] << "x" << count[1] << "x" << count[2]
            << " Nx*Ny = " << value.getNx() << "x" << value.getNy() << "x"
            << value.getNz() << " memStart = " << memStart[0] << "x" << memStart[1] << "x"
            << memStart[2] << " memCount = " << memCount[0] << "x" << memCount[1] << "x"
            << memCount[2] << std::endl;*/
  var.SetSelection({start, count});
  var.SetMemorySelection({memStart, memCount});
  stream.engine.Put<BoutReal>(var, &value(0, 0, 0));
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

  // Where the actual data starts in data pointer (to exclude ghost cells)
  adios2::Dims memStart = {static_cast<size_t>(mesh->MapLocalX),
                           static_cast<size_t>(mesh->MapLocalZ)};

  // The actual size of data pointer in memory (including ghost cells)
  adios2::Dims memCount = {static_cast<size_t>(value.getNx()),
                           static_cast<size_t>(value.getNz())};

  adios2::Variable<BoutReal> var = stream.GetArrayVariable<BoutReal>(varname, shape);
  /* std::cout << "PutVar FieldPerp rank " << BoutComm::rank() << " var = " << varname
            << " shape = " << shape[0] << "x" << shape[1] << " count = " << count[0]
            << "x" << count[1] << " Nx*Ny = " << value.getNx() << "x" << value.getNy()
            << " memStart = " << memStart[0] << "x" << memStart[1]
            << " memCount = " << memCount[0] << "x" << memCount[1] << std::endl; */
  var.SetSelection({start, count});
  var.SetMemorySelection({memStart, memCount});
  stream.engine.Put<BoutReal>(var, &value(0, 0));
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
      stream.io.SetEngine("BP5");
    }
  }

  /* Open file once and keep it open, close in stream desctructor
     or close after writing if singleWriteFile == true
    */
  if (!stream.engine) {
    adios2::Mode amode =
        (file_mode == FileMode::append ? adios2::Mode::Append : adios2::Mode::Write);
    stream.engine = stream.io.Open(filename, amode);
    if (!stream.engine) {
      throw BoutException("Could not open ADIOS file '{:s}' for writing", filename);
    }
  }

  /* Multiple write() calls allowed in a single adios step to output multiple
     Options objects in the same step. verifyTimesteps() will indicate the 
     completion of the step (and adios will publish the step).
  */
  if (!stream.isInStep) {
    stream.engine.BeginStep();
    stream.isInStep = true;
    stream.adiosStep = stream.engine.CurrentStep();
  }

  writeGroup(options, stream, "", time_dim);

  /* In singleWriteFile mode, we complete the step and close the file */
  if (singleWriteFile) {
    stream.engine.EndStep();
    stream.engine.Close();
  }
}

} // namespace bout

#endif // BOUT_HAS_ADIOS
