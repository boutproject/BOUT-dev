#include "bout/build_config.hxx"

#if BOUT_HAS_ADIOS

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

struct ADIOSStream {
  adios2::IO io;
  adios2::Engine engine;
  adios2::Variable<int> vStep;
  adios2::Variable<double> vTime;
  int adiosStep = 0;
};

/// Name of the attribute used to track individual variable's time indices
constexpr auto current_time_index_name = "current_time_index";

template <class T>
bool readVariable(adios2::Engine& reader, adios2::IO& io, const std::string& name,
                  const std::string& type, Options& result) {
  bool ret = false;
  std::vector<T> data;
  adios2::Variable<T> variable = io.InquireVariable<T>(name);

  if (variable.ShapeID() == adios2::ShapeID::GlobalValue) {
    T value;
    io.Get<T>(variable, &value, adios2::Mode::Sync);
    result[name] = value;
    return true;
  }

  if (variable.ShapeID() == adios2::ShapeID::LocalArray) {
    if (!rank)
      std::cout << "    LocalArray not supported" << std::endl;
    return ret;
  }

  auto dims = variable.Shape();
  auto ndim = shape.size();

  switch (ndims) {
  case 1: {
    Array<T> value(static_cast<int>(dims[0]));
    io.Get<T>(variable, value.data(), adios2::Mode::Sync);
    result[name] = value;
  }
  case 2: {
    Matrix<T> value(static_cast<int>(dims[0].getSize()),
                    static_cast<int>(dims[1].getSize()));
    io.Get<T>(variable, value.data(), adios2::Mode::Sync);
    result[name] = value;
  }
  case 3: {
    Tensor<T> value(static_cast<int>(dims[0].getSize()),
                    static_cast<int>(dims[1].getSize()),
                    static_cast<int>(dims[2].getSize()));
    io.Get<T>(variable, value.data(), adios2::Mode::Sync);
    result[name] = value;
  }
  }

  if (!rank)
    std::cout << "    array has " << variable.Steps() << " steps" << std::endl;

  /* Need to read the data here */
  result[name] = data.data();
  return ret;
}

bool readVariable(adios2::Engine& reader, adios2::IO& io, const std::string& name,
                  const std::string& type, Options& result) {
  bool ret;
#define declare_template_instantiation(T)                  \
  if (type == adios2::GetType<T>()) {                      \
    ret = readVariable<T>(reader, io, name, type, result); \
  }
  ADIOS2_FOREACH_ATTRIBUTE_PRIMITIVE_STDTYPE_1ARG(declare_template_instantiation)
#undef declare_template_instantiation
  return ret;
}

bool readAttribute(adios2::Engine& reader, adios2::IO& io, const std::string& name,
                   const std::string& type, Options& result) {
  bool ret;
#define declare_template_instantiation(T)                  \
  if (type == adios2::GetType<T>()) {                      \
    adios2::Attribute<T> a = io.InquireAtrribute<T>(name); \
    result[name] = a.Data().data();                        \
  }
  ADIOS2_FOREACH_ATTRIBUTE_PRIMITIVE_STDTYPE_1ARG(declare_template_instantiation)
#undef declare_template_instantiation
  return ret;
} // namespace bout

void readGroup(const std::string& filename, ADIOSStruct& group, Options& result) {}

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
      const auto& att = attpair.second;     // NcVarAtt object
      auto it = attpair.second.find("Type");
      const std::string& att_type = it->second;
      readAttribute(reader, io, att_name, att_type, result);
    }
  }

  return result;
}

void writeGroup(const Options& options, NcGroup group,
                const std::string& time_dimension) {

  for (const auto& childpair : options.getChildren()) {
    const auto& name = childpair.first;
    const auto& child = childpair.second;

    if (child.isValue()) {
      try {
        auto nctype = bout::utils::visit(NcTypeVisitor(), child.value);

        if (nctype.isNull()) {
          continue; // Skip this value
        }

        // Get spatial dimensions
        auto spatial_dims = bout::utils::visit(NcDimVisitor(group), child.value);

        // Vector of all dimensions, including time
        std::vector<NcDim> dims{spatial_dims};

        // Get the time dimension
        NcDim time_dim; ///< Time dimension (Null -> none)
        auto time_it = child.attributes.find("time_dimension");
        if (time_it != child.attributes.end()) {
          // Has a time dimension

          const auto& time_name = bout::utils::get<std::string>(time_it->second);

          // Only write time-varying values that match current time
          // dimension being written
          if (time_name != time_dimension) {
            continue;
          }

          time_dim = group.getDim(time_name, NcGroup::ParentsAndCurrent);
          if (time_dim.isNull()) {
            time_dim = group.addDim(time_name);
          }

          // prepend to vector of dimensions
          dims.insert(dims.begin(), time_dim);
        }

        // Check if the variable exists
        auto var = group.getVar(name);
        if (var.isNull()) {
          // Variable doesn't exist yet
          // Create variable
          // Temporary NcType as a workaround for bug in ADIOS 4.4.0 and
          // ADIOS-CXX4 4.2.0
          var = group.addVar(name, NcType{group, nctype.getId()}, dims);
          if (!time_dim.isNull()) {
            // Time evolving variable, so we'll need to keep track of its time index
            var.putAtt(current_time_index_name, ncInt, 0);
          }
        } else {
          // Variable does exist

          // Check types are the same
          if (var.getType() != nctype) {
            throw BoutException(
                "Changed type of variable '{:s}'. Was '{:s}', now writing '{:s}'", name,
                var.getType().getName(), nctype.getName());
          }

          // Check that the dimensions are correct
          auto var_dims = var.getDims();

          // Same number of dimensions?
          if (var_dims.size() != dims.size()) {
            throw BoutException(
                "Changed dimensions for variable '{:s}'\nIn file has {:d} "
                "dimensions, now writing {:d}\n",
                name, var_dims.size(), dims.size());
          }
          // Dimensions compatible?
          for (std::vector<netCDF::NcDim>::size_type i = 0; i < dims.size(); ++i) {
            if (var_dims[i] == dims[i]) {
              continue; // The same dimension -> ok
            }
            if (var_dims[i].isUnlimited() != dims[i].isUnlimited()) {
              throw BoutException("Unlimited dimension changed for variable '{:s}'",
                                  name);
            }
            if (var_dims[i].getSize() != dims[i].getSize()) {
              throw BoutException("Dimension size changed for variable '{:s}'", name);
            }
          }
          // All ok. Set dimensions to the variable's NcDims
          dims = var_dims;

          if (!time_dim.isNull()) {
            // A time dimension
            time_dim = dims[0];
          }
        }

        // Write the variable

        if (time_dim.isNull()) {
          // No time index

          // Put the data into the variable
          bout::utils::visit(NcPutVarVisitor(var), child.value);

        } else {
          // Has a time index, so need the record index

          const int current_time_index = getCurrentTimeIndex(var);

          std::vector<size_t> start_index; ///< Starting index where data will be inserted
          std::vector<size_t> count_index; ///< Size of each dimension

          // Dimensions, including time
          for (const auto& dim : dims) {
            start_index.push_back(0);
            count_index.push_back(dim.getSize());
          }
          // Time dimension
          start_index[0] = current_time_index;
          count_index[0] = 1; // Writing one record

          // Put the data into the variable
          bout::utils::visit(NcPutVarCountVisitor(var, start_index, count_index),
                             child.value);

          // We've just written a new time slice, so we need to update
          // the attribute to track it
          var.putAtt(current_time_index_name, ncInt, current_time_index + 1);
        }

        // Write attributes
        for (const auto& attribute : child.attributes) {
          const std::string& att_name = attribute.first;
          const auto& att = attribute.second;

          bout::utils::visit(NcPutAttVisitor(var, att_name), att);
        }

      } catch (const std::exception& e) {
        throw BoutException("Error while writing value '{:s}' : {:s}", name, e.what());
      }
    }

    if (child.isSection()) {
      // Check if the group exists
      TRACE("Writing group '{:s}'", name);

      auto subgroup = group.getGroup(name);
      if (subgroup.isNull()) {
        // Doesn't exist yet, so create it
        subgroup = group.addGroup(name);
      }

      writeGroup(child, subgroup, time_dimension);
    }
  }
}

/// Helper struct for returning errors from verifyTimesteps(NcGroup)
struct TimeDimensionError {
  std::string variable_name;
  std::string time_name;
  std::size_t expected_size;
  std::size_t current_size;
};

std::vector<TimeDimensionError> verifyTimesteps(const NcGroup& group) {

  // Map of dimension -> size
  std::map<NcDim, std::size_t> seen_time_dimensions;
  // Variables with mismatched dimension sizes. Note that this might
  // be a little odd: if the first variable we come across has the
  // "wrong" dimension size, we will actually list all the others as
  // being wrong!
  std::vector<TimeDimensionError> errors;

  // For each variable, check its time dimension against what we've
  // seen already. Note that this assumes a single time dimension per
  // variable whose is in the attribute "time_dimension"
  for (const auto& varpair : group.getVars()) {
    const auto& var_name = varpair.first; // Name of the variable
    const auto& var = varpair.second;     // The NcVar object

    // Get the name of the time dimension from the attribute
    const auto& attributes = var.getAtts();
    const auto time_it = attributes.find("time_dimension");
    if (time_it == attributes.end()) {
      // No "time_dimension" attribute so presumably not a
      // time-evolving variable
      continue;
    }

    // Use the attribute value to get the actual dimension
    std::string time_name;
    time_it->second.getValues(time_name);
    const auto time_dim = group.getDim(time_name, NcGroup::ParentsAndCurrent);

    // Check if we've already seen this dimension
    auto seen_it = seen_time_dimensions.find(time_dim);
    if (seen_it == seen_time_dimensions.end()) {
      // If we haven't, add it to the map with current time index
      seen_time_dimensions[time_dim] = time_dim.getSize();
      continue;
    }

    // If we have, check if the variable current time index matches time size
    const auto current_time = static_cast<std::size_t>(getCurrentTimeIndex(var));
    if (current_time == time_dim.getSize()) {
      continue;
    }
    // If not, add to list of errors
    errors.push_back({var_name, time_dim.getName(), time_dim.getSize(), current_time});
  }

  // Recurse down into subgroups, shoving any new errors into what
  // we've already got. Don't bother reserving the new size, this
  // shouldn't be big!
  for (const auto& child : group.getGroups()) {
    auto child_errors = verifyTimesteps(child.second);
    errors.insert(errors.end(), child_errors.begin(), child_errors.end());
  }

  return errors;
}

OptionsADIOS::OptionsADIOS() : data_file(nullptr) {}

OptionsADIOS::OptionsADIOS(std::string filename, FileMode mode)
    : filename(std::move(filename)), file_mode(mode), data_file(nullptr) {}

OptionsADIOS::~OptionsADIOS() = default;
OptionsADIOS::OptionsADIOS(OptionsADIOS&&) noexcept = default;
OptionsADIOS& OptionsADIOS::operator=(OptionsADIOS&&) noexcept = default;

void OptionsADIOS::verifyTimesteps() const {
  NcFile dataFile(filename, NcFile::read);
  auto errors = ::verifyTimesteps(dataFile);

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

/// Write options to file
void OptionsADIOS::write(const Options& options, const std::string& time_dim) {
  Timer timer("io");

  adios2::Mode mode =
      (file_mode == FileMode::replace) ? adios2::Mode::Write : adios2::Mode::Append;

  if (not data_file) {
    data_file = std::make_unique<ADIOSStream>();
  }

  // Open file
  ADIOSPtr adiosp = GetADIOSPtr();
  try {
    data_file->io = adiosp->AtIO(filename);
  } catch (const std::invalid_argument& e) {
    std::cerr << e.what() << '\n';
    data_file->io = adiosp->DeclareIO(filename);
  }

  data_file->engine = data_file->io.Open(filename, mode);
  if (!data_file->engine) {
    throw BoutException("Could not open ADIOS file '{:s}' for writing", filename);
  }

  writeGroup(options, *data_file, time_dim);

  data_file->sync();
}

} // namespace bout

#endif // BOUT_HAS_ADIOS
