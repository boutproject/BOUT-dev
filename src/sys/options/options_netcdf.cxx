#include "bout/build_config.hxx"

#if BOUT_HAS_NETCDF && !BOUT_HAS_LEGACY_NETCDF

#include "options_netcdf.hxx"

#include <exception>
#include <netcdf>
#include <vector>
#include <iostream>

using namespace netCDF;

namespace {
/// Name of the attribute used to track individual variable's time indices
constexpr auto current_time_index_name = "current_time_index";

/// NetCDF doesn't keep track of the current for each variable
/// (although the underlying HDF5 file does!), so we need to do it
/// ourselves. We'll use an attribute in the file to do so, which
/// means we don't need to keep track of it in the code
int getCurrentTimeIndex(const NcVar& var) {
  const auto atts_map = var.getAtts();
  const auto it = atts_map.find(current_time_index_name);
  if (it == atts_map.end()) {
    // Attribute doesn't exist, so let's start at zero. There
    // are various ways this might break, for example, if the
    // variable was added to the file by a different
    // program. But note, that if we use the size of the time
    // dimension here, this will increase every time we add a
    // new variable! So zero is probably the only sensible
    // choice for this
    return 0;
  }
  int current_time_index;
  it->second.getValues(&current_time_index);
  return current_time_index;
}

void readGroup(const std::string& filename, const NcGroup& group, Options& result) {

  // Iterate over all variables
  for (const auto& varpair : group.getVars()) {
    const auto& var_name = varpair.first; // Name of the variable
    const auto& var = varpair.second;     // The NcVar object
    
    auto var_type = var.getType(); // Variable type 
    auto ndims = var.getDimCount(); // Number of dimensions
    auto dims = var.getDims(); // Vector of dimensions
    
    switch (ndims) {
    case 0: {
      // Scalar variables
      
      if (var_type == ncDouble) {
        double value;
        var.getVar(&value);
        result[var_name] = value;
        result[var_name].attributes["source"] = filename;
      } else if (var_type == ncFloat) {
        float value;
        var.getVar(&value);
        result[var_name] = value;
        result[var_name].attributes["source"] = filename;
      } else if (var_type == ncInt) {
        int value;
        var.getVar(&value);
        result[var_name] = value;
        result[var_name].attributes["source"] = filename;
      } else if (var_type == ncString) {
        char* value;
        var.getVar(&value);
        result[var_name] = std::string(value);
        result[var_name].attributes["source"] = filename;
      }
      // Note: NetCDF does not support boolean atoms
      // else ignore
      break;
    }
    case 1: {
      if (var_type == ncDouble) {
        Array<double> value(dims[0].getSize());
        var.getVar(value.begin());
        result[var_name] = value;
        result[var_name].attributes["source"] = filename;
      }
      break;
    }
    case 2: {
      if (var_type == ncDouble) {
        Matrix<double> value(dims[0].getSize(), dims[1].getSize());
        var.getVar(value.begin());
        result[var_name] = value;
        result[var_name].attributes["source"] = filename;
      }
      break;
    }
    case 3: {
      if (var_type == ncDouble) {
        Tensor<double> value(dims[0].getSize(), dims[1].getSize(), dims[2].getSize());
        var.getVar(value.begin());
        result[var_name] = value;
        result[var_name].attributes["source"] = filename;
      }
    }
    }

    // Get variable attributes
    for (const auto& attpair : var.getAtts()) {
      const auto &att_name = attpair.first;   // Attribute name
      const auto &att = attpair.second;   // NcVarAtt object
      
      auto att_type = att.getType(); // Type of the attribute

      if (att_type == ncInt) {
        int value;
        att.getValues(&value);
        result[var_name].attributes[att_name] = value;
      } else if (att_type == ncFloat) {
        float value;
        att.getValues(&value);
        result[var_name].attributes[att_name] = value;
      } else if (att_type == ncDouble) {
        double value;
        att.getValues(&value);
        result[var_name].attributes[att_name] = value;
      } else if ((att_type == ncString) or (att_type == ncChar)) {
        std::string value;
        att.getValues(value);
        result[var_name].attributes[att_name] = value;
      }
      // Else ignore
    }
  }

  // Iterate over groups
  for (const auto& grouppair : group.getGroups()) {
    const auto& name = grouppair.first;
    const auto& subgroup = grouppair.second;

    readGroup(filename, subgroup, result[name]);
  }
}
} // namespace


namespace bout {
namespace experimental {

Options OptionsNetCDF::read() {
  // Open file
  NcFile dataFile(filename, NcFile::read);

  if (dataFile.isNull()) {
    throw BoutException("Could not open NetCDF file '{:s}'", filename);
  }

  Options result;
  readGroup(filename, dataFile, result);

  return result;
}

} // experimental
} // bout

namespace {

/// Convert variant into NcType
/// If the type is not recognised then NcType null object is returned
struct NcTypeVisitor {
  template <typename T>
  NcType operator()(const T& UNUSED(t)) {
    return {}; // Null object by default
  }
};

template <>
NcType NcTypeVisitor::operator()<bool>(const bool& UNUSED(t)) {
  return ncInt;
}

template <>
NcType NcTypeVisitor::operator()<int>(const int& UNUSED(t)) {
  return ncInt;
}

template <>
NcType NcTypeVisitor::operator()<double>(const double& UNUSED(t)) {
  return ncDouble;
}

template <>
MAYBE_UNUSED()
NcType NcTypeVisitor::operator()<float>(const float& UNUSED(t)) {
  return ncFloat;
}

template <>
NcType NcTypeVisitor::operator()<std::string>(const std::string& UNUSED(t)) {
  return ncString;
}

template <>
NcType NcTypeVisitor::operator()<Field2D>(const Field2D& UNUSED(t)) {
  return operator()<BoutReal>(0.0);
}

template <>
NcType NcTypeVisitor::operator()<Field3D>(const Field3D& UNUSED(t)) {
  return operator()<BoutReal>(0.0);
}

template <>
NcType NcTypeVisitor::operator()<FieldPerp>(const FieldPerp& UNUSED(t)) {
  return operator()<BoutReal>(0.0);
}

/// Visit a variant type, returning dimensions
struct NcDimVisitor {
  NcDimVisitor(NcGroup& group) : group(group) {}
  template <typename T>
  std::vector<NcDim> operator()(const T& UNUSED(t)) {
    return {};
  }

private:
  NcGroup& group;
};

NcDim findDimension(NcGroup& group, const std::string& name, unsigned int size) {
  // Get the dimension
  try {
    auto dim = group.getDim(name, NcGroup::ParentsAndCurrent);
    if (dim.isNull()) {
      // Dimension doesn't yet exist
      dim = group.addDim(name, size);
    } else {
      // Dimension exists, check it's the right size
      if (dim.getSize() != size) {
        // wrong size. Check this group
        dim = group.getDim(name, NcGroup::Current);
        if (!dim.isNull()) {
          // Already defined in this group
          return {}; // Return null object
        }
        // Define in this group
        dim = group.addDim(name, size);
      }
    }
    return dim;
  } catch (const std::exception& e) {
    throw BoutException("Error in findDimension('{:s}'): {:s}", name, e.what());
  }
}

template <>
std::vector<NcDim> NcDimVisitor::operator()<Field2D>(const Field2D& value) {
  auto xdim = findDimension(group, "x", value.getNx());
  ASSERT0(!xdim.isNull());

  auto ydim = findDimension(group, "y", value.getNy());
  ASSERT0(!ydim.isNull());

  return {xdim, ydim};
}

template <>
std::vector<NcDim> NcDimVisitor::operator()<Field3D>(const Field3D& value) {
  auto xdim = findDimension(group, "x", value.getNx());
  ASSERT0(!xdim.isNull());

  auto ydim = findDimension(group, "y", value.getNy());
  ASSERT0(!ydim.isNull());

  auto zdim = findDimension(group, "z", value.getNz());
  ASSERT0(!zdim.isNull());

  return {xdim, ydim, zdim};
}

template <>
std::vector<NcDim> NcDimVisitor::operator()<FieldPerp>(const FieldPerp& value) {
  auto xdim = findDimension(group, "x", value.getNx());
  ASSERT0(!xdim.isNull());

  auto zdim = findDimension(group, "z", value.getNz());
  ASSERT0(!zdim.isNull());

  return {xdim, zdim};
}

/// Visit a variant type, and put the data into a NcVar
struct NcPutVarVisitor {
  NcPutVarVisitor(NcVar& var) : var(var) {}
  template <typename T>
  void operator()(const T& value) {
    var.putVar(&value);
  }

private:
  NcVar& var;
};

template <>
void NcPutVarVisitor::operator()<bool>(const bool& value) {
  int int_val = value ? 1 : 0;
  var.putVar(&int_val);
}

template <>
void NcPutVarVisitor::operator()<std::string>(const std::string& value) {
  const char* cstr = value.c_str();
  var.putVar(&cstr);
}
  
/// In addition to writing the data, set the "cell_location" attribute
template <>
void NcPutVarVisitor::operator()<Field2D>(const Field2D& value) {
  // Pointer to data. Assumed to be contiguous array
  var.putVar(&value(0, 0));
  // Set cell location attribute
  var.putAtt("cell_location", toString(value.getLocation()));
  var.putAtt("direction_y", toString(value.getDirectionY()));
  var.putAtt("direction_z", toString(value.getDirectionZ()));
}
  
/// In addition to writing the data, set the "cell_location" attribute
template <>
void NcPutVarVisitor::operator()<Field3D>(const Field3D& value) {
  // Pointer to data. Assumed to be contiguous array
  var.putVar(&value(0, 0, 0));

  // Set cell location attribute
  var.putAtt("cell_location", toString(value.getLocation()));
  var.putAtt("direction_y", toString(value.getDirectionY()));
  var.putAtt("direction_z", toString(value.getDirectionZ()));
}

template <>
void NcPutVarVisitor::operator()<FieldPerp>(const FieldPerp& value) {
  // Pointer to data. Assumed to be contiguous array
  var.putVar(&value(0, 0));

  // Set cell location attribute
  var.putAtt("cell_location", toString(value.getLocation()));
  var.putAtt("direction_y", toString(value.getDirectionY()));
  var.putAtt("direction_z", toString(value.getDirectionZ()));
  var.putAtt("yindex_global", ncInt, value.getGlobalIndex());
}

/// Visit a variant type, and put the data into a NcVar
struct NcPutVarCountVisitor {
  NcPutVarCountVisitor(NcVar& var, const std::vector<size_t>& start,
                       const std::vector<size_t>& count)
      : var(var), start(start), count(count) {}
  template <typename T>
  void operator()(const T& value) {
    var.putVar(start, &value);
  }

private:
  NcVar& var;
  const std::vector<size_t>& start; ///< Starting (corner) index
  const std::vector<size_t>& count; ///< Index count in each dimension
};

template <>
void NcPutVarCountVisitor::operator()<std::string>(const std::string& value) {
  const char* cstr = value.c_str();
  var.putVar(start, &cstr);
}
template <>
void NcPutVarCountVisitor::operator()<Field2D>(const Field2D& value) {
  // Pointer to data. Assumed to be contiguous array
  var.putVar(start, count, &value(0, 0));
}
template <>
void NcPutVarCountVisitor::operator()<Field3D>(const Field3D& value) {
  // Pointer to data. Assumed to be contiguous array
  var.putVar(start, count, &value(0, 0, 0));
}
template <>
void NcPutVarCountVisitor::operator()<FieldPerp>(const FieldPerp& value) {
  // Pointer to data. Assumed to be contiguous array
  var.putVar(start, count, &value(0, 0));
}

/// Visit a variant type, and put the data into an attributute
struct NcPutAttVisitor {
  NcPutAttVisitor(NcVar& var, std::string name) : var(var), name(std::move(name)) {}
  template <typename T>
  void operator()(const T& UNUSED(value)) {
    // Default is to ignore if unhandled
  }
private:
  NcVar& var;
  std::string name;
};

template <>
void NcPutAttVisitor::operator()(const bool& value) {
  int ival = value ? 1 : 0;
  var.putAtt(name, ncInt, ival);
}
template <>
void NcPutAttVisitor::operator()(const int& value) {
  var.putAtt(name, ncInt, value);
}
template <>
void NcPutAttVisitor::operator()(const double& value) {
  var.putAtt(name, ncDouble, value);
}
template <>
MAYBE_UNUSED()
void NcPutAttVisitor::operator()(const float& value) {
  var.putAtt(name, ncFloat, value);
}
template <>
void NcPutAttVisitor::operator()(const std::string& value) {
  var.putAtt(name, value);
}
  
void writeGroup(const Options& options, NcGroup group) {

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
          // Temporary NcType as a workaround for bug in NetCDF 4.4.0 and NetCDF-CXX4 4.2.0
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
        for (const auto& it: child.attributes) {
          const std::string& att_name = it.first;
          const auto& att = it.second;

          bout::utils::visit(NcPutAttVisitor(var, att_name), att);
        }

      } catch (const std::exception &e) {
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

      writeGroup(child, subgroup);
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

} // namespace

namespace bout {
namespace experimental {

void OptionsNetCDF::verifyTimesteps() const {
  NcFile dataFile(filename, NcFile::read);
  auto errors = ::verifyTimesteps(dataFile);

  if (errors.empty()) {
    // No errors
    return;
  }

  std::string error_string = "";
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
void OptionsNetCDF::write(const Options& options) {
  // Check the file mode to use
  auto ncmode = NcFile::replace;
  if (file_mode == FileMode::append) {
    // NetCDF doesn't have a "read-write, create if exists" mode, so
    // we need to check ourselves if the file already exists; if it
    // doesn't, tell NetCDF to create it
    std::ifstream file(filename);
    ncmode = file.good() ? NcFile::FileMode::write : NcFile::FileMode::newFile;
  }

  NcFile dataFile(filename, ncmode);

  if (dataFile.isNull()) {
    throw BoutException("Could not open NetCDF file '{:s}' for writing", filename);
  }

  writeGroup(options, dataFile);
}

} // experimental
} // bout

#endif // BOUT_HAS_NETCDF
