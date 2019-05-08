
#ifdef NCDF4

#include "options_netcdf.hxx"

#include <exception>
#include <netcdf>
#include <vector>

using namespace netCDF;

namespace {
void readGroup(const std::string& filename, NcGroup group, Options& result) {

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
    throw BoutException("Could not open NetCDF file '%s'", filename.c_str());
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
    throw BoutException("Error in findDimension('%s'): %s", name.c_str(), e.what());
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
}
  
/// In addition to writing the data, set the "cell_location" attribute
template <>
void NcPutVarVisitor::operator()<Field3D>(const Field3D& value) {
  // Pointer to data. Assumed to be contiguous array
  var.putVar(&value(0, 0, 0));

  // Set cell location attribute
  var.putAtt("cell_location", toString(value.getLocation()));
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

/// Visit a variant type, and put the data into an attributute
struct NcPutAttVisitor {
  NcPutAttVisitor(NcVar& var, std::string name) : var(var), name(name) {}
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
void NcPutAttVisitor::operator()(const float& value) {
  var.putAtt(name, ncFloat, value);
}
template <>
void NcPutAttVisitor::operator()(const std::string& value) {
  var.putAtt(name, value);
}
  
void writeGroup(const Options& options, NcGroup group,
                std::map<int, size_t>& time_index) {

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

          auto time_name = bout::utils::get<std::string>(time_it->second);
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
        } else {
          // Variable does exist

          // Check types are the same
          if (var.getType() != nctype) {
            throw BoutException(
                "Changed type of variable '%s'. Was '%s', now writing '%s'", name.c_str(),
                var.getType().getName().c_str(), nctype.getName().c_str());
          }

          // Check that the dimensions are correct
          auto var_dims = var.getDims();

          // Same number of dimensions?
          if (var_dims.size() != dims.size()) {
            throw BoutException("Changed dimensions for variable '%s'\nIn file has %zu "
                                "dimensions, now writing %zu\n",
                                name.c_str(), var_dims.size(), dims.size());
          }
          // Dimensions compatible?
          for (std::vector<netCDF::NcDim>::size_type i = 0; i < dims.size(); ++i) {
            if (var_dims[i] == dims[i]) {
              continue; // The same dimension -> ok
            }
            if (var_dims[i].isUnlimited() != dims[i].isUnlimited()) {
              throw BoutException("Unlimited dimension changed for variable '%s'",
                                  name.c_str());
            }
            if (var_dims[i].getSize() != dims[i].getSize()) {
              throw BoutException("Dimension size changed for variable '%s'",
                                  name.c_str());
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

          // Get the index from the map storing the current index
          // This is needed because NetCDF doesn't provide a way to get
          // the size of this variable along an unlimited dimension.
          // Instead the dimension is shared between variables.

          auto time_index_it = time_index.find(time_dim.getId());
          if (time_index_it == time_index.end()) {
            // Haven't seen this index before
            time_index[time_dim.getId()] = time_dim.getSize();
          }

          std::vector<size_t> start_index; ///< Starting index where data will be inserted
          std::vector<size_t> count_index; ///< Size of each dimension

          

          // Dimensions, including time
          for (const auto& dim : dims) {
            start_index.push_back(0);
            count_index.push_back(dim.getSize());
          }
          // Time dimension
          start_index[0] = time_index[time_dim.getId()];
          count_index[0] = 1; // Writing one record

          // Put the data into the variable
          bout::utils::visit(NcPutVarCountVisitor(var, start_index, count_index),
                             child.value);
        }

        // Write attributes
        for (const auto& it: child.attributes) {
          const std::string& att_name = it.first;
          const auto& att = it.second;
          
          bout::utils::visit(NcPutAttVisitor(var, att_name), att);
        }
        
      } catch (const std::exception &e) {
        throw BoutException("Error while writing value '%s' : %s", name.c_str(), e.what());
      }
    }

    if (child.isSection()) {
      // Check if the group exists
      TRACE("Writing group '%s'", name.c_str());

      auto subgroup = group.getGroup(name);
      if (subgroup.isNull()) {
        // Doesn't exist yet, so create it
        subgroup = group.addGroup(name);
      }
      
      writeGroup(child, subgroup, time_index);
    }
  }
}

} // namespace

namespace bout {
namespace experimental {

/// Write options to file
void OptionsNetCDF::write(const Options& options) {
  // Check the file mode to use
  auto ncmode = NcFile::replace;
  if (file_mode == FileMode::append) {
    ncmode = NcFile::write;
  }

  NcFile dataFile(filename, ncmode);

  if (dataFile.isNull()) {
    throw BoutException("Could not open NetCDF file '%s' for writing", filename.c_str());
  }

  writeGroup(options, dataFile, time_index);
}

} // experimental
} // bout

#endif // NCDF4
