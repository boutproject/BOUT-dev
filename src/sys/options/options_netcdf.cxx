
#ifdef NCDF4

#include "options_netcdf.hxx"

#include <vector>
#include <netcdf>

using namespace netCDF;

namespace {
void readGroup(const std::string &filename, NcGroup group, Options& result) {
  
  // Iterate over all variables
  for (const auto& varpair : group.getVars()) {
    const auto& var_name = varpair.first; // Name of the variable
    const auto& var = varpair.second;     // The NcVar object

    if (var.getDimCount() == 0) {
      // Scalar variables

      auto var_type = var.getType();

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
      // else ignore
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

namespace {

/// Convert variant into NcType
/// If the type is not recognised then NcType null object is returned
struct NcTypeVisitor {
  template <typename T>
  NcType operator()(const T& UNUSED(t)) {
    return {}; // Null object by defaul
  }
};

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
  void operator()(const T& UNUSED(t)) {}

private:
  NcVar& var;
};

template <>
void NcPutVarVisitor::operator()<int>(const int& value) {
  var.putVar(&value);
}
template <>
void NcPutVarVisitor::operator()<double>(const double& value) {
  var.putVar(&value);
}
template <>
void NcPutVarVisitor::operator()<float>(const float& value) {
  var.putVar(&value);
}
template <>
void NcPutVarVisitor::operator()<std::string>(const std::string& value) {
  const char* cstr = value.c_str();
  var.putVar(&cstr);
}
template <>
void NcPutVarVisitor::operator()<Field2D>(const Field2D& value) {
  // Pointer to data. Assumed to be contiguous array
  var.putVar(&value(0,0));
}
template <>
void NcPutVarVisitor::operator()<Field3D>(const Field3D& value) {
  // Pointer to data. Assumed to be contiguous array
  var.putVar(&value(0,0,0));
}
  
void writeGroup(const Options& options, NcGroup group) {

  for (const auto& childpair : options.getChildren()) {
    const auto& name = childpair.first;
    const auto& child = childpair.second;

    if (child.isValue()) {
      auto nctype = bout::utils::visit(NcTypeVisitor(), child.value);

      if (nctype.isNull()) {
        continue; // Skip this value
      }

      auto dims = bout::utils::visit(NcDimVisitor(group), child.value);
      
      auto var = group.addVar(name, nctype, dims);
      
      // Put the data into the variable
      bout::utils::visit(NcPutVarVisitor(var), child.value);
    }

    if (child.isSection()) {
      writeGroup(child, group.addGroup(name));
    }
  }
}
  
} // namespace

/// Write options to file
void OptionsNetCDF::write(const Options& options) {
  NcFile dataFile(filename, NcFile::replace);

  if (dataFile.isNull()) {
    throw BoutException("Could not open NetCDF file '%s' for writing", filename.c_str());
  }

  writeGroup(options, dataFile);
}

#endif // NCDF4
