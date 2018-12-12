
#ifdef NCDF4

#include "options_netcdf.hxx"

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

void writeGroup(const Options& options, NcGroup group) {

  for (const auto& childpair : options.getChildren()) {
    const auto& name = childpair.first;
    const auto& child = childpair.second;

    if (child.isValue()) {
      auto nctype = bout::utils::visit(NcTypeVisitor(), child.value);

      if (nctype.isNull()) {
        continue; // Skip this value
      }

      auto var = group.addVar(name, nctype);

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
