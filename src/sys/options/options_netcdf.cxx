
#ifdef NCDF4

#include "options_netcdf.hxx"

#include <netcdf>

using namespace netCDF;

Options OptionsNetCDF::read() {
  // Open file
  NcFile dataFile(filename, NcFile::read);

  if (dataFile.isNull()) {
    throw BoutException("Could not open NetCDF file '%s'", filename.c_str()); 
  }
  
  Options result;

  // Iterate over all variables
  for (const auto& varpair : dataFile.getVars()) {
    const auto& var_name = varpair.first; // Name of the variable
    const auto& var = varpair.second; // The NcVar object
    
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
  return result;
}

/// Write options to file
void OptionsNetCDF::write(const Options &options) {
  
}


#endif // NCDF4
