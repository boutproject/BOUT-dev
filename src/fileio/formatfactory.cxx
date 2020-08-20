#include "bout/build_config.hxx"

#include <globals.hxx>

#include "formatfactory.hxx"

#include "impls/emptyformat.hxx"

#include "impls/netcdf4/ncxx4.hxx"
#include "impls/netcdf/nc_format.hxx"
#include "impls/hdf5/h5_format.hxx"
#include "impls/pnetcdf/pnetcdf.hxx"

#include <boutexception.hxx>
#include <output.hxx>
#include <utils.hxx>
#include <cstring>

FormatFactory *FormatFactory::instance = nullptr;

FormatFactory* FormatFactory::getInstance() {
  if (instance == nullptr) {
    // Create the singleton object
    instance = new FormatFactory();
  }
  return instance;
}

// Work out which data format to use for given filename
std::unique_ptr<DataFormat> FormatFactory::createDataFormat(const char *filename,
                                                            bool parallel,
                                                            Mesh* mesh_in) {
  if ((filename == nullptr) || (strcasecmp(filename, "default") == 0)) {
    // Return default file format
    

    if (parallel) {
#ifdef PNCDF
      return bout::utils::make_unique<PncFormat>(mesh_in);
#else
    }

#if BOUT_HAS_NETCDF
    return bout::utils::make_unique<Ncxx4>(mesh_in);
#else

#if BOUT_HAS_LEGACY_NETCDF
    return bout::utils::make_unique<NcFormat>(mesh_in);
#else

#if BOUT_HAS_HDF5
    return bout::utils::make_unique<H5Format>(mesh_in);
#else

#error No file format available; aborting.

#endif // BOUT_HAS_HDF5
#endif // BOUT_HAS_LEGACY_NETCDF
#endif // BOUT_HAS_NETCDF
#endif // PNCDF
    throw BoutException("Parallel I/O disabled, no serial library found");
  }

  // Extract the file extension

  int len = strlen(filename);

  int ind = len-1;  
  while((ind != -1) && (filename[ind] != '.')) {
    ind--;
  }
  
  const char *s = filename + ind+1;

  // Match strings
  
#ifdef PNCDF
  if(parallel) {
    const char *pncdf_match[] = {"cdl", "nc", "ncdf"};
    if(matchString(s, 3, pncdf_match) != -1) {
      output.write("\tUsing Parallel NetCDF format for file '{:s}'\n", filename);
      return bout::utils::make_unique<PncFormat>();
    }
  }
#endif

#if BOUT_HAS_NETCDF
  const char *ncdf_match[] = {"cdl", "nc", "ncdf"};
  if(matchString(s, 3, ncdf_match) != -1) {
    output.write("\tUsing NetCDF4 format for file '{:s}'\n", filename);
    return bout::utils::make_unique<Ncxx4>();
  }
#endif

#if BOUT_HAS_LEGACY_NETCDF
  const char *ncdf_match[] = {"cdl", "nc", "ncdf"};
  if(matchString(s, 3, ncdf_match) != -1) {
    output.write("\tUsing NetCDF format for file '{:s}'\n", filename);
    return bout::utils::make_unique<NcFormat>();
  }
#endif

#if BOUT_HAS_HDF5
  const char *hdf5_match[] = {"h5","hdf","hdf5"};
  if(matchString(s, 3, hdf5_match) != -1) {
    output.write("\tUsing HDF5 format for file '{:s}'\n", filename);
#ifdef PHDF5
    return bout::utils::make_unique<H5Format>(parallel);
#else
    return bout::utils::make_unique<H5Format>();
#endif
  }
#endif

  throw BoutException("\tFile extension not recognised for '{:s}'\n", filename);
  return nullptr;
}

////////////////////// Private functions /////////////////////////////

int FormatFactory::matchString(const char *str, int n, const char **match) {
  for(int i=0;i<n;i++) {
    if(strcasecmp(str, match[i]) == 0) {
      return i;
    }
  }
  return -1;
}

////////////////////// Depreciated function ///////////////////////////

std::unique_ptr<DataFormat> data_format(const char *filename) {
  return FormatFactory::getInstance()->createDataFormat(filename);
}
