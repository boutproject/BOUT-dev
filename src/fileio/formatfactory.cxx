
#include <globals.hxx>

#include "formatfactory.hxx"

#include "impls/emptyformat.hxx"

#include "impls/netcdf4/ncxx4.hxx"
#include "impls/netcdf/nc_format.hxx"
#include "impls/hdf5/h5_format.hxx"
#include "impls/pnetcdf/pnetcdf.hxx"

#include <boutexception.hxx>
#include <output.hxx>
#include <string.h>

FormatFactory* FormatFactory::instance = NULL;

FormatFactory* FormatFactory::getInstance() {
  if(instance == NULL) {
    // Create the singleton object
    instance = new FormatFactory();
  }
  return instance;
}

// Work out which data format to use for given filename
DataFormat* FormatFactory::createDataFormat(const char *filename, bool parallel) {
  if((filename == NULL) || (strcasecmp(filename, "default") == 0)) {
    // Return default file format
    

#ifdef PNCDF
    if(parallel)
      return new PncFormat;
#else

#ifdef NCDF4
    return new Ncxx4;
#else

#ifdef NCDF
    //output.write("\tUsing default format (NetCDF)\n");
    return new NcFormat;
#else

#ifdef HDF5
    return new H5Format;
#else

#error No file format available; aborting.

#endif // HDF5
#endif // NCDF
#endif // NCDF4
#endif // PNCDF
    throw new BoutException("Parallel I/O disabled, no serial library found");
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
      output.write("\tUsing Parallel NetCDF format for file '%s'\n", filename);
    return new PncFormat;
    }
  }
#endif

#ifdef NCDF4
  const char *ncdf_match[] = {"cdl", "nc", "ncdf"};
  if(matchString(s, 3, ncdf_match) != -1) {
    output.write("\tUsing NetCDF4 format for file '%s'\n", filename);
    return new Ncxx4;
  }
#endif

#ifdef NCDF
  const char *ncdf_match[] = {"cdl", "nc", "ncdf"};
  if(matchString(s, 3, ncdf_match) != -1) {
    output.write("\tUsing NetCDF format for file '%s'\n", filename);
    return new NcFormat;
  }
#endif

#ifdef HDF5
  const char *hdf5_match[] = {"h5","hdf","hdf5"};
  if(matchString(s, 3, hdf5_match) != -1) {
    output.write("\tUsing HDF5 format for file '%s'\n", filename);
#ifdef PHDF5
    return new H5Format(parallel);
#else
    return new H5Format();
#endif
  }
#endif

  throw BoutException("\tFile extension not recognised for '%s'\n", filename);
  return NULL;
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

DataFormat* data_format(const char *filename) {
  return FormatFactory::getInstance()->createDataFormat(filename);
}
