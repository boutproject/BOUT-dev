
#include <globals.hxx>

#include "formatfactory.hxx"

#include "impls/emptyformat.hxx"

#include "impls/pdb/pdb_format.hxx"
#include "impls/netcdf/nc_format.hxx"

#include <boutexception.hxx>

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
DataFormat* FormatFactory::createDataFormat(const char *filename) {
  if((filename == NULL) || (strcasecmp(filename, "default") == 0)) {
    // Return default file format
    
#ifdef PDBF
    //output.write("\tUsing default format (PDB)\n");
    return new PdbFormat;
#else

#ifdef NCDF
    //output.write("\tUsing default format (NetCDF)\n");
    return new NcFormat;
#else

#error No file format available; aborting.

#endif // NCDF
#endif // PDBF
  }

  // Extract the file extension

  int len = strlen(filename);

  int ind = len-1;  
  while((ind != -1) && (filename[ind] != '.')) {
    ind--;
  }
  
  const char *s = filename + ind+1;

  // Match strings
  
#ifdef PDBF
  const char *pdb_match[] = {"pdb"};
  if(matchString(s, 1, pdb_match) != -1) {
    output.write("\tUsing PDB format for file '%s'\n", filename);
    return new PdbFormat;
  }
#endif

#ifdef NCDF
  const char *ncdf_match[] = {"cdl", "nc", "ncdf"};
  if(matchString(s, 3, ncdf_match) != -1) {
    output.write("\tUsing NetCDF format for file '%s'\n", filename);
    return new NcFormat;
  }
#endif

  output.write("\tFile extension not recognised for '%s'\n", filename);
  // Set to the default
  return createDataFormat();
}

////////////////////// Private functions /////////////////////////////

int FormatFactory::matchString(const char *str, int n, const char **match) {
  for(int i=0;i<n;i++)
    if(strcasecmp(str, match[i]) == 0)
      return i;
  return -1;
}

////////////////////// Depreciated function ///////////////////////////

DataFormat* data_format(const char *filename) {
  return FormatFactory::getInstance()->createDataFormat(filename);
}
