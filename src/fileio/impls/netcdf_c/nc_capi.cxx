
#ifdef NCDFCAPI

#include "nc_capi.hxx"

#include <netcdf.h>

#include <boutexception.hxx>
#include <utils.hxx>
#include <output.hxx>
#include <globals.hxx>

NcdfCapi::NcdfCapi() : datafile(0) {
  
}

NcdfCapi::NcdfCapi(const std::string &name) : datafile(0) {
  if(!openr(name))
    throw BoutException("Couldn't open file %s", name.c_str());
}

NcdfCapi::~NcdfCapi() {
  close();
}

bool NcdfCapi::openr(const std::string &name) {
  close();
  
  // Open existing file for reading
  if( nc_open(name.c_str(), NC_NOWRITE, &datafile) != NC_NOERR ) {
    // error
  }

  // Get dimensions
  if(nc_inq_dimid(datafile, "x", &xDim) != NC_NOERR) {
    output.write("WARNING: NetCDF file should have an 'x' dimension\n");
  }
  if(nc_inq_dimid(datafile, "y", &yDim) != NC_NOERR) {
    output.write("WARNING: NetCDF file should have an 'y' dimension\n");
  }
  if(nc_inq_dimid(datafile, "z", &zDim) != NC_NOERR) {
    // Z dimension optional, and could be any size
  }
  if(nc_inq_dimid(datafile, "t", &tDim) != NC_NOERR) {
    // t dimension optional
  }
  
  recDimList[0] = tDim;
  recDimList[1] = xDim;
  recDimList[2] = yDim;
  recDimList[3] = zDim;
  
  dimList = recDimList+1;
  
  return true;
}

bool NcdfCapi::openw(const std::string &name, bool append) {
  close();

  if(!append) {
    // Creating new file
    if( nc_create(name.c_str(), NC_CLASSIC_MODEL, &datafile) != NC_NOERR ) {
      // error
      output.write("ERROR: NetCDF C API\n");
      return false;
    }
    
    // Add dimensions
    if( nc_def_dim(datafile, "x", mesh->ngx, &xDim) != NC_NOERR) {
      // error
    }
    if( nc_def_dim(datafile, "y", mesh->ngy, &yDim) != NC_NOERR) {
      // error
    }
    if( nc_def_dim(datafile, "z", mesh->ngz, &zDim) != NC_NOERR) {
      // error
    }
    if( nc_def_dim(datafile, "t", NC_UNLIMITED, &tDim) != NC_NOERR) {
      // error
    }
    
    // Finish define mode
    if( nc_enddef(datafile) != NC_NOERR ) {
      // error
    }
  }else {
    // Appending to an existing file
    if( nc_open(name.c_str(), NC_WRITE, &datafile) != NC_NOERR ) {
      // error
    }
    
    // Get dimensions from the file
    if(nc_inq_dimid(datafile, "x", &xDim) != NC_NOERR) {
      output.write("ERROR: NetCDF file should have an 'x' dimension\n");
      close();
      return false;
    }
    if(nc_inq_dimid(datafile, "y", &yDim) != NC_NOERR) {
      output.write("ERROR: NetCDF file should have an 'y' dimension\n");
      close();
      return false;
    }
    if(nc_inq_dimid(datafile, "z", &zDim) != NC_NOERR) {
      output.write("ERROR: NetCDF file should have an 'z' dimension\n");
      close();
      return false;
    }
    if(nc_inq_dimid(datafile, "t", &tDim) != NC_NOERR) {
      output.write("ERROR: NetCDF file should have an 't' dimension\n");
      close();
      return false;
    }
    
    // Test they're the right size (and t is unlimited)
    size_t len;
    if( nc_inq_dimlen(datafile, xDim, &len) != NC_NOERR) {
      output.write("ERROR: NetCDF");
      close();
      return false;
    }
    if(len != mesh->ngx) {
      output.write("ERROR: NetCDF x dimension has length %d. Expected %d\n", len, mesh->ngx);
      close();
      return false;
    }
    
    if( nc_inq_dimlen(datafile, yDim, &len) != NC_NOERR) {
      // error
    }
    if(len != mesh->ngy) {
      output.write("ERROR: NetCDF y dimension has length %d. Expected %d\n", len, mesh->ngy);
      close();
      return false;
    }
    
    if( nc_inq_dimlen(datafile, zDim, &len) != NC_NOERR) {
      // error
    }
    if(len != mesh->ngz) {
      output.write("ERROR: NetCDF z dimension has length %d. Expected %d\n", len, mesh->ngz);
      close();
      return false;
    }
    
  }
  
  recDimList[0] = tDim;
  recDimList[1] = xDim;
  recDimList[2] = yDim;
  recDimList[3] = zDim;
  
  dimList = recDimList+1;
  
  return true;
}

bool NcdfCapi::is_valid() {
  return datafile != 0;
}

void NcdfCapi::close() {
  if(datafile == 0)
    return;
  nc_close(datafile);
  datafile = 0;
}
  
void NcdfCapi::flush() {
  if(datafile == 0)
    return;
  nc_sync(datafile);
}

const vector<int> NcdfCapi::getSize(const string &var) {
  vector<int> size;
  
  if(!is_valid())
    return size;
  
  // Find the variable ID
  int varid;
  if(nc_inq_varid(datafile, var.c_str(), &varid) != NC_NOERR) {
    // error: no variable
    return size;
  }

  // Get number of dimensions
  int ndims;
  if(nc_inq_varndims(datafile, varid, &ndims) != NC_NOERR) {
    // error
    throw BoutException("NetCDF error in nc_inq_varndims");
  }

  if(ndims == 0) {
    // Just a scalar
    size.push_back(1);
    return size;
  }

  size.resize(ndims);
  
  // Get the dimension IDs
  if(nc_inq_vardimid(datafile, varid, &size[0]) != NC_NOERR) {
    // error
    throw BoutException("NetCDF error in nc_inq_vardimid");
  }
  
  // Convert dimension ID to a length, overwriting ID
  for(int i=0;i<ndims;i++) {
    size_t len;
    if(nc_inq_dimlen(datafile, size[i], &len) != NC_NOERR) {
      // error
      throw BoutException("NetCDF error in nc_inq_dimlen");
    }
    size[i] = len;
  }
  
  return size;
  
}

bool NcdfCapi::setGlobalOrigin(int x, int y, int z) {
  x0 = x;
  y0 = y;
  z0 = z;
  
  return true;
}

bool NcdfCapi::setRecord(int t) {
  t0 = t;

  return true;
}

bool NcdfCapi::read(int *var, const std::string &name, int lx, int ly, int lz) {
  if(!is_valid())
    return false;
  
  if((lx < 0) || (ly < 0) || (lz < 0))
    return false;
  
  int varid;
  if(nc_inq_varid(datafile, name.c_str(), &varid) != NC_NOERR) {
    return false; // Variable not found
  }
  
  size_t start[3], counts[3];
  start[0] = x0; start[1] = y0; start[2] = z0;
  counts[0] = lx; counts[1] = ly; counts[2] = lz;
  
  if(nc_get_vara_int(datafile, varid, start, counts, var) != NC_NOERR) {
    return false;
  }

  return true;
}

bool NcdfCapi::read(BoutReal *var, const std::string &name, int lx, int ly, int lz) {
  if(!is_valid())
    return false;
  
  if((lx < 0) || (ly < 0) || (lz < 0))
    return false;
  
  int varid;
  if(nc_inq_varid(datafile, name.c_str(), &varid) != NC_NOERR) {
    return false; // Variable not found
  }
  
  size_t start[3], counts[3];
  start[0] = x0; start[1] = y0; start[2] = z0;
  counts[0] = lx; counts[1] = ly; counts[2] = lz;
  
  if(nc_get_vara_double(datafile, varid, start, counts, var) != NC_NOERR) {
    return false;
  }
  
  return true;
}

bool NcdfCapi::write(int *var, const std::string &name, int lx, int ly, int lz) {
  if(!is_valid())
    return false;
  
  if((lx < 0) || (ly < 0) || (lz < 0))
    return false;
  
  int varid;
  if(nc_inq_varid(datafile, name.c_str(), &varid) != NC_NOERR) {
    // Variable not in file, so add it
    
    if(nc_redef(datafile) != NC_NOERR) {
      // error
      return false;
    }

    int nd = 0; // Number of dimensions
    if(lx != 0) nd = 1;
    if(ly != 0) nd = 2;
    if(lz != 0) nd = 3;
     
    if(nc_def_var(datafile, name.c_str(), NC_INT, nd, dimList, &varid) != NC_NOERR) {
      output.write("ERROR: NetCDF could not add int '%s'\n", name.c_str());
      return false;
    }
    
    if( nc_enddef(datafile) != NC_NOERR ) {
      // error
    }
  }
  
  size_t start[3], counts[3];
  start[0] = x0; start[1] = y0; start[2] = z0;
  counts[0] = lx; counts[1] = ly; counts[2] = lz;
  
  if(nc_put_vara_int(datafile, varid, start, counts, var) != NC_NOERR) {
    return false;
  }
  return true;
}

bool NcdfCapi::write(BoutReal *var, const std::string &name, int lx, int ly, int lz) {
  if(!is_valid())
    return false;
  
  int varid;
  if(nc_inq_varid(datafile, name.c_str(), &varid) != NC_NOERR) {
    // Variable not in file, so add it
    
    if(nc_redef(datafile) != NC_NOERR) {
      // error
      return false;
    }

    int nd = 0; // Number of dimensions
    if(lx != 0) nd = 1;
    if(ly != 0) nd = 2;
    if(lz != 0) nd = 3;
    
    nc_type type = NC_DOUBLE;
    if(lowPrecision)
      type = NC_FLOAT;
    
    if(nc_def_var(datafile, name.c_str(), type, nd, dimList, &varid) != NC_NOERR) {
      output.write("ERROR: NetCDF could not add BoutReal '%s'\n", name.c_str());
      return false;
    }
    
    if( nc_enddef(datafile) != NC_NOERR ) {
      // error
    }
  }
  
  size_t start[3], counts[3];
  start[0] = x0; start[1] = y0; start[2] = z0;
  counts[0] = lx; counts[1] = ly; counts[2] = lz;
  
  if(nc_put_vara_double(datafile, varid, start, counts, var) != NC_NOERR) {
    return false;
  }
  return true;
}

bool NcdfCapi::read_rec(int *var, const std::string &name, int lx, int ly, int lz) {
  if(!is_valid())
    return false;

  if((lx < 0) || (ly < 0) || (lz < 0))
    return false;
  
  int varid;
  if(nc_inq_varid(datafile, name.c_str(), &varid) != NC_NOERR) {
    return false;
  }
  
  size_t start[4], counts[4];
  start[0] = t0; start[1] = x0; start[2] = y0; start[3] = z0;
  counts[0] = 1; counts[1] = lx; counts[2] = ly; counts[3] = lz;
  
  if(nc_get_vara_int(datafile, varid, start, counts, var) != NC_NOERR) {
    return false;
  }

  return true;
}

bool NcdfCapi::read_rec(BoutReal *var, const std::string &name, int lx, int ly, int lz) {
  
  if(!is_valid())
    return false;

  if((lx < 0) || (ly < 0) || (lz < 0))
    return false;
  
  int varid;
  if(nc_inq_varid(datafile, name.c_str(), &varid) != NC_NOERR) {
    return false;
  }
  
  size_t start[4], counts[4];
  start[0] = t0; start[1] = x0; start[2] = y0; start[3] = z0;
  counts[0] = 1; counts[1] = lx; counts[2] = ly; counts[3] = lz;
  
  if(nc_get_vara_double(datafile, varid, start, counts, var) != NC_NOERR) {
    return false;
  }

  return true;
}

bool NcdfCapi::write_rec(int *var, const std::string &name, int lx, int ly, int lz) {
  if(!is_valid())
    return false;
  
  if((lx < 0) || (ly < 0) || (lz < 0))
    return false;
  
  int varid;
  if(nc_inq_varid(datafile, name.c_str(), &varid) != NC_NOERR) {
    // Variable not in file, so add it
    
    if(nc_redef(datafile) != NC_NOERR) {
      // error
      return false;
    }

    int nd = 1; // Number of dimensions
    if(lx != 0) nd = 2;
    if(ly != 0) nd = 3;
    if(lz != 0) nd = 4;
     
    if(nc_def_var(datafile, name.c_str(), NC_INT, nd, recDimList, &varid) != NC_NOERR) {
      output.write("ERROR: NetCDF could not add int '%s'\n", name.c_str());
      return false;
    }
    
    rec_nr[name] = default_rec; // Starting record
    
    if( nc_enddef(datafile) != NC_NOERR ) {
      // error
    }
  }else {
    // Get record number
    if(rec_nr.find(name) == rec_nr.end()) {
      // Add to map
      rec_nr[name] = default_rec;
    }
  }
  size_t start[4], counts[4];
  start[0] = rec_nr[name]; start[1] = x0; start[2] = y0; start[3] = z0;
  counts[0] = 1; counts[1] = lx; counts[2] = ly; counts[3] = lz;
  
  if(nc_put_vara_int(datafile, varid, start, counts, var) != NC_NOERR) {
    return false;
  }

  // Increment record number
  rec_nr[name] += 1;

  return true;
}

bool NcdfCapi::write_rec(BoutReal *var, const std::string &name, int lx, int ly, int lz) {
  if(!is_valid())
    return false;
  
  if((lx < 0) || (ly < 0) || (lz < 0))
    return false;
  
  int varid;
  if(nc_inq_varid(datafile, name.c_str(), &varid) != NC_NOERR) {
    // Variable not in file, so add it
    
    if(nc_redef(datafile) != NC_NOERR) {
      // error
      return false;
    }

    int nd = 1; // Number of dimensions
    if(lx != 0) nd = 2;
    if(ly != 0) nd = 3;
    if(lz != 0) nd = 4;
     
    nc_type vartype = NC_DOUBLE;
    if(lowPrecision)
      vartype = NC_FLOAT;

    if(nc_def_var(datafile, name.c_str(), vartype, nd, recDimList, &varid) != NC_NOERR) {
      output.write("ERROR: NetCDF could not add int '%s'\n", name.c_str());
      return false;
    }
    
    rec_nr[name] = default_rec; // Starting record
    
    if( nc_enddef(datafile) != NC_NOERR ) {
      // error
    }
  }else {
    // Get record number
    if(rec_nr.find(name) == rec_nr.end()) {
      // Add to map
      rec_nr[name] = default_rec;
    }
  }

  if(lowPrecision) {
    // An out of range value can make the conversion
    // corrupt the whole dataset. Make sure everything
    // is in the range of a float
    
    for(int i=0;i<lx*ly*lz;i++) {
      if(var[i] > 1e20)
	var[i] = 1e20;
      if(var[i] < -1e20)
	var[i] = -1e20;
    }
  }
  
  for(int i=0;i<lx*ly*lz;i++) {
    if(!finite(var[i]))
      var[i] = 0.0;
  }

  size_t start[4], counts[4];
  start[0] = rec_nr[name]; start[1] = x0; start[2] = y0; start[3] = z0;
  counts[0] = 1; counts[1] = lx; counts[2] = ly; counts[3] = lz;

  if(nc_put_vara_double(datafile, varid, start, counts, var) != NC_NOERR) {
    return false;
  }

  // Increment record number
  rec_nr[name] += 1;

  return true;
}

#endif // NCDFCAPI
