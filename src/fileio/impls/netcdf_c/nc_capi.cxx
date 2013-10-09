
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
  
  
}

bool NcdfCapi::write(int *var, const std::string &name, int lx, int ly, int lz) {
  if(!is_valid())
    return false;
  
}

bool NcdfCapi::write(BoutReal *var, const std::string &name, int lx, int ly, int lz) {
  if(!is_valid())
    return false;
  
  
}

bool NcdfCapi::read_rec(int *var, const std::string &name, int lx, int ly, int lz) {

}

bool NcdfCapi::read_rec(BoutReal *var, const std::string &name, int lx, int ly, int lz) {

}

bool NcdfCapi::write_rec(int *var, const std::string &name, int lx, int ly, int lz) {
  if(!is_valid())
    return false;
  
}

bool NcdfCapi::write_rec(BoutReal *var, const std::string &name, int lx, int ly, int lz) {
  if(!is_valid())
    return false;
  
  
}

#endif // NCDFCAPI
