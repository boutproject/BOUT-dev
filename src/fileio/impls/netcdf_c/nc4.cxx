
#include "nc4.hxx"

#include <netcdf.h>

Nc4::Nc4() : datafile(0) {
  
}

Nc4::Nc4(const std::string &name) : datafile(0) {
  if(!openr(name))
    throw BoutException("Couldn't open file %s", name.c_str());
}

Nc4::~Nc4() {
  close();
}

bool Nc4::openr(const std::string &name) {
  close();
  
  // Open existing file for reading
  fname = copy_string(name.c_str());
  if( nc_open(fname, NC_NOWRITE, &datafile) != NC_NOERR ) {
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

bool Nc4::openw(const std::string &name, bool append=false) {
  close();
  
  fname = copy_string(name.c_str());

  if(!append) {
    // Creating new file
    if( nc_create(fname, NC_CLASSIC_MODEL, &datafile) != NC_NOERR ) {
      // error
      
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
    if( nc_open(fname, NC_WRITE, &datafile) != NC_NOERR ) {
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
    int len;
    if( nc_inq_dimlen(datafile, xDim, &len) != NC_NOERR) {
      // error
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

bool Nc4::is_valid() {
  return datafile != 0;
}

void Nc4::close() {
  if(datafile == 0)
    return;
  nc_close(datafile);
  free(fname);
  fname = NULL;
  datafile = 0;
}
  
void Nc4::flush() {
  if(datafile == 0)
    return;
  nc_sync(dataset);
}

const vector<int> Nc4::getSize(const string &var) {
  vector<int> size;
  if(!is_valid())
    return size;
  
  // Find the variable ID
  int varid;
  if(nc_inq_varid(datafile, var.c_str(), &varid) != NC_NOERR) {
    // error: no variable
  }
  
  // Get number of dimensions
  int ndims;
  if(nc_inq_varndims(datafile, varid, &ndims) != NC_NOERR) {
    // error
  }
  size.resize(ndims);
  
  // Get the dimension IDs
  if(nc_inq_vardimid(datafile, varid, &size[0]) != NC_NOERR) {
    // error
  }
  
  // Convert dimension ID to a length, overwriting ID
  for(int i=0;i<ndims;i++) {
    size_t len;
    if(nc_inq_dimlen(datafile, size[i], &len) != NC_NOERR) {
      // error
    }
    size[i] = len;
  }
  
  return size;
}

bool Nc4::setGlobalOrigin(int x, int y, int z) {
  x0 = x;
  y0 = y;
  z0 = z;
  
  return true;
}

bool Nc4::setRecord(int t) {
  t0 = t;

  return true;
}

bool Nc4::read(int *var, const std::string &name, int lx = 1, int ly = 0, int lz = 0) {
  if(!is_valid())
    return false;
  
  if((lx < 0) || (ly < 0) || (lz < 0))
    return false;
  
  int varid;
  if(nc_inq_varid(datafile, name.c_str(), &varid) != NC_NOERR) {
    return false; // Variable not found
  }
  
  
}

bool Nc4::read(BoutReal *var, const std::string &name, int lx = 1, int ly = 0, int lz = 0) {
  if(!is_valid())
    return false;
  
  if((lx < 0) || (ly < 0) || (lz < 0))
    return false;
  
  int varid;
  if(nc_inq_varid(datafile, name.c_str(), &varid) != NC_NOERR) {
    return false; // Variable not found
  }
  
  
}

bool Nc4::write(int *var, const std::string &name, int lx = 0, int ly = 0, int lz = 0) {

}

bool Nc4::write(BoutReal *var, const std::string &name, int lx = 0, int ly = 0, int lz = 0) {

}

