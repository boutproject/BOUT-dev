/**************************************************************************
 * Copyright 2010 B.D.Dudson, S.Farley, M.V.Umansky, X.Q.Xu
 *
 * Contact: Ben Dudson, bd512@york.ac.uk
 * 
 * This file is part of BOUT++.
 *
 * BOUT++ is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * BOUT++ is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with BOUT++.  If not, see <http://www.gnu.org/licenses/>.
 *
 **************************************************************************/

#include <globals.hxx>
#include <boutcomm.hxx>

#include "pnetcdf.hxx"


#ifdef PNCDF

#include <pnetcdf.h> // Parallel NetCDF library

#include <utils.hxx>
#include <cmath>

#include <output.hxx>
#include <msg_stack.hxx>

#define CHKERR(ret)                                                                  \
  {                                                                                  \
    if (ret != NC_NOERR)                                                             \
      throw BoutException("pnetcdf line {:d}: {:s}", __LINE__, ncmpi_strerror(ret)); \
  }

// Define this to see loads of info messages
//#define NCDF_VERBOSE

PncFormat::PncFormat(Mesh* mesh_in) : DataFormat(mesh_in) {
  x0 = y0 = z0 = t0 = 0;
  lowPrecision = false;
  dimList = recDimList+1;
  default_rec = 0;
  rec_nr.clear();

  fname = nullptr;
}

PncFormat::PncFormat(const char *name, Mesh* mesh_in) : DataFormat(mesh_in) {
  x0 = y0 = z0 = t0 = 0;
  lowPrecision = false;
  dimList = recDimList+1;
  default_rec = 0;
  rec_nr.clear();

  openr(name);
}

PncFormat::~PncFormat() {
  close();
  rec_nr.clear();
}

bool PncFormat::openr(const char *name) {
  TRACE("PncFormat::openr");

  if(is_valid()) // Already open. Close then re-open
    close(); 

  int ret;
  MPI_Comm comm = BoutComm::get();
  
  /// Open file for reading
  
  ret = ncmpi_open(comm, name, NC_NOWRITE, MPI_INFO_NULL, &ncfile); CHKERR(ret);
  
  /// Get the dimensions from the file
  ret = ncmpi_inq_dimid(ncfile, "x", &xDim);
  if(ret != NC_NOERR) {
    output_warn.write("WARNING: NetCDF file should have an 'x' dimension\n");
  }
  
  ret = ncmpi_inq_dimid(ncfile, "y", &yDim);
  if(ret != NC_NOERR) {
    output_warn.write("WARNING: NetCDF file should have an 'y' dimension\n");
  }
  
  // z and t dimensions less crucial
  ret = ncmpi_inq_dimid(ncfile, "z", &zDim);
  if(ret != NC_NOERR) {
#ifdef NCDF_VERBOSE
    output_info.write("INFO: NetCDF file has no 'z' coordinate\n");
#endif
  }
  ret = ncmpi_inq_dimid(ncfile, "t", &tDim);
  if(ret != NC_NOERR) {
#ifdef NCDF_VERBOSE
    output_info.write("INFO: NetCDF file has no 'z' coordinate\n");
#endif
  }
  
  recDimList[0] = tDim;
  recDimList[1] = xDim;
  recDimList[2] = yDim;
  recDimList[3] = zDim;

  fname = copy_string(name);

  return true;
}

bool PncFormat::openw(const char *name, bool append) {
  TRACE("PncFormat::openw");
  
  if(is_valid()) // Already open. Close then re-open
    close(); 

  int ret; // Return status
  MPI_Comm comm = BoutComm::get();

  if(append) {
    ret = ncmpi_open(comm, name, NC_WRITE, MPI_INFO_NULL, &ncfile); CHKERR(ret);
    
    /// Get the dimensions from the file
    
    ret = ncmpi_inq_dimid(ncfile, "x", &xDim); CHKERR(ret);
    ret = ncmpi_inq_dimid(ncfile, "y", &yDim); CHKERR(ret);
    ret = ncmpi_inq_dimid(ncfile, "z", &zDim); CHKERR(ret);
    ret = ncmpi_inq_dimid(ncfile, "t", &tDim); CHKERR(ret);
    
    /// Test they're the right size (and t is unlimited)
    
    MPI_Offset len;
    ret = ncmpi_inq_dimlen(ncfile, xDim, &len); CHKERR(ret);
    if(len != mesh->GlobalNx)
      throw BoutException(
          "ERROR: x dimension length ({:d}) incompatible with mesh size ({:d})", (int)len,
          mesh->GlobalNx);

    ret = ncmpi_inq_dimlen(ncfile, yDim, &len); CHKERR(ret);
    if(len != mesh->GlobalNy)
      throw BoutException(
          "ERROR: y dimension length ({:d}) incompatible with mesh size ({:d})", (int)len,
          mesh->GlobalNy);

    ret = ncmpi_inq_dimlen(ncfile, zDim, &len); CHKERR(ret);
    if(len != mesh->GlobalNz)
      throw BoutException(
          "ERROR: z dimension length ({:d}) incompatible with mesh size ({:d})", (int)len,
          mesh->GlobalNz);

    // Get the size of the 't' dimension for records
    ret = ncmpi_inq_dimlen(ncfile, tDim, &len); CHKERR(ret);
    default_rec = (int) len;
    
  }else {
    /// Create file using collective operation
    ret = ncmpi_create(comm,
                       name,
                       NC_WRITE|NC_64BIT_OFFSET,
                       MPI_INFO_NULL,
                       &ncfile); CHKERR(ret);
    
    /// Add the dimensions
    
    ret = ncmpi_def_dim(ncfile, "x", mesh->GlobalNx, &xDim); CHKERR(ret);
    
    ret = ncmpi_def_dim(ncfile, "y", mesh->GlobalNy, &yDim); CHKERR(ret);
    
    ret = ncmpi_def_dim(ncfile, "z", mesh->GlobalNz, &zDim); CHKERR(ret);
    
    ret = ncmpi_def_dim(ncfile, "t", NC_UNLIMITED, &tDim); CHKERR(ret); // Unlimited 
    
    /// Finish define mode
    ret = ncmpi_enddef(ncfile); CHKERR(ret);

    default_rec = 0; // Starting at record 0
  }
  
  recDimList[0] = tDim;
  recDimList[1] = xDim;
  recDimList[2] = yDim;
  recDimList[3] = zDim;

  fname = copy_string(name);

  return true;
}

void PncFormat::close() {
  TRACE("PncFormat::close");
  
  if(!is_valid())
    return; // Already closed

  int ret = ncmpi_close(ncfile);
  free(fname);
  fname = nullptr;
}

void PncFormat::flush() {
  if(!is_valid())
    return;
  int ret = ncmpi_sync(ncfile);
}

const vector<int> PncFormat::getSize(const char *name) {
  TRACE("PncFormat::getSize");

  vector<int> size;

  if(!is_valid())
    return size;

  int ret, var;
  if(ret = ncmpi_inq_varid(ncfile, name, &var)) {
    // Variable not in file
    return size;
  }
  
  // Number of dimensions
  int nd;
  if(ret = ncmpi_inq_varndims (ncfile, var, &nd)) 
    return size; // Not valid

  if(nd == 0) {
    size.push_back(1);
    return size;
  }
  
  // Get the dimension IDs
  int dimid[nd];
  
  ret = ncmpi_inq_vardimid(ncfile, var, dimid); CHKERR(ret);
  
  for(int i=0;i<nd;i++) {
    // Get length of each dimension
    MPI_Offset len;
    ret = ncmpi_inq_dimlen(ncfile, dimid[i], &len); CHKERR(ret);
    size.push_back((int) len);
  }

  return size;
}

bool PncFormat::setGlobalOrigin(int x, int y, int z) {
  x0 = x;
  y0 = y;
  z0 = z;
  return true;
}

bool PncFormat::setRecord(int t) {
  t0 = t;

  return true;
}

bool PncFormat::read(int *data, const char *name, int lx, int ly, int lz) {
  TRACE("PncFormat::read(int)");

  if(!is_valid())
    return false;

  if((lx < 0) || (ly < 0) || (lz < 0))
    return false;

  int ret;
  int var;
  if(ret = ncmpi_inq_varid(ncfile, name, &var)) {
    // Variable not in file
#ifdef NCDF_VERBOSE
    output_info.write("INFO: Parallel NetCDF variable '{:s}' not found\n", name);
#endif
    return false;
  }
  
  // Number of dimensions
  int nd;
  ret = ncmpi_inq_varndims(ncfile, var, &nd); CHKERR(ret);
  
  if(nd == 0) {
    ret = ncmpi_get_var_int_all(ncfile, var, data); CHKERR(ret);
    return true;
  }
  
  MPI_Offset start[3], count[3];
  start[0] = x0; start[1] = y0; start[2] = z0;
  count[0] = lx; count[1] = ly; count[2] = lz;
  
  ret = ncmpi_get_vara_int_all(ncfile, var, start, count, data);
  
  return true;
}

bool PncFormat::read(int *var, const string &name, int lx, int ly, int lz) {
  return read(var, name.c_str(), lx, ly, lz);
}

// Helper functions to select the correct ncmpi function
int pnc_get_var_all(int ncfile, int var, double* data) {
  return ncmpi_get_var_double_all(ncfile, var, data);
}

int pnc_get_var_all(int ncfile, int var, float* data) {
  return ncmpi_get_var_float_all(ncfile, var, data);
}

int pnc_get_vara_all(int ncfile, int var, MPI_Offset* start, MPI_Offset* count, double* data) {
  return ncmpi_get_vara_double_all(ncfile, var, start, count, data);
}

int pnc_get_vara_all(int ncfile, int var, MPI_Offset* start, MPI_Offset* count, float* data) {
  return ncmpi_get_vara_float_all(ncfile, var, start, count, data);
}

bool PncFormat::read(BoutReal *data, const char *name, int lx, int ly, int lz) {
  TRACE("PncFormat::read(BoutReal)");

  if(!is_valid())
    return false;

  if((lx < 0) || (ly < 0) || (lz < 0))
    return false;
  
  int ret;
  int var;
  if(ret = ncmpi_inq_varid(ncfile, name, &var)) {
    // Variable not in file
#ifdef NCDF_VERBOSE
    output_info.write("INFO: Parallel NetCDF variable '{:s}' not found\n", name);
#endif
    return false;
  }
  
  // Number of dimensions
  int nd;
  ret = ncmpi_inq_varndims(ncfile, var, &nd); CHKERR(ret);
  
  if(nd == 0) {
    ret = pnc_get_var_all(ncfile, var, data); CHKERR(ret);
    return true;
  }
  
  
  MPI_Offset start[3], count[3];
  start[0] = x0; start[1] = y0; start[2] = z0;
  count[0] = lx; count[1] = ly; count[2] = lz;

  ret = pnc_get_vara_all(ncfile, var, start, count, data); CHKERR(ret);

  return true;
}

bool PncFormat::read(BoutReal *var, const string &name, int lx, int ly, int lz) {
  return read(var, name.c_str(), lx, ly, lz);
}

bool PncFormat::write(int *data, const char *name, int lx, int ly, int lz) {
  TRACE("PncFormat::write(int)");

  if(!is_valid())
    return false;

  if((lx < 0) || (ly < 0) || (lz < 0))
    return false;
  
  int nd = 0; // Number of dimensions
  if(lx != 0) nd = 1;
  if(ly != 0) nd = 2;
  if(lz != 0) nd = 3;

  int ret;
  int var;
  if(ret = ncmpi_inq_varid(ncfile, name, &var)) {
    // Variable not in file
    
    // Put into define mode
    ret = ncmpi_redef(ncfile); CHKERR(ret); 
    // Define variable
    ret = ncmpi_def_var(ncfile, name, NC_INT, nd, dimList, &var); CHKERR(ret);
    // Back out of define mode
    ret = ncmpi_enddef(ncfile); CHKERR(ret);
  }
  
  if(nd == 0) {
    // Writing a scalar
    ret = ncmpi_put_var_int_all(ncfile, var, data); CHKERR(ret);
    return true;
  }
  
  // An array of values
  
  MPI_Offset start[3], count[3];
  start[0] = x0; start[1] = y0; start[2] = z0;
  count[0] = lx; count[1] = ly; count[2] = lz;
  
  ret = ncmpi_put_vara_int_all(ncfile, var, start, count, data); CHKERR(ret);

  return true;
}

bool PncFormat::write(int *var, const string &name, int lx, int ly, int lz) {
  return write(var, name.c_str(), lx, ly, lz);
}

// Helper functions to select the correct ncmpi function
int pnc_put_var_all(int ncfile, int var, double* data) {
  return ncmpi_put_var_double_all(ncfile, var, data);
}

int pnc_put_var_all(int ncfile, int var, float* data) {
  return ncmpi_put_var_float_all(ncfile, var, data);
}

// Helper functions to select the correct ncmpi function
int pnc_put_vara_all(int ncfile, int var, MPI_Offset* start, MPI_Offset* count, double* data) {
  return ncmpi_put_vara_double_all(ncfile, var, start, count, data);
}

int pnc_put_vara_all(int ncfile, int var, MPI_Offset* start, MPI_Offset* count, float* data) {
  return ncmpi_put_vara_float_all(ncfile, var, start, count, data);
}

bool PncFormat::write(BoutReal *data, const char *name, int lx, int ly, int lz) {
  TRACE("PncFormat::write(BoutReal)");

  if(!is_valid())
    return false;

  if((lx < 0) || (ly < 0) || (lz < 0))
    return false;

  int nd = 0; // Number of dimensions
  if(lx != 0) nd = 1;
  if(ly != 0) nd = 2;
  if(lz != 0) nd = 3;
  
  int ret;
  int var;
  if(ret = ncmpi_inq_varid(ncfile, name, &var)) {
    // Variable not in file
    
    // Put into define mode
    ret = ncmpi_redef(ncfile); CHKERR(ret); 

    nc_type type = (lowPrecision) ? NC_FLOAT : NC_DOUBLE;

    // Define variable
    ret = ncmpi_def_var(ncfile, name, type, nd, dimList, &var); CHKERR(ret);
    // Back out of define mode
    ret = ncmpi_enddef(ncfile); CHKERR(ret);
  }
  
  if(nd == 0) {
    // Writing a scalar
    ret = pnc_put_var_all(ncfile, var, data); CHKERR(ret);
    return true;
  }
  
  if(lowPrecision) {
    // An out of range value can make the conversion
    // corrupt the whole dataset. Make sure everything
    // is in the range of a float
    
    for(int i=0;i<lx*ly*lz;i++) {
      if(data[i] > 1e20)
	data[i] = 1e20;
      if(data[i] < -1e20)
	data[i] = -1e20;
    }
  }
  
  // Similarly, non-finite numbers can have nasty effects
  for(int i=0;i<lx*ly*lz;i++) {
    if(!finite(data[i]))
      data[i] = 0.0;
  }

  MPI_Offset start[3], count[3];
  start[0] = x0; start[1] = y0; start[2] = z0;
  count[0] = lx; count[1] = ly; count[2] = lz;
  
  ret = pnc_put_vara_all(ncfile, var, start, count, data); CHKERR(ret);

  return true;
}

bool PncFormat::write(BoutReal *var, const string &name, int lx, int ly, int lz) {
  return write(var, name.c_str(), lx, ly, lz);
}

/***************************************************************************
 * Record-based (time-dependent) data
 ***************************************************************************/

bool PncFormat::read_rec(int *data, const char *name, int lx, int ly, int lz) {
  if(!is_valid())
    return false;

  if((lx < 0) || (ly < 0) || (lz < 0))
    return false;
  
  int ret;
  int var;
  if(ret = ncmpi_inq_varid(ncfile, name, &var)) {
    // Variable not in file
#ifdef NCDF_VERBOSE
    output_info.write("INFO: Parallel NetCDF variable '{:s}' not found\n", name);
#endif
    return false;
  }
  
  // NOTE: Probably should do something here to check t0

  MPI_Offset start[4], count[4];
  start[0] = t0; start[1] = x0; start[2] = y0; start[3] = z0;
  count[0] = 1;  count[1] = lx; count[2] = ly; count[3] = lz;
  
  ret = ncmpi_get_vara_int_all(ncfile, var, start, count, data);
  
  return true;
}

bool PncFormat::read_rec(int *var, const string &name, int lx, int ly, int lz) {
  return read_rec(var, name.c_str(), lx, ly, lz);
}

bool PncFormat::read_rec(BoutReal *data, const char *name, int lx, int ly, int lz) {
  if(!is_valid())
    return false;

  if((lx < 0) || (ly < 0) || (lz < 0))
    return false;
  
  int ret;
  int var;
  if(ret = ncmpi_inq_varid(ncfile, name, &var)) {
    // Variable not in file
#ifdef NCDF_VERBOSE
    output_info.write("INFO: Parallel NetCDF variable '{:s}' not found\n", name);
#endif
    return false;
  }
  
  // NOTE: Probably should do something here to check t0

  MPI_Offset start[4], count[4];
  start[0] = t0; start[1] = x0; start[2] = y0; start[3] = z0;
  count[0] = 1;  count[1] = lx; count[2] = ly; count[3] = lz;

  ret = pnc_get_vara_all(ncfile, var, start, count, data); CHKERR(ret);
  
  return true;
}

bool PncFormat::read_rec(BoutReal *var, const string &name, int lx, int ly, int lz) {
  return read_rec(var, name.c_str(), lx, ly, lz);
}

bool PncFormat::write_rec(int *data, const char *name, int lx, int ly, int lz) {
  if(!is_valid())
    return false;

  if((lx < 0) || (ly < 0) || (lz < 0))
    return false;

  int nd = 1; // Number of dimensions
  if(lx != 0) nd = 2;
  if(ly != 0) nd = 3;
  if(lz != 0) nd = 4;
  
  int ret;
  int var;
  if(ret = ncmpi_inq_varid(ncfile, name, &var)) {
    // Variable not in file
    
    // Put into define mode
    ret = ncmpi_redef(ncfile); CHKERR(ret); 
    // Define variable
    ret = ncmpi_def_var(ncfile, name, NC_INT, nd, recDimList, &var); CHKERR(ret);
    // Back out of define mode
    ret = ncmpi_enddef(ncfile); CHKERR(ret);
    
    rec_nr[name] = default_rec; // Starting record
  }else {
    // Get record number
    if(rec_nr.find(name) == rec_nr.end()) {
      // Add to map
      rec_nr[name] = default_rec;
    }
  }
  
  MPI_Offset start[4], count[4];
  start[0] = rec_nr[name]; start[1] = x0; start[2] = y0; start[3] = z0;
  count[0] = 1;  count[1] = lx; count[2] = ly; count[3] = lz;
  
  ret = ncmpi_put_vara_int_all(ncfile, var, start, count, data); CHKERR(ret);
  
  // Increment record number
  rec_nr[name] += 1;

  return true;
}

bool PncFormat::write_rec(int *var, const string &name, int lx, int ly, int lz) {
  return write_rec(var, name.c_str(), lx, ly, lz);
}

bool PncFormat::write_rec(BoutReal *data, const char *name, int lx, int ly, int lz) {
  TRACE("PncFormat::write_rec(BoutReal)");

  if(!is_valid())
    return false;

  if((lx < 0) || (ly < 0) || (lz < 0))
    return false;

  int nd = 1; // Number of dimensions
  if(lx != 0) nd = 2;
  if(ly != 0) nd = 3;
  if(lz != 0) nd = 4;

  int ret;
  int var;
  if(ret = ncmpi_inq_varid(ncfile, name, &var)) {
    // Variable not in file
    
    // Put into define mode
    ret = ncmpi_redef(ncfile); CHKERR(ret); 
    // Define variable
    nc_type type = (lowPrecision) ? NC_FLOAT : NC_DOUBLE;
    ret = ncmpi_def_var(ncfile, name, type, nd, recDimList, &var); CHKERR(ret);
    // Back out of define mode
    ret = ncmpi_enddef(ncfile); CHKERR(ret);
    
    rec_nr[name] = default_rec; // Starting record
  }else {
    // Get record number
    if(rec_nr.find(name) == rec_nr.end()) {
      // Add to map
      rec_nr[name] = default_rec;
    }
  }

#ifdef NCDF_VERBOSE
  output_info.write("INFO: NetCDF writing record {:d} of '{:s}' in '{:s}'\n",t, name, fname); 
#endif

  if(lowPrecision) {
    // An out of range value can make the conversion
    // corrupt the whole dataset. Make sure everything
    // is in the range of a float
    
    for(int i=0;i<lx*ly*lz;i++) {
      if(data[i] > 1e20)
	data[i] = 1e20;
      if(data[i] < -1e20)
	data[i] = -1e20;
    }
  }
  
  for(int i=0;i<lx*ly*lz;i++) {
    if(!finite(data[i]))
      data[i] = 0.0;
  }

  MPI_Offset start[4], count[4];
  start[0] = rec_nr[name]; start[1] = x0; start[2] = y0; start[3] = z0;
  count[0] = 1;  count[1] = lx; count[2] = ly; count[3] = lz;

  // Add the record

  ret = pnc_put_vara_all(ncfile, var, start, count, data); CHKERR(ret);
  
  // Increment record number
  rec_nr[name] += 1;

  return true;
}

bool PncFormat::write_rec(BoutReal *var, const string &name, int lx, int ly, int lz) {
  return write_rec(var, name.c_str(), lx, ly, lz);
}

/***************************************************************************
 * Private functions
 ***************************************************************************/

#endif // PNCDF

