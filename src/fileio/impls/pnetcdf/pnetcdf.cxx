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

#include <utils.hxx>
#include <cmath>

#include <output.hxx>
#include <msg_stack.hxx>

#define CHKERR(ret) { if( ret != NC_NOERR ) throw new BoutException("pnetcdf line %d: %s", __LINE__, ncmpi_strerror(ret)); }

// Define this to see loads of info messages
//#define NCDF_VERBOSE

PncFormat::PncFormat() {
  x0 = y0 = z0 = t0 = 0;
  lowPrecision = false;
  dimList = recDimList+1;
  default_rec = 0;
  rec_nr.clear();

  fname = NULL;
}

PncFormat::PncFormat(const char *name) {
  x0 = y0 = z0 = t0 = 0;
  lowPrecision = false;
  dimList = recDimList+1;
  default_rec = 0;
  rec_nr.clear();

  openr(name);
}

PncFormat::PncFormat(const string &name) {
  dataFile = NULL;
  x0 = y0 = z0 = t0 = 0;
  recDimList = new const NcDim*[4];
  dimList = recDimList+1;
  lowPrecision = false;

  default_rec = 0;
  rec_nr.clear();

  openr(name);
}

PncFormat::~PncFormat() {
  close();
  rec_nr.clear();
}

bool PncFormat::openr(const char *name) {
#ifdef CHECK
  msg_stack.push("PncFormat::openr");
#endif

  int ret;
  MPI_Comm comm = BoutComm::get();
  
  /// Open file for reading
  ret = ncmpi_open(comm, name, NC_NOWRITE, MPI_INFO_NULL, &ncfile); CHKERR(ret);

  if(dataFile != NULL) // Already open. Close then re-open
    close(); 


  /// Get the dimensions from the file
  ret = ncmpi_inq_dimid(ncid, "x", &xDim);
  if(ret != NC_NOERR) {
    output.write("WARNING: NetCDF file should have an 'x' dimension\n");
  }
  
  ret = ncmpi_inq_dimid(ncid, "y", &yDim);
  if(ret != NC_NOERR) {
    output.write("WARNING: NetCDF file should have an 'y' dimension\n");
  }
  
  // z and t dimensions less crucial
  ret = ncmpi_inq_dimid(ncid, "z", &zDim);
  if(ret != NC_NOERR) {
#ifdef NCDF_VERBOSE
    output.write("INFO: NetCDF file has no 'z' coordinate\n");
#endif
  }
  ret = ncmpi_inq_dimid(ncid, "t", &tDim);
  if(ret != NC_NOERR) {
#ifdef NCDF_VERBOSE
    output.write("INFO: NetCDF file has no 'z' coordinate\n");
#endif
  }
  
  recDimList[0] = tDim;
  recDimList[1] = xDim;
  recDimList[2] = yDim;
  recDimList[3] = zDim;

  fname = copy_string(name);

#ifdef CHECK
  msg_stack.pop();
#endif

  return true;
}

bool PncFormat::openw(const char *name, bool append) {
#ifdef CHECK
  msg_stack.push("PncFormat::openw");
#endif
  
  if(dataFile != NULL) // Already open. Close then re-open
    close(); 

  int ret; // Return status
  MPI_Comm comm = BoutComm::get();

  if(append) {
    ret = ncmpi_open(MPI_COMM_WORLD, "foo.nc", NC_WRITE, MPI_INFO_NULL, &ncid); CHKERR(ret);
    
    /// Get the dimensions from the file
    
    ret = ncmpi_inq_dimid(ncid, "x", &xDim); CHKERR(ret);
    ret = ncmpi_inq_dimid(ncid, "y", &yDim); CHKERR(ret);
    ret = ncmpi_inq_dimid(ncid, "z", &zDim); CHKERR(ret);
    ret = ncmpi_inq_dimid(ncid, "t", &tDim); CHKERR(ret);
    
    /// Test they're the right size (and t is unlimited)
    
    MPI_Offset len;
    ret = ncmpi_inq_dimlen(ncfile, xDim, &len); CHKERR(ret);
    if(len != mesh->GlobalNx)
      throw new BoutException("ERROR: x dimension length (%d) incompatible with mesh size (%d)", (int) len, mesh->GlobalNx);

    ret = ncmpi_inq_dimlen(ncfile, yDim, &len); CHKERR(ret);
    if(len != mesh->GlobalNy)
      throw new BoutException("ERROR: y dimension length (%d) incompatible with mesh size (%d)", (int) len, mesh->GlobalNy);
    
    ret = ncmpi_inq_dimlen(ncfile, zDim, &len); CHKERR(ret);
    if(len != mesh->GlobalNz)
      throw new BoutException("ERROR: z dimension length (%d) incompatible with mesh size (%d)", (int) len, mesh->GlobalNz);
    
    // Get the size of the 't' dimension for records
    ret = ncmpi_inq_dimlen(ncfile, tDim, &len); CHKERR(ret);
    default_rec = (int) len;
    
  }else {
    /// Create file using collective operation
    ret = ncmpi_create(comm,
                       name,
                       NC_WRITE|NC_64BIT_OFFSET,
                       MPI_INFO_NULL,
                       &ncfile);
    if( ret != NC_NOERR ) throw new BoutException(ncmpi_strerror(err));
    
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

#ifdef CHECK
  msg_stack.pop();
#endif

  return true;
}

bool PncFormat::is_valid()
{
  if(dataFile == NULL)
    return false;
  return dataFile->is_valid();
}

void PncFormat::close() {
#ifdef CHECK
  msg_stack.push("PncFormat::close");
#endif
  
  int ret = ncmpi_close(ncfile);
  free(fname);
  fname = NULL;

#ifdef CHECK
  msg_stack.pop();
#endif
}

const vector<int> PncFormat::getSize(const char *name)
{
  vector<int> size;

#ifdef NCDF_VERBOSE
  NcError err(NcError::verbose_nonfatal);
#else
  NcError err(NcError::silent_nonfatal);
#endif

  if(!is_valid())
    return size;

#ifdef CHECK
  msg_stack.push("PncFormat::getSize");
#endif

  NcVar *var;
  
  if(!(var = dataFile->get_var(name)))
    return size;
  
  if(!var->is_valid())
    return size;

  int nd = var->num_dims(); // Number of dimensions

  if(nd == 0) {
    size.push_back(1);
    return size;
  }
  
  long *ls = var->edges();
  for(int i=0;i<nd;i++)
    size.push_back((int) ls[i]);

  delete[] ls;
  
#ifdef CHECK
  msg_stack.pop();
#endif

  return size;
}

const vector<int> PncFormat::getSize(const string &var)
{
  return getSize(var.c_str());
}

bool PncFormat::setOrigin(int x, int y, int z)
{
  x0 = x;
  y0 = y;
  z0 = z;
  
  return true;
}

bool PncFormat::setRecord(int t)
{
  t0 = t;

  return true;
}

bool PncFormat::read(int *data, const char *name, int lx, int ly, int lz)
{
  if(!is_valid())
    return false;

  if((lx < 0) || (ly < 0) || (lz < 0))
    return false;
  
#ifdef CHECK
  msg_stack.push("PncFormat::read(int)");
#endif

  NcVar *var;
  
  if(!(var = dataFile->get_var(name))) {
#ifdef NCDF_VERBOSE
    output.write("INFO: NetCDF variable '%s' not found\n", name);
#endif
#ifdef CHECK
    msg_stack.pop();
#endif
    return false;
  }
  
  long cur[3], counts[3];
  cur[0] = x0;    cur[1] = y0;    cur[2] = z0;
  counts[0] = lx; counts[1] = ly; counts[2] = lz;

  if(!(var->set_cur(cur))) {
#ifdef NCDF_VERBOSE
    output.write("INFO: NetCDF Could not set cur(%d,%d,%d) for variable '%s'\n", 
		 x0,y0,z0, name);
#endif
#ifdef CHECK
    msg_stack.pop();
#endif
    return false;
  }

  if(!(var->get(data, counts))) {
#ifdef NCDF_VERBOSE
    output.write("INFO: NetCDF could not read data for '%s'\n", name);
#endif
#ifdef CHECK
    msg_stack.pop();
#endif
    return false;
  }
  
#ifdef CHECK
  msg_stack.pop();
#endif

  return true;
}

bool PncFormat::read(int *var, const string &name, int lx, int ly, int lz)
{
  return read(var, name.c_str(), lx, ly, lz);
}

bool PncFormat::read(BoutReal *data, const char *name, int lx, int ly, int lz)
{
  if(!is_valid())
    return false;

  if((lx < 0) || (ly < 0) || (lz < 0))
    return false;

#ifdef CHECK
  msg_stack.push("PncFormat::read(BoutReal)");
#endif

  // Create an error object so netCDF doesn't exit
#ifdef NCDF_VERBOSE
  NcError err(NcError::verbose_nonfatal);
#else
  NcError err(NcError::silent_nonfatal);
#endif

  NcVar *var;
  
  if(!(var = dataFile->get_var(name))) {
#ifdef CHECK
    msg_stack.pop();
#endif
    return false;
  }
  
  long cur[3], counts[3];
  cur[0] = x0;    cur[1] = y0;    cur[2] = z0;
  counts[0] = lx; counts[1] = ly; counts[2] = lz;

  if(!(var->set_cur(cur))) {
#ifdef CHECK
    msg_stack.pop();
#endif
    return false;
  }

  if(!(var->get(data, counts))) {
#ifdef CHECK
    msg_stack.pop();
#endif
    return false;
  }
  
#ifdef CHECK
  msg_stack.pop();
#endif

  return true;
}

bool PncFormat::read(BoutReal *var, const string &name, int lx, int ly, int lz) {
  return read(var, name.c_str(), lx, ly, lz);
}

bool PncFormat::write(int *data, const char *name, int lx, int ly, int lz) {
  if(!is_valid())
    return false;

  if((lx < 0) || (ly < 0) || (lz < 0))
    return false;
  
  int nd = 0; // Number of dimensions
  if(lx != 0) nd = 1;
  if(ly != 0) nd = 2;
  if(lz != 0) nd = 3;

#ifdef CHECK
  msg_stack.push("PncFormat::write(int)");
#endif

  int ret;
  int var;
  if(!(var = dataFile->get_var(name))) {
    // Variable not in file
    
    // Put into define mode
    ret = ncmpi_redef(ncfile); CHKERR(ret); 
    // Define variable
    ret = ncmpi_def_var(ncfile, name, NC_INT, nd, dimList, &var); CHKERR(ret);
  }
  
  long cur[3], counts[3];
  cur[0] = x0;    cur[1] = y0;    cur[2] = z0;
  counts[0] = lx; counts[1] = ly; counts[2] = lz;

  if(!(var->set_cur(cur)))
    return false;
  
  if(!(var->put(data, counts)))
    return false;

#ifdef CHECK
  msg_stack.pop();
#endif

  return true;
}

bool PncFormat::write(int *var, const string &name, int lx, int ly, int lz) {
  return write(var, name.c_str(), lx, ly, lz);
}

bool PncFormat::write(BoutReal *data, const char *name, int lx, int ly, int lz)
{
  if(!is_valid())
    return false;

  if((lx < 0) || (ly < 0) || (lz < 0))
    return false;
  
#ifdef CHECK
  msg_stack.push("PncFormat::write(BoutReal)");
#endif

  int nd = 0; // Number of dimensions
  if(lx != 0) nd = 1;
  if(ly != 0) nd = 2;
  if(lz != 0) nd = 3;
  
#ifdef NCDF_VERBOSE
  NcError err(NcError::verbose_nonfatal);
#else
  NcError err(NcError::silent_nonfatal);
#endif

  NcVar *var;
  if(!(var = dataFile->get_var(name))) {
    // Variable not in file, so add it.
    if(lowPrecision) {
      var = dataFile->add_var(name, ncFloat, nd, dimList);
    }else
      var = dataFile->add_var(name, ncDouble, nd, dimList);

    if(!var->is_valid()) {
      output.write("ERROR: NetCDF could not add BoutReal '%s' to file '%s'\n", name, fname);
      return false;
    }
  }  

  long cur[3], counts[3];
  cur[0] = x0;    cur[1] = y0;    cur[2] = z0;
  counts[0] = lx; counts[1] = ly; counts[2] = lz;

  if(!(var->set_cur(cur)))
    return false;

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

  if(!(var->put(data, counts)))
    return false;

#ifdef CHECK
  msg_stack.pop();
#endif

  return true;
}

bool PncFormat::write(BoutReal *var, const string &name, int lx, int ly, int lz)
{
  return write(var, name.c_str(), lx, ly, lz);
}

/***************************************************************************
 * Record-based (time-dependent) data
 ***************************************************************************/

bool PncFormat::read_rec(int *data, const char *name, int lx, int ly, int lz)
{
  if(!is_valid())
    return false;

  if((lx < 0) || (ly < 0) || (lz < 0))
    return false;

  // Create an error object so netCDF doesn't exit
#ifdef NCDF_VERBOSE
  NcError err(NcError::verbose_nonfatal);
#else
  NcError err(NcError::silent_nonfatal);
#endif

  NcVar *var;
  
  if(!(var = dataFile->get_var(name)))
    return false;
  
  // NOTE: Probably should do something here to check t0

  long cur[4], counts[4];
  cur[0] = t0; cur[1] = x0;    cur[2] = y0;    cur[3] = z0;
  counts[0] = 1; counts[1] = lx; counts[2] = ly; counts[3] = lz;
  
  if(!(var->set_cur(cur)))
    return false;
  
  if(!(var->get(data, counts)))
    return false;
  
  return true;
}

bool PncFormat::read_rec(int *var, const string &name, int lx, int ly, int lz)
{
  return read_rec(var, name.c_str(), lx, ly, lz);
}

bool PncFormat::read_rec(BoutReal *data, const char *name, int lx, int ly, int lz)
{
  if(!is_valid())
    return false;

  if((lx < 0) || (ly < 0) || (lz < 0))
    return false;
  
  // Create an error object so netCDF doesn't exit
#ifdef NCDF_VERBOSE
  NcError err(NcError::verbose_nonfatal);
#else
  NcError err(NcError::silent_nonfatal);
#endif

  NcVar *var;
  
  if(!(var = dataFile->get_var(name)))
    return false;
  
  // NOTE: Probably should do something here to check t0

  long cur[4], counts[4];
  cur[0] = t0; cur[1] = x0;    cur[2] = y0;    cur[3] = z0;
  counts[0] = 1; counts[1] = lx; counts[2] = ly; counts[3] = lz;
  
  if(!(var->set_cur(cur)))
    return false;
  
  if(!(var->get(data, counts)))
    return false;
  
  return true;
}

bool PncFormat::read_rec(BoutReal *var, const string &name, int lx, int ly, int lz)
{
  return read_rec(var, name.c_str(), lx, ly, lz);
}

bool PncFormat::write_rec(int *data, const char *name, int lx, int ly, int lz)
{
  if(!is_valid())
    return false;

  if((lx < 0) || (ly < 0) || (lz < 0))
    return false;

  int nd = 1; // Number of dimensions
  if(lx != 0) nd = 2;
  if(ly != 0) nd = 3;
  if(lz != 0) nd = 4;

#ifdef NCDF_VERBOSE
  NcError err(NcError::verbose_nonfatal);
#else
  NcError err(NcError::silent_nonfatal);
#endif

  NcVar *var;
  
  // Try to find variable
  if(!(var = dataFile->get_var(name))) {
    // Need to add to file

    var = dataFile->add_var(name, ncInt, nd, recDimList);

    rec_nr[name] = default_rec; // Starting record

    if(!var->is_valid()) {
#ifdef NCDF_VERBOSE
      output.write("ERROR: NetCDF Could not add variable '%s' to file '%s'\n", name, fname);
#endif
      return false;
    }
  }else {
    // Get record number
    if(rec_nr.find(name) == rec_nr.end()) {
      // Add to map
      rec_nr[name] = default_rec;
    }
  }
  
  if(!var->put_rec(data, rec_nr[name]))
    return false;

  var->sync();
  
  // Increment record number
  rec_nr[name] = rec_nr[name] + 1;

  return true;
}

bool PncFormat::write_rec(int *var, const string &name, int lx, int ly, int lz)
{
  return write_rec(var, name.c_str(), lx, ly, lz);
}

bool PncFormat::write_rec(BoutReal *data, const char *name, int lx, int ly, int lz)
{
  if(!is_valid())
    return false;

  if((lx < 0) || (ly < 0) || (lz < 0))
    return false;

#ifdef CHECK
  msg_stack.push("PncFormat::write_rec(BoutReal)");
#endif

  int nd = 1; // Number of dimensions
  if(lx != 0) nd = 2;
  if(ly != 0) nd = 3;
  if(lz != 0) nd = 4;

#ifdef NCDF_VERBOSE
  NcError err(NcError::verbose_nonfatal);
#else
  NcError err(NcError::silent_nonfatal);
#endif

  NcVar *var;

  // Try to find variable
  if(!(var = dataFile->get_var(name))) {
    // Need to add to file
    
    NcType vartype = ncDouble;
    if(lowPrecision)
      vartype = ncFloat;
  
    var = dataFile->add_var(name, vartype, nd, recDimList);
    
    rec_nr[name] = default_rec; // Starting record

    if(!var->is_valid()) {
#ifdef NCDF_VERBOSE
      output.write("ERROR: NetCDF Could not add variable '%s' to file '%s'\n", name, fname);
#endif
#ifdef CHECK
  msg_stack.pop();
#endif
      return false;
    }
  }else {
    // Get record number
    if(rec_nr.find(name) == rec_nr.end()) {
      // Add to map
      rec_nr[name] = default_rec;
    }
  }
  
  int t = rec_nr[name];

#ifdef NCDF_VERBOSE
  output.write("INFO: NetCDF writing record %d of '%s' in '%s'\n",t, name, fname); 
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

  // Add the record
  if(!var->put_rec(data, t))
    return false;

  var->sync();
  
  // Increment record number
  rec_nr[name] = rec_nr[name] + 1;

#ifdef CHECK
  msg_stack.pop();
#endif

  return true;
}

bool PncFormat::write_rec(BoutReal *var, const string &name, int lx, int ly, int lz) {
  return write_rec(var, name.c_str(), lx, ly, lz);
}

/***************************************************************************
 * Private functions
 ***************************************************************************/

#endif // PNCDF

