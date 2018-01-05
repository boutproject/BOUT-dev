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
#include "nc_format.hxx"

#ifdef NCDF

#include <utils.hxx>
#include <cmath>

#include <output.hxx>
#include <msg_stack.hxx>

// Define this to see loads of info messages
//#define NCDF_VERBOSE

NcFormat::NcFormat() {
  dataFile = NULL;
  x0 = y0 = z0 = t0 = 0;
  recDimList = new const NcDim*[4];
  dimList = recDimList+1;
  lowPrecision = false;

  default_rec = 0;
  rec_nr.clear();

  fname = NULL;
}

NcFormat::NcFormat(const char *name) {
  dataFile = NULL;
  x0 = y0 = z0 = t0 = 0;
  recDimList = new const NcDim*[4];
  dimList = recDimList+1;
  lowPrecision = false;

  default_rec = 0;
  rec_nr.clear();

  openr(name);
}

NcFormat::~NcFormat() {
  delete[] recDimList;
  close();
  rec_nr.clear();
}

bool NcFormat::openr(const char *name) {
  TRACE("NcFormat::openr");

  if(dataFile != NULL) // Already open. Close then re-open
    close(); 

  // Create an error object so netCDF doesn't exit
#ifdef NCDF_VERBOSE
  NcError err(NcError::verbose_nonfatal);
#else
  NcError err(NcError::silent_nonfatal);
#endif

  dataFile = new NcFile(name, NcFile::ReadOnly);

  if(!dataFile->is_valid()) {
    delete dataFile;
    dataFile = NULL;
    return false;
  }

  /// Get the dimensions from the file

  if(!(xDim = dataFile->get_dim("x"))) {
    output_warn.write("WARNING: NetCDF file should have an 'x' dimension\n");
    /*
    delete dataFile;
    dataFile = NULL;
    return false;
    */
    xDim = NULL;
  }else if(mesh != NULL) {
    // Check that the dimension size is correct
    if(xDim->size() != mesh->LocalNx) {
      throw BoutException("X dimension incorrect. Expected %d, got %d", mesh->LocalNx, xDim->size());
    }
  }
  
  if(!(yDim = dataFile->get_dim("y"))) {
    output_warn.write("WARNING: NetCDF file should have a 'y' dimension\n");
    /*
    delete dataFile;
    dataFile = NULL;
    return false;
    */
    yDim = NULL;
  }else if(mesh != NULL) {
    // Check that the dimension size is correct
    if(yDim->size() != mesh->LocalNy) {
      throw BoutException("Y dimension incorrect. Expected %d, got %d", mesh->LocalNy, yDim->size());
    }
  }
  
  if(!(zDim = dataFile->get_dim("z"))) {
    // Z dimension optional, and could be any size (Fourier harmonics)
#ifdef NCDF_VERBOSE
    output_info.write("INFO: NetCDF file has no 'z' coordinate\n");
#endif
    zDim = NULL;
  }else if(mesh != NULL) {
    // Check that the dimension size is correct
    if(zDim->size() != mesh->LocalNz) {
      throw BoutException("Z dimension incorrect. Expected %d, got %d", mesh->LocalNz, zDim->size());
    }
  }
  
  if(!(tDim = dataFile->get_dim("t"))) {
    // T dimension optional
#ifdef NCDF_VERBOSE
    output_info.write("INFO: NetCDF file has no 't' coordinate\n");
#endif
    tDim = NULL;
  }
  
  recDimList[0] = tDim;
  recDimList[1] = xDim;
  recDimList[2] = yDim;
  recDimList[3] = zDim;

  fname = copy_string(name);

  return true;
}

bool NcFormat::openw(const char *name, bool append) {

  TRACE("NcFormat::openw");

  // Create an error object so netCDF doesn't exit
#ifdef NCDF_VERBOSE
  NcError err(NcError::verbose_nonfatal);
#else
  NcError err(NcError::silent_nonfatal);
#endif

  if(dataFile != NULL) // Already open. Close then re-open
    close(); 

  if(append) {
    dataFile = new NcFile(name, NcFile::Write);

    if(!dataFile->is_valid()) {
      delete dataFile;
      dataFile = NULL;
      return false;
    }

    /// Get the dimensions from the file

    if(!(xDim = dataFile->get_dim("x"))) {
      output_error.write("ERROR: NetCDF file should have an 'x' dimension\n");
      delete dataFile;
      dataFile = NULL;
      return false;
    }

    if(!(yDim = dataFile->get_dim("y"))) {
      output_error.write("ERROR: NetCDF file should have a 'y' dimension\n");
      delete dataFile;
      dataFile = NULL;
      return false;
    }

    if(!(zDim = dataFile->get_dim("z"))) {
      output_error.write("ERROR: NetCDF file should have a 'z' dimension\n");
      delete dataFile;
      dataFile = NULL;
      return false;
    }

    if(!(tDim = dataFile->get_dim("t"))) {
      output_error.write("ERROR: NetCDF file should have a 't' dimension\n");
      delete dataFile;
      dataFile = NULL;
      return false;
    }

    /// Test they're the right size (and t is unlimited)
    
    if((xDim->size() != mesh->LocalNx) || (yDim->size() != mesh->LocalNy) || (zDim->size() != mesh->LocalNz)
       || (!tDim->is_unlimited()) ) {
      delete dataFile;
      dataFile = NULL;
      return false;
    }

    // Get the size of the 't' dimension for records
    default_rec = tDim->size();
    
  }else {
    dataFile = new NcFile(name, NcFile::Replace);
    
    if(!dataFile->is_valid()) {
      delete dataFile;
      dataFile = NULL;
      return false;
    }

    /// Add the dimensions
    
    if(!(xDim = dataFile->add_dim("x", mesh->LocalNx))) {
      delete dataFile;
      dataFile = NULL;
      return false;
    }
  
    if(!(yDim = dataFile->add_dim("y", mesh->LocalNy))) {
      delete dataFile;
      dataFile = NULL;
      return false;
    }
    
    if(!(zDim = dataFile->add_dim("z", mesh->LocalNz))) {
      delete dataFile;
      dataFile = NULL;
      return false;
    }
    
    if(!(tDim = dataFile->add_dim("t"))) { // unlimited dimension
      delete dataFile;
      dataFile = NULL;
      return false;
    }

    default_rec = 0; // Starting at record 0
  }

  recDimList[0] = tDim;
  recDimList[1] = xDim;
  recDimList[2] = yDim;
  recDimList[3] = zDim;

  fname = copy_string(name);

  return true;
}

bool NcFormat::is_valid() {
  if(dataFile == NULL)
    return false;
  return dataFile->is_valid();
}

void NcFormat::close() {
  if(dataFile == NULL)
    return;
  
  TRACE("NcFormat::close");

#ifdef NCDF_VERBOSE
  NcError err(NcError::verbose_nonfatal);
#else
  NcError err(NcError::silent_nonfatal);
#endif

  dataFile->close();
  delete dataFile;
  dataFile = NULL;
  
  free(fname);
  fname = NULL;
}

void NcFormat::flush() {
  if(!is_valid())
    return;
  
  dataFile->sync();
}

const vector<int> NcFormat::getSize(const char *name) {
  vector<int> size;

#ifdef NCDF_VERBOSE
  NcError err(NcError::verbose_nonfatal);
#else
  NcError err(NcError::silent_nonfatal);
#endif

  if(!is_valid())
    return size;

  TRACE("NcFormat::getSize");

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
    size.push_back(static_cast<int>(ls[i]));

  delete[] ls;

  return size;
}

const vector<int> NcFormat::getSize(const string &var) {
  return getSize(var.c_str());
}

bool NcFormat::setGlobalOrigin(int x, int y, int z) {
  x0 = x;
  y0 = y;
  z0 = z;
  
  return true;
}

bool NcFormat::setRecord(int t) {
  t0 = t;

  return true;
}

bool NcFormat::read(int *data, const char *name, int lx, int ly, int lz) {
  if(!is_valid())
    return false;

  if((lx < 0) || (ly < 0) || (lz < 0))
    return false;
  
  // Check for valid name
  checkName(name);

  TRACE("NcFormat::read(int)");

  // Create an error object so netCDF doesn't exit
#ifdef NCDF_VERBOSE
  NcError err(NcError::verbose_nonfatal);
#else
  NcError err(NcError::silent_nonfatal);
#endif

  NcVar *var;
  
  if(!(var = dataFile->get_var(name))) {
#ifdef NCDF_VERBOSE
    output_info.write("INFO: NetCDF variable '%s' not found\n", name);
#endif
    return false;
  }
  
  long cur[3], counts[3];
  cur[0] = x0;    cur[1] = y0;    cur[2] = z0;
  counts[0] = lx; counts[1] = ly; counts[2] = lz;

  if(!(var->set_cur(cur))) {
#ifdef NCDF_VERBOSE
    output_info.write("INFO: NetCDF Could not set cur(%d,%d,%d) for variable '%s'\n", 
		 x0,y0,z0, name);
#endif
    return false;
  }

  if(!(var->get(data, counts))) {
#ifdef NCDF_VERBOSE
    output_info.write("INFO: NetCDF could not read data for '%s'\n", name);
#endif
    return false;
  }

  return true;
}

bool NcFormat::read(int *var, const string &name, int lx, int ly, int lz) {
  return read(var, name.c_str(), lx, ly, lz);
}

bool NcFormat::read(BoutReal *data, const char *name, int lx, int ly, int lz) {
  if(!is_valid())
    return false;

  if((lx < 0) || (ly < 0) || (lz < 0))
    return false;

  TRACE("NcFormat::read(BoutReal)");

  // Create an error object so netCDF doesn't exit
#ifdef NCDF_VERBOSE
  NcError err(NcError::verbose_nonfatal);
#else
  NcError err(NcError::silent_nonfatal);
#endif

  NcVar *var;
  
  if(!(var = dataFile->get_var(name))) {
    return false;
  }
  
  long cur[3], counts[3];
  cur[0] = x0;    cur[1] = y0;    cur[2] = z0;
  counts[0] = lx; counts[1] = ly; counts[2] = lz;

  if(!(var->set_cur(cur))) {
    return false;
  }

  if(!(var->get(data, counts))) {
    return false;
  }

  return true;
}

bool NcFormat::read(BoutReal *var, const string &name, int lx, int ly, int lz) {
  return read(var, name.c_str(), lx, ly, lz);
}

bool NcFormat::write(int *data, const char *name, int lx, int ly, int lz) {
  if(!is_valid())
    return false;

  if((lx < 0) || (ly < 0) || (lz < 0))
    return false;
  
  // Check for valid name
  checkName(name);

  int nd = 0; // Number of dimensions
  if(lx != 0) nd = 1;
  if(ly != 0) nd = 2;
  if(lz != 0) nd = 3;

  TRACE("NcFormat::write(int)");

#ifdef NCDF_VERBOSE
  NcError err(NcError::verbose_nonfatal);
#else
  NcError err(NcError::silent_nonfatal);
#endif
  
  NcVar *var;
  if(!(var = dataFile->get_var(name))) {
    // Variable not in file, so add it.
    
    var = dataFile->add_var(name, ncInt, nd, dimList);
    if(var == NULL) {
      output_error.write("ERROR: NetCDF could not add int '%s' to file '%s'\n", name, fname);
      return false;
    }
    if(!var->is_valid()) {
      output_error.write("ERROR: NetCDF could not add int '%s' to file '%s'\n", name, fname);
      return false;
    }
  }
  
  long cur[3], counts[3];
  cur[0] = x0;    cur[1] = y0;    cur[2] = z0;
  counts[0] = lx; counts[1] = ly; counts[2] = lz;

  if(!(var->set_cur(cur)))
    return false;
  
  if(!(var->put(data, counts)))
    return false;

  return true;
}

bool NcFormat::write(int *var, const string &name, int lx, int ly, int lz) {
  return write(var, name.c_str(), lx, ly, lz);
}

bool NcFormat::write(BoutReal *data, const char *name, int lx, int ly, int lz) {
  if(!is_valid())
    return false;

  if((lx < 0) || (ly < 0) || (lz < 0))
    return false;

  // Check for valid name
  checkName(name);
  
  TRACE("NcFormat::write(BoutReal)");

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
      output_error.write("ERROR: NetCDF could not add BoutReal '%s' to file '%s'\n", name, fname);
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
    int i_max=1;
    if (lx>0) i_max*=lx;
    if (ly>0) i_max*=ly;
    if (lz>0) i_max*=lz;
    for(int i=0;i<i_max;i++) {
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

  return true;
}

bool NcFormat::write(BoutReal *var, const string &name, int lx, int ly, int lz) {
  return write(var, name.c_str(), lx, ly, lz);
}

/***************************************************************************
 * Record-based (time-dependent) data
 ***************************************************************************/

bool NcFormat::read_rec(int *data, const char *name, int lx, int ly, int lz) {
  if(!is_valid())
    return false;

  if((lx < 0) || (ly < 0) || (lz < 0))
    return false;

  // Check for valid name
  checkName(name);

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

bool NcFormat::read_rec(int *var, const string &name, int lx, int ly, int lz) {
  return read_rec(var, name.c_str(), lx, ly, lz);
}

bool NcFormat::read_rec(BoutReal *data, const char *name, int lx, int ly, int lz) {
  if(!is_valid())
    return false;

  if((lx < 0) || (ly < 0) || (lz < 0))
    return false;
  
  // Check for valid name
  checkName(name);

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

bool NcFormat::read_rec(BoutReal *var, const string &name, int lx, int ly, int lz) {
  return read_rec(var, name.c_str(), lx, ly, lz);
}

bool NcFormat::write_rec(int *data, const char *name, int lx, int ly, int lz) {
  if(!is_valid())
    return false;

  if((lx < 0) || (ly < 0) || (lz < 0))
    return false;

  // Check for valid name
  checkName(name);

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
      output_error.write("ERROR: NetCDF Could not add variable '%s' to file '%s'\n", name, fname);
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

bool NcFormat::write_rec(int *var, const string &name, int lx, int ly, int lz) {
  return write_rec(var, name.c_str(), lx, ly, lz);
}

bool NcFormat::write_rec(BoutReal *data, const char *name, int lx, int ly, int lz) {
  if(!is_valid())
    return false;

  if((lx < 0) || (ly < 0) || (lz < 0))
    return false;

  // Check the name
  checkName(name);

  TRACE("NcFormat::write_rec(BoutReal*)");

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
    ASSERT1(var != 0);
    
    rec_nr[name] = default_rec; // Starting record
    
    if(!var->is_valid()) {
#ifdef NCDF_VERBOSE
      output_error.write("ERROR: NetCDF Could not add variable '%s' to file '%s'\n", name, fname);
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
  output_info.write("INFO: NetCDF writing record %d of '%s' in '%s'\n",t, name, fname); 
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
  int i_max=1;
  if (lx>0) i_max*=lx;
  if (ly>0) i_max*=ly;
  if (lz>0) i_max*=lz;
  for(int i=0;i<i_max;i++) {
    if(!finite(data[i]))
      data[i] = 0.0;
  }

  // Add the record
  if(!var->put_rec(data, t))
    return false;

  var->sync();
  
  // Increment record number
  rec_nr[name] = rec_nr[name] + 1;
  
  return true;
}

bool NcFormat::write_rec(BoutReal *var, const string &name, int lx, int ly, int lz) {
  return write_rec(var, name.c_str(), lx, ly, lz);
}

/***************************************************************************
 * Private functions
 ***************************************************************************/

void NcFormat::checkName(const char* name) {
  // Check if this name contains an invalid character
  
  const char* c = name;
  while(*c != 0) {
    if(*c == '*')
      throw BoutException("Invalid character (*) in NetCDF variable name '%s'", name);
    c++;
  }
}

#endif // NCDF

