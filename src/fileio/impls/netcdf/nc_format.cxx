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

NcFormat::NcFormat()
{
  dataFile = NULL;
  x0 = y0 = z0 = t0 = 0;
  recDimList = new const NcDim*[4];
  dimList = recDimList+1;
  lowPrecision = false;

  default_rec = 0;
  rec_nr.clear();

  fname = NULL;
}

NcFormat::NcFormat(const char *name)
{
  dataFile = NULL;
  x0 = y0 = z0 = t0 = 0;
  recDimList = new const NcDim*[4];
  dimList = recDimList+1;
  lowPrecision = false;

  default_rec = 0;
  rec_nr.clear();

  openr(name);
}

NcFormat::NcFormat(const string &name)
{
  dataFile = NULL;
  x0 = y0 = z0 = t0 = 0;
  recDimList = new const NcDim*[4];
  dimList = recDimList+1;
  lowPrecision = false;

  default_rec = 0;
  rec_nr.clear();

  openr(name);
}

NcFormat::~NcFormat()
{
  delete[] recDimList;
  close();
  rec_nr.clear();
}

bool NcFormat::openr(const string &name)
{
  return openr(name.c_str());
}

bool NcFormat::openr(const char *name)
{
#ifdef CHECK
  msg_stack.push("NcFormat::openr");
#endif

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
    output.write("WARNING: NetCDF file should have an 'x' dimension\n");
    /*
    delete dataFile;
    dataFile = NULL;
    return false;
    */
    xDim = NULL;
  }
  
  if(!(yDim = dataFile->get_dim("y"))) {
    output.write("WARNING: NetCDF file should have a 'y' dimension\n");
    /*
    delete dataFile;
    dataFile = NULL;
    return false;
    */
    yDim = NULL;
  }
  
  if(!(zDim = dataFile->get_dim("z"))) {
    // Z dimension optional, and could be any size (Fourier harmonics)
#ifdef NCDF_VERBOSE
    output.write("INFO: NetCDF file has no 'z' coordinate\n");
#endif
    zDim = NULL;
  }
  
  if(!(tDim = dataFile->get_dim("t"))) {
    // T dimension optional
#ifdef NCDF_VERBOSE
    output.write("INFO: NetCDF file has no 't' coordinate\n");
#endif
    tDim = NULL;
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

bool NcFormat::openw(const string &name, bool append)
{
  return openw(name.c_str(), append);
}

bool NcFormat::openw(const char *name, bool append)
{
#ifdef CHECK
  msg_stack.push("NcFormat::openw");
#endif

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
      output.write("ERROR: NetCDF file should have an 'x' dimension\n");
      delete dataFile;
      dataFile = NULL;
      return false;
    }

    if(!(yDim = dataFile->get_dim("y"))) {
      output.write("ERROR: NetCDF file should have a 'y' dimension\n");
      delete dataFile;
      dataFile = NULL;
      return false;
    }

    if(!(zDim = dataFile->get_dim("z"))) {
      output.write("ERROR: NetCDF file should have a 'z' dimension\n");
      delete dataFile;
      dataFile = NULL;
      return false;
    }

    if(!(tDim = dataFile->get_dim("t"))) {
      output.write("ERROR: NetCDF file should have a 't' dimension\n");
      delete dataFile;
      dataFile = NULL;
      return false;
    }

    /// Test they're the right size (and t is unlimited)
    
    if((xDim->size() != mesh->ngx) || (yDim->size() != mesh->ngy) || (zDim->size() != mesh->ngz)
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
    
    if(!(xDim = dataFile->add_dim("x", mesh->ngx))) {
      delete dataFile;
      dataFile = NULL;
      return false;
    }
  
    if(!(yDim = dataFile->add_dim("y", mesh->ngy))) {
      delete dataFile;
      dataFile = NULL;
      return false;
    }
    
    if(!(zDim = dataFile->add_dim("z", mesh->ngz))) {
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

#ifdef CHECK
  msg_stack.pop();
#endif

  return true;
}

bool NcFormat::is_valid()
{
  if(dataFile == NULL)
    return false;
  return dataFile->is_valid();
}

void NcFormat::close()
{
  if(dataFile == NULL)
    return;

#ifdef CHECK
  msg_stack.push("NcFormat::close");
#endif

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

#ifdef CHECK
  msg_stack.pop();
#endif
}

void NcFormat::flush() {
  if(!is_valid())
    return;
  
  dataFile->sync();
}

const vector<int> NcFormat::getSize(const char *name)
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
  msg_stack.push("NcFormat::getSize");
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

const vector<int> NcFormat::getSize(const string &var)
{
  return getSize(var.c_str());
}

bool NcFormat::setGlobalOrigin(int x, int y, int z)
{
  x0 = x;
  y0 = y;
  z0 = z;
  
  return true;
}

bool NcFormat::setRecord(int t)
{
  t0 = t;

  return true;
}

bool NcFormat::read(int *data, const char *name, int lx, int ly, int lz)
{
  if(!is_valid())
    return false;

  if((lx < 0) || (ly < 0) || (lz < 0))
    return false;
  
#ifdef CHECK
  msg_stack.push("NcFormat::read(int)");
#endif

  // Create an error object so netCDF doesn't exit
#ifdef NCDF_VERBOSE
  NcError err(NcError::verbose_nonfatal);
#else
  NcError err(NcError::silent_nonfatal);
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

bool NcFormat::read(int *var, const string &name, int lx, int ly, int lz)
{
  return read(var, name.c_str(), lx, ly, lz);
}

bool NcFormat::read(BoutReal *data, const char *name, int lx, int ly, int lz)
{
  if(!is_valid())
    return false;

  if((lx < 0) || (ly < 0) || (lz < 0))
    return false;

#ifdef CHECK
  msg_stack.push("NcFormat::read(BoutReal)");
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

bool NcFormat::read(BoutReal *var, const string &name, int lx, int ly, int lz)
{
  return read(var, name.c_str(), lx, ly, lz);
}

bool NcFormat::write(int *data, const char *name, int lx, int ly, int lz)
{
  if(!is_valid())
    return false;

  if((lx < 0) || (ly < 0) || (lz < 0))
    return false;
  
  int nd = 0; // Number of dimensions
  if(lx != 0) nd = 1;
  if(ly != 0) nd = 2;
  if(lz != 0) nd = 3;

#ifdef CHECK
  msg_stack.push("NcFormat::write(int)");
#endif

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
      output.write("ERROR: NetCDF could not add int '%s' to file '%s'\n", name, fname);
      return false;
    }
    if(!var->is_valid()) {
      output.write("ERROR: NetCDF could not add int '%s' to file '%s'\n", name, fname);
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

#ifdef CHECK
  msg_stack.pop();
#endif

  return true;
}

bool NcFormat::write(int *var, const string &name, int lx, int ly, int lz)
{
  return write(var, name.c_str(), lx, ly, lz);
}

bool NcFormat::write(BoutReal *data, const char *name, int lx, int ly, int lz)
{
  if(!is_valid())
    return false;

  if((lx < 0) || (ly < 0) || (lz < 0))
    return false;
  
#ifdef CHECK
  msg_stack.push("NcFormat::write(BoutReal)");
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

bool NcFormat::write(BoutReal *var, const string &name, int lx, int ly, int lz)
{
  return write(var, name.c_str(), lx, ly, lz);
}

/***************************************************************************
 * Record-based (time-dependent) data
 ***************************************************************************/

bool NcFormat::read_rec(int *data, const char *name, int lx, int ly, int lz)
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

bool NcFormat::read_rec(int *var, const string &name, int lx, int ly, int lz)
{
  return read_rec(var, name.c_str(), lx, ly, lz);
}

bool NcFormat::read_rec(BoutReal *data, const char *name, int lx, int ly, int lz)
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

bool NcFormat::read_rec(BoutReal *var, const string &name, int lx, int ly, int lz)
{
  return read_rec(var, name.c_str(), lx, ly, lz);
}

bool NcFormat::write_rec(int *data, const char *name, int lx, int ly, int lz)
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

bool NcFormat::write_rec(int *var, const string &name, int lx, int ly, int lz)
{
  return write_rec(var, name.c_str(), lx, ly, lz);
}

bool NcFormat::write_rec(BoutReal *data, const char *name, int lx, int ly, int lz)
{
  if(!is_valid())
    return false;

  if((lx < 0) || (ly < 0) || (lz < 0))
    return false;

#ifdef CHECK
  msg_stack.push("NcFormat::write_rec(BoutReal)");
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

bool NcFormat::write_rec(BoutReal *var, const string &name, int lx, int ly, int lz)
{
  return write_rec(var, name.c_str(), lx, ly, lz);
}

/***************************************************************************
 * Private functions
 ***************************************************************************/

#endif // NCDF

