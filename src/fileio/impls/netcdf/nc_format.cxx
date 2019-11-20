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

#include <bout/mesh.hxx>
#include <output.hxx>
#include <msg_stack.hxx>

using std::string;
using std::vector;

// Define this to see loads of info messages
//#define NCDF_VERBOSE

NcFormat::NcFormat(Mesh* mesh_in) : DataFormat(mesh_in) {
  dataFile = nullptr;
  x0 = y0 = z0 = t0 = 0;
  recDimList = new const NcDim*[4];
  dimList = recDimList+1;
  lowPrecision = false;

  default_rec = 0;
  rec_nr.clear();

  fname = nullptr;
}

NcFormat::NcFormat(const char *name, Mesh* mesh_in) : DataFormat(mesh_in) {
  dataFile = nullptr;
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

  if (dataFile != nullptr) // Already open. Close then re-open
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
    dataFile = nullptr;
    return false;
  }

  /// Get the dimensions from the file

  if(!(xDim = dataFile->get_dim("x"))) {
    output_warn.write("WARNING: NetCDF file should have an 'x' dimension\n");
    /*
    delete dataFile;
    dataFile = nullptr;
    return false;
    */
    xDim = nullptr;
  } else if (mesh != nullptr) {
    // Check that the dimension size is correct
    if (xDim->size() != mesh->LocalNx) {
      throw BoutException("X dimension incorrect. Expected {:d}, got {:d}", mesh->LocalNx,
                          xDim->size());
    }
  }
  
  if(!(yDim = dataFile->get_dim("y"))) {
    output_warn.write("WARNING: NetCDF file should have a 'y' dimension\n");
    /*
    delete dataFile;
    dataFile = nullptr;
    return false;
    */
    yDim = nullptr;
  } else if (mesh != nullptr) {
    // Check that the dimension size is correct
    if(yDim->size() != mesh->LocalNy) {
      throw BoutException("Y dimension incorrect. Expected {:d}, got {:d}", mesh->LocalNy,
                          yDim->size());
    }
  }
  
  if(!(zDim = dataFile->get_dim("z"))) {
    // Z dimension optional, and could be any size (Fourier harmonics)
#ifdef NCDF_VERBOSE
    output_info.write("INFO: NetCDF file has no 'z' coordinate\n");
#endif
    zDim = nullptr;
  } else if (mesh != nullptr) {
    // Check that the dimension size is correct
    if(zDim->size() != mesh->LocalNz) {
      throw BoutException("Z dimension incorrect. Expected {:d}, got {:d}", mesh->LocalNz,
                          zDim->size());
    }
  }
  
  if(!(tDim = dataFile->get_dim("t"))) {
    // T dimension optional
#ifdef NCDF_VERBOSE
    output_info.write("INFO: NetCDF file has no 't' coordinate\n");
#endif
    tDim = nullptr;
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

  if (dataFile != nullptr) // Already open. Close then re-open
    close(); 

  if(append) {
    dataFile = new NcFile(name, NcFile::Write);

    if(!dataFile->is_valid()) {
      delete dataFile;
      dataFile = nullptr;
      return false;
    }

    /// Get the dimensions from the file

    if(!(xDim = dataFile->get_dim("x"))) {
      output_error.write("ERROR: NetCDF file should have an 'x' dimension\n");
      delete dataFile;
      dataFile = nullptr;
      return false;
    }

    if(!(yDim = dataFile->get_dim("y"))) {
      output_error.write("ERROR: NetCDF file should have a 'y' dimension\n");
      delete dataFile;
      dataFile = nullptr;
      return false;
    }

    if(!(zDim = dataFile->get_dim("z"))) {
      output_error.write("ERROR: NetCDF file should have a 'z' dimension\n");
      delete dataFile;
      dataFile = nullptr;
      return false;
    }

    if(!(tDim = dataFile->get_dim("t"))) {
      output_error.write("ERROR: NetCDF file should have a 't' dimension\n");
      delete dataFile;
      dataFile = nullptr;
      return false;
    }

    /// Test they're the right size (and t is unlimited)
    
    if((xDim->size() != mesh->LocalNx) || (yDim->size() != mesh->LocalNy) || (zDim->size() != mesh->LocalNz)
       || (!tDim->is_unlimited()) ) {
      delete dataFile;
      dataFile = nullptr;
      return false;
    }

    // Get the size of the 't' dimension for records
    default_rec = tDim->size();
    
  }else {
    dataFile = new NcFile(name, NcFile::Replace);
    
    if(!dataFile->is_valid()) {
      delete dataFile;
      dataFile = nullptr;
      return false;
    }

    /// Add the dimensions
    
    if(!(xDim = dataFile->add_dim("x", mesh->LocalNx))) {
      delete dataFile;
      dataFile = nullptr;
      return false;
    }
  
    if(!(yDim = dataFile->add_dim("y", mesh->LocalNy))) {
      delete dataFile;
      dataFile = nullptr;
      return false;
    }
    
    if(!(zDim = dataFile->add_dim("z", mesh->LocalNz))) {
      delete dataFile;
      dataFile = nullptr;
      return false;
    }
    
    if(!(tDim = dataFile->add_dim("t"))) { // unlimited dimension
      delete dataFile;
      dataFile = nullptr;
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
  if (dataFile == nullptr)
    return false;
  return dataFile->is_valid();
}

void NcFormat::close() {
  if (dataFile == nullptr)
    return;
  
  TRACE("NcFormat::close");

#ifdef NCDF_VERBOSE
  NcError err(NcError::verbose_nonfatal);
#else
  NcError err(NcError::silent_nonfatal);
#endif

  dataFile->close();
  delete dataFile;
  dataFile = nullptr;

  free(fname);
  fname = nullptr;
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

// Add a variable to the file
bool NcFormat::addVarInt(const string &name, bool repeat) {
  if(!is_valid())
    return false;

  // Create an error object so netCDF doesn't exit
#ifdef NCDF_VERBOSE
  NcError err(NcError::verbose_nonfatal);
#else
  NcError err(NcError::silent_nonfatal);
#endif

  NcVar* var;
  if (!(var = dataFile->get_var(name.c_str()))) {
    // Variable not in file, so add it.
    if (repeat)
      var = dataFile->add_var(name.c_str(), ncInt, 1, recDimList);
    else
      var = dataFile->add_var(name.c_str(), ncInt, 0, dimList);

    if(!var->is_valid()) {
      output_error.write("ERROR: NetCDF could not add int '%s' to file '%s'\n", name.c_str(), fname);
      return false;
    }
  }
  return true;
}

bool NcFormat::addVarBoutReal(const string &name, bool repeat) {
  if(!is_valid())
    return false;

  // Create an error object so netCDF doesn't exit
#ifdef NCDF_VERBOSE
  NcError err(NcError::verbose_nonfatal);
#else
  NcError err(NcError::silent_nonfatal);
#endif

  NcVar* var;
  if (!(var = dataFile->get_var(name.c_str()))) {
    // Variable not in file, so add it.
    auto nc_float_type = lowPrecision ? ncFloat : ncDouble;
    if (repeat)
      var = dataFile->add_var(name.c_str(), nc_float_type, 1, recDimList);
    else
      var = dataFile->add_var(name.c_str(), nc_float_type, 0, dimList);

    if(!var->is_valid()) {
      output_error.write("ERROR: NetCDF could not add BoutReal '%s' to file '%s'\n", name.c_str(), fname);
      return false;
    }
  }
  return true;
}

bool NcFormat::addVarField2D(const string &name, bool repeat) {
  if(!is_valid())
    return false;

  // Create an error object so netCDF doesn't exit
#ifdef NCDF_VERBOSE
  NcError err(NcError::verbose_nonfatal);
#else
  NcError err(NcError::silent_nonfatal);
#endif

  NcVar* var;
  if (!(var = dataFile->get_var(name.c_str()))) {
    // Variable not in file, so add it.
    auto nc_float_type = lowPrecision ? ncFloat : ncDouble;
    if (repeat)
      var = dataFile->add_var(name.c_str(), nc_float_type, 3, recDimList);
    else
      var = dataFile->add_var(name.c_str(), nc_float_type, 2, dimList);

    if(!var->is_valid()) {
      output_error.write("ERROR: NetCDF could not add Field2D '%s' to file '%s'\n", name.c_str(), fname);
      return false;
    }
  }
  return true;
}

bool NcFormat::addVarField3D(const string &name, bool repeat) {
  if(!is_valid())
    return false;

  // Create an error object so netCDF doesn't exit
#ifdef NCDF_VERBOSE
  NcError err(NcError::verbose_nonfatal);
#else
  NcError err(NcError::silent_nonfatal);
#endif

  NcVar* var;
  if (!(var = dataFile->get_var(name.c_str()))) {
    // Variable not in file, so add it.
    auto nc_float_type = lowPrecision ? ncFloat : ncDouble;
    if (repeat)
      var = dataFile->add_var(name.c_str(), nc_float_type, 4, recDimList);
    else
      var = dataFile->add_var(name.c_str(), nc_float_type, 3, dimList);

    if(!var->is_valid()) {
      output_error.write("ERROR: NetCDF could not add Field3D '%s' to file '%s'\n", name.c_str(), fname);
      return false;
    }
  }
  return true;
}

bool NcFormat::addVarFieldPerp(const string &name, bool repeat) {
  if(!is_valid())
    return false;

  // Create an error object so netCDF doesn't exit
#ifdef NCDF_VERBOSE
  NcError err(NcError::verbose_nonfatal);
#else
  NcError err(NcError::silent_nonfatal);
#endif

  NcVar* var;
  if (!(var = dataFile->get_var(name.c_str()))) {
    // Variable not in file, so add it.
    auto nc_float_type = lowPrecision ? ncFloat : ncDouble;
    if (repeat){
      const NcDim * dims[3] = {tDim, xDim, zDim};
      var = dataFile->add_var(name.c_str(), nc_float_type, 3, dims);
    } else {
      const NcDim * dims[2] = {xDim, zDim};
      var = dataFile->add_var(name.c_str(), nc_float_type, 2, dims);
    }

    if(!var->is_valid()) {
      output_error.write("ERROR: NetCDF could not add FieldPerp '%s' to file '%s'\n", name.c_str(), fname);
      return false;
    }
  }
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

bool NcFormat::read_perp(BoutReal *data, const std::string& name, int lx, int lz) {
  if(!is_valid())
    return false;

  if((lx < 0) || (lz < 0))
    return false;

  TRACE("NcFormat::read_perp(BoutReal)");

  // Create an error object so netCDF doesn't exit
#ifdef NCDF_VERBOSE
  NcError err(NcError::verbose_nonfatal);
#else
  NcError err(NcError::silent_nonfatal);
#endif

  NcVar *var;

  if(!(var = dataFile->get_var(name.c_str()))) {
    return false;
  }

  long cur[2], counts[2];
  cur[0] = x0;    cur[1] = z0;
  counts[0] = lx; counts[1] = lz;

  if(!(var->set_cur(cur))) {
    return false;
  }

  if(!(var->get(data, counts))) {
    return false;
  }

  return true;
}

bool NcFormat::write(int *data, const char *name, int lx, int ly, int lz) {
  if(!is_valid())
    return false;

  if((lx < 0) || (ly < 0) || (lz < 0))
    return false;
  
  // Check for valid name
  checkName(name);

  TRACE("NcFormat::write(int)");

#ifdef NCDF_VERBOSE
  NcError err(NcError::verbose_nonfatal);
#else
  NcError err(NcError::silent_nonfatal);
#endif
  
  NcVar *var;
  if(!(var = dataFile->get_var(name))) {
    output_error.write("ERROR: NetCDF int variable '%s' has not been added to file '%s'\n", name, fname);
    return false;
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

#ifdef NCDF_VERBOSE
  NcError err(NcError::verbose_nonfatal);
#else
  NcError err(NcError::silent_nonfatal);
#endif

  NcVar *var;
  if(!(var = dataFile->get_var(name))) {
    output_error.write("ERROR: NetCDF BoutReal variable '%s' has not been added to file '%s'\n", name, fname);
    return false;
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

bool NcFormat::write_perp(BoutReal *data, const std::string& name, int lx, int lz) {
  if(!is_valid())
    return false;

  if((lx < 0) || (lz < 0))
    return false;

  // Check for valid name
  checkName(name.c_str());

  TRACE("NcFormat::write_perp(BoutReal)");

#ifdef NCDF_VERBOSE
  NcError err(NcError::verbose_nonfatal);
#else
  NcError err(NcError::silent_nonfatal);
#endif

  NcVar *var;
  if(!(var = dataFile->get_var(name.c_str()))) {
    output_error.write("ERROR: NetCDF BoutReal variable '%s' has not been added to file '%s'\n", name.c_str(), fname);
    return false;
  }

  long cur[2], counts[2];
  cur[0] = x0;    cur[1] = z0;
  counts[0] = lx; counts[1] = lz;

  if(!(var->set_cur(cur)))
    return false;

  if(lowPrecision) {
    // An out of range value can make the conversion
    // corrupt the whole dataset. Make sure everything
    // is in the range of a float
    int i_max=1;
    if (lx>0) i_max*=lx;
    if (lz>0) i_max*=lz;
    for(int i=0;i<i_max;i++) {
      if(data[i] > 1e20)
        data[i] = 1e20;
      if(data[i] < -1e20)
        data[i] = -1e20;
    }
  }

  for(int i=0;i<lx*lz;i++) {
    if(!finite(data[i]))
      data[i] = 0.0;
  }

  if(!(var->put(data, counts)))
    return false;

  return true;
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

bool NcFormat::read_rec_perp(BoutReal *data, const std::string& name, int lx, int lz) {
  if(!is_valid())
    return false;

  if((lx < 0) || (lz < 0))
    return false;

  // Check for valid name
  checkName(name.c_str());

  // Create an error object so netCDF doesn't exit
#ifdef NCDF_VERBOSE
  NcError err(NcError::verbose_nonfatal);
#else
  NcError err(NcError::silent_nonfatal);
#endif

  NcVar *var;

  if(!(var = dataFile->get_var(name.c_str())))
    return false;

  // NOTE: Probably should do something here to check t0

  long cur[3], counts[3];
  cur[0] = t0; cur[1] = x0; cur[2] = z0;
  counts[0] = 1; counts[1] = lx; counts[2] = lz;

  if(!(var->set_cur(cur)))
    return false;

  if(!(var->get(data, counts)))
    return false;

  return true;
}

bool NcFormat::write_rec(int *data, const char *name, int lx, int ly, int lz) {
  if(!is_valid())
    return false;

  if((lx < 0) || (ly < 0) || (lz < 0))
    return false;

  // Check for valid name
  checkName(name);

#ifdef NCDF_VERBOSE
  NcError err(NcError::verbose_nonfatal);
#else
  NcError err(NcError::silent_nonfatal);
#endif

  NcVar *var;
  
  // Try to find variable
  if(!(var = dataFile->get_var(name))) {
    output_error.write("ERROR: NetCDF int variable '%s' has not been added to file '%s'\n", name, fname);
    return false;
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

#ifdef NCDF_VERBOSE
  NcError err(NcError::verbose_nonfatal);
#else
  NcError err(NcError::silent_nonfatal);
#endif

  NcVar *var;

  // Try to find variable
  if(!(var = dataFile->get_var(name))) {
    output_error.write("ERROR: NetCDF BoutReal variable '%s' has not been added to file '%s'\n", name, fname);
    return false;
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

bool NcFormat::write_rec_perp(BoutReal *data, const std::string& name, int lx, int lz) {
  if(!is_valid())
    return false;

  if((lx < 0) || (lz < 0))
    return false;

  // Check the name
  checkName(name.c_str());

  TRACE("NcFormat::write_rec_perp(BoutReal*)");

#ifdef NCDF_VERBOSE
  NcError err(NcError::verbose_nonfatal);
#else
  NcError err(NcError::silent_nonfatal);
#endif

  NcVar *var;

  // Try to find variable
  if(!(var = dataFile->get_var(name.c_str()))) {
    output_error.write("ERROR: NetCDF BoutReal variable '%s' has not been added to file '%s'\n", name.c_str(), fname);
    return false;
  }else {
    // Get record number
    if(rec_nr.find(name.c_str()) == rec_nr.end()) {
      // Add to map
      rec_nr[name] = default_rec;
    }
  }

  int t = rec_nr[name];

#ifdef NCDF_VERBOSE
  output_info.write("INFO: NetCDF writing record %d of '%s' in '%s'\n",t, name.c_str(), fname);
#endif

  if(lowPrecision) {
    // An out of range value can make the conversion
    // corrupt the whole dataset. Make sure everything
    // is in the range of a float

    for(int i=0;i<lx*lz;i++) {
      if(data[i] > 1e20)
        data[i] = 1e20;
      if(data[i] < -1e20)
        data[i] = -1e20;
    }
  }
  int i_max=1;
  if (lx>0) i_max*=lx;
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

/***************************************************************************
 * Attributes
 ***************************************************************************/

void NcFormat::setAttribute(const std::string &varname, const std::string &attrname,
                            const std::string &text) {
  TRACE("NcFormat::setAttribute(string)");

#ifdef NCDF_VERBOSE
  NcError err(NcError::verbose_nonfatal);
#else
  NcError err(NcError::silent_nonfatal);
#endif

  std::string existing_att;
  if (getAttribute(varname, attrname, existing_att)) {
    if (text != existing_att) {
      output_warn.write("Overwriting attribute '%s' of variable '%s' with '%s', was previously '%s'",
          attrname.c_str(), varname.c_str(), text.c_str(), existing_att.c_str());
    }
  }
  // else: attribute does not exist, so just write it

  if (varname == "" ) {
    // file attribute
    dataFile->add_att(attrname.c_str(), text.c_str());
  } else {
    // variable attribute
    NcVar* var = dataFile->get_var(varname.c_str());
    if (var == nullptr or !var->is_valid()) {
      throw BoutException("Variable '{:s}' not in NetCDF file", varname);
    }

    var->add_att(attrname.c_str(), text.c_str());
  }
}

void NcFormat::setAttribute(const std::string &varname, const std::string &attrname,
                            int value) {
  TRACE("NcFormat::setAttribute(int)");

#ifdef NCDF_VERBOSE
  NcError err(NcError::verbose_nonfatal);
#else
  NcError err(NcError::silent_nonfatal);
#endif

  int existing_att;
  if (getAttribute(varname, attrname, existing_att)) {
    if (value != existing_att) {
      output_warn.write("Overwriting attribute '%s' of variable '%s' with '%i', was previously '%i'",
          attrname.c_str(), varname.c_str(), value, existing_att);
    }
  }
  // else: attribute does not exist, so just write it

  if (varname == "") {
    // attribute of file
    dataFile->add_att(attrname.c_str(), value);
  } else {
    // attribute of variable
    NcVar* var = dataFile->get_var(varname.c_str());
    if (var == nullptr or !var->is_valid()) {
      throw BoutException("Variable '{:s}' not in NetCDF file", varname);
    }

    var->add_att(attrname.c_str(), value);
  }
}

void NcFormat::setAttribute(const std::string &varname, const std::string &attrname,
                            BoutReal value) {
  TRACE("NcFormat::setAttribute(BoutReal)");

#ifdef NCDF_VERBOSE
  NcError err(NcError::verbose_nonfatal);
#else
  NcError err(NcError::silent_nonfatal);
#endif

  int existing_att;
  if (getAttribute(varname, attrname, existing_att)) {
    if (value != existing_att) {
      output_warn.write("Overwriting attribute '%s' of variable '%s' with '%f', was previously '%d'",
			attrname.c_str(), varname.c_str(), value, existing_att);
    }
  }
  // else: attribute does not exist, so just write it

  if (varname == "") {
    // attribute of file
    dataFile->add_att(attrname.c_str(), value);
  } else {
    // attribute of variable
    NcVar* var = dataFile->get_var(varname.c_str());
    if (var == nullptr or !var->is_valid()) {
      throw BoutException("Variable '{:s}' not in NetCDF file", varname);
    }

    var->add_att(attrname.c_str(), value);
  }
}

bool NcFormat::getAttribute(const std::string &varname, const std::string &attrname, std::string &text) {
  TRACE("NcFormat::getAttribute(string)");

#ifdef NCDF_VERBOSE
  NcError err(NcError::verbose_nonfatal);
#else
  NcError err(NcError::silent_nonfatal);
#endif

  if (varname == "") {
    // attribute of file
    NcAtt* fileAtt;
    if (!(fileAtt = dataFile->get_att(attrname.c_str()))) {
      return false;
    }

    auto values = fileAtt->values();
    if (values == nullptr)
      return false;

    text = values->as_string(0);

    return true;
  } else {
    NcVar* var = dataFile->get_var(varname.c_str());
    if (var == nullptr or !var->is_valid()) {
      throw BoutException("Variable '{:s}' not in NetCDF file", varname);
    }

    NcAtt* varAtt;
    if (!(varAtt = var->get_att(attrname.c_str()))) {
      return false;
    }

    auto values = varAtt->values();
    if (values == nullptr)
      return false;

    text = values->as_string(0);

    return true;
  }
}

bool NcFormat::getAttribute(const std::string &varname, const std::string &attrname, int &value) {
  TRACE("NcFormat::getAttribute(int)");

#ifdef NCDF_VERBOSE
  NcError err(NcError::verbose_nonfatal);
#else
  NcError err(NcError::silent_nonfatal);
#endif

  if (varname == "") {
    // attribute of file
    NcAtt* fileAtt;
    if (!(fileAtt = dataFile->get_att(attrname.c_str()))) {
      return false;
    }

    auto values = fileAtt->values();
    if (values == nullptr)
      return false;

    value = values->as_int(0);

    return true;
  } else {
    // attribute of variable
    NcVar* var;
    if (!(var = dataFile->get_var(varname.c_str()))) {
      throw BoutException("Variable '{:s}' not in NetCDF file", varname);
    }

    NcAtt* varAtt;
    if (!(varAtt = var->get_att(attrname.c_str())))
      return false;

    auto values = varAtt->values();
    if (values == nullptr)
      return false;

    value = values->as_int(0);

    return true;
  }
}

bool NcFormat::getAttribute(const std::string &varname, const std::string &attrname, BoutReal &value) {
  TRACE("NcFormat::getAttribute(BoutReal)");

#ifdef NCDF_VERBOSE
  NcError err(NcError::verbose_nonfatal);
#else
  NcError err(NcError::silent_nonfatal);
#endif

  if (varname == "") {
    // attribute of file
    NcAtt* fileAtt;
    if (!(fileAtt = dataFile->get_att(attrname.c_str()))) {
      return false;
    }

    auto values = fileAtt->values();
    if (values == nullptr)
      return false;

    value = values->as_double(0);

    return true;
  } else {
    // attribute of variable
    NcVar* var;
    if (!(var = dataFile->get_var(varname.c_str()))) {
      throw BoutException("Variable '{:s}' not in NetCDF file", varname);
    }

    NcAtt* varAtt;
    if (!(varAtt = var->get_att(attrname.c_str())))
      return false;

    auto values = varAtt->values();
    if (values == nullptr)
      return false;

    value = values->as_double(0);

    return true;
  }
}


/***************************************************************************
 * Private functions
 ***************************************************************************/

void NcFormat::checkName(const char* name) {
  // Check if this name contains an invalid character
  
  const char* c = name;
  while(*c != 0) {
    if(*c == '*')
      throw BoutException("Invalid character (*) in NetCDF variable name '{:s}'", name);
    c++;
  }
}

#endif // NCDF

