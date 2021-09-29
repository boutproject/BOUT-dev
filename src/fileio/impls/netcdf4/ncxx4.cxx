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

#include "bout/build_config.hxx"

#include "ncxx4.hxx"

#if BOUT_HAS_NETCDF && !BOUT_HAS_LEGACY_NETCDF

#include <bout/mesh.hxx>
#include <globals.hxx>
#include <utils.hxx>
#include <cmath>

#include <output.hxx>
#include <msg_stack.hxx>

using std::string;

using namespace netCDF;

// Define this to see loads of info messages
//#define NCDF_VERBOSE

Ncxx4::Ncxx4(Mesh* mesh_in) : DataFormat(mesh_in) {
  dataFile = nullptr;
  x0 = y0 = z0 = t0 = 0;
  recDimList = new const NcDim*[4];
  dimList = recDimList+1;
  lowPrecision = false;

  default_rec = 0;
  rec_nr.clear();

  fname = nullptr;
}

Ncxx4::Ncxx4(const char *name, Mesh* mesh_in) : DataFormat(mesh_in) {
  dataFile = nullptr;
  x0 = y0 = z0 = t0 = 0;
  recDimList = new const NcDim*[4];
  dimList = recDimList+1;
  lowPrecision = false;

  default_rec = 0;
  rec_nr.clear();

  Ncxx4::openr(name);
}

Ncxx4::~Ncxx4() {
  delete[] recDimList;
  Ncxx4::close();
  rec_nr.clear();
}

bool Ncxx4::openr(const char *name) {
  TRACE("Ncxx4::openr");

#ifdef NCDF_VERBOSE
  output.write("Ncxx4:: openr({:s})\n", name);
#endif

  if (dataFile != nullptr) {
    // Already open. Close then re-open
    Ncxx4::close();
  }

  try {
    dataFile = new NcFile(name, NcFile::read);
  } catch (netCDF::exceptions::NcException&) {
    return false;
  }

  if(dataFile->isNull()) {
    delete dataFile;
    dataFile = nullptr;
    return false;
  }

  /// Get the dimensions from the file

  xDim = dataFile->getDim("x");
  if(xDim.isNull())
    output_warn.write("WARNING: NetCDF file should have an 'x' dimension\n");

  yDim = dataFile->getDim("y");
  if(yDim.isNull())
    output_warn.write("WARNING: NetCDF file should have a 'y' dimension\n");

  zDim = dataFile->getDim("z");
  if(zDim.isNull()) {
    // Z dimension optional, and could be any size (Fourier harmonics)
#ifdef NCDF_VERBOSE
    output_info.write("INFO: NetCDF file has no 'z' coordinate\n");
#endif
  }

  tDim = dataFile->getDim("t");
  if(tDim.isNull()) {
    // T dimension optional
#ifdef NCDF_VERBOSE
    output_info.write("INFO: NetCDF file has no 't' coordinate\n");
#endif
  }

  recDimList[0] = &tDim;
  recDimList[1] = &xDim;
  recDimList[2] = &yDim;
  recDimList[3] = &zDim;

  fname = copy_string(name);

  return true;
}

bool Ncxx4::openw(const char *name, bool append) {
  TRACE("Ncxx4::openw");
  
#ifdef NCDF_VERBOSE
  output.write("Ncxx4:: openw({:s}, {})\n", name, append);
#endif

  if (dataFile != nullptr) // Already open. Close then re-open
    close(); 

  if(append) {
    try {
      dataFile = new NcFile(name, NcFile::write);
    } catch (netCDF::exceptions::NcException&) {
      return false;
    }

    if(dataFile->isNull()) {
      delete dataFile;
      dataFile = nullptr;
      return false;
    }

    /// Get the dimensions from the file

    xDim = dataFile->getDim("x");
    if(xDim.isNull()) {
      output_error.write("ERROR: NetCDF file should have an 'x' dimension\n");
      delete dataFile;
      dataFile = nullptr;
      return false;
    }

    yDim = dataFile->getDim("y");
    if(yDim.isNull()) {
      output_error.write("ERROR: NetCDF file should have a 'y' dimension\n");
      delete dataFile;
      dataFile = nullptr;
      return false;
    }

    zDim = dataFile->getDim("z");
    if(zDim.isNull()) {
      output_error.write("ERROR: NetCDF file should have a 'z' dimension\n");
      delete dataFile;
      dataFile = nullptr;
      return false;
    }

    tDim = dataFile->getDim("t");
    if(tDim.isNull()) {
      output_error.write("ERROR: NetCDF file should have a 't' dimension\n");
      delete dataFile;
      dataFile = nullptr;
      return false;
    }

    /// Test they're the right size (and t is unlimited)
    if ((xDim.getSize() != static_cast<size_t>(mesh->LocalNx)) ||
        (yDim.getSize() != static_cast<size_t>(mesh->LocalNy)) ||
        (zDim.getSize() != static_cast<size_t>(mesh->LocalNz)) || (!tDim.isUnlimited())) {
      delete dataFile;
      dataFile = nullptr;
      return false;
    }

    // Get the size of the 't' dimension for records
    default_rec = tDim.getSize();
    
  }else {
    try {
      dataFile = new NcFile(name, NcFile::replace);
    } catch (netCDF::exceptions::NcException&) {
      return false;
    }
    
    if(dataFile->isNull()) {
      delete dataFile;
      dataFile = nullptr;
      return false;
    }

    /// Add the dimensions
    
    xDim = dataFile->addDim("x", mesh->LocalNx);
    if(xDim.isNull()) {
      delete dataFile;
      dataFile = nullptr;
      return false;
    }
  
    yDim = dataFile->addDim("y", mesh->LocalNy);
    if(yDim.isNull()) {
      delete dataFile;
      dataFile = nullptr;
      return false;
    }
    
    zDim = dataFile->addDim("z", mesh->LocalNz);
    if(zDim.isNull()) {
      delete dataFile;
      dataFile = nullptr;
      return false;
    }
    
    tDim = dataFile->addDim("t");
    if(tDim.isNull()) { // unlimited dimension
      delete dataFile;
      dataFile = nullptr;
      return false;
    }

    default_rec = 0; // Starting at record 0
  }

  recDimList[0] = &tDim;
  recDimList[1] = &xDim;
  recDimList[2] = &yDim;
  recDimList[3] = &zDim;

  fname = copy_string(name);

  return true;
}

bool Ncxx4::is_valid() {
  if (dataFile == nullptr)
    return false;
  return !dataFile->isNull();
}

void Ncxx4::close() {
  TRACE("Ncxx4::close");

#ifdef NCDF_VERBOSE
  output.write("Ncxx4:: close()\n"); 
#endif

  if (dataFile == nullptr)
    return;
  
  delete dataFile;
  dataFile = nullptr;

  free(fname);
  fname = nullptr;
}

void Ncxx4::flush() {
}

const std::vector<int> Ncxx4::getSize(const char *name) {
  TRACE("Ncxx4::getSize");

#ifdef NCDF_VERBOSE
  output.write("Ncxx4:: getSize({:s})\n", name); 
#endif
  std::vector<int> size;

  if(!is_valid())
    return size;

  NcVar var;
  
  var = dataFile->getVar(name);
  if(var.isNull())
    return size;
  
  if(var.getDimCount() == 0) {
    size.push_back(1);
    return size;
  }
  
  for(const auto& dim: var.getDims()) {
    size.push_back(dim.getSize());
  }

  return size;
}

const std::vector<int> Ncxx4::getSize(const std::string &var) {
  return getSize(var.c_str());
}

bool Ncxx4::setGlobalOrigin(int x, int y, int z) {
  x0 = x;
  y0 = y;
  z0 = z;
  
  return true;
}

bool Ncxx4::setRecord(int t) {
  t0 = t;

  return true;
}

// Add a variable to the file
bool Ncxx4::addVarInt(const string &name, bool repeat) {
  if(!is_valid())
    return false;

  NcVar var = dataFile->getVar(name);
  if(var.isNull()) {
    // Variable not in file, so add it.
    if (repeat)
      var = dataFile->addVar(name, ncInt, getRecDimVec(1));
    else
      var = dataFile->addVar(name, ncInt, getDimVec(0));

    if(var.isNull()) {
      output_error.write("ERROR: NetCDF could not add int '{:s}' to file '{:s}'\n", name, fname);
      return false;
    }
  }
  return true;
}

bool Ncxx4::addVarIntVec(const string &name, bool repeat, size_t size) {
  if(!is_valid())
    return false;

  NcVar var = dataFile->getVar(name);
  const auto vec_dim_name = "vec" + std::to_string(size);
  auto vec_dim = dataFile->getDim(vec_dim_name);
  if (vec_dim.isNull()) {
    vec_dim = dataFile->addDim(vec_dim_name, size);
  }
  if(var.isNull()) {
    // Variable not in file, so add it.
    if (repeat) {
      auto dims = getRecDimVec(2);
      dims[1] = vec_dim;
      var = dataFile->addVar(name, ncInt, dims);
    } else {
      auto dims = getDimVec(1);
      dims[0] = vec_dim;
      var = dataFile->addVar(name, ncInt, dims);
    }

    if(var.isNull()) {
      output_error.write("ERROR: NetCDF could not add int '{:s}' to file '{:s}'\n", name, fname);
      return false;
    }
  } else {
    // Check the existing variable is consistent with what's being added
    if (repeat) {
      ASSERT0(var.getDimCount() == 2);
      if (var.getDim(1).getSize() != size) {
        throw BoutException("Found existing variable '{:s}' with size {}. Trying to add "
                            "with size {}.", name, var.getDim(1).getSize(), size);
      }
    } else {
      ASSERT0(var.getDimCount() == 1);
      if (var.getDim(0).getSize() != size) {
        throw BoutException("Found existing variable '{:s}' with size {}. Trying to add "
                            "with size {}.", name, var.getDim(0).getSize(), size);
      }
    }
  }
  return true;
}

bool Ncxx4::addVarString(const string &name, bool repeat, size_t size) {
  if (!is_valid()) {
    return false;
  }

  NcVar var = dataFile->getVar(name);
  const auto vec_dim_name = "char" + std::to_string(size);
  auto vec_dim = dataFile->getDim(vec_dim_name);
  if (vec_dim.isNull()) {
    vec_dim = dataFile->addDim(vec_dim_name, size);
  }
  if(var.isNull()) {
    // Variable not in file, so add it.
    if (repeat) {
      auto dims = getRecDimVec(2);
      dims[1] = vec_dim;
      var = dataFile->addVar(name, ncChar, dims);
    } else {
      auto dims = getDimVec(1);
      dims[0] = vec_dim;
      var = dataFile->addVar(name, ncChar, dims);
    }

    if(var.isNull()) {
      output_error.write("ERROR: NetCDF could not add char '{:s}' to file '{:s}'\n", name, fname);
      return false;
    }
  } else {
    // Check the existing variable is consistent with what's being added
    if (repeat) {
      ASSERT0(var.getDimCount() == 2);
      if (var.getDim(1).getSize() != size) {
        throw BoutException("Found existing variable '{:s}' with size {}. Trying to add "
                            "with size {}.", name, var.getDim(1).getSize(), size);
      }
    } else {
      ASSERT0(var.getDimCount() == 1);
      if (var.getDim(0).getSize() != size) {
        throw BoutException("Found existing variable '{:s}' with size {}. Trying to add "
                            "with size {}.", name, var.getDim(0).getSize(), size);
      }
    }
  }
  return true;
}

bool Ncxx4::addVarBoutReal(const string &name, bool repeat) {
  if(!is_valid())
    return false;

  NcVar var = dataFile->getVar(name);
  if(var.isNull()) {
    // Variable not in file, so add it.
    if(lowPrecision) {
      if (repeat)
        var = dataFile->addVar(name, ncFloat, getRecDimVec(1));
      else
        var = dataFile->addVar(name, ncFloat, getDimVec(0));
    } else {
      if (repeat)
        var = dataFile->addVar(name, ncDouble, getRecDimVec(1));
      else
        var = dataFile->addVar(name, ncDouble, getDimVec(0));
    }

    if(var.isNull()) {
      output_error.write("ERROR: NetCDF could not add BoutReal '{:s}' to file '{:s}'\n", name, fname);
      return false;
    }
  }
  return true;
}

bool Ncxx4::addVarField2D(const string &name, bool repeat) {
  if(!is_valid())
    return false;

  NcVar var = dataFile->getVar(name);
  if(var.isNull()) {
    // Variable not in file, so add it.
    if(lowPrecision) {
      if (repeat)
        var = dataFile->addVar(name, ncFloat, getRecDimVec(3));
      else
        var = dataFile->addVar(name, ncFloat, getDimVec(2));
    } else {
      if (repeat)
        var = dataFile->addVar(name, ncDouble, getRecDimVec(3));
      else
        var = dataFile->addVar(name, ncDouble, getDimVec(2));
    }

    if(var.isNull()) {
      output_error.write("ERROR: NetCDF could not add Field2D '{:s}' to file '{:s}'\n", name, fname);
      return false;
    }
  }
  return true;
}

bool Ncxx4::addVarField3D(const string &name, bool repeat) {
  if(!is_valid())
    return false;

  NcVar var = dataFile->getVar(name);
  if(var.isNull()) {
    // Variable not in file, so add it.
    if(lowPrecision) {
      if (repeat)
        var = dataFile->addVar(name, ncFloat, getRecDimVec(4));
      else
        var = dataFile->addVar(name, ncFloat, getDimVec(3));
    } else {
      if (repeat)
        var = dataFile->addVar(name, ncDouble, getRecDimVec(4));
      else
        var = dataFile->addVar(name, ncDouble, getDimVec(3));
    }

    if(var.isNull()) {
      output_error.write("ERROR: NetCDF could not add Field3D '{:s}' to file '{:s}'\n", name, fname);
      return false;
    }
  }
  return true;
}

bool Ncxx4::addVarFieldPerp(const string &name, bool repeat) {
  if(!is_valid())
    return false;

  NcVar var = dataFile->getVar(name);
  if(var.isNull()) {
    // Variable not in file, so add it.
    if (repeat) {
      std::vector<NcDim> vec = {*recDimList[0], *recDimList[1], *recDimList[3]};
      if(lowPrecision) {
        var = dataFile->addVar(name, ncFloat, vec);
      } else {
        var = dataFile->addVar(name, ncDouble, vec);
      }
    } else {
      std::vector<NcDim> vec = {*dimList[0], *dimList[2]};
      if(lowPrecision) {
        var = dataFile->addVar(name, ncFloat, vec);
      } else {
        var = dataFile->addVar(name, ncDouble, vec);
      }
    }

    if(var.isNull()) {
      output_error.write("ERROR: NetCDF could not add FieldPerp '{:s}' to file '{:s}'\n", name, fname);
      return false;
    }
  }
  return true;
}

bool Ncxx4::read(int *data, const char *name, int lx, int ly, int lz) {
  TRACE("Ncxx4::read(int)");

#ifdef NCDF_VERBOSE
  output.write("Ncxx4:: read(int, {:s})\n", name);
#endif

  if(!is_valid())
    return false;

  if((lx < 0) || (ly < 0) || (lz < 0))
    return false;
    
  NcVar var = dataFile->getVar(name);
  if(var.isNull()) {
#ifdef NCDF_VERBOSE
    output_info.write("INFO: NetCDF variable '{:s}' not found\n", name);
#endif
    return false;
  }
  
  std::vector<size_t> start(3);
  start[0] = x0; start[1] = y0; start[2] = z0;
  std::vector<size_t> counts(3);
  counts[0] = lx; counts[1] = ly; counts[2] = lz;
  
  var.getVar(start, counts, data);

  return true;
}

bool Ncxx4::read(int *var, const std::string &name, int lx, int ly, int lz) {
  return read(var, name.c_str(), lx, ly, lz);
}

bool Ncxx4::read(char *data, const char *name, int n) {
  TRACE("Ncxx4::read(char)");

#ifdef NCDF_VERBOSE
  output.write("Ncxx4:: read(char, {:s})\n", name);
#endif

  if (!is_valid()) {
    return false;
  }

  if (n < 0) {
    return false;
  }

  NcVar var = dataFile->getVar(name);
  if (var.isNull()) {
#ifdef NCDF_VERBOSE
    output_info.write("INFO: NetCDF variable '{:s}' not found\n", name);
#endif
    return false;
  }

  std::vector<size_t> start = {0};
  std::vector<size_t> counts = {size_t(n)};

  var.getVar(start, counts, data);

  return true;
}

bool Ncxx4::read(char *var, const std::string &name, int n) {
  return read(var, name.c_str(), n);
}

bool Ncxx4::read(BoutReal *data, const char *name, int lx, int ly, int lz) {
  TRACE("Ncxx4::read(BoutReal)");

#ifdef NCDF_VERBOSE
  output.write("Ncxx4:: read(BoutReal, {:s})\n", name); 
#endif
  if(!is_valid())
    return false;

  if((lx < 0) || (ly < 0) || (lz < 0))
    return false;
  
  NcVar var = dataFile->getVar(name);
  
  if(var.isNull()) {
    return false;
  }

  std::vector<size_t> start(3);
  start[0] = x0; start[1] = y0; start[2] = z0;
  std::vector<size_t> counts(3);
  counts[0] = lx; counts[1] = ly; counts[2] = lz;
  
  var.getVar(start, counts, data);
  
  return true;
}

bool Ncxx4::read(BoutReal *var, const std::string &name, int lx, int ly, int lz) {
  return read(var, name.c_str(), lx, ly, lz);
}

bool Ncxx4::read_perp(BoutReal *data, const std::string& name, int lx, int lz) {
  TRACE("Ncxx4::read_perp(BoutReal)");

#ifdef NCDF_VERBOSE
  output.write("Ncxx4:: read(BoutReal, {:s})\n", name);
#endif
  if(!is_valid())
    return false;

  if((lx < 0) || (lz < 0))
    return false;

  NcVar var = dataFile->getVar(name);

  if(var.isNull()) {
    return false;
  }

  std::vector<size_t> start(2);
  start[0] = x0; start[1] = z0;
  std::vector<size_t> counts(2);
  counts[0] = lx; counts[1] = lz;

  var.getVar(start, counts, data);

  return true;
}

bool Ncxx4::write(int *data, const char *name, int lx, int ly, int lz) {
  TRACE("Ncxx4::write(int)");

#ifdef NCDF_VERBOSE
  output.write("Ncxx4:: write(int, {:s})\n", name);
#endif
  if(!is_valid())
    return false;

  if((lx < 0) || (ly < 0) || (lz < 0))
    return false;

  NcVar var = dataFile->getVar(name);
  if(var.isNull()) {
    output_error.write("ERROR: NetCDF int variable '{:s}' has not been added to file '{:s}'\n", name, fname);
    return false;
  }

#ifdef NCDF_VERBOSE
  output.write("Ncxx4:: write { Writing Variable } \n");
#endif

  std::vector<size_t> start(3);
  start[0] = x0; start[1] = y0; start[2] = z0;
  std::vector<size_t> counts(3);
  counts[0] = lx; counts[1] = ly; counts[2] = lz;

  var.putVar(start, counts, data);

#ifdef NCDF_VERBOSE
  output.write("Ncxx4:: write { Done } \n");
#endif

  return true;
}

bool Ncxx4::write(int *var, const std::string &name, int lx, int ly, int lz) {
  return write(var, name.c_str(), lx, ly, lz);
}

bool Ncxx4::write(char *data, const char *name, int n) {
  TRACE("Ncxx4::write(char)");

#ifdef NCDF_VERBOSE
  output.write("Ncxx4:: write(char, {:s})\n", name);
#endif
  if (!is_valid()) {
    return false;
  }

  if (n < 0) {
    return false;
  }

  NcVar var = dataFile->getVar(name);
  if (var.isNull()) {
    output_error.write("ERROR: NetCDF char variable '{:s}' has not been added to file '{:s}'\n", name, fname);
    return false;
  }

#ifdef NCDF_VERBOSE
  output.write("Ncxx4:: write { Writing Variable } \n");
#endif

  std::vector<size_t> start = {0};
  std::vector<size_t> counts = {size_t(n)};

  var.putVar(start, counts, data);

#ifdef NCDF_VERBOSE
  output.write("Ncxx4:: write { Done } \n");
#endif

  return true;
}

bool Ncxx4::write(char *var, const std::string &name, int n) {
  return write(var, name.c_str(), n);
}

bool Ncxx4::write(BoutReal *data, const char *name, int lx, int ly, int lz) {
  TRACE("Ncxx4::write(BoutReal)");

#ifdef NCDF_VERBOSE
  output.write("Ncxx4:: write(BoutReal, {:s})\n", name); 
#endif
  if(!is_valid())
    return false;

  if((lx < 0) || (ly < 0) || (lz < 0))
    return false;
  
  NcVar var = dataFile->getVar(name);
  if(var.isNull()) {
    output_error.write("ERROR: NetCDF BoutReal variable '{:s}' has not been added to file '{:s}'\n", name, fname);
    return false;
  }

  std::vector<size_t> start(3);
  start[0] = x0; start[1] = y0; start[2] = z0;
  std::vector<size_t> counts(3);
  counts[0] = lx; counts[1] = ly; counts[2] = lz;

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

  var.putVar(start, counts, data);

  return true;
}

bool Ncxx4::write(BoutReal *var, const std::string &name, int lx, int ly, int lz) {
  return write(var, name.c_str(), lx, ly, lz);
}

bool Ncxx4::write_perp(BoutReal *data, const std::string& name, int lx, int lz) {
  TRACE("Ncxx4::write_perp(BoutReal)");

#ifdef NCDF_VERBOSE
  output.write("Ncxx4:: write_perp(BoutReal, {:s})\n", name);
#endif
  if(!is_valid())
    return false;

  if((lx < 0) || (lz < 0))
    return false;

  NcVar var = dataFile->getVar(name);
  if(var.isNull()) {
    output_error.write(
        "ERROR: NetCDF BoutReal variable '{:s}' has not been added to file '{:s}'\n", name, fname);
    return false;
  }

  std::vector<size_t> start(2);
  start[0] = x0; start[1] = z0;
  std::vector<size_t> counts(2);
  counts[0] = lx; counts[1] = lz;

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

  for(int i=0;i<lx*lz;i++) {
    if(!finite(data[i]))
      data[i] = 0.0;
  }

  var.putVar(start, counts, data);

  return true;
}

/***************************************************************************
 * Record-based (time-dependent) data
 ***************************************************************************/

bool Ncxx4::read_rec(int *data, const char *name, int lx, int ly, int lz) {
#ifdef NCDF_VERBOSE
  output.write("Ncxx4:: read_rec(int, {:s})\n", name); 
#endif
  if(!is_valid())
    return false;

  if((lx < 0) || (ly < 0) || (lz < 0))
    return false;

  NcVar var = dataFile->getVar(name);
  
  if(var.isNull())
    return false;
  
  // NOTE: Probably should do something here to check t0

  std::vector<size_t> start(4);
  start[0] = t0; start[1] = x0; start[2] = y0; start[3] = z0;
  std::vector<size_t> counts(4);
  counts[0] = 1; counts[1] = lx; counts[2] = ly; counts[3] = lz;
  
  var.getVar(start, counts, data);
  
  return true;
}

bool Ncxx4::read_rec(int *var, const std::string &name, int lx, int ly, int lz) {
  return read_rec(var, name.c_str(), lx, ly, lz);
}

bool Ncxx4::read_rec(char *data, const char *name, int n) {
#ifdef NCDF_VERBOSE
  output.write("Ncxx4:: read_rec(char, {:s})\n", name);
#endif
  if (!is_valid()) {
    return false;
  }

  if (n < 0) {
    return false;
  }

  NcVar var = dataFile->getVar(name);

  if (var.isNull()) {
    return false;
  }

  // NOTE: Probably should do something here to check t0

  std::vector<size_t> start = {size_t(t0), 0};
  std::vector<size_t> counts = {1, size_t(n)};

  var.getVar(start, counts, data);

  return true;
}

bool Ncxx4::read_rec(char *var, const std::string &name, int n) {
  return read_rec(var, name.c_str(), n);
}

bool Ncxx4::read_rec(BoutReal *data, const char *name, int lx, int ly, int lz) {
#ifdef NCDF_VERBOSE
  output.write("Ncxx4:: read_rec(BoutReal, {:s})\n", name); 
#endif
  if(!is_valid())
    return false;

  if((lx < 0) || (ly < 0) || (lz < 0))
    return false;

  NcVar var = dataFile->getVar(name);
  
  if(var.isNull())
    return false;
  
  // NOTE: Probably should do something here to check t0

  std::vector<size_t> start(4);
  start[0] = t0; start[1] = x0; start[2] = y0; start[3] = z0;
  std::vector<size_t> counts(4);
  counts[0] = 1; counts[1] = lx; counts[2] = ly; counts[3] = lz;
  
  var.getVar(start, counts, data);
  
  return true;
}

bool Ncxx4::read_rec(BoutReal *var, const std::string &name, int lx, int ly, int lz) {
  return read_rec(var, name.c_str(), lx, ly, lz);
}

bool Ncxx4::read_rec_perp(BoutReal *data, const std::string& name, int lx, int lz) {
#ifdef NCDF_VERBOSE
  output.write("Ncxx4:: read_rec_perp(BoutReal, {:s})\n", name);
#endif
  if(!is_valid())
    return false;

  if((lx < 0) || (lz < 0))
    return false;

  NcVar var = dataFile->getVar(name);

  if(var.isNull())
    return false;

  // NOTE: Probably should do something here to check t0

  std::vector<size_t> start(3);
  start[0] = t0; start[1] = x0; start[2] = z0;
  std::vector<size_t> counts(4);
  counts[0] = 1; counts[1] = lx; counts[2] = lz;

  var.getVar(start, counts, data);

  return true;
}

bool Ncxx4::write_rec(int *data, const char *name, int lx, int ly, int lz) {
#ifdef NCDF_VERBOSE
  output.write("Ncxx4:: write_rec(int, {:s})\n", name); 
#endif
  if(!is_valid())
    return false;

  if((lx < 0) || (ly < 0) || (lz < 0))
    return false;
  
  // Try to find variable
  NcVar var = dataFile->getVar(name);
  if(var.isNull()) {
    output_error.write("ERROR: NetCDF int variable '{:s}' has not been added to file '{:s}'\n", name, fname);
    return false;
  }else {
    // Get record number
    if(rec_nr.find(name) == rec_nr.end()) {
      // Add to map
      rec_nr[name] = default_rec;
    }
  }
  
  std::vector<size_t> start(4);
  start[0] = rec_nr[name]; start[1] = x0; start[2] = y0; start[3] = z0;
  std::vector<size_t> counts(4);
  counts[0] = 1; counts[1] = lx; counts[2] = ly; counts[3] = lz;

#ifdef NCDF_VERBOSE
  output.write("Ncxx4:: write_rec { Writing variable } \n");
#endif

  var.putVar(start, counts, data);
  
  // Increment record number
  rec_nr[name] = rec_nr[name] + 1;

  return true;
}

bool Ncxx4::write_rec(int *var, const std::string &name, int lx, int ly, int lz) {
  return write_rec(var, name.c_str(), lx, ly, lz);
}

bool Ncxx4::write_rec(char *data, const char *name, int n) {
#ifdef NCDF_VERBOSE
  output.write("Ncxx4:: write_rec(char, {:s})\n", name);
#endif
  if (!is_valid()) {
    return false;
  }

  if (n < 0) {
    return false;
  }

  // Try to find variable
  NcVar var = dataFile->getVar(name);
  if (var.isNull()) {
    output_error.write("ERROR: NetCDF char variable '{:s}' has not been added to file '{:s}'\n", name, fname);
    return false;
  } else {
    // Get record number
    if (rec_nr.find(name) == rec_nr.end()) {
      // Add to map
      rec_nr[name] = default_rec;
    }
  }

  std::vector<size_t> start = {size_t(rec_nr[name]), 0};
  std::vector<size_t> counts = {1, size_t(n)};

#ifdef NCDF_VERBOSE
  output.write("Ncxx4:: write_rec { Writing variable } \n");
#endif

  var.putVar(start, counts, data);

  // Increment record number
  rec_nr[name] = rec_nr[name] + 1;

  return true;
}

bool Ncxx4::write_rec(char *var, const std::string &name, int n) {
  return write_rec(var, name.c_str(), n);
}

bool Ncxx4::write_rec(BoutReal *data, const char *name, int lx, int ly, int lz) {
  TRACE("Ncxx4::write_rec(BoutReal)");

#ifdef NCDF_VERBOSE
  output.write("Ncxx4::write_rec(BoutReal, {:s})\n", name); 
#endif
  if(!is_valid())
    return false;

  if((lx < 0) || (ly < 0) || (lz < 0))
    return false;

  // Try to find variable
  NcVar var = dataFile->getVar(name);
  if(var.isNull()) {
    output_error.write("ERROR: NetCDF BoutReal variable '{:s}' has not been added to file '{:s}'\n", name, fname);
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

  std::vector<size_t> start(4);
  start[0] = t; start[1] = x0; start[2] = y0; start[3] = z0;
  std::vector<size_t> counts(4);
  counts[0] = 1; counts[1] = lx; counts[2] = ly; counts[3] = lz;

  // Add the record
  var.putVar(start, counts, data);
  
  // Increment record number
  rec_nr[name] = rec_nr[name] + 1;

  return true;
}

bool Ncxx4::write_rec(BoutReal *var, const std::string &name, int lx, int ly, int lz) {
  return write_rec(var, name.c_str(), lx, ly, lz);
}

bool Ncxx4::write_rec_perp(BoutReal *data, const std::string& name, int lx, int lz) {
  TRACE("Ncxx4::write_rec_perp(BoutReal)");

#ifdef NCDF_VERBOSE
  output.write("Ncxx4::write_rec_perp(BoutReal, {:s})\n", name);
#endif
  if(!is_valid())
    return false;

  if((lx < 0) || (lz < 0))
    return false;

  // Try to find variable
  NcVar var = dataFile->getVar(name);
  if(var.isNull()) {
    output_error.write(
        "ERROR: NetCDF BoutReal variable '{:s}' has not been added to file '{:s}'\n", name, fname);
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
  output_info.write("INFO: NetCDF writing record {:d} of '{:s}' in '{:s}'\n",t, name, fname);
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

  for(int i=0;i<lx*lz;i++) {
    if(!finite(data[i]))
      data[i] = 0.0;
  }

  std::vector<size_t> start(3);
  start[0] = t; start[1] = x0; start[2] = z0;
  std::vector<size_t> counts(3);
  counts[0] = 1; counts[1] = lx; counts[2] = lz;

  // Add the record
  var.putVar(start, counts, data);

  // Increment record number
  rec_nr[name] = rec_nr[name] + 1;

  return true;
}

/***************************************************************************
 * Attributes
 ***************************************************************************/

void Ncxx4::setAttribute(const std::string &varname, const std::string &attrname,
                         const std::string &text) {
  TRACE("Ncxx4::setAttribute(string)");

  std::string existing_att;
  if (getAttribute(varname, attrname, existing_att)) {
    if (text != existing_att) {
      output_warn.write("Overwriting attribute '{:s}' of variable '{:s}' with '{:s}', "
                        "was previously '{:s}'",
                        attrname, varname, text, existing_att);
    }
  }
  // else: attribute does not exist, so just write it

  if (varname == "") {
    // write attribute of file
    dataFile->putAtt(attrname, text);
  } else {
    // write attribute of variable
    NcVar var = dataFile->getVar(varname);
    if (var.isNull()) {
      throw BoutException("Variable '{:s}' not in NetCDF file", varname);
    }

    var.putAtt(attrname, text);
  }
}

void Ncxx4::setAttribute(const std::string &varname, const std::string &attrname,
                         int value) {
  TRACE("Ncxx4::setAttribute(int)");

  int existing_att;
  if (getAttribute(varname, attrname, existing_att)) {
    if (value != existing_att) {
      output_warn.write("Overwriting attribute '{:s}' of variable '{:s}' with '{:d}', "
                        "was previously '{:d}'",
                        attrname, varname, value, existing_att);
    }
  }
  // else: attribute does not exist, so just write it

  if (varname == "") {
    // write attribute of file
    dataFile->putAtt(attrname, NcType::nc_INT, value);
  } else {
    // write attribute of variable
    NcVar var = dataFile->getVar(varname);
    if (var.isNull()) {
      throw BoutException("Variable '{:s}' not in NetCDF file", varname);
    }

    var.putAtt(attrname, NcType::nc_INT, value);
  }
}

void Ncxx4::setAttribute(const std::string &varname, const std::string &attrname,
                         BoutReal value) {
  TRACE("Ncxx4::setAttribute(BoutReal)");

  BoutReal existing_att;
  if (getAttribute(varname, attrname, existing_att)) {
    if (value != existing_att) {
      output_warn.write("Overwriting attribute '{:s}' of variable '{:s}' with '{:f}', "
                        "was previously '{:f}'",
                        attrname, varname, value, existing_att);
    }
  }
  // else: attribute does not exist, so just write it

  if (varname == "") {
    // write attribute of file
    dataFile->putAtt(attrname, NcType::nc_DOUBLE, value);
  } else {
    // write attribute of variable
    NcVar var = dataFile->getVar(varname);
    if (var.isNull()) {
      throw BoutException("Variable '{:s}' not in NetCDF file", varname);
    }

    var.putAtt(attrname, NcType::nc_DOUBLE, value);
  }
}

bool Ncxx4::getAttribute(const std::string &varname, const std::string &attrname, std::string &text) {
  TRACE("Ncxx4::getAttribute(string)");

  if (varname == "") {
    // attribute of file
    // Check if attribute exists without throwing exception when it doesn't
    std::multimap<std::string, NcGroupAtt> fileAtts_list = dataFile->getAtts();
    if (fileAtts_list.find(attrname) == fileAtts_list.end()) {
      return false;
    } else {
      NcGroupAtt fileAtt = dataFile->getAtt(attrname);
      fileAtt.getValues(text);

      return true;
    }
  } else {
    // attribute of variable
    NcVar var = dataFile->getVar(varname);
    if (var.isNull()) {
      throw BoutException("Variable '{:s}' not in NetCDF file", varname);
    }

    // Check if attribute exists without throwing exception when it doesn't
    auto varAtts_list = var.getAtts();
    if (varAtts_list.find(attrname) == varAtts_list.end()) {
      return false;
    } else {
      NcVarAtt varAtt = var.getAtt(attrname);
      varAtt.getValues(text);

      return true;
    }
  }
}

bool Ncxx4::getAttribute(const std::string &varname, const std::string &attrname, int &value) {
  TRACE("Ncxx4::getAttribute(int)");

  if (varname == "") {
    // attribute of file
    // Check if attribute exists without throwing exception when it doesn't
    std::multimap<string, NcGroupAtt> fileAtts_list = dataFile->getAtts();
    if (fileAtts_list.find(attrname) == fileAtts_list.end()) {
      return false;
    } else {
      NcGroupAtt fileAtt = dataFile->getAtt(attrname);
      fileAtt.getValues(&value);

      return true;
    }
  } else {
    NcVar var = dataFile->getVar(varname);
    if (var.isNull()) {
      throw BoutException("Variable '{:s}' not in NetCDF file", varname);
    }

    // Check if attribute exists without throwing exception when it doesn't
    auto varAtts_list = var.getAtts();
    if (varAtts_list.find(attrname) == varAtts_list.end()) {
      return false;
    } else {
      NcVarAtt varAtt = var.getAtt(attrname);
      varAtt.getValues(&value);

      return true;
    }
  }
}

bool Ncxx4::getAttribute(const std::string &varname, const std::string &attrname, BoutReal &value) {
  TRACE("Ncxx4::getAttribute(BoutReal)");

  if (varname == "") {
    // attribute of file
    // Check if attribute exists without throwing exception when it doesn't
    std::multimap<std::string, NcGroupAtt> fileAtts_list = dataFile->getAtts();
    if (fileAtts_list.find(attrname) == fileAtts_list.end()) {
      return false;
    } else {
      NcGroupAtt fileAtt = dataFile->getAtt(attrname);
      fileAtt.getValues(&value);

      return true;
    }
  } else {
    NcVar var = dataFile->getVar(varname);
    if (var.isNull()) {
      throw BoutException("Variable '{:s}' not in NetCDF file", varname);
    }

    // Check if attribute exists without throwing exception when it doesn't
    auto varAtts_list = var.getAtts();
    if (varAtts_list.find(attrname) == varAtts_list.end()) {
      return false;
    } else {
      NcVarAtt varAtt = var.getAtt(attrname);
      varAtt.getValues(&value);

      return true;
    }
  }
}


/***************************************************************************
 * Private functions
 ***************************************************************************/

std::vector<NcDim> Ncxx4::getDimVec(int nd) {
  std::vector<NcDim> vec(nd);
  for(int i=0;i<nd;i++)
    vec[i] = *dimList[i];
  return vec;
}

std::vector<NcDim> Ncxx4::getRecDimVec(int nd) {
  std::vector<NcDim> vec(nd);
  for(int i=0;i<nd;i++)
    vec[i] = *recDimList[i];
  return vec;
}

#endif // BOUT_HAS_NETCDF
