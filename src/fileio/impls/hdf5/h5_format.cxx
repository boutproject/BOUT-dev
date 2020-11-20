/**************************************************************************
 * Copyright 2015 J.T.Omotani, B.D.Dudson, S.Farley, M.V.Umansky, X.Q.Xu
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
#include "h5_format.hxx"

#ifdef HDF5

#include <utils.hxx>
#include <cmath>
#include <string>
#include <mpi.h>

#include <bout/mesh.hxx>
#include <output.hxx>
#include <msg_stack.hxx>
#include <boutcomm.hxx>

H5Format::H5Format(bool parallel_in, Mesh* mesh_in) : DataFormat(mesh_in) {
  parallel = parallel_in;
  x0 = y0 = z0 = t0 = 0;
  lowPrecision = false;
  fname = nullptr;
  dataFile = -1;
  chunk_length = 10; // could change this to try to optimize IO performance (i.e. allocate new chunks of disk space less often)
  
  dataFile_plist = H5Pcreate(H5P_FILE_ACCESS);
  if (dataFile_plist < 0)
    throw BoutException("Failed to create dataFile_plist");

#ifdef PHDF5
  if (parallel)
    if (H5Pset_fapl_mpio(dataFile_plist, BoutComm::get(), MPI_INFO_NULL) < 0)
      throw BoutException("Failed to set dataFile_plist");
#endif

  dataSet_plist = H5Pcreate(H5P_DATASET_XFER);
  if (dataSet_plist < 0)
    throw BoutException("Failed to create dataSet_plist");

#ifdef PHDF5
  if (parallel)
    if (H5Pset_dxpl_mpio(dataSet_plist, H5FD_MPIO_INDEPENDENT) < 0) // Default, independent writes
//     if (H5Pset_dxpl_mpio(dataSet_plist, H5FD_MPIO_COLLECTIVE) < 0) // Alternative, collective writes
      throw BoutException("Failed to set dataSet_plist");
#endif

  // Disable automatic printing of error messages so that we can catch
  // errors without printing error messages to stdout
  if (H5Eset_auto(H5E_DEFAULT, nullptr, nullptr) < 0)
    throw BoutException("Failed to set error stack to not print errors");
}

H5Format::H5Format(const char *name, bool parallel_in, Mesh* mesh_in)
  : DataFormat(mesh_in) {
  parallel = parallel_in;
  x0 = y0 = z0 = t0 = 0;
  lowPrecision = false;
  fname = nullptr;
  dataFile = -1;
  chunk_length = 10; // could change this to try to optimize IO performance (i.e. allocate new chunks of disk space less often)
  
  dataFile_plist = H5Pcreate(H5P_FILE_ACCESS);
  if (dataFile_plist < 0)
    throw BoutException("Failed to create dataFile_plist");

#ifdef PHDF5
  if (parallel)
    if (H5Pset_fapl_mpio(dataFile_plist, BoutComm::get(), MPI_INFO_NULL) < 0)
      throw BoutException("Failed to set dataFile_plist");
#endif  

  dataSet_plist = H5Pcreate(H5P_DATASET_XFER);
  if (dataSet_plist < 0)
    throw BoutException("Failed to create dataSet_plist");

#ifdef PHDF5
  if (parallel)
    if (H5Pset_dxpl_mpio(dataSet_plist, H5FD_MPIO_COLLECTIVE) < 0)
      throw BoutException("Failed to set dataSet_plist");
#endif

  // Disable automatic printing of error messages so that we can catch
  // errors without printing error messages to stdout
  if (H5Eset_auto(H5E_DEFAULT, nullptr, nullptr) < 0)
    throw BoutException("Failed to set error stack to not print errors");

  H5Format::openr(name);
}

H5Format::~H5Format() {
  H5Format::close();
  H5Pclose(dataFile_plist);
}

bool H5Format::openr(const char *name) {
  TRACE("H5Format::openr");

  if(dataFile > 0) // Already open. Close then re-open
    close();

  dataFile = H5Fopen(name, H5F_ACC_RDONLY, dataFile_plist);
  if( dataFile < 0 ) {
    throw BoutException("Failed to open dataFile");
  }

  return true;
}

bool H5Format::openw(const char *name, bool append) {
  TRACE("H5Format::openw");
  
  if(dataFile > 0) // Already open. Close then re-open
    close(); 

  if(append) {
    dataFile = H5Fopen(name, H5F_ACC_RDWR, dataFile_plist);
  }
  else {
    dataFile = H5Fcreate(name, H5F_ACC_TRUNC, H5P_DEFAULT, dataFile_plist);
  }
  if( dataFile < 0 ) {
    throw BoutException("Failed to open dataFile");
  }

  fname = copy_string(name);

  return true;
}

bool H5Format::is_valid() { return dataFile >= 0; }

void H5Format::close() {
  TRACE("H5Format::close");
  
  if (H5Format::is_valid()) {
    H5Fclose(dataFile);
    dataFile = -1;
  }
}

void H5Format::flush() {
  if(!is_valid())
    return;
  
  if (H5Fflush(dataFile,H5F_SCOPE_LOCAL) < 0) {
    throw BoutException("Failed to flush dataFile");
  }
}

const std::vector<int> H5Format::getSize(const char *name) {
  TRACE("H5Format::getSize");

  std::vector<int> size;

  if(!is_valid())
    return size;
  
  hid_t dataSet = H5Dopen(dataFile, name, H5P_DEFAULT);
  if (dataSet < 0) {
    // variable does not exist in file, so return empty size
    return size;
  }
  
  hid_t dataSpace = H5Dget_space(dataSet);
  if (dataSpace < 0)
    throw BoutException("Failed to create dataSpace");
  
  int nd = H5Sget_simple_extent_ndims(dataSpace);
  if (nd < 0)
    throw BoutException("Failed to get dataSpace ndims");
  
  if (nd==0) {
    if (H5Sclose(dataSpace) < 0)
      throw BoutException("Failed to close dataSpace");
    if (H5Dclose(dataSet) < 0)
      throw BoutException("Failed to close dataSet");
    
    size.push_back(1);
    return size;
  } else {
    std::vector<hsize_t> dims(nd);
    int error = H5Sget_simple_extent_dims(dataSpace, dims.data(), nullptr);
    if (error < 0)
      throw BoutException("Failed to get dimensions of dataSpace");
    
    if (H5Sclose(dataSpace) < 0)
      throw BoutException("Failed to close dataSpace");
    if (H5Dclose(dataSet) < 0)
      throw BoutException("Failed to close dataSet");
    
    std::copy(begin(dims), end(dims), std::back_inserter(size));
  }

  return size;
}

const std::vector<int> H5Format::getSize(const std::string &var) {
  return getSize(var.c_str());
}

bool H5Format::setGlobalOrigin(int x, int y, int z) {
  x0 = x;
  y0 = y;
  z0 = z;
  x0_local = 0;
  y0_local = 0;
  z0_local = 0;
  
  return true;
}

bool H5Format::setLocalOrigin(int x, int y, int z, int offset_x, int offset_y, int offset_z) {
  
  if(!setGlobalOrigin(x + mesh->OffsetX, y + mesh->OffsetY, z + mesh->OffsetZ))
    return false;
  
  x0_local = offset_x;
  y0_local = offset_y;
  z0_local = offset_z;
  
  return true;
}

bool H5Format::setRecord(int t) {
  t0 = t;

  return true;
}

// Add a variable to the file
bool H5Format::addVar(const std::string &name, bool repeat, hid_t write_hdf5_type,
    std::string datatype, int lx, int ly, int lz) {
  hid_t dataSet = H5Dopen(dataFile, name.c_str(), H5P_DEFAULT);
  if (dataSet >= 0) { // >=0 means variable already exists, so return.
    if (H5Dclose(dataSet) < 0)
      throw BoutException("Failed to close dataSet");
    return true;
  }

  int nd = 0;
  if (datatype == "scalar") nd = 0;
  else if (datatype == "vector") nd = 1;
  else if (datatype == "FieldX") nd = 1;
  else if (datatype == "Field2D") nd = 2;
  else if (datatype == "FieldPerp") nd = 2;
  else if (datatype == "Field3D") nd = 3;
  else throw BoutException("Unrecognized datatype '"+datatype+"'");

  if (repeat) {
    // add time dimension
    datatype += "_t";
    nd += 1;

    hsize_t init_size[4];
    if (parallel) {
      init_size[0]=0;
      init_size[1] = lx == 0 ? mesh->GlobalNx-2*mesh->xstart : lx;
      if (datatype == "FieldPerp_t") {
        init_size[2] = lz == 0 ? mesh->GlobalNz : lz;
      } else {
        init_size[2] = ly == 0 ? mesh->GlobalNy - 2*mesh->ystart : ly;
      }
      init_size[3] = lz == 0 ? mesh->GlobalNz : lz;
    }
    else {
      init_size[0]=0;
      init_size[1] = lx == 0 ? mesh->LocalNx : lx;
      if (datatype == "FieldPerp_t") {
        init_size[2] = lz == 0 ? mesh->LocalNz : lz;
      } else {
        init_size[2] = ly == 0 ? mesh->LocalNy : ly;
      }
      init_size[3] = lz == 0 ? mesh->LocalNz : lz;
    }

    // Modify dataset creation properties, i.e. enable chunking.
    hid_t propertyList = H5Pcreate(H5P_DATASET_CREATE);
    if (propertyList < 0)
      throw BoutException("Failed to create propertyList");
    hsize_t chunk_dims[4],max_dims[4];
    max_dims[0] = H5S_UNLIMITED; max_dims[1]=init_size[1]; max_dims[2]=init_size[2]; max_dims[3]=init_size[3];
    chunk_dims[0] = chunk_length; chunk_dims[1]=init_size[1]; chunk_dims[2]=init_size[2]; chunk_dims[3]=init_size[3];
    if (H5Pset_chunk(propertyList, nd, chunk_dims) < 0)
      throw BoutException("Failed to set chunk property");

    hid_t init_space = H5Screate_simple(nd, init_size, max_dims);
    if (init_space < 0)
      throw BoutException("Failed to create init_space");
    dataSet = H5Dcreate(dataFile, name.c_str(), write_hdf5_type, init_space, H5P_DEFAULT, propertyList, H5P_DEFAULT);
    if (dataSet < 0)
      throw BoutException("Failed to create dataSet");

    // Add attribute to say what kind of field this is

    // Create new dataspace for attribute
    hid_t attribute_dataspace = H5Screate(H5S_SCALAR);
    if (attribute_dataspace < 0)
      throw BoutException("Failed to create attribute_dataspace");

    // Create new string datatype for attribute
    hid_t variable_length_string_type = H5Tcopy(H5T_C_S1);
    if (variable_length_string_type < 0)
      throw BoutException("Failed to create variable_length_string_type");
    if (H5Tset_size(variable_length_string_type, H5T_VARIABLE) < 0)
      throw BoutException("Failed to create string type");

    // Create attribute and write to it
    hid_t myatt_in = H5Acreate(dataSet, "bout_type", variable_length_string_type, attribute_dataspace, H5P_DEFAULT, H5P_DEFAULT);
    if (myatt_in < 0)
      throw BoutException("Failed to create attribute");
    if (H5Awrite(myatt_in, variable_length_string_type, &datatype) < 0)
      throw BoutException("Failed to write attribute");

    if (H5Pclose(propertyList) < 0)
      throw BoutException("Failed to close propertyList");
    if (H5Sclose(init_space) < 0)
      throw BoutException("Failed to close init_space");
    if (H5Sclose(attribute_dataspace) < 0)
      throw BoutException("Failed to close attribute_dataspace");
    if (H5Tclose(variable_length_string_type) < 0)
      throw BoutException("Failed to close variable_length_string_type");
    if (H5Aclose(myatt_in) < 0)
      throw BoutException("Failed to close myatt_in");
  } else {
    dataSet = H5Dopen(dataFile, name.c_str(), H5P_DEFAULT);
    if (dataSet < 0) {
      // Negative value indicates error, i.e. file does not exist, so create:
      hsize_t init_size[3];
      if (parallel) {
        init_size[0] = lx == 0 ? mesh->GlobalNx - 2 * mesh->xstart : lx;
        if (datatype == "FieldPerp") {
          init_size[1] = lz == 0 ? mesh->GlobalNz : lz;
        } else {
          init_size[1] = ly == 0 ? mesh->GlobalNy - 2 * mesh->ystart : ly;
        }
        init_size[2] = lz == 0 ? mesh->GlobalNz : lz;
      } else {
        init_size[0] = lx == 0 ? mesh->LocalNx : lx;
        if (datatype == "FieldPerp") {
          init_size[1] = lz == 0 ? mesh->LocalNz : lz;
        } else {
          init_size[1] = ly == 0 ? mesh->LocalNy : ly;
        }
        init_size[2] = lz == 0 ? mesh->LocalNz : lz;
      }

      // Create value for attribute to say what kind of field this is

      if (nd==0) {
        // Need to write a scalar, not a 0-d array
        nd = 1;
        init_size[0] = 1;
      }

      hid_t init_space = H5Screate_simple(nd, init_size, init_size);
      if (init_space < 0)
        throw BoutException("Failed to create init_space");
      dataSet = H5Dcreate(dataFile, name.c_str(), write_hdf5_type, init_space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      if (dataSet < 0)
        throw BoutException("Failed to create dataSet");

      // Add attribute to say what kind of field this is
      setAttribute(dataSet, "bout_type", datatype);
    }
  }

  if (H5Dclose(dataSet) < 0)
    throw BoutException("Failed to close dataSet");
  return true;
}

bool H5Format::addVarInt(const std::string &name, bool repeat) {
  return addVar(name, repeat, H5T_NATIVE_INT, "scalar");
}

bool H5Format::addVarIntVec(const std::string &name, bool repeat, size_t size) {
  return addVar(name, repeat, H5T_NATIVE_INT, "vector", size);
}

bool H5Format::addVarBoutReal(const std::string &name, bool repeat) {
  auto h5_float_type = lowPrecision ? H5T_NATIVE_FLOAT : H5T_NATIVE_DOUBLE;
  return addVar(name, repeat, h5_float_type, "scalar");
}

bool H5Format::addVarField2D(const std::string &name, bool repeat) {
  auto h5_float_type = lowPrecision ? H5T_NATIVE_FLOAT : H5T_NATIVE_DOUBLE;
  return addVar(name, repeat, h5_float_type, "Field2D");
}

bool H5Format::addVarField3D(const std::string &name, bool repeat) {
  auto h5_float_type = lowPrecision ? H5T_NATIVE_FLOAT : H5T_NATIVE_DOUBLE;
  return addVar(name, repeat, h5_float_type, "Field3D");
}

bool H5Format::addVarFieldPerp(const std::string &name, bool repeat) {
  auto h5_float_type = lowPrecision ? H5T_NATIVE_FLOAT : H5T_NATIVE_DOUBLE;
  return addVar(name, repeat, h5_float_type, "FieldPerp");
}

bool H5Format::read(int *data, const char *name, int lx, int ly, int lz) {
  return read(data, H5T_NATIVE_INT, name, lx, ly, lz);
}

bool H5Format::read(int *var, const std::string &name, int lx, int ly, int lz) {
  return read(var, name.c_str(), lx, ly, lz);
}

bool H5Format::read(BoutReal *data, const char *name, int lx, int ly, int lz) {
  return read(data, H5T_NATIVE_DOUBLE, name, lx, ly, lz);
}

bool H5Format::read(BoutReal *var, const std::string &name, int lx, int ly, int lz) {
  return read(var, name.c_str(), lx, ly, lz);
}

bool H5Format::read(void *data, hid_t hdf5_type, const char *name, int lx, int ly, int lz) {
  TRACE("H5Format::read(void)");

  if(!is_valid())
    return false;

  if((lx < 0) || (ly < 0) || (lz < 0))
    return false;
  
  int nd = 0; // Number of dimensions
  if(lx != 0) nd = 1;
  if(ly != 0) nd = 2;
  if(lz != 0) nd = 3;
  hsize_t counts[3],offset[3],offset_local[3],init_size_local[3];
  counts[0]=lx; counts[1]=ly; counts[2]=lz;
  offset[0]=x0; offset[1]=y0; offset[2]=z0;
  offset_local[0]=x0_local;
  offset_local[1]=y0_local;
  offset_local[2]=z0_local;

  // Want to be able to use without needing mesh to be initialised; makes hyperslab selection redundant
  init_size_local[0]=offset_local[0]+counts[0];
  init_size_local[1]=offset_local[1]+counts[1];
  init_size_local[2]=offset_local[2]+counts[2];
  
  hid_t mem_space = H5Screate_simple(nd, init_size_local, init_size_local);
  if (mem_space < 0)
    throw BoutException("Failed to create mem_space");

  hid_t dataSet = H5Dopen(dataFile, name, H5P_DEFAULT);
  if (dataSet < 0) {
    return false;
  }
  
  hid_t dataSpace = H5Dget_space(dataSet);
  if (dataSpace < 0)
    throw BoutException("Failed to create dataSpace");
  if (nd > 0 && !(nd==1 && lx==1))
    if (H5Sselect_hyperslab(dataSpace, H5S_SELECT_SET, offset, /*stride=*/nullptr, counts,
                            /*block=*/nullptr) < 0)
      throw BoutException("Failed to select hyperslab");
  
  if (H5Dread(dataSet, hdf5_type, mem_space, dataSpace, H5P_DEFAULT, data) < 0)
    throw BoutException("Failed to read data");
  
  if (H5Sclose(mem_space) < 0)
    throw BoutException("Failed to close mem_space");
  if (H5Sclose(dataSpace) < 0)
    throw BoutException("Failed to close dataSpace");
  if (H5Dclose(dataSet) < 0)
    throw BoutException("Failed to close dataSet");

  return true;
}

bool H5Format::read_perp(BoutReal *data, const std::string& name, int lx, int lz) {
  TRACE("H5Format::read(void)");

  hid_t hdf5_type = H5T_NATIVE_DOUBLE;

  if(!is_valid())
    return false;

  if((lx < 0) || (lz < 0))
    return false;

  int nd = 0; // Number of dimensions
  if(lx != 0) nd = 1;
  if(lz != 0) nd = 2;
  hsize_t counts[2],offset[2],offset_local[2],init_size_local[2];
  counts[0]=lx; counts[1]=lz;
  offset[0]=x0; offset[1]=z0;
  offset_local[0]=x0_local;
  offset_local[1]=z0_local;

  // Want to be able to use without needing mesh to be initialised; makes hyperslab selection redundant
  init_size_local[0]=offset_local[0]+counts[0];
  init_size_local[1]=offset_local[1]+counts[1];

  hid_t mem_space = H5Screate_simple(nd, init_size_local, init_size_local);
  if (mem_space < 0)
    throw BoutException("Failed to create mem_space");

  hid_t dataSet = H5Dopen(dataFile, name.c_str(), H5P_DEFAULT);
  if (dataSet < 0) {
    return false;
  }

  hid_t dataSpace = H5Dget_space(dataSet);
  if (dataSpace < 0)
    throw BoutException("Failed to create dataSpace");
  if (nd > 0 && !(nd==1 && lx==1))
    if (H5Sselect_hyperslab(dataSpace, H5S_SELECT_SET, offset, /*stride=*/nullptr, counts,
                            /*block=*/nullptr) < 0)
      throw BoutException("Failed to select hyperslab");

  if (H5Dread(dataSet, hdf5_type, mem_space, dataSpace, H5P_DEFAULT, data) < 0)
    throw BoutException("Failed to read data");

  if (H5Sclose(mem_space) < 0)
    throw BoutException("Failed to close mem_space");
  if (H5Sclose(dataSpace) < 0)
    throw BoutException("Failed to close dataSpace");
  if (H5Dclose(dataSet) < 0)
    throw BoutException("Failed to close dataSet");

  return true;
}

bool H5Format::write(int *data, const char *name, int lx, int ly, int lz) {
  return write(data, H5T_NATIVE_INT, name, lx, ly, lz);
}

bool H5Format::write(int *var, const std::string &name, int lx, int ly, int lz) {
  return write(var, name.c_str(), lx, ly, lz);
}

bool H5Format::write(BoutReal *data, const char *name, int lx, int ly, int lz) {
  
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
    
    return write(data, H5T_NATIVE_DOUBLE, name, lx, ly, lz);
  }
  else {
    return write(data, H5T_NATIVE_DOUBLE, name, lx, ly, lz);
  }
  
}

bool H5Format::write(BoutReal *var, const std::string &name, int lx, int ly, int lz) {
  return write(var, name.c_str(), lx, ly, lz);
}

bool H5Format::write(void *data, hid_t mem_hdf5_type, const char *name, int lx, int ly,
                     int lz) {
  TRACE("H5Format::write(void)");

  if(!is_valid())
    return false;

  if((lx < 0) || (ly < 0) || (lz < 0))
    return false;

  int nd = 0; // Number of dimensions
  if(lx != 0) nd = 1;
  if(ly != 0) nd = 2;
  if(lz != 0) nd = 3;
  hsize_t counts[3], offset[3], offset_local[3], init_size_local[3];
  counts[0] = lx;
  counts[1] = ly;
  counts[2] = lz;
  offset[0] = x0;
  offset[1] = y0;
  offset[2] = z0;
  offset_local[0] = x0_local;
  offset_local[1] = y0_local;
  offset_local[2] = z0_local;
  init_size_local[0] = mesh->LocalNx;
  init_size_local[1] = mesh->LocalNy;
  init_size_local[2] = mesh->LocalNz;

  if (nd==0) {
    // Need to write a scalar, not a 0-d array
    nd = 1;
    counts[0] = 1;
    offset[0] = 0;
    offset_local[0] = 0;
    init_size_local[0] = 1;
  }
  
  hid_t mem_space = H5Screate_simple(nd, init_size_local, init_size_local);
  if (mem_space < 0)
    throw BoutException("Failed to create mem_space");
  if (H5Sselect_hyperslab(mem_space, H5S_SELECT_SET, offset_local, /*stride=*/nullptr,
                          counts, /*block=*/nullptr) < 0)
    throw BoutException("Failed to select hyperslab");
  
  hid_t dataSet = H5Dopen(dataFile, name, H5P_DEFAULT);
  if (dataSet < 0) {
    output_error.write("ERROR: HDF5 variable '%s' has not been added to file '%s'\n", name, fname);
    return false;
  }
  
  hid_t dataSpace = H5Dget_space(dataSet);
  if (dataSpace < 0)
    throw BoutException("Failed to create dataSpace");
  if (H5Sselect_hyperslab(dataSpace, H5S_SELECT_SET, offset, /*stride=*/nullptr, counts,
                          /*block=*/nullptr) < 0)
    throw BoutException("Failed to select hyperslab");
  
  if (H5Dwrite(dataSet, mem_hdf5_type, mem_space, dataSpace, dataSet_plist, data) < 0)
    throw BoutException("Failed to write data");
  
  if (H5Sclose(mem_space) < 0)
    throw BoutException("Failed to close mem_space");
  if (H5Sclose(dataSpace) < 0)
    throw BoutException("Failed to close dataSpace");
  if (H5Dclose(dataSet) < 0)
    throw BoutException("Failed to close dataSet");

  return true;
}

bool H5Format::write_perp(BoutReal *data, const std::string& name, int lx, int lz) {
  TRACE("H5Format::write_perp(void)");

  hid_t mem_hdf5_type = H5T_NATIVE_DOUBLE;

  if(!is_valid())
    return false;

  if((lx < 0) || (lz < 0))
    return false;

  int nd = 0; // Number of dimensions
  if(lx != 0) nd = 1;
  if(lz != 0) nd = 2;
  hsize_t counts[2], offset[2], offset_local[2], init_size_local[2];
  counts[0] = lx;
  counts[1] = lz;
  offset[0] = x0;
  offset[1] = z0;
  offset_local[0] = x0_local;
  offset_local[1] = z0_local;
  init_size_local[0] = mesh->LocalNx;
  init_size_local[1] = mesh->LocalNz;

  if (nd==0) {
    // Need to write a scalar, not a 0-d array
    nd = 1;
    counts[0] = 1;
    offset[0] = 0;
    offset_local[0] = 0;
    init_size_local[0] = 1;
  }

  hid_t mem_space = H5Screate_simple(nd, init_size_local, init_size_local);
  if (mem_space < 0)
    throw BoutException("Failed to create mem_space");
  if (H5Sselect_hyperslab(mem_space, H5S_SELECT_SET, offset_local, /*stride=*/nullptr,
                          counts, /*block=*/nullptr) < 0)
    throw BoutException("Failed to select hyperslab");

  hid_t dataSet = H5Dopen(dataFile, name.c_str(), H5P_DEFAULT);
  if (dataSet < 0) {
    output_error.write("ERROR: HDF5 variable '%s' has not been added to file '%s'\n", name.c_str(), fname);
    return false;
  }

  hid_t dataSpace = H5Dget_space(dataSet);
  if (dataSpace < 0)
    throw BoutException("Failed to create dataSpace");
  if (H5Sselect_hyperslab(dataSpace, H5S_SELECT_SET, offset, /*stride=*/nullptr, counts,
                          /*block=*/nullptr) < 0)
    throw BoutException("Failed to select hyperslab");

  if (H5Dwrite(dataSet, mem_hdf5_type, mem_space, dataSpace, dataSet_plist, data) < 0)
    throw BoutException("Failed to write data");

  if (H5Sclose(mem_space) < 0)
    throw BoutException("Failed to close mem_space");
  if (H5Sclose(dataSpace) < 0)
    throw BoutException("Failed to close dataSpace");
  if (H5Dclose(dataSet) < 0)
    throw BoutException("Failed to close dataSet");

  return true;
}

/***************************************************************************
 * Record-based (time-dependent) data
 ***************************************************************************/

bool H5Format::read_rec(int *data, const char *name, int lx, int ly, int lz) {
  return read_rec(data, H5T_NATIVE_INT, name, lx, ly, lz);
}

bool H5Format::read_rec(int *var, const std::string &name, int lx, int ly, int lz) {
  return read_rec(var, name.c_str(), lx, ly, lz);
}

bool H5Format::read_rec(BoutReal *data, const char *name, int lx, int ly, int lz) {
  
  return read_rec(data, H5T_NATIVE_DOUBLE, name, lx, ly, lz);
  
}

bool H5Format::read_rec(BoutReal *var, const std::string &name, int lx, int ly, int lz) {
  return read_rec(var, name.c_str(), lx, ly, lz);
}

bool H5Format::read_rec(void *data, hid_t hdf5_type, const char *name, int lx, int ly,
                        int lz) {
  if (!is_valid()) {
    return false;
  }

  if ((lx < 0) || (ly < 0) || (lz < 0)) {
    return false;
  }

  int nd = 1; // Number of dimensions
  if (lx != 0) {
    nd = 2;
  }
  if (ly != 0) {
    nd = 3;
  }
  if (lz != 0) {
    nd = 4;
  }
  hsize_t counts[4], offset[4];
  hsize_t offset_local[3], init_size_local[3];
  counts[0] = 1;
  counts[1] = lx;
  counts[2] = ly;
  counts[3] = lz;
  offset[0] = t0;
  offset[1] = x0;
  offset[2] = y0;
  offset[3] = z0;
  offset_local[0] = x0_local;
  offset_local[1] = y0_local;
  offset_local[2] = z0_local;
  init_size_local[0] = mesh->LocalNx;
  init_size_local[1] = mesh->LocalNy;
  init_size_local[2] = mesh->LocalNz;

  if (nd == 1) {
    // Need to write a time-series of scalars
    nd = 1;
    counts[1] = 1;
    offset[1] = 0;
    init_size_local[0] = 1;
  }

  hid_t mem_space = H5Screate_simple(nd, init_size_local, init_size_local);
  if (mem_space < 0)
    throw BoutException("Failed to create mem_space");
  if (H5Sselect_hyperslab(mem_space, H5S_SELECT_SET, offset_local, /*stride=*/nullptr,
                          counts, /*block=*/nullptr) < 0)
    throw BoutException("Failed to select hyperslab");

  hid_t dataSet = H5Dopen(dataFile, name, H5P_DEFAULT);
  if (dataSet < 0)
    throw BoutException("Failed to open dataSet");

  hid_t dataSpace = H5Dget_space(dataSet);
  if (dataSpace < 0)
    throw BoutException("Failed to create dataSpace");
  if (H5Sselect_hyperslab(dataSpace, H5S_SELECT_SET, offset, /*stride=*/nullptr, counts,
                          /*block=*/nullptr) < 0)
    throw BoutException("Failed to select hyperslab");

  if (H5Dread(dataSet, hdf5_type, mem_space, dataSpace, H5P_DEFAULT, data) < 0)
    throw BoutException("Failed to read data");

  if (H5Sclose(mem_space) < 0)
    throw BoutException("Failed to close mem_space");
  if (H5Sclose(dataSpace) < 0)
    throw BoutException("Failed to close dataSpace");
  if (H5Dclose(dataSet) < 0)
    throw BoutException("Failed to close dataSet");

  return true;
}

bool H5Format::read_rec_perp(BoutReal *data, const std::string& name, int lx, int lz) {
  if (!is_valid()) {
    return false;
  }

  hid_t hdf5_type = H5T_NATIVE_DOUBLE;

  if ((lx < 0) || (lz < 0)) {
    return false;
  }

  int nd = 1; // Number of dimensions
  if (lx != 0) {
    nd = 2;
  }
  if (lz != 0) {
    nd = 3;
  }
  hsize_t counts[3], offset[3];
  hsize_t offset_local[2], init_size_local[2];
  counts[0] = 1;
  counts[1] = lx;
  counts[2] = lz;
  offset[0] = t0;
  offset[1] = x0;
  offset[2] = z0;
  offset_local[0] = x0_local;
  offset_local[1] = z0_local;
  init_size_local[0] = mesh->LocalNx;
  init_size_local[1] = mesh->LocalNz;

  if (nd == 1) {
    // Need to write a time-series of scalars
    nd = 1;
    counts[1] = 1;
    offset[1] = 0;
    init_size_local[0] = 1;
  }

  hid_t mem_space = H5Screate_simple(nd, init_size_local, init_size_local);
  if (mem_space < 0)
    throw BoutException("Failed to create mem_space");
  if (H5Sselect_hyperslab(mem_space, H5S_SELECT_SET, offset_local, /*stride=*/nullptr,
                          counts, /*block=*/nullptr) < 0)
    throw BoutException("Failed to select hyperslab");

  hid_t dataSet = H5Dopen(dataFile, name.c_str(), H5P_DEFAULT);
  if (dataSet < 0)
    throw BoutException("Failed to open dataSet");

  hid_t dataSpace = H5Dget_space(dataSet);
  if (dataSpace < 0)
    throw BoutException("Failed to create dataSpace");
  if (H5Sselect_hyperslab(dataSpace, H5S_SELECT_SET, offset, /*stride=*/nullptr, counts,
                          /*block=*/nullptr) < 0)
    throw BoutException("Failed to select hyperslab");

  if (H5Dread(dataSet, hdf5_type, mem_space, dataSpace, H5P_DEFAULT, data) < 0)
    throw BoutException("Failed to read data");

  if (H5Sclose(mem_space) < 0)
    throw BoutException("Failed to close mem_space");
  if (H5Sclose(dataSpace) < 0)
    throw BoutException("Failed to close dataSpace");
  if (H5Dclose(dataSet) < 0)
    throw BoutException("Failed to close dataSet");

  return true;
}

bool H5Format::write_rec(int *data, const char *name, int lx, int ly, int lz) {
  return write_rec(data, H5T_NATIVE_INT, name, lx, ly, lz);
}

bool H5Format::write_rec(int *var, const std::string &name, int lx, int ly, int lz) {
  return write_rec(var, name.c_str(), lx, ly, lz);
}

bool H5Format::write_rec(BoutReal *data, const char *name, int lx, int ly, int lz) {
  
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
    return write_rec(data, H5T_NATIVE_DOUBLE, name, lx, ly, lz);
  }
  
  return write_rec(data, H5T_NATIVE_DOUBLE, name, lx, ly, lz);
}

bool H5Format::write_rec(BoutReal *var, const std::string &name, int lx, int ly, int lz) {
  return write_rec(var, name.c_str(), lx, ly, lz);
}

bool H5Format::write_rec(void *data, hid_t mem_hdf5_type, const char *name, int lx, int ly, int lz) {
  if(!is_valid())
    return false;

  if((lx < 0) || (ly < 0) || (lz < 0))
    return false;

  int nd = 1; // Number of dimensions
  if(lx != 0) nd = 2;
  if(ly != 0) nd = 3;
  if(lz != 0) nd = 4;
  int nd_local = nd-1;
  hsize_t counts[4], offset[4];
  hsize_t counts_local[3], offset_local[3], init_size_local[3];
  counts[0] = 1;
  counts[1] = lx;
  counts[2] = ly;
  counts[3] = lz;
  counts_local[0] = lx;
  counts_local[1] = ly;
  counts_local[2] = lz;
  // Do this later, after setting t0//  offset[0]=t0;
  offset[1] = x0;
  offset[2] = y0;
  offset[3] = z0;
  offset_local[0] = x0_local;
  offset_local[1] = y0_local;
  offset_local[2] = z0_local;
  init_size_local[0] = mesh->LocalNx;
  init_size_local[1] = mesh->LocalNy;
  init_size_local[2] = mesh->LocalNz;

  if (nd_local == 0) {
    nd_local = 1;
    // Need to write a time-series of scalars
    counts_local[0] = 1;
    offset_local[0] = 0;
    init_size_local[0] = 1;
  }
  
  hid_t mem_space = H5Screate_simple(nd_local, init_size_local, init_size_local);
  if (mem_space < 0)
    throw BoutException("Failed to create mem_space");
  if (H5Sselect_hyperslab(mem_space, H5S_SELECT_SET, offset_local, /*stride=*/nullptr,
                          counts_local, /*block=*/nullptr) < 0)
    throw BoutException("Failed to select hyperslab");
  
  hid_t dataSet = H5Dopen(dataFile, name, H5P_DEFAULT);
  if (dataSet >= 0) { // >=0 means file exists, so open. Else error.
    
    hsize_t dims[4] = {};
    hid_t dataSpace = H5Dget_space(dataSet);
    if (dataSpace < 0)
      throw BoutException("Failed to create dataSpace");
    if (H5Sget_simple_extent_dims(dataSpace, dims, /*maxdims=*/nullptr) < 0)
      throw BoutException("Failed to get dims");
    dims[0]+=1;
    if (t0 == -1) {
      // Want t0 to be last record
      t0 = dims[0]-1;
    }
    
    if (H5Dset_extent(dataSet, dims) < 0)
      throw BoutException("Failed to extend dataSet");
    
    if (H5Sclose(dataSpace) < 0)
      throw BoutException("Failed to close dataSpace");
    
  }
  else {
    output_error.write("ERROR: HDF5 variable '%s' has not been added to file '%s'\n", name, fname);
    return false;
  }

  offset[0]=t0;
  
  hid_t dataSpace = H5Dget_space(dataSet);
  if (dataSpace < 0)
    throw BoutException("Failed to create dataSpace");
  if (H5Sselect_hyperslab(dataSpace, H5S_SELECT_SET, offset, /*stride=*/nullptr, counts,
                          /*block=*/nullptr) < 0)
    throw BoutException("Failed to select hyperslab");
  
  if (H5Dwrite(dataSet, mem_hdf5_type, mem_space, dataSpace, dataSet_plist, data) < 0)
    throw BoutException("Failed to write data");
  
  if (H5Sclose(mem_space) < 0)
    throw BoutException("Failed to close mem_space");
  if (H5Sclose(dataSpace) < 0)
    throw BoutException("Failed to close dataSpace");
  if (H5Dclose(dataSet) < 0)
    throw BoutException("Failed to close dataSet");
  
  return true;
}

bool H5Format::write_rec_perp(BoutReal *data, const std::string& name, int lx, int lz) {
  if(!is_valid())
    return false;

  hid_t mem_hdf5_type = H5T_NATIVE_DOUBLE;

  if((lx < 0) || (lz < 0))
    return false;

  int nd = 1; // Number of dimensions
  if(lx != 0) nd = 2;
  if(lz != 0) nd = 3;
  int nd_local = nd-1;
  hsize_t counts[3], offset[3];
  hsize_t counts_local[2], offset_local[2], init_size_local[2];
  counts[0] = 1;
  counts[1] = lx;
  counts[2] = lz;
  counts_local[0] = lx;
  counts_local[1] = lz;
  // Do this later, after setting t0//  offset[0]=t0;
  offset[1] = x0;
  offset[2] = z0;
  offset_local[0] = x0_local;
  offset_local[1] = z0_local;
  init_size_local[0] = mesh->LocalNx;
  init_size_local[1] = mesh->LocalNz;

  if (nd_local == 0) {
    nd_local = 1;
    // Need to write a time-series of scalars
    counts_local[0] = 1;
    offset_local[0] = 0;
    init_size_local[0] = 1;
  }

  hid_t mem_space = H5Screate_simple(nd_local, init_size_local, init_size_local);
  if (mem_space < 0)
    throw BoutException("Failed to create mem_space");
  if (H5Sselect_hyperslab(mem_space, H5S_SELECT_SET, offset_local, /*stride=*/nullptr,
                          counts_local, /*block=*/nullptr) < 0)
    throw BoutException("Failed to select hyperslab");

  hid_t dataSet = H5Dopen(dataFile, name.c_str(), H5P_DEFAULT);
  if (dataSet >= 0) { // >=0 means file exists, so open. Else error.

    hsize_t dims[3] = {};
    hid_t dataSpace = H5Dget_space(dataSet);
    if (dataSpace < 0)
      throw BoutException("Failed to create dataSpace");
    if (H5Sget_simple_extent_dims(dataSpace, dims, /*maxdims=*/nullptr) < 0)
      throw BoutException("Failed to get dims");
    dims[0]+=1;
    if (t0 == -1) {
      // Want t0 to be last record
      t0 = dims[0]-1;
    }

    if (H5Dset_extent(dataSet, dims) < 0)
      throw BoutException("Failed to extend dataSet");

    if (H5Sclose(dataSpace) < 0)
      throw BoutException("Failed to close dataSpace");

  }
  else {
    output_error.write("ERROR: HDF5 variable '%s' has not been added to file '%s'\n", name.c_str(), fname);
    return false;
  }

  offset[0]=t0;

  hid_t dataSpace = H5Dget_space(dataSet);
  if (dataSpace < 0)
    throw BoutException("Failed to create dataSpace");
  if (H5Sselect_hyperslab(dataSpace, H5S_SELECT_SET, offset, /*stride=*/nullptr, counts,
                          /*block=*/nullptr) < 0)
    throw BoutException("Failed to select hyperslab");

  if (H5Dwrite(dataSet, mem_hdf5_type, mem_space, dataSpace, dataSet_plist, data) < 0)
    throw BoutException("Failed to write data");

  if (H5Sclose(mem_space) < 0)
    throw BoutException("Failed to close mem_space");
  if (H5Sclose(dataSpace) < 0)
    throw BoutException("Failed to close dataSpace");
  if (H5Dclose(dataSet) < 0)
    throw BoutException("Failed to close dataSet");

  return true;
}


/***************************************************************************
 * Attributes
 ***************************************************************************/

void H5Format::setAttribute(const std::string &varname, const std::string &attrname,
                         const std::string &text) {
  TRACE("H5Format::setAttribute(varname, attrname, string)");

  std::string existing_att;
  if (getAttribute(varname, attrname, existing_att)) {
    if (text != existing_att) {
      output_warn.write("Overwriting attribute '%s' of variable '%s' with '%s', was previously '%s'",
          attrname.c_str(), varname.c_str(), text.c_str(), existing_att.c_str());
    }
  }
  // else: attribute does not exist, so just write it

  if (varname == "") {
    // attribute of file
    setAttribute(dataFile, attrname, text);
  } else {
    // attribute of variable
    hid_t dataSet = H5Dopen(dataFile, varname.c_str(), H5P_DEFAULT);
    if (dataSet < 0) {
      // Negative value indicates error, i.e. variable does not exist
      throw BoutException("Trying to create attribute for variable that does not exist");
    }

    setAttribute(dataSet, attrname, text);

    if (H5Dclose(dataSet) < 0) {
      throw BoutException("Failed to close dataSet");
    }
  }
}

void H5Format::setAttribute(const std::string &varname, const std::string &attrname,
                         int value) {
  TRACE("H5Format::setAttribute(varname, attrname, int)");

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
    setAttribute(dataFile, attrname, value);
  } else {
    // attribute of variable
    hid_t dataSet = H5Dopen(dataFile, varname.c_str(), H5P_DEFAULT);
    if (dataSet < 0) {
      // Negative value indicates error, i.e. variable does not exist
      throw BoutException("Trying to create attribute for variable that does not exist");
    }

    setAttribute(dataSet, attrname, value);

    if (H5Dclose(dataSet) < 0) {
      throw BoutException("Failed to close dataSet");
    }
  }
}

void H5Format::setAttribute(const std::string &varname, const std::string &attrname,
                         BoutReal value) {
  TRACE("H5Format::setAttribute(varname, attrname, BoutReal)");

  BoutReal existing_att;
  if (getAttribute(varname, attrname, existing_att)) {
    if (value != existing_att) {
      output_warn.write("Overwriting attribute '%s' of variable '%s' with '%f', was previously '%f'",
          attrname.c_str(), varname.c_str(), value, existing_att);
    }
  }
  // else: attribute does not exist, so just write it

  if (varname == "") {
    // attribute of file
    setAttribute(dataFile, attrname, value);
  } else {
    // attribute of variable
    hid_t dataSet = H5Dopen(dataFile, varname.c_str(), H5P_DEFAULT);
    if (dataSet < 0) {
      // Negative value indicates error, i.e. variable does not exist
      throw BoutException("Trying to create attribute for variable that does not exist");
    }

    setAttribute(dataSet, attrname, value);

    if (H5Dclose(dataSet) < 0) {
      throw BoutException("Failed to close dataSet");
    }
  }
}

void H5Format::setAttribute(const hid_t &dataSet, const std::string &attrname,
                         const std::string &text) {
  TRACE("H5Format::setAttribute(dataSet, attrname, string)");

  // Create new dataspace for attribute
  hid_t attribute_dataspace = H5Screate(H5S_SCALAR);
  if (attribute_dataspace < 0)
    throw BoutException("Failed to create attribute_dataspace");

  // Create new string datatype for attribute
  hid_t variable_length_string_type = H5Tcopy(H5T_C_S1);
  if (variable_length_string_type < 0)
    throw BoutException("Failed to create variable_length_string_type");
  if (H5Tset_size(variable_length_string_type, H5T_VARIABLE) < 0)
    throw BoutException("Failed to create string type");

  // Create attribute if it does not exist and write to it
  hid_t myatt_in = H5Aopen(dataSet, attrname.c_str(), H5P_DEFAULT);
  if (myatt_in < 0) {
    // Need to create attribute
    myatt_in = H5Acreate(dataSet, attrname.c_str(), variable_length_string_type, attribute_dataspace, H5P_DEFAULT, H5P_DEFAULT);
    if (myatt_in < 0)
      throw BoutException("Failed to create attribute");
  }
  // Need to pass `const char**` to HDF5 for reasons, and
  // `&text.c_str()` isn't valid (can't take address of an
  // r-value/temporary), so we need an intermediate variable
  const char* c_text = text.c_str();
  if (H5Awrite(myatt_in, variable_length_string_type, &c_text) < 0)
    throw BoutException("Failed to write attribute");

  if (H5Sclose(attribute_dataspace) < 0)
    throw BoutException("Failed to close attribute_dataspace");
  if (H5Tclose(variable_length_string_type) < 0)
    throw BoutException("Failed to close variable_length_string_type");
  if (H5Aclose(myatt_in) < 0)
    throw BoutException("Failed to close myatt_in");
}

void H5Format::setAttribute(const hid_t &dataSet, const std::string &attrname,
                         int value) {
  TRACE("H5Format::setAttribute(dataSet, attrname, int)");

  // Create new dataspace for attribute
  hid_t attribute_dataspace = H5Screate(H5S_SCALAR);
  if (attribute_dataspace < 0)
    throw BoutException("Failed to create attribute_dataspace");

  // Create attribute and write to it
  hid_t myatt_in = H5Acreate(dataSet, attrname.c_str(), H5T_NATIVE_INT, attribute_dataspace, H5P_DEFAULT, H5P_DEFAULT);
  if (myatt_in < 0)
    throw BoutException("Failed to create attribute");
  if (H5Awrite(myatt_in, H5T_NATIVE_INT, &value) < 0)
    throw BoutException("Failed to write attribute");

  if (H5Sclose(attribute_dataspace) < 0)
    throw BoutException("Failed to close attribute_dataspace");
  if (H5Aclose(myatt_in) < 0)
    throw BoutException("Failed to close myatt_in");
}

void H5Format::setAttribute(const hid_t &dataSet, const std::string &attrname,
                         BoutReal value) {
  TRACE("H5Format::setAttribute(dataSet, attrname, BoutReal)");

  // Create new dataspace for attribute
  hid_t attribute_dataspace = H5Screate(H5S_SCALAR);
  if (attribute_dataspace < 0)
    throw BoutException("Failed to create attribute_dataspace");

  // Create attribute and write to it
  hid_t myatt_in = H5Acreate(dataSet, attrname.c_str(), H5T_NATIVE_DOUBLE, attribute_dataspace, H5P_DEFAULT, H5P_DEFAULT);
  if (myatt_in < 0)
    throw BoutException("Failed to create attribute");
  if (H5Awrite(myatt_in, H5T_NATIVE_DOUBLE, &value) < 0)
    throw BoutException("Failed to write attribute");

  if (H5Sclose(attribute_dataspace) < 0)
    throw BoutException("Failed to close attribute_dataspace");
  if (H5Aclose(myatt_in) < 0)
    throw BoutException("Failed to close myatt_in");
}

bool H5Format::getAttribute(const std::string &varname, const std::string &attrname, std::string &text) {
  TRACE("H5Format::getAttribute(varname, attrname, string)");

  if (varname == "") {
    // attribute of file
    return getAttribute(dataFile, attrname, text);
  } else {
    // attribute of variable
    hid_t dataSet = H5Dopen(dataFile, varname.c_str(), H5P_DEFAULT);
    if (dataSet < 0) {
      // Negative value indicates error, i.e. variable does not exist
      throw BoutException("Trying to read attribute for variable that does not exist");
    }

    bool result = getAttribute(dataSet, attrname, text);

    if (H5Dclose(dataSet) < 0)
      throw BoutException("Failed to close dataSet");

    return result;
  }
}

bool H5Format::getAttribute(const std::string &varname, const std::string &attrname, int &value) {
  TRACE("H5Format::getAttribute(varname, attrname, int)");

  if (varname == "") {
    // attribute of file
    return getAttribute(dataFile, attrname, value);
  } else {
    // attribute of variable
    hid_t dataSet = H5Dopen(dataFile, varname.c_str(), H5P_DEFAULT);
    if (dataSet < 0) {
      // Negative value indicates error, i.e. variable does not exist
      throw BoutException("Trying to read attribute for variable that does not exist");
    }

    bool result = getAttribute(dataSet, attrname, value);

    if (H5Dclose(dataSet) < 0)
      throw BoutException("Failed to close dataSet");

    return result;
  }
}

bool H5Format::getAttribute(const std::string &varname, const std::string &attrname, BoutReal &value) {
  TRACE("H5Format::getAttribute(varname, attrname, BoutReal)");

  if (varname == "") {
    // attribute of file
    return getAttribute(dataFile, attrname, value);
  } else {
    // attribute of variable
    hid_t dataSet = H5Dopen(dataFile, varname.c_str(), H5P_DEFAULT);
    if (dataSet < 0) {
      // Negative value indicates error, i.e. variable does not exist
      throw BoutException("Trying to read attribute for variable that does not exist");
    }

    bool result = getAttribute(dataSet, attrname, value);

    if (H5Dclose(dataSet) < 0)
      throw BoutException("Failed to close dataSet");

    return result;
  }
}

bool H5Format::getAttribute(const hid_t &dataSet, const std::string &attrname, std::string &text) {
  TRACE("H5Format::getAttribute(hid_t, attrname, string)");

  // Open attribute
  hid_t myatt = H5Aopen(dataSet, attrname.c_str(), H5P_DEFAULT);
  if (myatt < 0) {
    return false;
  }

  // Create new string datatype for attribute
  hid_t variable_length_string_type = H5Tcopy(H5T_C_S1);
  if (variable_length_string_type < 0)
    throw BoutException("Failed to create variable_length_string_type");
  if (H5Tset_size(variable_length_string_type, H5T_VARIABLE) < 0)
    throw BoutException("Failed to create string type");

  // Read attribute: Need to pass `char**` to HDF5 for reasons, but
  // luckliy it will allocate c_text for us
  char* c_text;
  if (H5Aread(myatt, variable_length_string_type, &c_text) < 0)
    throw BoutException("Failed to read attribute");
  text = c_text;
  // Release resources allocated by HDF5
  free(c_text);

  if (H5Tclose(variable_length_string_type) < 0)
    throw BoutException("Failed to close variable_length_string_type");
  if (H5Aclose(myatt) < 0)
    throw BoutException("Failed to close myatt_in");

  return true;
}

bool H5Format::getAttribute(const hid_t &dataSet, const std::string &attrname, int &value) {
  TRACE("H5Format::getAttribute(hid_t, attrname, int)");

  // Open attribute
  hid_t myatt = H5Aopen(dataSet, attrname.c_str(), H5P_DEFAULT);
  if (myatt < 0) {
    return false;
  }

  // Read attribute
  if (H5Aread(myatt, H5T_NATIVE_INT, &value) < 0)
    throw BoutException("Failed to read attribute");

  if (H5Aclose(myatt) < 0)
    throw BoutException("Failed to close myatt_in");

  return true;
}

bool H5Format::getAttribute(const hid_t &dataSet, const std::string &attrname, BoutReal &value) {
  TRACE("H5Format::getAttribute(hid_t, attrname, BoutReal)");

  // Open attribute
  hid_t myatt = H5Aopen(dataSet, attrname.c_str(), H5P_DEFAULT);
  if (myatt < 0) {
    return false;
  }

  // Read attribute
  if (H5Aread(myatt, H5T_NATIVE_DOUBLE, &value) < 0)
    throw BoutException("Failed to read attribute");

  if (H5Aclose(myatt) < 0)
    throw BoutException("Failed to close myatt_in");

  return true;
}

#endif // HDF5

