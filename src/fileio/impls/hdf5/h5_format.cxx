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

#include <bout/globals.hxx>
#include "h5_format.hxx"

#ifdef HDF5

#include <bout/utils.hxx>
#include <cmath>
#include <string>
#include <mpi.h>

#include <bout/output.hxx>
#include <bout/msg_stack.hxx>
#include <bout/boutcomm.hxx>

H5Format::H5Format(bool parallel_in) {
  parallel = parallel_in;
  x0 = y0 = z0 = t0 = 0;
  lowPrecision = false;
  fname = NULL;
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

  if (H5Eset_auto(H5E_DEFAULT, NULL, NULL) < 0) // Disable automatic printing of error messages so that we can catch errors without printing error messages to stdout
    throw BoutException("Failed to set error stack to not print errors");
}

H5Format::H5Format(const char *name, bool parallel_in) {
  parallel = parallel_in;
  x0 = y0 = z0 = t0 = 0;
  lowPrecision = false;
  fname = NULL;
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

  if (H5Eset_auto(H5E_DEFAULT, NULL, NULL) < 0) // Disable automatic printing of error messages so that we can catch errors without printing error messages to stdout
    throw BoutException("Failed to set error stack to not print errors");

  openr(name);
}

H5Format::~H5Format() {
  close();
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

bool H5Format::is_valid() {
  if(dataFile<0)
    return false;
  return true;
}

void H5Format::close() {
  TRACE("H5Format::close");
  
  if (is_valid()) {
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

const vector<int> H5Format::getSize(const char *name) {
  TRACE("H5Format::getSize");

  vector<int> size;

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
  }
  else {
    hsize_t* dims = new hsize_t[nd];
    int error = H5Sget_simple_extent_dims(dataSpace, dims, NULL);
    if (error < 0)
      throw BoutException("Failed to get dimensions of dataSpace");
    
    if (H5Sclose(dataSpace) < 0)
      throw BoutException("Failed to close dataSpace");
    if (H5Dclose(dataSet) < 0)
      throw BoutException("Failed to close dataSet");
    
    for (int i=0; i<nd; i++)
      size.push_back(dims[i]);

    delete[] dims;
  }

  return size;
}

const vector<int> H5Format::getSize(const string &var) {
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

bool H5Format::read(int *data, const char *name, int lx, int ly, int lz) {
  return read(data, H5T_NATIVE_INT, name, lx, ly, lz);
}

bool H5Format::read(int *var, const string &name, int lx, int ly, int lz) {
  return read(var, name.c_str(), lx, ly, lz);
}

bool H5Format::read(BoutReal *data, const char *name, int lx, int ly, int lz) {
  return read(data, H5T_NATIVE_DOUBLE, name, lx, ly, lz);
}

bool H5Format::read(BoutReal *var, const string &name, int lx, int ly, int lz) {
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
  offset_local[0]=x0_local;offset_local[1]=y0_local;offset_local[2]=z0_local;
  init_size_local[0]=offset_local[0]+counts[0]; init_size_local[1]=offset_local[1]+counts[1]; init_size_local[2]=offset_local[2]+counts[2]; // Want to be able to use without needing mesh to be initialised; makes hyperslab selection redundant
  
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
    if (H5Sselect_hyperslab(dataSpace, H5S_SELECT_SET, offset, /*stride=*/NULL, counts, /*block=*/NULL) < 0)
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
  return write(data, H5T_NATIVE_INT, H5T_NATIVE_INT, name, lx, ly, lz);
}

bool H5Format::write(int *var, const string &name, int lx, int ly, int lz) {
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
    
    return write(data, H5T_NATIVE_DOUBLE, H5T_NATIVE_FLOAT, name, lx, ly, lz);
  }
  else {
    return write(data, H5T_NATIVE_DOUBLE, H5T_NATIVE_DOUBLE, name, lx, ly, lz);
  }
  
}

bool H5Format::write(BoutReal *var, const string &name, int lx, int ly, int lz) {
  return write(var, name.c_str(), lx, ly, lz);
}

bool H5Format::write(void *data, hid_t mem_hdf5_type, hid_t write_hdf5_type, const char *name, int lx, int ly, int lz) {
  TRACE("H5Format::write(void)");

  if(!is_valid())
    return false;

  if((lx < 0) || (ly < 0) || (lz < 0))
    return false;

  int nd = 0; // Number of dimensions
  if(lx != 0) nd = 1;
  if(ly != 0) nd = 2;
  if(lz != 0) nd = 3;
  hsize_t counts[3],offset[3],offset_local[3],init_size[3],init_size_local[3];
  counts[0]=lx; counts[1]=ly; counts[2]=lz;
  offset[0]=x0; offset[1]=y0; offset[2]=z0;
  offset_local[0]=x0_local;offset_local[1]=y0_local;offset_local[2]=z0_local;
  if (parallel) {
    init_size[0]=mesh->GlobalNx-2*mesh->xstart; init_size[1]=mesh->GlobalNy-2*mesh->ystart; init_size[2]=mesh->GlobalNz;
  }
  else {
    init_size[0]=mesh->LocalNx; init_size[1]=mesh->LocalNy; init_size[2]=mesh->LocalNz;
  }
  init_size_local[0]=mesh->LocalNx; init_size_local[1]=mesh->LocalNy; init_size_local[2]=mesh->LocalNz;
  
  if (nd==0) {
    // Need to write a scalar, not a 0-d array
    nd = 1;
    counts[0] = 1;
    offset[0] = 0;
    offset_local[0] = 0;
    init_size[0] = 1;
    init_size_local[0] = 1;
  }
  
  hid_t mem_space = H5Screate_simple(nd, init_size_local, init_size_local);
  if (mem_space < 0)
    throw BoutException("Failed to create mem_space");
  if (H5Sselect_hyperslab(mem_space, H5S_SELECT_SET, offset_local, /*stride=*/NULL, counts, /*block=*/NULL) < 0)
    throw BoutException("Failed to select hyperslab");
  
  hid_t dataSet = H5Dopen(dataFile, name, H5P_DEFAULT);
  if (dataSet < 0) {
    // Negative value indicates error, i.e. file does not exist, so create:
    hid_t init_space = H5Screate_simple(nd, init_size, init_size);
    if (init_space < 0)
      throw BoutException("Failed to create init_space");
    dataSet = H5Dcreate(dataFile, name, write_hdf5_type, init_space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (dataSet < 0)
      throw BoutException("Failed to create dataSet");
    
    // Add attribute to say what kind of field this is
    std::string datatype = "scalar";
    if(lx != 0) datatype = "FieldX";
    if(ly != 0) datatype = "Field2D";
    if(lz != 0) datatype = "Field3D";
    
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
    hid_t myatt_in = H5Acreate(dataSet, "type", variable_length_string_type, attribute_dataspace, H5P_DEFAULT, H5P_DEFAULT);
    if (myatt_in < 0)
      throw BoutException("Failed to create attribute");
    if (H5Awrite(myatt_in, variable_length_string_type, &datatype) < 0)
      throw BoutException("Failed to write attribute");
    
    if (H5Sclose(init_space) < 0)
      throw BoutException("Failed to close init_space");
    if (H5Sclose(attribute_dataspace) < 0)
      throw BoutException("Failed to close attribute_dataspace");
    if (H5Tclose(variable_length_string_type) < 0)
      throw BoutException("Failed to close variable_length_string_type");
    if (H5Aclose(myatt_in) < 0)
      throw BoutException("Failed to close myatt_in");
    
  }
  
  hid_t dataSpace = H5Dget_space(dataSet);
  if (dataSpace < 0)
    throw BoutException("Failed to create dataSpace");
  if (H5Sselect_hyperslab(dataSpace, H5S_SELECT_SET, offset, /*stride=*/NULL, counts, /*block=*/NULL) < 0)
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
}/***************************************************************************
 * Record-based (time-dependent) data
 ***************************************************************************/

bool H5Format::read_rec(int *data, const char *name, int lx, int ly, int lz) {
  return read_rec(data, H5T_NATIVE_INT, name, lx, ly, lz);
}

bool H5Format::read_rec(int *var, const string &name, int lx, int ly, int lz) {
  return read_rec(var, name.c_str(), lx, ly, lz);
}

bool H5Format::read_rec(BoutReal *data, const char *name, int lx, int ly, int lz) {
  
  return read_rec(data, H5T_NATIVE_DOUBLE, name, lx, ly, lz);
  
}

bool H5Format::read_rec(BoutReal *var, const string &name, int lx, int ly, int lz) {
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
  if (H5Sselect_hyperslab(mem_space, H5S_SELECT_SET, offset_local, /*stride=*/NULL,
                          counts, /*block=*/NULL) < 0)
    throw BoutException("Failed to select hyperslab");

  hid_t dataSet = H5Dopen(dataFile, name, H5P_DEFAULT);
  if (dataSet < 0)
    throw BoutException("Failed to open dataSet");

  hid_t dataSpace = H5Dget_space(dataSet);
  if (dataSpace < 0)
    throw BoutException("Failed to create dataSpace");
  if (H5Sselect_hyperslab(dataSpace, H5S_SELECT_SET, offset, /*stride=*/NULL, counts,
                          /*block=*/NULL) < 0)
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
  return write_rec(data, H5T_NATIVE_INT, H5T_NATIVE_INT, name, lx, ly, lz);
}

bool H5Format::write_rec(int *var, const string &name, int lx, int ly, int lz) {
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
    return write_rec(data, H5T_NATIVE_DOUBLE, H5T_NATIVE_FLOAT, name, lx, ly, lz);
  }
  
  return write_rec(data, H5T_NATIVE_DOUBLE, H5T_NATIVE_DOUBLE, name, lx, ly, lz);
}

bool H5Format::write_rec(BoutReal *var, const string &name, int lx, int ly, int lz) {
  return write_rec(var, name.c_str(), lx, ly, lz);
}

bool H5Format::write_rec(void *data, hid_t mem_hdf5_type, hid_t write_hdf5_type, const char *name, int lx, int ly, int lz) {
  if(!is_valid())
    return false;

  if((lx < 0) || (ly < 0) || (lz < 0))
    return false;

  int nd = 1; // Number of dimensions
  if(lx != 0) nd = 2;
  if(ly != 0) nd = 3;
  if(lz != 0) nd = 4;
  int nd_local = nd-1;
  hsize_t counts[4],offset[4],init_size[4];
  hsize_t counts_local[3],offset_local[3],init_size_local[3];
  counts[0]=1; counts[1]=lx; counts[2]=ly; counts[3]=lz;
  counts_local[0]=lx; counts_local[1]=ly; counts_local[2]=lz;
// Do this later, after setting t0//  offset[0]=t0;
  offset[1]=x0; offset[2]=y0; offset[3]=z0;
  offset_local[0]=x0_local;offset_local[1]=y0_local;offset_local[2]=z0_local;
  if (parallel) {
    init_size[0]=1;init_size[1]=mesh->GlobalNx-2*mesh->xstart; init_size[2]=mesh->GlobalNy-2*mesh->ystart; init_size[3]=mesh->GlobalNz;
  }
  else {
    init_size[0]=1;init_size[1]=mesh->LocalNx; init_size[2]=mesh->LocalNy; init_size[3]=mesh->LocalNz;
  }
  init_size_local[0]=mesh->LocalNx; init_size_local[1]=mesh->LocalNy; init_size_local[2]=mesh->LocalNz;
  
  if (nd_local==0) {
    nd_local = 1;
    // Need to write a time-series of scalars
    counts_local[0] = 1;
    offset_local[0] = 0;
    init_size_local[0] = 1;
  }
  
  hid_t mem_space = H5Screate_simple(nd_local, init_size_local, init_size_local);
  if (mem_space < 0)
    throw BoutException("Failed to create mem_space");
  if (H5Sselect_hyperslab(mem_space, H5S_SELECT_SET, offset_local, /*stride=*/NULL, counts_local, /*block=*/NULL) < 0)
    throw BoutException("Failed to select hyperslab");
  
  hid_t dataSet = H5Dopen(dataFile, name, H5P_DEFAULT);
  if (dataSet >= 0) { // >=0 means file exists, so open. Else create it.
    
    hsize_t dims[4] = {};
    hid_t dataSpace = H5Dget_space(dataSet);
    if (dataSpace < 0)
      throw BoutException("Failed to create dataSpace");
    if (H5Sget_simple_extent_dims(dataSpace, dims, /*maxdims=*/NULL) < 0)
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
    dataSet = H5Dcreate(dataFile, name, write_hdf5_type, init_space, H5P_DEFAULT, propertyList, H5P_DEFAULT);
    if (dataSet < 0)
      throw BoutException("Failed to create dataSet");
    
    // Add attribute to say what kind of field this is
    std::string datatype = "scalar_t";
    if(lx != 0) datatype = "FieldX_t";
    if(ly != 0) datatype = "Field2D_t";
    if(lz != 0) datatype = "Field3D_t";
    
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
    hid_t myatt_in = H5Acreate(dataSet, "type", variable_length_string_type, attribute_dataspace, H5P_DEFAULT, H5P_DEFAULT);
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
    
    if (t0 == -1) {
      // Want t0 to be last record, we are creating the first record, so must have t0=0
      t0 = 0;
    }
    
  }

  offset[0]=t0;
  
  hid_t dataSpace = H5Dget_space(dataSet);
  if (dataSpace < 0)
    throw BoutException("Failed to create dataSpace");
  if (H5Sselect_hyperslab(dataSpace, H5S_SELECT_SET, offset, /*stride=*/NULL, counts, /*block=*/NULL) < 0)
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

#endif // HDF5

