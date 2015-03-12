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

#include <output.hxx>
#include <msg_stack.hxx>

H5Format::H5Format() {
  x0 = y0 = z0 = t0 = 0;
  lowPrecision = false;
  fname = NULL;
  dataFile = NULL;
  chunk_length = 10; // could change this to try to optimize IO performance (i.e. allocate new chunks of disk space less often)
  H5::Exception::dontPrint(); // Disable automatic printing of error messages so that we can catch exceptions without printing error messages to stdout
}

H5Format::H5Format(const char *name) {
  x0 = y0 = z0 = t0 = 0;
  lowPrecision = false;
  dataFile = NULL;
  chunk_length = 10; // could change this to try to optimize IO performance (i.e. allocate new chunks of disk space less often)
  H5::Exception::dontPrint(); // Disable automatic printing of error messages so that we can catch exceptions without printing error messages to stdout
  openr(name);
}

H5Format::H5Format(const string &name) {
  x0 = y0 = z0 = t0 = 0;
  lowPrecision = false;
  dataFile = NULL;
  chunk_length = 10; // could change this to try to optimize IO performance (i.e. allocate new chunks of disk space less often)
  H5::Exception::dontPrint(); // Disable automatic printing of error messages so that we can catch exceptions without printing error messages to stdout
  openr(name);
}

H5Format::~H5Format() {
  close();
}

bool H5Format::openr(const string &name) {
  return openr(name.c_str());
}

bool H5Format::openr(const char *name) {
#ifdef CHECK
  msg_stack.push("H5Format::openr");
#endif

  if(dataFile != NULL) // Already open. Close then re-open
    close();

  if(!is_valid()) {
    return false;
  }

  dataFile->openFile(name,H5F_ACC_RDONLY);
  
#ifdef CHECK
  msg_stack.pop();
#endif

  return true;
}

bool H5Format::openw(const string &name, bool append) {
  return openw(name.c_str(), append);
}

bool H5Format::openw(const char *name, bool append) {
#ifdef CHECK
  msg_stack.push("H5Format::openw");
#endif
  
  if (parallel) {
    // Need to do stuff here???
    throw BoutException("HDF5 does not support parallel file access in the C++ interface. Boo!");
  }

  if(dataFile != NULL) // Already open. Close then re-open
    close(); 

  if(append) {
    dataFile = new H5::H5File();
    dataFile->openFile(name, H5F_ACC_RDWR);
  }
  else {
    dataFile = new H5::H5File(name, H5F_ACC_TRUNC);
  }

  fname = copy_string(name);

#ifdef CHECK
  msg_stack.pop();
#endif

  return true;
}

bool H5Format::is_valid() {
  if(dataFile == NULL)
    return false;
// Should do some check on whether dataFile is a valid H5::H5File here, maybe using getID()???
  return true;
}

void H5Format::close() {
  #ifdef CHECK
  msg_stack.push("H5Format::close");
#endif
  
  if (is_valid()) {
    dataFile->close();
    delete dataFile;
    dataFile = NULL;
  }
  
#ifdef CHECK
  msg_stack.pop();
#endif
}

void H5Format::flush() {
  if(!is_valid())
    return;
  
  dataFile->flush(H5F_SCOPE_LOCAL);
}

const vector<int> H5Format::getSize(const char *name) {
  throw BoutException("H5Format::getSize not implemented");
  vector<int> size;
// 
// #ifdef NCDF_VERBOSE
//   NcError err(NcError::verbose_nonfatal);
// #else
//   NcError err(NcError::silent_nonfatal);
// #endif
// 
//   if(!is_valid())
//     return size;
// 
// #ifdef CHECK
//   msg_stack.push("H5Format::getSize");
// #endif
// 
//   NcVar *var;
//   
//   if(!(var = dataFile->get_var(name)))
//     return size;
//   
//   if(!var->is_valid())
//     return size;
// 
//   int nd = var->num_dims(); // Number of dimensions
// 
//   if(nd == 0) {
//     size.push_back(1);
//     return size;
//   }
//   
//   long *ls = var->edges();
//   for(int i=0;i<nd;i++)
//     size.push_back((int) ls[i]);
// 
//   delete[] ls;
//   
#ifdef CHECK
  msg_stack.pop();
#endif

  return size;
}

const vector<int> H5Format::getSize(const string &var) {
  return getSize(var.c_str());
}

bool H5Format::setGlobalOrigin(int x, int y, int z) {
  x0 = x;
  y0 = y;
  z0 = z;
  
  return true;
}

bool H5Format::setLocalOrigin(int x, int y, int z, int offset_x, int offset_y, int offset_z) {
  x0_local = offset_x;
  y0_local = offset_y;
  z0_local = offset_z;
  return setGlobalOrigin(x + mesh->OffsetX, y + mesh->OffsetY, z + mesh->OffsetZ);
}

bool H5Format::setRecord(int t) {
  t0 = t;

  return true;
}

bool H5Format::read(int *data, const char *name, int lx, int ly, int lz) {
  if(!is_valid())
    return false;

  if((lx < 0) || (ly < 0) || (lz < 0))
    return false;
  
//   // Check for valid name
//   checkName(name);

#ifdef CHECK
  msg_stack.push("H5Format::read(int)");
#endif
  
  int nd = 0; // Number of dimensions
  if(lx != 0) nd = 1;
  if(ly != 0) nd = 2;
  if(lz != 0) nd = 3;
  hsize_t counts[3],offset[3];
  counts[0]=lx; counts[1]=ly; counts[2]=lz;
  offset[0]=x0; offset[1]=y0; offset[2]=z0;
  
  H5::DataSpace mem_space(nd, counts);
  
  H5::DataSet dataSet = dataFile->openDataSet(name);
  
  H5::DataSpace dataSpace = dataSet.getSpace();
  dataSpace.selectHyperslab(H5S_SELECT_SET, counts, offset);
  
  dataSet.read(data, H5::PredType::NATIVE_INT, mem_space, dataSpace);
  
#ifdef CHECK
  msg_stack.pop();
#endif

  return true;
}

bool H5Format::read(int *var, const string &name, int lx, int ly, int lz) {
  return read(var, name.c_str(), lx, ly, lz);
}

bool H5Format::read(BoutReal *data, const char *name, int lx, int ly, int lz) {
  if(!is_valid())
    return false;

  if((lx < 0) || (ly < 0) || (lz < 0))
    return false;

#ifdef CHECK
  msg_stack.push("H5Format::read(BoutReal)");
#endif
  
  int nd = 0; // Number of dimensions
  if(lx != 0) nd = 1;
  if(ly != 0) nd = 2;
  if(lz != 0) nd = 3;
  hsize_t counts[3],offset[3];
  counts[0]=lx; counts[1]=ly; counts[2]=lz;
  offset[0]=x0; offset[1]=y0; offset[2]=z0;
  
  H5::DataSpace mem_space(nd, counts);
  
  H5::DataSet dataSet = dataFile->openDataSet(name);
  
  H5::DataSpace dataSpace = dataSet.getSpace();
  dataSpace.selectHyperslab(H5S_SELECT_SET, counts, offset);
  
  dataSet.read(data, H5::PredType::NATIVE_DOUBLE, mem_space, dataSpace);
  
#ifdef CHECK
  msg_stack.pop();
#endif

  return true;
}

bool H5Format::read(BoutReal *var, const string &name, int lx, int ly, int lz) {
  return read(var, name.c_str(), lx, ly, lz);
}

bool H5Format::write(int *data, const char *name, int lx, int ly, int lz) {
  if(!is_valid())
    return false;

  if((lx < 0) || (ly < 0) || (lz < 0))
    return false;
  
//   // Check for valid name
//   checkName(name);
  
#ifdef CHECK
  msg_stack.push("H5Format::write(int)");
#endif
  
  int nd = 0; // Number of dimensions
  if(lx != 0) nd = 1;
  if(ly != 0) nd = 2;
  if(lz != 0) nd = 3;
  hsize_t counts[3],offset[3], offset_local[3],init_size[3],init_size_local[3];
  counts[0]=lx; counts[1]=ly; counts[2]=lz;
  offset[0]=x0; offset[1]=y0; offset[2]=z0;
  offset_local[0]=x0_local;offset_local[1]=y0_local;offset_local[2]=z0_local;
  if (parallel) {
    init_size[0]=mesh->GlobalNx-2*mesh->mxg; init_size[1]=mesh->GlobalNy-2*mesh->myg; init_size[2]=mesh->GlobalNz;
  }
  else {
    init_size[0]=mesh->ngx; init_size[1]=mesh->ngy; init_size[2]=mesh->ngz;
  }
  init_size_local[0]=mesh->ngx; init_size_local[1]=mesh->ngy; init_size_local[2]=mesh->ngz;
  
  if (nd==0) {
    // Need to write a scalar, not a 0-d array
    nd = 1;
    counts[0] = 1;
    offset[0] = 0;
    init_size[0] = 1;
    init_size_local[0] = 1;
  }
  
  H5::DataSpace mem_space(nd, init_size_local);
  mem_space.selectHyperslab(H5S_SELECT_SET, counts, offset_local);
  
  H5::DataSet dataSet;
  try {
    dataSet = dataFile->openDataSet(name);
  }
  catch (H5::FileIException exception) {
    H5::DataSpace init_space(nd, init_size);
    dataSet = dataFile->createDataSet(name, H5::PredType::NATIVE_INT, init_space);
    
    // Add attribute to say what kind of field this is
    // Determine type of data
//     char* datatype; // type of the data
//     datatype = copy_string("scalar");
//     if(lx != 0) datatype = copy_string("FieldX");
//     if(ly != 0) copy_string("Field2D");
//     if(lz != 0) copy_string("Field3D");
    
    // Create attribute
//     /*H5::StrType char9(H5::PredType::NATIVE_CHAR, 9);
// //     hsize_t attribute_size[1];
// //     attribute_size[0] = 9;
// //     H5::DataSpace attribute_space = H5::DataSpace(1, attribute_size);
// //     H5::Attribute attribute = dataSet.createAttribute("type", H5::PredType::NATIVE_CHAR, attribute_space);
//     hsize_t attribute_size[1];
//     attribute_size[0] = 1;
//     H5::DataSpace attribute_space = H5::DataSpace(1, attribute_size);
//     H5::Attribute attribute = dataSet.createAttribute("type", char9, attribute_space);
//     
//     // Write attribute
// //     attribute.write(H5::PredType::NATIVE_CHAR, datatype);
//     attribute.write(char9, datatype);*/
//     delete [] datatype;
    
    H5std_string datatype; // type of the data
    datatype = "scalar";
    if(lx != 0) datatype = "FieldX";
    if(ly != 0) datatype = "Field2D";
    if(lz != 0) datatype = "Field3D";
    
    // Create new dataspace for attribute
    H5::DataSpace attribute_dataspace = H5::DataSpace(H5S_SCALAR);
    
    // Create new string datatype for attribute
    H5::StrType stringdatatype(H5::PredType::C_S1, 9); // of length 9 characters
    
    // Create attribute and write to it
    H5::Attribute myatt_in = dataSet.createAttribute("type", stringdatatype, attribute_dataspace);
    myatt_in.write(stringdatatype, datatype);
    
  }
  
  H5::DataSpace dataSpace = dataSet.getSpace();
  dataSpace.selectHyperslab(H5S_SELECT_SET, counts, offset);
  
  dataSet.write(data, H5::PredType::NATIVE_INT, mem_space, dataSpace);

#ifdef CHECK
  msg_stack.pop();
#endif

  return true;
}

bool H5Format::write(int *var, const string &name, int lx, int ly, int lz) {
  return write(var, name.c_str(), lx, ly, lz);
}

bool H5Format::write(BoutReal *data, const char *name, int lx, int ly, int lz) {
  if(!is_valid())
    return false;

  if((lx < 0) || (ly < 0) || (lz < 0))
    return false;

//   // Check for valid name
//   checkName(name);
  
#ifdef CHECK
  msg_stack.push("H5Format::write(BoutReal)");
#endif
  
  int nd = 0; // Number of dimensions
  if(lx != 0) nd = 1;
  if(ly != 0) nd = 2;
  if(lz != 0) nd = 3;
  hsize_t counts[3],offset[3],init_size[3],init_size_local[3];
  counts[0]=lx; counts[1]=ly; counts[2]=lz;
  offset[0]=x0; offset[1]=y0; offset[2]=z0;
  init_size[0]=mesh->ngx; init_size[1]=mesh->ngy; init_size[2]=mesh->ngz;
  init_size_local[0]=mesh->ngx; init_size_local[1]=mesh->ngy; init_size_local[2]=mesh->ngz;
  
  if (nd==0) {
    // Need to write a scalar, not a 0-d array
    nd = 1;
    counts[0] = 1;
    offset[0] = 0;
    init_size[0] = 1;
    init_size_local[0] = 1;
  }
  
  H5::DataSpace mem_space(nd, counts);
  
  
  H5::DataSet dataSet;
  try {
    dataSet = dataFile->openDataSet(name);
  }
  catch (H5::FileIException exception) {
    H5::DataSpace init_space(nd, init_size);
    if (lowPrecision) {
      dataSet = dataFile->createDataSet(name, H5::PredType::NATIVE_FLOAT, init_space);
      
      // Add attribute to say what kind of field this is
      // Determine type of data
//       char* datatype; // type of the data
//       datatype = copy_string("scalar");
//       if(lx != 0) datatype = copy_string("FieldX");
//       if(ly != 0) datatype = copy_string("Field2D");
//       if(lz != 0) datatype = copy_string("Field3D");
//       
//       // Create attribute
//       H5::StrType char9(H5::PredType::NATIVE_CHAR, 9);
// //       hsize_t attribute_size[1];
// //       attribute_size[0] = 9;
// //       H5::DataSpace attribute_space = H5::DataSpace(1, attribute_size);
// //       H5::Attribute attribute = dataSet.createAttribute("type", H5::PredType::NATIVE_CHAR, attribute_space);
//       hsize_t attribute_size[1];
//       attribute_size[0] = 1;
//       H5::DataSpace attribute_space = H5::DataSpace(1, attribute_size);
//       H5::Attribute attribute = dataSet.createAttribute("type", char9, attribute_space);
//       
//       // Write attribute
// //       attribute.write(H5::PredType::NATIVE_CHAR, datatype);
//       attribute.write(char9, datatype);
//       
//       delete [] datatype;
      
      H5std_string datatype; // type of the data
      datatype = "scalar";
      if(lx != 0) datatype = "FieldX";
      if(ly != 0) datatype = "Field2D";
      if(lz != 0) datatype = "Field3D";
      
      // Create new dataspace for attribute
      H5::DataSpace attribute_dataspace = H5::DataSpace(H5S_SCALAR);
      
      // Create new string datatype for attribute
      H5::StrType stringdatatype(H5::PredType::C_S1, 9); // of length 9 characters
      
      // Create attribute and write to it
      H5::Attribute myatt_in = dataSet.createAttribute("type", stringdatatype, attribute_dataspace);
      myatt_in.write(stringdatatype, datatype);
    }
    else {
      dataSet = dataFile->createDataSet(name, H5::PredType::NATIVE_DOUBLE, init_space);
      
      // Add attribute to say what kind of field this is
      // Determine type of data
//       char* datatype; // type of the data
//       datatype = copy_string("scalar");
//       if(lx != 0) datatype = copy_string("FieldX");
//       if(ly != 0) copy_string("Field2D");
//       if(lz != 0) copy_string("Field3D");
//       
      // Create attribute
//       H5::StrType char9(H5::PredType::NATIVE_CHAR, 9);
// //       hsize_t attribute_size[1];
// //       attribute_size[0] = 9;
// //       H5::DataSpace attribute_space = H5::DataSpace(1, attribute_size);
// //       H5::Attribute attribute = dataSet.createAttribute("type", H5::PredType::NATIVE_CHAR, attribute_space);
//       hsize_t attribute_size[1];
//       attribute_size[0] = 1;
//       H5::DataSpace attribute_space = H5::DataSpace(1, attribute_size);
//       H5::Attribute attribute = dataSet.createAttribute("type", char9, attribute_space);
//       
//       // Write attribute
// //       attribute.write(H5::PredType::NATIVE_CHAR, datatype);
//       attribute.write(char9, datatype);
//       
//       delete [] datatype;
      H5std_string datatype; // type of the data
      datatype = "scalar";
      if(lx != 0) datatype = "FieldX";
      if(ly != 0) datatype = "Field2D";
      if(lz != 0) datatype = "Field3D";
      
      // Create new dataspace for attribute
      H5::DataSpace attribute_dataspace = H5::DataSpace(H5S_SCALAR);
      
      // Create new string datatype for attribute
      H5::StrType stringdatatype(H5::PredType::C_S1, 9); // of length 9 characters
      
      // Create attribute and write to it
      H5::Attribute myatt_in = dataSet.createAttribute("type", stringdatatype, attribute_dataspace);
      myatt_in.write(stringdatatype, datatype);
    }
  }
  
  H5::DataSpace dataSpace = dataSet.getSpace();
  dataSpace.selectHyperslab(H5S_SELECT_SET, counts, offset);
  
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
  
  dataSet.write(data, H5::PredType::NATIVE_DOUBLE, mem_space, dataSpace);

#ifdef CHECK
  msg_stack.pop();
#endif

  return true;
}

bool H5Format::write(BoutReal *var, const string &name, int lx, int ly, int lz) {
  return write(var, name.c_str(), lx, ly, lz);
}

/***************************************************************************
 * Record-based (time-dependent) data
 ***************************************************************************/

bool H5Format::read_rec(int *data, const char *name, int lx, int ly, int lz) {
  if(!is_valid())
    return false;

  if((lx < 0) || (ly < 0) || (lz < 0))
    return false;

//   // Check for valid name
//   checkName(name);
  
  int nd = 1; // Number of dimensions
  if(lx != 0) nd = 2;
  if(ly != 0) nd = 3;
  if(lz != 0) nd = 4;
  hsize_t counts[4],offset[4], offset_local[4];
  hsize_t init_size[3],init_size_local[3];
  counts[0]=1; counts[1]=lx; counts[2]=ly; counts[3]=lz;
  offset[0]=t0; offset[1]=x0; offset[2]=y0; offset[3]=z0;
  offset_local[0]=t0;offset_local[1]=x0_local;offset_local[2]=y0_local;offset_local[3]=z0_local;
  if (parallel) {
    init_size[0]=mesh->GlobalNx-2*mesh->mxg; init_size[1]=mesh->GlobalNy-2*mesh->myg; init_size[2]=mesh->GlobalNz;
  }
  else {
    init_size[0]=mesh->ngx; init_size[1]=mesh->ngy; init_size[2]=mesh->ngz;
  }
  init_size_local[0]=mesh->ngx; init_size_local[1]=mesh->ngy; init_size_local[2]=mesh->ngz;
  
  if (nd==1) {
    // Need to write a time-series of scalars
    nd = 1;
    counts[1] = 1;
    offset[1] = 0;
    init_size[0] = 1;
    init_size_local[0] = 1;
  }
  
  H5::DataSpace mem_space(nd, counts);
  
  H5::DataSet dataSet = dataFile->openDataSet(name);
  
  H5::DataSpace dataSpace = dataSet.getSpace();
  dataSpace.selectHyperslab(H5S_SELECT_SET, counts, offset);
  
  dataSet.read(data, H5::PredType::NATIVE_INT, mem_space, dataSpace);
  
  return true;
}

bool H5Format::read_rec(int *var, const string &name, int lx, int ly, int lz) {
  return read_rec(var, name.c_str(), lx, ly, lz);
}

bool H5Format::read_rec(BoutReal *data, const char *name, int lx, int ly, int lz) {
  if(!is_valid())
    return false;

  if((lx < 0) || (ly < 0) || (lz < 0))
    return false;
  
//   // Check for valid name
//   checkName(name);
  
  int nd = 1; // Number of dimensions
  if(lx != 0) nd = 2;
  if(ly != 0) nd = 3;
  if(lz != 0) nd = 4;
  hsize_t counts[4],offset[4];
  counts[0]=1; counts[1]=lx; counts[2]=ly; counts[3]=lz;
  offset[0]=t0; offset[1]=x0; offset[2]=y0; offset[3]=z0;
  
  H5::DataSpace mem_space(nd, counts);
  
  H5::DataSet dataSet = dataFile->openDataSet(name);
  
  H5::DataSpace dataSpace = dataSet.getSpace();
  dataSpace.selectHyperslab(H5S_SELECT_SET, counts, offset);
  
  dataSet.write(data, H5::PredType::NATIVE_DOUBLE, mem_space, dataSpace);
  
  return true;
}

bool H5Format::read_rec(BoutReal *var, const string &name, int lx, int ly, int lz) {
  return read_rec(var, name.c_str(), lx, ly, lz);
}

bool H5Format::write_rec(int *data, const char *name, int lx, int ly, int lz) {
  if(!is_valid())
    return false;

  if((lx < 0) || (ly < 0) || (lz < 0))
    return false;

//   // Check for valid name
//   checkName(name);
  
  int nd = 1; // Number of dimensions
  if(lx != 0) nd = 2;
  if(ly != 0) nd = 3;
  if(lz != 0) nd = 4;
  hsize_t counts[4],offset[4], offset_local[4];
  hsize_t init_size[3],init_size_local[3];
  counts[0]=1; counts[1]=lx; counts[2]=ly; counts[3]=lz;
  offset[0]=t0; offset[1]=x0; offset[2]=y0; offset[3]=z0;
  offset_local[0]=t0;offset_local[1]=x0_local;offset_local[2]=y0_local;offset_local[3]=z0_local;
  if (parallel) {
    init_size[0]=mesh->GlobalNx-2*mesh->mxg; init_size[1]=mesh->GlobalNy-2*mesh->myg; init_size[2]=mesh->GlobalNz;
  }
  else {
    init_size[0]=mesh->ngx; init_size[1]=mesh->ngy; init_size[2]=mesh->ngz;
  }
  init_size_local[0]=mesh->ngx; init_size_local[1]=mesh->ngy; init_size_local[2]=mesh->ngz;
  
  if (nd==1) {
    // Need to write a time-series of scalars
    nd = 1;
    counts[1] = 1;
    offset[1] = 0;
    init_size[0] = 1;
    init_size_local[0] = 1;
  }
  
  H5::DataSpace mem_space(nd, counts);
  
  H5::DataSet dataSet;
  try {
    
    dataSet = dataFile->openDataSet(name);
    
    hsize_t dims[4] = {};
    H5::DataSpace currentSpace = dataSet.getSpace();
    currentSpace.getSimpleExtentDims(dims);
    dims[0]+=1;
    if (t0 == -1) {
      // Want t0 to be last record
      t0 = dims[0]-1;
    }
    
    dataSet.extend(dims);
    
  }
  catch (H5::FileIException exception) {
    
    hsize_t init_size[4];
    init_size[0]=1; init_size[1]=mesh->ngx; init_size[2]=mesh->ngy; init_size[3]=mesh->ngz;
    
    // Modify dataset creation properties, i.e. enable chunking.
    H5::DSetCreatPropList propertyList;
    hsize_t chunk_dims[4],max_dims[4];
    max_dims[0] = H5S_UNLIMITED; max_dims[1]=mesh->ngx; max_dims[2]=mesh->ngy; max_dims[3]=mesh->ngz;
    chunk_dims[0] = chunk_length; chunk_dims[1]=mesh->ngx; chunk_dims[2]=mesh->ngy; chunk_dims[3]=mesh->ngz;
    propertyList.setChunk(nd, chunk_dims);
    
    H5::DataSpace init_space(nd, init_size, max_dims);
    dataSet = dataFile->createDataSet(name, H5::PredType::NATIVE_INT, init_space, propertyList);
    
    // Add attribute to say what kind of field this is
    // Determine type of data
//     char* datatype; // type of the data
//     datatype = copy_string("scalar_t");
//     if(lx != 0) datatype = copy_string("FieldX_t");
//     if(ly != 0) datatype = copy_string("Field2D_t");
//     if(lz != 0) datatype = copy_string("Field3D_t");
//     
//     // Create attribute
//     H5::StrType char9(H5::PredType::NATIVE_CHAR, 9);
// //       hsize_t attribute_size[1];
// //       attribute_size[0] = 9;
// //       H5::DataSpace attribute_space = H5::DataSpace(1, attribute_size);
// //       H5::Attribute attribute = dataSet.createAttribute("type", H5::PredType::NATIVE_CHAR, attribute_space);
//     hsize_t attribute_size[1];
//     attribute_size[0] = 1;
//     H5::DataSpace attribute_space = H5::DataSpace(1, attribute_size);
//     H5::Attribute attribute = dataSet.createAttribute("type", char9, attribute_space);
//     
//     // Write attribute
// //       attribute.write(H5::PredType::NATIVE_CHAR, datatype);
//     attribute.write(char9, datatype);
//     
//     delete [] datatype;
    H5std_string datatype; // type of the data
    datatype = "scalar_t";
    if(lx != 0) datatype = "FieldX_t";
    if(ly != 0) datatype = "Field2D_t";
    if(lz != 0) datatype = "Field3D_t";
    
    // Create new dataspace for attribute
    H5::DataSpace attribute_dataspace = H5::DataSpace(H5S_SCALAR);
    
    // Create new string datatype for attribute
    H5::StrType stringdatatype(H5::PredType::C_S1, 9); // of length 9 characters
    
    // Create attribute and write to it
    H5::Attribute myatt_in = dataSet.createAttribute("type", stringdatatype, attribute_dataspace);
    myatt_in.write(stringdatatype, datatype);
    
    if (t0 == -1) {
      // Want t0 to be last record, we are creating the first record, so must have t0=0
      t0 = 0;
    }
    
  }
//   catch (H5::GroupIException exception) {
//     
//     hsize_t init_size[4];
//     init_size[0]=1; init_size[1]=mesh->ngx; init_size[2]=mesh->ngy; init_size[3]=mesh->ngz;
//     
//     // Modify dataset creation properties, i.e. enable chunking.
//     H5::DSetCreatPropList propertyList;
//     hsize_t chunk_dims[4],max_dims[4];
//     chunk_dims[0] = H5S_UNLIMITED; max_dims[1]=mesh->ngx; max_dims[2]=mesh->ngy; max_dims[3]=mesh->ngz;
//     chunk_dims[0] = chunk_length; chunk_dims[1]=mesh->ngx; chunk_dims[2]=mesh->ngy; chunk_dims[3]=mesh->ngz;
//     propertyList.setChunk(nd, chunk_dims);
//     
//     H5::DataSpace init_space(nd, init_size, max_dims);
//     dataSet = dataFile->createDataSet(name, H5::PredType::NATIVE_INT, init_space, propertyList);
//     
//   }
  
  offset[0]=t0; offset[1]=x0; offset[2]=y0; offset[3]=z0;
  
  H5::DataSpace dataSpace = dataSet.getSpace();
  dataSpace.selectHyperslab(H5S_SELECT_SET, counts, offset);
  
  dataSet.write(data, H5::PredType::NATIVE_INT, mem_space, dataSpace);
  
  return true;
}

bool H5Format::write_rec(int *var, const string &name, int lx, int ly, int lz) {
  return write_rec(var, name.c_str(), lx, ly, lz);
}

bool H5Format::write_rec(BoutReal *data, const char *name, int lx, int ly, int lz) {
  if(!is_valid())
    return false;

  if((lx < 0) || (ly < 0) || (lz < 0))
    return false;

//   // Check for valid name
//   checkName(name);
  
  int nd = 1; // Number of dimensions
  if(lx != 0) nd = 2;
  if(ly != 0) nd = 3;
  if(lz != 0) nd = 4;
  hsize_t counts[4],offset[4];
  counts[0]=1; counts[1]=lx; counts[2]=ly; counts[3]=lz;
  
  H5::DataSpace mem_space(nd, counts);
  
  H5::DataSet dataSet;
  try {
    
    dataSet = dataFile->openDataSet(name);
    
    hsize_t dims[4] = {};
    H5::DataSpace currentSpace = dataSet.getSpace();
    currentSpace.getSimpleExtentDims(dims);
    dims[0]+=1;
    if (t0 == -1) {
      // Want t0 to be last record
      t0 = dims[0]-1;
    }
    
    dataSet.extend(dims);
    
  }
  catch (H5::FileIException exception) {
    
    hsize_t init_size[4];
    init_size[0]=1; init_size[1]=mesh->ngx; init_size[2]=mesh->ngy; init_size[3]=mesh->ngz;
    
    // Modify dataset creation properties, i.e. enable chunking.
    H5::DSetCreatPropList propertyList;
    hsize_t chunk_dims[4],max_dims[4];
    max_dims[0] = H5S_UNLIMITED; max_dims[1]=mesh->ngx; max_dims[2]=mesh->ngy; max_dims[3]=mesh->ngz;
    chunk_dims[0] = chunk_length; chunk_dims[1]=mesh->ngx; chunk_dims[2]=mesh->ngy; chunk_dims[3]=mesh->ngz;
    propertyList.setChunk(nd, chunk_dims);
    
    H5::DataSpace init_space(nd, init_size, max_dims);
    if (lowPrecision) {
      dataSet = dataFile->createDataSet(name, H5::PredType::NATIVE_FLOAT, init_space, propertyList);
      
      // Add attribute to say what kind of field this is
      // Determine type of data
//       char* datatype; // type of the data
//       datatype = copy_string("scalar_t");
//       if(lx != 0) datatype = copy_string("FieldX_t");
//       if(ly != 0) datatype = copy_string("Field2D_t");
//       if(lz != 0) datatype = copy_string("Field3D_t");
//       
//       // Create attribute
//       H5::StrType char9(H5::PredType::NATIVE_CHAR, 9);
// //       hsize_t attribute_size[1];
// //       attribute_size[0] = 9;
// //       H5::DataSpace attribute_space = H5::DataSpace(1, attribute_size);
// //       H5::Attribute attribute = dataSet.createAttribute("type", H5::PredType::NATIVE_CHAR, attribute_space);
//       hsize_t attribute_size[1];
//       attribute_size[0] = 1;
//       H5::DataSpace attribute_space = H5::DataSpace(1, attribute_size);
//       H5::Attribute attribute = dataSet.createAttribute("type", char9, attribute_space);
//       
//       // Write attribute
// //       attribute.write(H5::PredType::NATIVE_CHAR, datatype);
//       attribute.write(char9, datatype);
// 
//       
//       delete [] datatype;
      H5std_string datatype; // type of the data
      datatype = "scalar_t";
      if(lx != 0) datatype = "FieldX_t";
      if(ly != 0) datatype = "Field2D_t";
      if(lz != 0) datatype = "Field3D_t";
      
      // Create new dataspace for attribute
      H5::DataSpace attribute_dataspace = H5::DataSpace(H5S_SCALAR);
      
      // Create new string datatype for attribute
      H5::StrType stringdatatype(H5::PredType::C_S1, 9); // of length 9 characters
      
      // Create attribute and write to it
      H5::Attribute myatt_in = dataSet.createAttribute("type", stringdatatype, attribute_dataspace);
      myatt_in.write(stringdatatype, datatype);
      
    }
    else {
      dataSet = dataFile->createDataSet(name, H5::PredType::NATIVE_DOUBLE, init_space, propertyList);
      
      // Add attribute to say what kind of field this is
      // Determine type of data
//       char* datatype; // type of the data
//       datatype = copy_string("scalar_t");
//       if(lx != 0) datatype = copy_string("FieldX_t");
//       if(ly != 0) datatype = copy_string("Field2D_t");
//       if(lz != 0) datatype = copy_string("Field3D_t");
//       
//       // Create attribute
//       H5::StrType char9(H5::PredType::NATIVE_CHAR, 9);
// //       hsize_t attribute_size[1];
// //       attribute_size[0] = 9;
// //       H5::DataSpace attribute_space = H5::DataSpace(1, attribute_size);
// //       H5::Attribute attribute = dataSet.createAttribute("type", H5::PredType::NATIVE_CHAR, attribute_space);
//       hsize_t attribute_size[1];
//       attribute_size[0] = 1;
//       H5::DataSpace attribute_space = H5::DataSpace(1, attribute_size);
//       H5::Attribute attribute = dataSet.createAttribute("type", char9, attribute_space);
//       
//       // Write attribute
// //       attribute.write(H5::PredType::NATIVE_CHAR, datatype);
//       attribute.write(char9, datatype);
//       
//       delete [] datatype;
      H5std_string datatype; // type of the data
      datatype = "scalar_t";
      if(lx != 0) datatype = "FieldX_t";
      if(ly != 0) datatype = "Field2D_t";
      if(lz != 0) datatype = "Field3D_t";
      
      // Create new dataspace for attribute
      H5::DataSpace attribute_dataspace = H5::DataSpace(H5S_SCALAR);
      
      // Create new string datatype for attribute
      H5::StrType stringdatatype(H5::PredType::C_S1, 9); // of length 9 characters
      
      // Create attribute and write to it
      H5::Attribute myatt_in = dataSet.createAttribute("type", stringdatatype, attribute_dataspace);
      myatt_in.write(stringdatatype, datatype);
      
    }
    
    if (t0 == -1) {
      // Want t0 to be last record, we are creating the first record, so must have t0=0
      t0 = 0;
    }
    
  }
//   catch (H5::GroupIException exception) {
//     
//     hsize_t init_size[4];
//     init_size[0]=1; init_size[1]=mesh->ngx; init_size[2]=mesh->ngy; init_size[3]=mesh->ngz;
//     
//     // Modify dataset creation properties, i.e. enable chunking.
//     H5::DSetCreatPropList propertyList;
//     hsize_t chunk_dims[4],max_dims[4];
//     chunk_dims[0] = H5S_UNLIMITED; max_dims[1]=mesh->ngx; max_dims[2]=mesh->ngy; max_dims[3]=mesh->ngz;
//     chunk_dims[0] = chunk_length; chunk_dims[1]=mesh->ngx; chunk_dims[2]=mesh->ngy; chunk_dims[3]=mesh->ngz;
//     propertyList.setChunk(nd, chunk_dims);
//     
//     H5::DataSpace init_space(nd, init_size, max_dims);
//     dataSet = dataFile->createDataSet(name, H5::PredType::NATIVE_INT, init_space, propertyList);
//     
//   }
  
  offset[0]=t0; offset[1]=x0; offset[2]=y0; offset[3]=z0;
  
  H5::DataSpace dataSpace = dataSet.getSpace();
  dataSpace.selectHyperslab(H5S_SELECT_SET, counts, offset);
  
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
  
  dataSet.write(data, H5::PredType::NATIVE_DOUBLE, mem_space, dataSpace);
  
  return true;
}

bool H5Format::write_rec(BoutReal *var, const string &name, int lx, int ly, int lz) {
  return write_rec(var, name.c_str(), lx, ly, lz);
}

/***************************************************************************
 * Private functions
 ***************************************************************************/

// void H5Format::checkName(const char* name) {
//   // Check if this name contains an invalid character
//   
//   const char* c = name;
//   while(*c != 0) {
//     if(*c == '*')
//       throw BoutException("Invalid character (*) in NetCDF variable name '%s'", name);
//     c++;
//   }
// }

#endif // HDF5

