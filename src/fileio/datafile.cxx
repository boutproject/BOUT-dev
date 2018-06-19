/*********************************************************
 * Data file class
 * 
 * Changelog:
 *
 * 2009-08 Ben Dudson <bd512@york.ac.uk>
 *    * Added enable and wtime global properties
 *
 * 2009-07 Ben Dudson <bd512@york.ac.uk>
 *    * Adapted to handle PDB and/or NetCDF files
 *
 * 2009-04 Ben Dudson <bd512@york.ac.uk>
 *    * Initial version
 *
 **************************************************************************
 * Copyright 2010 B.D.Dudson, S.Farley, M.V.Umansky, X.Q.Xu
 *
 * Contact Ben Dudson, bd512@york.ac.uk
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
 *********************************************************/
//#include "mpi.h" // For MPI_Wtime()

#include <globals.hxx>
#include <bout/sys/timer.hxx>
#include <datafile.hxx>
#include <boutexception.hxx>
#include <output.hxx>
#include <boutcomm.hxx>
#include <utils.hxx>
#include <msg_stack.hxx>
#include <cstring>
#include "formatfactory.hxx"

Datafile::Datafile(Options *opt) : parallel(false), flush(true), guards(true), floats(false), openclose(true), enabled(true), shiftOutput(false), shiftInput(false), flushFrequencyCounter(0), flushFrequency(1), file(nullptr) {
  filenamelen=FILENAMELEN;
  filename=new char[filenamelen];
  filename[0] = 0; // Terminate the string
  
  if(opt == NULL)
    return; // To allow static initialisation
  // Read options
  
  OPTION(opt, parallel, false); // By default no parallel formats for now
  OPTION(opt, flush, true);     // Safer. Disable explicitly if required
  OPTION(opt, guards, true);    // Compatible with old behavior
  OPTION(opt, floats, false); // High precision by default
  OPTION(opt, openclose, true); // Open and close every write or read
  OPTION(opt, enabled, true);
  OPTION(opt, init_missing, false); // Initialise missing variables?
  OPTION(opt, shiftOutput, false); // Do we want to write 3D fields in shifted space?
  OPTION(opt, shiftInput, false); // Do we want to read 3D fields in shifted space?
  OPTION(opt, flushFrequency, 1); // How frequently do we flush the file
  
}

Datafile::Datafile(Datafile &&other) noexcept
    : parallel(other.parallel), flush(other.flush), guards(other.guards),
      floats(other.floats), openclose(other.openclose), Lx(other.Lx), Ly(other.Ly),
      Lz(other.Lz), enabled(other.enabled), shiftOutput(other.shiftOutput),
      shiftInput(other.shiftInput), flushFrequencyCounter(other.flushFrequencyCounter),
      flushFrequency(other.flushFrequency), file(std::move(other.file)),
      int_arr(std::move(other.int_arr)), BoutReal_arr(std::move(other.BoutReal_arr)),
      f2d_arr(std::move(other.f2d_arr)), f3d_arr(std::move(other.f3d_arr)),
      v2d_arr(std::move(other.v2d_arr)), v3d_arr(std::move(other.v3d_arr)) {
  filenamelen = other.filenamelen;
  filename = other.filename;
  other.filenamelen = 0;
  other.filename = nullptr;
  other.file = nullptr;
}

Datafile::Datafile(const Datafile &other) :
  parallel(other.parallel), flush(other.flush), guards(other.guards),
  floats(other.floats), openclose(other.openclose), Lx(other.Lx), Ly(other.Ly), Lz(other.Lz),
  enabled(other.enabled), shiftOutput(other.shiftOutput), shiftInput(other.shiftInput), flushFrequencyCounter(other.flushFrequencyCounter), flushFrequency(other.flushFrequency), 
  file(nullptr), int_arr(other.int_arr),
  BoutReal_arr(other.BoutReal_arr), f2d_arr(other.f2d_arr),
  f3d_arr(other.f3d_arr), v2d_arr(other.v2d_arr), v3d_arr(other.v3d_arr) {
  filenamelen=other.filenamelen;
  filename=new char[filenamelen];
  strncpy(filename,other.filename,filenamelen);
  // Same added variables, but the file not the same 
}

Datafile& Datafile::operator=(Datafile &&rhs) noexcept {
  parallel     = rhs.parallel;
  flush        = rhs.flush;
  guards       = rhs.guards;
  floats       = rhs.floats;
  openclose    = rhs.openclose;
  enabled      = rhs.enabled;
  init_missing = rhs.init_missing;
  shiftOutput  = rhs.shiftOutput;
  shiftInput   = rhs.shiftInput;
  flushFrequencyCounter = 0;
  flushFrequency = rhs.flushFrequency;
  file         = std::move(rhs.file);
  int_arr      = std::move(rhs.int_arr);
  BoutReal_arr = std::move(rhs.BoutReal_arr);
  f2d_arr      = std::move(rhs.f2d_arr);
  f3d_arr      = std::move(rhs.f3d_arr);
  v2d_arr      = std::move(rhs.v2d_arr);
  v3d_arr      = std::move(rhs.v3d_arr);
  if (filenamelen < rhs.filenamelen){
    delete[] filename;
    filenamelen=rhs.filenamelen;
    filename=new char[filenamelen];
  }
  strncpy(filename,rhs.filename,filenamelen);
  return *this;
}

Datafile::~Datafile() {
  if (filename != nullptr){
    delete[] filename;
    filename=nullptr;
    filenamelen=0;
  }
}

bool Datafile::openr(const char *format, ...) {
  if(format == (const char*) NULL) 
    throw BoutException("Datafile::open: No argument given for opening file!");

  bout_vsnprintf(filename,filenamelen, format);
  
  // Get the data format
  file = FormatFactory::getInstance()->createDataFormat(filename, parallel);
  
  if(!file)
    throw BoutException("Datafile::open: Factory failed to create a DataFormat!");
  
  // If parallel do not want to write ghost points, and it is easier then to ignore the boundary guard cells as well
  if (parallel) {
    file->setLocalOrigin(0, 0, 0, mesh->xstart, mesh->ystart, 0);
  }
  else {
    file->setGlobalOrigin(0,0,0);
  }
  
  if(!openclose) {
    // Open the file now. Otherwise defer until later
    int MYPE;
    MPI_Comm_rank(BoutComm::get(), &MYPE);
    if(!file->openr(filename, MYPE))
      throw BoutException("Datafile::open: Failed to open file!");
  }
  
  return true;
}

bool Datafile::openw(const char *format, ...) {
  if(!enabled)
    return true;
  
  if(format == (const char*) NULL)
    throw BoutException("Datafile::open: No argument given for opening file!");

  bout_vsnprintf(filename, filenamelen, format);
  
  // Get the data format
  file = FormatFactory::getInstance()->createDataFormat(filename, parallel);
  
  if(!file)
    throw BoutException("Datafile::open: Factory failed to create a DataFormat!");
  
  // If parallel do not want to write ghost points, and it is easier then to ignore the boundary guard cells as well
  if (parallel) {
    file->setLocalOrigin(0, 0, 0, mesh->xstart, mesh->ystart, 0);
    Lx = mesh->LocalNx-2*mesh->xstart;
    Ly = mesh->LocalNy-2*mesh->ystart;
    Lz = mesh->LocalNz;
  }
  else {
    file->setGlobalOrigin(0,0,0);
    Lx = mesh->LocalNx;
    Ly = mesh->LocalNy;
    Lz = mesh->LocalNz;
  }
  
  appending = false;
  if(!openclose) {
    // Open the file
    int MYPE;
    MPI_Comm_rank(BoutComm::get(), &MYPE);
    if(!file->openw(filename, MYPE, appending))
      throw BoutException("Datafile::open: Failed to open file!");
  }
  
  return true;
}

bool Datafile::opena(const char *format, ...) {
  if(!enabled)
    return true;
  
  if(format == (const char*) NULL)
    throw BoutException("Datafile::open: No argument given for opening file!");

  bout_vsnprintf(filename, filenamelen, format);

  // Get the data format
  file = FormatFactory::getInstance()->createDataFormat(filename, parallel);
  
  if(!file)
    throw BoutException("Datafile::open: Factory failed to create a DataFormat!");

  // If parallel do not want to write ghost points, and it is easier then to ignore the boundary guard cells as well
  if (parallel) {
    file->setLocalOrigin(0, 0, 0, mesh->xstart, mesh->ystart, 0);
    Lx = mesh->LocalNx-2*mesh->xstart;
    Ly = mesh->LocalNy-2*mesh->ystart;
    Lz = mesh->LocalNz;
  }
  else {
    file->setGlobalOrigin(0,0,0);
    Lx = mesh->LocalNx;
    Ly = mesh->LocalNy;
    Lz = mesh->LocalNz;
  }
  
  appending = true;
  if(!openclose) {
    // Open the file
    int MYPE;
    MPI_Comm_rank(BoutComm::get(), &MYPE);
    if(!file->openw(filename, MYPE, true))
      throw BoutException("Datafile::open: Failed to open file!");
  }
  return true;
}

bool Datafile::isValid() {
  if(!enabled)
    return true; // Pretend to be valid
  
  if(!file)
    return false;
  
  return file->is_valid();
}

void Datafile::close() {
  if(!file)
    return;
  if(!openclose)
    file->close();
  // free:
  file = nullptr;
}

void Datafile::setLowPrecision() {
  if(!enabled)
    return;
  floats = true;
  file->setLowPrecision();
}

void Datafile::add(int &i, const char *name, bool save_repeat) {
  TRACE("DataFile::add(int)");
  if (varAdded(string(name))) {
    // Check if it's the same variable
    if (&i == varPtr(string(name))) {
      output_warn.write("WARNING: variable '%s' added again to Datafile\n", name);
    } else {
      throw BoutException("Variable '%s' already added to Datafile", name);
    }
  }

  VarStr<int> d;

  d.ptr = &i;
  d.name = string(name);
  d.save_repeat = save_repeat;
  d.covar = false;
  
  int_arr.push_back(d);
}

void Datafile::add(BoutReal &r, const char *name, bool save_repeat) {
  TRACE("DataFile::add(BoutReal)");
  if (varAdded(string(name))) {
    // Check if it's the same variable
    if (&r == varPtr(string(name))) {
      output_warn.write("WARNING: variable '%s' added again to Datafile\n", name);
    } else {
      throw BoutException("Variable '%s' already added to Datafile", name);
    }
  }

  VarStr<BoutReal> d;

  d.ptr = &r;
  d.name = string(name);
  d.save_repeat = save_repeat;
  d.covar = false;
  
  BoutReal_arr.push_back(d);
}

void Datafile::add(Field2D &f, const char *name, bool save_repeat) {
  TRACE("DataFile::add(Field2D)");
  if (varAdded(string(name))) {
    // Check if it's the same variable
    if (&f == varPtr(string(name))) {
      output_warn.write("WARNING: variable '%s' added again to Datafile", name);
    } else {
      throw BoutException("Variable '%s' already added to Datafile", name);
    }
  }

  VarStr<Field2D> d;

  d.ptr = &f;
  d.name = string(name);
  d.save_repeat = save_repeat;
  d.covar = false;
  
  f2d_arr.push_back(d);
}

void Datafile::add(Field3D &f, const char *name, bool save_repeat) {
  TRACE("DataFile::add(Field3D)");
  if (varAdded(string(name))) {
    // Check if it's the same variable
    if (&f == varPtr(string(name))) {
      output_warn.write("WARNING: variable '%s' added again to Datafile\n", name);
    } else {
      throw BoutException("Variable '%s' already added to Datafile", name);
    }
  }

  VarStr<Field3D> d;

  d.ptr = &f;
  d.name = string(name);
  d.save_repeat = save_repeat;
  d.covar = false;
  
  f3d_arr.push_back(d);
}

void Datafile::add(Vector2D &f, const char *name, bool save_repeat) {
  TRACE("DataFile::add(Vector2D)");
  if (varAdded(string(name))) {
    // Check if it's the same variable
    if (&f == varPtr(string(name))) {
      output_warn.write("WARNING: variable '%s' added again to Datafile\n", name);
    } else {
      throw BoutException("Variable '%s' already added to Datafile", name);
    }
  }

  VarStr<Vector2D> d;

  d.ptr = &f;
  d.name = string(name);
  d.save_repeat = save_repeat;
  d.covar = f.covariant;

  v2d_arr.push_back(d);
}

void Datafile::add(Vector3D &f, const char *name, bool save_repeat) {
  TRACE("DataFile::add(Vector3D)");
  if (varAdded(string(name))) {
    // Check if it's the same variable
    if (&f == varPtr(string(name))) {
      output_warn.write("WARNING: variable '%s' added again to Datafile\n", name);
    } else {
      throw BoutException("Variable '%s' already added to Datafile", name);
    }
  }

  VarStr<Vector3D> d;

  d.ptr = &f;
  d.name = string(name);
  d.save_repeat = save_repeat;
  d.covar = f.covariant;

  v3d_arr.push_back(d);
}

bool Datafile::read() {
  Timer timer("io");  ///< Start timer. Stops when goes out of scope

  if(openclose) {
    // Open the file
    int MYPE;
    MPI_Comm_rank(BoutComm::get(), &MYPE);
    if(!file->openr(filename, MYPE))
      throw BoutException("Datafile::read: Failed to open file!");
  }
  
  if(!file->is_valid())
    throw BoutException("Datafile::read: File is not valid!");

  file->setRecord(-1); // Read the latest record

  // Read integers
  for(const auto& var : int_arr) {
    if(var.save_repeat) {
      if(!file->read_rec(var.ptr, var.name.c_str())) {
        if(!init_missing) {
          throw BoutException("Missing data for %s in input. Set init_missing=true to set to zero.", var.name.c_str());
        }
        output_warn.write("\tWARNING: Could not read integer %s. Setting to zero\n", var.name.c_str());
        *(var.ptr) = 0;
        continue;
      }
    } else {
      if(!file->read(var.ptr, var.name.c_str())) {
        if(!init_missing) {
          throw BoutException("Missing data for %s in input. Set init_missing=true to set to zero.", var.name.c_str());
        }
        output_warn.write("\tWARNING: Could not read integer %s. Setting to zero\n", var.name.c_str());
        *(var.ptr) = 0;
        continue;
      }
    }
  }

  // Read BoutReals
  for(const auto& var : BoutReal_arr) {
    if(var.save_repeat) {
      if(!file->read_rec(var.ptr, var.name)) {
        if(!init_missing) {
          throw BoutException("Missing data for %s in input. Set init_missing=true to set to zero.", var.name.c_str());
        }
        output_warn.write("\tWARNING: Could not read BoutReal %s. Setting to zero\n", var.name.c_str());
        *(var.ptr) = 0;
        continue;
      }
    } else {
      if(!file->read(var.ptr, var.name)) {
        if(!init_missing) {
          throw BoutException("Missing data for %s in input. Set init_missing=true to set to zero.", var.name.c_str());
        }
        output_warn.write("\tWARNING: Could not read BoutReal %s. Setting to zero\n", var.name.c_str());
        *(var.ptr) = 0;
        continue;
      }
    }
  }
  
  // Read 2D fields
  for(const auto& var : f2d_arr) {
    read_f2d(var.name, var.ptr, var.save_repeat);
  }

  // Read 3D fields
  for(const auto& var : f3d_arr) {
    read_f3d(var.name, var.ptr, var.save_repeat);
  }

  // 2D vectors
  for(const auto& var : v2d_arr) {
    if(var.covar) {
      // Reading covariant vector
      read_f2d(var.name+string("_x"), &(var.ptr->x), var.save_repeat);
      read_f2d(var.name+string("_y"), &(var.ptr->y), var.save_repeat);
      read_f2d(var.name+string("_z"), &(var.ptr->z), var.save_repeat);
    } else {
      read_f2d(var.name+string("x"), &(var.ptr->x), var.save_repeat);
      read_f2d(var.name+string("y"), &(var.ptr->y), var.save_repeat);
      read_f2d(var.name+string("z"), &(var.ptr->z), var.save_repeat);
    }

    var.ptr->covariant = var.covar;
  }

  // 3D vectors
  for(const auto& var : v3d_arr) {
    if(var.covar) {
      // Reading covariant vector
      read_f3d(var.name+string("_x"), &(var.ptr->x), var.save_repeat);
      read_f3d(var.name+string("_y"), &(var.ptr->y), var.save_repeat);
      read_f3d(var.name+string("_z"), &(var.ptr->z), var.save_repeat);
    } else {
      read_f3d(var.name+string("x"), &(var.ptr->x), var.save_repeat);
      read_f3d(var.name+string("y"), &(var.ptr->y), var.save_repeat);
      read_f3d(var.name+string("z"), &(var.ptr->z), var.save_repeat);
    }

    var.ptr->covariant = var.covar;
  }

  if(openclose) {
    // Close the file
    file->close();
  }

  return true;
}

bool Datafile::write() {
  if(!enabled)
    return true; // Just pretend it worked

  TRACE("Datafile::write()");

  if(!file)
    throw BoutException("Datafile::write: File is not valid!");

  if(openclose && (flushFrequencyCounter % flushFrequency == 0)) {
    // Open the file
    int MYPE;
    MPI_Comm_rank(BoutComm::get(), &MYPE);
    if(!file->openw(filename, MYPE, appending))
      throw BoutException("Datafile::write: Failed to open file!");
    appending = true;
    flushFrequencyCounter = 0;
  }
  
  if(!file->is_valid())
    throw BoutException("Datafile::open: File is not valid!");

  if(floats)
    file->setLowPrecision();
  
  Timer timer("io");
  
  file->setRecord(-1); // Latest record

  // Write integers
  for(const auto& var : int_arr) {
    write_int(var.name, var.ptr, var.save_repeat);
  }
  
  // Write BoutReals
  for(const auto& var : BoutReal_arr) {
    write_real(var.name, var.ptr, var.save_repeat);
  }

  // Write 2D fields
  for (const auto& var : f2d_arr) {
    write_f2d(var.name, var.ptr, var.save_repeat);

    // Add cell location
    file->setAttribute(var.name, "cell_location", CELL_LOC_STRING(var.ptr->getLocation()));
    
    // Add string attributes
    {
      auto it = attrib_string.find(var.name);
      if (it != attrib_string.end()) {
        // Some string attributes are set
        for (const auto &keyval : it->second) {
          // keyval->first is the attribute name; second is the value
          file->setAttribute(var.name, keyval.first, keyval.second);
        }
      }
    }
    // Add integer attributes
    {
      auto it = attrib_int.find(var.name);
      if (it != attrib_int.end()) {
        // Some string attributes are set
        for (const auto &keyval : it->second) {
          // keyval->first is the attribute name; second is the value
          file->setAttribute(var.name, keyval.first, keyval.second);
        }
      }
    }
  }

  // Write 3D fields
  for (const auto& var : f3d_arr) {
    write_f3d(var.name, var.ptr, var.save_repeat);

    // Add cell location
    file->setAttribute(var.name, "cell_location", CELL_LOC_STRING(var.ptr->getLocation()));
    
    // Add string attributes
    {
      auto it = attrib_string.find(var.name);
      if (it != attrib_string.end()) {
        // Some string attributes are set
        for (const auto& keyval : it->second) {
          // keyval->first is the attribute name; second is the value
          file->setAttribute(var.name, keyval.first, keyval.second);
        }
      }
    }
    // Add integer attributes
    {
      auto it = attrib_int.find(var.name);
      if (it != attrib_int.end()) {
        // Some string attributes are set
        for (const auto &keyval : it->second) {
          // keyval->first is the attribute name; second is the value
          file->setAttribute(var.name, keyval.first, keyval.second);
        }
      }
    }
  }
  
  // 2D vectors
  for(const auto& var : v2d_arr) {
    if(var.covar) {
      // Writing covariant vector
      Vector2D v  = *(var.ptr);
      v.toCovariant();
      
      write_f2d(var.name+string("_x"), &(v.x), var.save_repeat);
      write_f2d(var.name+string("_y"), &(v.y), var.save_repeat);
      write_f2d(var.name+string("_z"), &(v.z), var.save_repeat);

      // Add cell location
      file->setAttribute(var.name+string("_x"), "cell_location",
                         CELL_LOC_STRING(v.x.getLocation()));
      file->setAttribute(var.name+string("_y"), "cell_location",
                         CELL_LOC_STRING(v.y.getLocation()));
      file->setAttribute(var.name+string("_z"), "cell_location",
                         CELL_LOC_STRING(v.z.getLocation()));
    } else {
      // Writing contravariant vector
      Vector2D v  = *(var.ptr);
      v.toContravariant();
      
      write_f2d(var.name+string("x"), &(v.x), var.save_repeat);
      write_f2d(var.name+string("y"), &(v.y), var.save_repeat);
      write_f2d(var.name+string("z"), &(v.z), var.save_repeat);

      // Add cell location
      file->setAttribute(var.name+string("_x"), "cell_location",
                         CELL_LOC_STRING(v.x.getLocation()));
      file->setAttribute(var.name+string("_y"), "cell_location",
                         CELL_LOC_STRING(v.y.getLocation()));
      file->setAttribute(var.name+string("_z"), "cell_location",
                         CELL_LOC_STRING(v.z.getLocation()));
    }
  }

  // 3D vectors
  for(const auto& var : v3d_arr) {
    if(var.covar) {
      // Writing covariant vector
      Vector3D v  = *(var.ptr);
      v.toCovariant();
      
      write_f3d(var.name+string("_x"), &(v.x), var.save_repeat);
      write_f3d(var.name+string("_y"), &(v.y), var.save_repeat);
      write_f3d(var.name+string("_z"), &(v.z), var.save_repeat);

      // Add cell location
      file->setAttribute(var.name+string("_x"), "cell_location",
                         CELL_LOC_STRING(v.x.getLocation()));
      file->setAttribute(var.name+string("_y"), "cell_location",
                         CELL_LOC_STRING(v.y.getLocation()));
      file->setAttribute(var.name+string("_z"), "cell_location",
                         CELL_LOC_STRING(v.z.getLocation()));
    } else {
      // Writing contravariant vector
      Vector3D v  = *(var.ptr);
      v.toContravariant();
      
      write_f3d(var.name+string("x"), &(v.x), var.save_repeat);
      write_f3d(var.name+string("y"), &(v.y), var.save_repeat);
      write_f3d(var.name+string("z"), &(v.z), var.save_repeat);

      // Add cell location
      file->setAttribute(var.name+string("_x"), "cell_location",
                         CELL_LOC_STRING(v.x.getLocation()));
      file->setAttribute(var.name+string("_y"), "cell_location",
                         CELL_LOC_STRING(v.y.getLocation()));
      file->setAttribute(var.name+string("_z"), "cell_location",
                         CELL_LOC_STRING(v.z.getLocation()));
    }
  }
  
  if(openclose  && (flushFrequencyCounter+1 % flushFrequency == 0)){
    file->close();
  }
  flushFrequencyCounter++;
  return true;
}

bool Datafile::write(const char *format, ...) const {
  if(!enabled)
    return true;

  if(format == (const char*) NULL)
    throw BoutException("Datafile::write: No argument given!");

  int filenamelen=FILENAMELEN;
  char * filename=new char[filenamelen];

  bout_vsnprintf(filename, filenamelen, format);

  // Create a new datafile
  Datafile tmp(*this);

  tmp.openw(filename);
  bool ret = tmp.write();
  tmp.close();

  delete[] filename;

  return ret;
}

bool Datafile::writeVar(const int &i, const char *name) {
  // Should do this a better way...
  int *i2 = new int;
  *i2 = i;
  add(*i2, name);
  return true;
}

bool Datafile::writeVar(BoutReal r, const char *name) {
  BoutReal *r2 = new BoutReal;
  *r2 = r;
  add(*r2, name);
  return true;
}

/////////////////////////////////////////////////////////////

bool Datafile::read_f2d(const string &name, Field2D *f, bool save_repeat) {
  f->allocate();
  
  if(save_repeat) {
    if(!file->read_rec(&((*f)(0,0)), name, mesh->LocalNx, mesh->LocalNy)) {
      if(init_missing) {
        output_warn.write("\tWARNING: Could not read 2D field %s. Setting to zero\n", name.c_str());
        *f = 0.0;
      } else {
        throw BoutException("Missing 2D evolving field %s in input. Set init_missing=true to set to zero.", name.c_str());
      }
      return false;
    }
  }else {
    if(!file->read(&((*f)(0,0)), name, mesh->LocalNx, mesh->LocalNy)) {
      if(init_missing) {
        output_warn.write("\tWARNING: Could not read 2D field %s. Setting to zero\n", name.c_str());
        *f = 0.0;
      } else {
        throw BoutException("Missing 2D field %s in input. Set init_missing=true to set to zero.", name.c_str());
      }
      return false;
    }
  }
  
  return true;
}

bool Datafile::read_f3d(const string &name, Field3D *f, bool save_repeat) {
  f->allocate();
  
  if(save_repeat) {
    if(!file->read_rec(&((*f)(0,0,0)), name, mesh->LocalNx, mesh->LocalNy, mesh->LocalNz)) {
      if(init_missing) {
        output_warn.write("\tWARNING: Could not read 3D field %s. Setting to zero\n", name.c_str());
        *f = 0.0;
      }else {
        throw BoutException("Missing 3D evolving field %s in input. Set init_missing=true to set to zero.", name.c_str());
      }
      return false;
    }
  }else {
    if(!file->read(&((*f)(0,0,0)), name, mesh->LocalNx, mesh->LocalNy, mesh->LocalNz)) {
      if(init_missing) {
        output_warn.write("\tWARNING: Could not read 3D field %s. Setting to zero\n", name.c_str());
        *f = 0.0;
      }else {
        throw BoutException("Missing 3D field %s in input. Set init_missing=true to set to zero.", name.c_str());
      }
      return false;
    }
  }

  
  if (shiftInput) {
    // Input file is in field-aligned coordinates e.g. BOUT++ 3.x restart file
    *f = mesh->fromFieldAligned(*f);
  }
  
  return true;
}

bool Datafile::write_int(const string &name, int *f, bool save_repeat) {
  if(save_repeat) {
    return file->write_rec(f, name);
  }else {
    return file->write(f, name);
  }
}

bool Datafile::write_real(const string &name, BoutReal *f, bool save_repeat) {
  if(save_repeat) {
    return file->write_rec(f, name);
  }else {
    return file->write(f, name);
  }
}

bool Datafile::write_f2d(const string &name, Field2D *f, bool save_repeat) {
  if (!f->isAllocated()) {
    throw BoutException("Datafile::write_f2d: Field2D '%s' is not allocated!", name.c_str());
  }
  if (save_repeat) {
    if (!file->write_rec(&((*f)(0, 0)), name, mesh->LocalNx, mesh->LocalNy)) {
      throw BoutException("Datafile::write_f2d: Failed to write %s!", name.c_str());
    }
  } else {
    if (!file->write(&((*f)(0, 0)), name, mesh->LocalNx, mesh->LocalNy)) {
      throw BoutException("Datafile::write_f2d: Failed to write %s!", name.c_str());
    }
  }
  return true;
}

bool Datafile::write_f3d(const string &name, Field3D *f, bool save_repeat) {
  if (!f->isAllocated()) {
    throw BoutException("Datafile::write_f3d: Field3D '%s' is not allocated!", name.c_str());
  }

  //Deal with shifting the output
  Field3D f_out(f->getMesh());
  if(shiftOutput) {
    f_out = mesh->toFieldAligned(*f);
  }else {
    f_out = *f;
  }

  if(save_repeat) {
    return file->write_rec(&(f_out(0,0,0)), name, mesh->LocalNx, mesh->LocalNy, mesh->LocalNz);
  }else {
    return file->write(&(f_out(0,0,0)), name, mesh->LocalNx, mesh->LocalNy, mesh->LocalNz);
  }
}

bool Datafile::varAdded(const string &name) {
  for(const auto& var : int_arr ) {
    if(name == var.name)
      return true;
  }

  for(const auto& var : BoutReal_arr ) {
    if(name == var.name)
      return true;
  }

  for(const auto& var : f2d_arr ) {
    if(name == var.name)
      return true;
  }
  
  for(const auto& var : f3d_arr ) {
    if(name == var.name)
      return true;
  }
  
  for(const auto& var : v2d_arr ) {
    if(name == var.name)
      return true;
  }

  for(const auto& var : v3d_arr ) {
    if(name == var.name)
      return true;
  }
  return false;
}

void *Datafile::varPtr(const string &name) {
  for (const auto &var : int_arr) {
    if (name == var.name) {
      return static_cast<void *>(var.ptr);
    }
  }

  for (const auto &var : BoutReal_arr) {
    if (name == var.name) {
      return static_cast<void *>(var.ptr);
    }
  }

  for (const auto &var : f2d_arr) {
    if (name == var.name) {
      return static_cast<void *>(var.ptr);
    }
  }

  for (const auto &var : f3d_arr) {
    if (name == var.name) {
      return static_cast<void *>(var.ptr);
    }
  }

  for (const auto &var : v2d_arr) {
    if (name == var.name) {
      return static_cast<void *>(var.ptr);
    }
  }

  for (const auto &var : v3d_arr) {
    if (name == var.name) {
      return static_cast<void *>(var.ptr);
    }
  }
  return nullptr;
}
