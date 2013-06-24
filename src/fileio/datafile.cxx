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

#include "formatfactory.hxx"

Datafile::Datafile(Options *opt) : parallel(false), flush(true), guards(true), floats(false), openclose(true), enabled(true), file(NULL) {
  if(opt == NULL)
    return; // To allow static initialisation
  
  // Read options
  
  OPTION(opt, parallel, false); // By default no parallel formats for now
  OPTION(opt, flush, true);     // Safer. Disable explicitly if required
  OPTION(opt, guards, true);    // Compatible with old behavior
  OPTION(opt, floats, false); // High precision by default
  OPTION(opt, openclose, true); // Open and close every write or read
  OPTION(opt, enabled, true);
}

Datafile::Datafile(const Datafile &other) : parallel(other.parallel), flush(other.flush), guards(other.guards), 
                                            floats(other.floats), openclose(other.openclose), 
                                            enabled(other.enabled), file(NULL), int_arr(other.int_arr), 
                                            BoutReal_arr(other.BoutReal_arr), f2d_arr(other.f2d_arr), 
                                            f3d_arr(other.f3d_arr), v2d_arr(other.v2d_arr), v3d_arr(other.v3d_arr) {
  
  // Same added variables, but the file not the same 
}

Datafile& Datafile::operator=(const Datafile &rhs) {
  parallel     = rhs.parallel;
  flush        = rhs.flush;
  guards       = rhs.guards;
  floats     = rhs.floats;
  openclose    = rhs.openclose;
  enabled      = rhs.enabled;
  file         = NULL; // All values copied except this
  int_arr      = rhs.int_arr;
  BoutReal_arr = rhs.BoutReal_arr;
  f2d_arr      = rhs.f2d_arr;
  f3d_arr      = rhs.f3d_arr;
  v2d_arr      = rhs.v2d_arr;
  v3d_arr      = rhs.v3d_arr;
  return *this;
}

Datafile::~Datafile() {
  if(file != NULL)
    delete file;
}

bool Datafile::openr(const char *format, ...) {
  va_list ap;  // List of arguments
  if(format == (const char*) NULL)
    return 1;
  va_start(ap, format);
    vsprintf(filename, format, ap);
  va_end(ap);
  
  // Get the data format
  file = FormatFactory::getInstance()->createDataFormat(filename, parallel);
  
  if(!file)
    return false;
  
  if(!openclose) {
    // Open the file now. Otherwise defer until later
    int MYPE;
    MPI_Comm_rank(BoutComm::get(), &MYPE);
    if(!file->openr(filename, MYPE))
      return false;
  }
  
  return true;
}

bool Datafile::openw(const char *format, ...) {
  if(!enabled)
    return true;
  
  va_list ap;  // List of arguments
  if(format == (const char*) NULL)
    return 1;
  va_start(ap, format);
  vsprintf(filename, format, ap);
  va_end(ap);
  
  // Get the data format
  file = FormatFactory::getInstance()->createDataFormat(filename, parallel);
  
  if(!file)
    return false;
  
  appending = false;
  if(!openclose) {
    // Open the file
    int MYPE;
    MPI_Comm_rank(BoutComm::get(), &MYPE);
    if(!file->openw(filename, MYPE))
      return false;
  }
  
  return true;
}

bool Datafile::opena(const char *format, ...) {
  if(!enabled)
    return true;
  
  va_list ap;  // List of arguments
  if(format == (const char*) NULL)
    return 1;
  va_start(ap, format);
  vsprintf(filename, format, ap);
  va_end(ap);
  
  // Get the data format
  file = FormatFactory::getInstance()->createDataFormat(filename, parallel);
  
  if(!file)
    return false;

  appending = true;
  if(!openclose) {
    // Open the file
    int MYPE;
    MPI_Comm_rank(BoutComm::get(), &MYPE);
    if(!file->openw(filename, MYPE, true))
      return false;
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
}

void Datafile::setLowPrecision() {
  if(!enabled)
    return;
  floats = true;
  file->setLowPrecision();
}

void Datafile::add(int &i, const char *name, int grow) {
  if(varAdded(string(name)))
    throw BoutException("Variable '%s' already added to Datafile", name);

  VarStr<int> d;

  d.ptr = &i;
  d.name = string(name);
  d.grow = (grow > 0) ? true : false;
  
  int_arr.push_back(d);
}

void Datafile::add(BoutReal &r, const char *name, int grow) {
  if(varAdded(string(name)))
    throw BoutException("Variable '%s' already added to Datafile", name);
  
  VarStr<BoutReal> d;

  d.ptr = &r;
  d.name = string(name);
  d.grow = (grow > 0) ? true : false;
  
  BoutReal_arr.push_back(d);
}

void Datafile::add(Field2D &f, const char *name, int grow) {
  if(varAdded(string(name)))
    throw BoutException("Variable '%s' already added to Datafile", name);
  
  VarStr<Field2D> d;

  d.ptr = &f;
  d.name = string(name);
  d.grow = (grow > 0) ? true : false;
  
  f2d_arr.push_back(d);
}

void Datafile::add(Field3D &f, const char *name, int grow) {
  if(varAdded(string(name)))
    throw BoutException("Variable '%s' already added to Datafile", name);
  
  VarStr<Field3D> d;

  d.ptr = &f;
  d.name = string(name);
  d.grow = (grow > 0) ? true : false;
  
  f3d_arr.push_back(d);
}

void Datafile::add(Vector2D &f, const char *name, int grow) {
  if(varAdded(string(name)))
    throw BoutException("Variable '%s' already added to Datafile", name);
  
  VarStr<Vector2D> d;

  d.ptr = &f;
  d.name = string(name);
  d.grow = (grow > 0) ? true : false;
  d.covar = f.covariant;
  
  v2d_arr.push_back(d);
}

void Datafile::add(Vector3D &f, const char *name, int grow) {
  if(varAdded(string(name)))
    throw BoutException("Variable '%s' already added to Datafile", name);
  
  VarStr<Vector3D> d;

  d.ptr = &f;
  d.name = string(name);
  d.grow = (grow > 0) ? true : false;
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
      return false;
  }
  
  if(!file->is_valid())
    return false;  

  file->setRecord(-1); // Read the latest record

  // Read integers

  for(std::vector< VarStr<int> >::iterator it = int_arr.begin(); it != int_arr.end(); it++) {
    if(it->grow) {
      if(!file->read_rec(it->ptr, it->name.c_str())) {
	output.write("\tWARNING: Could not read integer %s. Setting to zero\n", it->name.c_str());
	*(it->ptr) = 0;
	continue;
      }
    }else {
      if(!file->read(it->ptr, it->name.c_str())) {
	output.write("\tWARNING: Could not read integer %s. Setting to zero\n", it->name.c_str());
	*(it->ptr) = 0;
	continue;
      }
    }
  }

  // Read BoutReals

  for(std::vector< VarStr<BoutReal> >::iterator it = BoutReal_arr.begin(); it != BoutReal_arr.end(); it++) {
    if(it->grow) {
      if(!file->read_rec(it->ptr, it->name)) {
	output.write("\tWARNING: Could not read BoutReal %s. Setting to zero\n", it->name.c_str());
	*(it->ptr) = 0;
	continue;
      }
    }else {
      if(!file->read(it->ptr, it->name)) {
	output.write("\tWARNING: Could not read BoutReal %s. Setting to zero\n", it->name.c_str());
	*(it->ptr) = 0;
	continue;
      }
    }
  }
  
  // Read 2D fields
  
  for(std::vector< VarStr<Field2D> >::iterator it = f2d_arr.begin(); it != f2d_arr.end(); it++) {
    read_f2d(it->name, it->ptr, it->grow);
  }

  // Read 3D fields
  
  for(std::vector< VarStr<Field3D> >::iterator it = f3d_arr.begin(); it != f3d_arr.end(); it++) {
    read_f3d(it->name, it->ptr, it->grow);
  }

  // 2D vectors
  
  for(std::vector< VarStr<Vector2D> >::iterator it = v2d_arr.begin(); it != v2d_arr.end(); it++) {
    if(it->covar) {
      // Reading covariant vector
      read_f2d(it->name+string("_x"), &(it->ptr->x), it->grow);
      read_f2d(it->name+string("_y"), &(it->ptr->y), it->grow);
      read_f2d(it->name+string("_z"), &(it->ptr->z), it->grow);
    }else {
      read_f2d(it->name+string("x"), &(it->ptr->x), it->grow);
      read_f2d(it->name+string("y"), &(it->ptr->y), it->grow);
      read_f2d(it->name+string("z"), &(it->ptr->z), it->grow);
    }

    it->ptr->covariant = it->covar;
  }

  // 3D vectors
  
  for(std::vector< VarStr<Vector3D> >::iterator it = v3d_arr.begin(); it != v3d_arr.end(); it++) {
    if(it->covar) {
      // Reading covariant vector
      read_f3d(it->name+string("_x"), &(it->ptr->x), it->grow);
      read_f3d(it->name+string("_y"), &(it->ptr->y), it->grow);
      read_f3d(it->name+string("_z"), &(it->ptr->z), it->grow);
    }else {
      read_f3d(it->name+string("x"), &(it->ptr->x), it->grow);
      read_f3d(it->name+string("y"), &(it->ptr->y), it->grow);
      read_f3d(it->name+string("z"), &(it->ptr->z), it->grow);
    }

    it->ptr->covariant = it->covar;
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
  
  if(!file)
    return false;
  
  if(openclose) {
    // Open the file
    int MYPE;
    MPI_Comm_rank(BoutComm::get(), &MYPE);
    if(!file->openw(filename, MYPE, appending))
      return false;
    appending = true;
  }

  
  if(!file->is_valid())
    return false;

  if(floats)
    file->setLowPrecision();
  
  Timer timer("io");
  
  file->setRecord(-1); // Latest record

  // Write integers
  for(std::vector< VarStr<int> >::iterator it = int_arr.begin(); it != int_arr.end(); it++) {
    if(it->grow) {
      file->write_rec(it->ptr, it->name);
    }else {
      file->write(it->ptr, it->name);
    }
  }
  
  // Write BoutReals
  for(std::vector< VarStr<BoutReal> >::iterator it = BoutReal_arr.begin(); it != BoutReal_arr.end(); it++) {
    if(it->grow) {
      file->write_rec(it->ptr, it->name);
    }else {
      file->write(it->ptr, it->name);
    }
  }

  // Write 2D fields
  
  for(std::vector< VarStr<Field2D> >::iterator it = f2d_arr.begin(); it != f2d_arr.end(); it++) {
    write_f2d(it->name, it->ptr, it->grow);
  }

  // Write 3D fields
  
  for(std::vector< VarStr<Field3D> >::iterator it = f3d_arr.begin(); it != f3d_arr.end(); it++) {
    write_f3d(it->name, it->ptr, it->grow);
  }
  
  // 2D vectors
  
  for(std::vector< VarStr<Vector2D> >::iterator it = v2d_arr.begin(); it != v2d_arr.end(); it++) {
    if(it->covar) {
      // Writing covariant vector
      Vector2D v  = *(it->ptr);
      v.toCovariant();
      
      write_f2d(it->name+string("_x"), &(v.x), it->grow);
      write_f2d(it->name+string("_y"), &(v.y), it->grow);
      write_f2d(it->name+string("_z"), &(v.z), it->grow);
    }else {
      // Writing contravariant vector
      Vector2D v  = *(it->ptr);
      v.toContravariant();
      
      write_f2d(it->name+string("x"), &(v.x), it->grow);
      write_f2d(it->name+string("y"), &(v.y), it->grow);
      write_f2d(it->name+string("z"), &(v.z), it->grow);
    }
  }

  // 3D vectors
  
  for(std::vector< VarStr<Vector3D> >::iterator it = v3d_arr.begin(); it != v3d_arr.end(); it++) {
    if(it->covar) {
      // Writing covariant vector
      Vector3D v  = *(it->ptr);
      v.toCovariant();
      
      write_f3d(it->name+string("_x"), &(v.x), it->grow);
      write_f3d(it->name+string("_y"), &(v.y), it->grow);
      write_f3d(it->name+string("_z"), &(v.z), it->grow);
    }else {
      // Writing contravariant vector
      Vector3D v  = *(it->ptr);
      v.toContravariant();
      
      write_f3d(it->name+string("x"), &(v.x), it->grow);
      write_f3d(it->name+string("y"), &(v.y), it->grow);
      write_f3d(it->name+string("z"), &(v.z), it->grow);
    }
  }
  
  if(openclose)
    file->close();

  return true;
}

bool Datafile::write(const char *format, ...) const {
  if(!enabled)
    return true;
  
  va_list ap;  // List of arguments
  if(format == (const char*) NULL)
    return false;
  char filename[512];
  va_start(ap, format);
  vsprintf(filename, format, ap);
  va_end(ap);

  // Create a new datafile
  Datafile tmp(*this);
  
  tmp.openw(filename);
  bool ret = tmp.write();
  tmp.close();
  
  return ret;
}

bool Datafile::writeVar(const int &i, const char *name) {
  // Should do this a better way...
  int *i2 = new int;
  *i2 = i;
  add(*i2, name);
  return true;
}

bool Datafile::writeVar(const BoutReal &r, const char *name) {
  BoutReal *r2 = new BoutReal;
  *r2 = r;
  add(*r2, name);
  return true;
}

/////////////////////////////////////////////////////////////

bool Datafile::read_f2d(const string &name, Field2D *f, bool grow) {
  f->allocate();
  
  if(grow) {
    if(!file->read_rec(*(f->getData()), name, mesh->ngx, mesh->ngy)) {
      output.write("\tWARNING: Could not read 2D field %s. Setting to zero\n", name.c_str());
      *f = 0.0;
      return false;
    }
  }else {
    if(!file->read(*(f->getData()), name, mesh->ngx, mesh->ngy)) {
      output.write("\tWARNING: Could not read 2D field %s. Setting to zero\n", name.c_str());
      *f = 0.0;
      return false;
    }
  }
  return true;
}

bool Datafile::read_f3d(const string &name, Field3D *f, bool grow) {
  f->allocate();
  
  if(grow) {
    if(!file->read_rec(**(f->getData()), name, mesh->ngx, mesh->ngy, mesh->ngz)) {
      output.write("\tWARNING: Could not read 3D field %s. Setting to zero\n", name.c_str());
      *f = 0.0;
      return false;
    }
  }else {
    if(!file->read(**(f->getData()), name, mesh->ngx, mesh->ngy, mesh->ngz)) {
      output.write("\tWARNING: Could not read 3D field %s. Setting to zero\n", name.c_str());
      *f = 0.0;
      return false;
    }
  }
  return true;
}

bool Datafile::write_f2d(const string &name, Field2D *f, bool grow) {
  if(!f->isAllocated())
    return false; // No data allocated
  
  if(grow) {
    return file->write_rec(*(f->getData()), name, mesh->ngx, mesh->ngy);
  }else {
    return file->write(*(f->getData()), name, mesh->ngx, mesh->ngy);
  }
}

bool Datafile::write_f3d(const string &name, Field3D *f, bool grow) {
  if(!f->isAllocated()) {
    //output << "Datafile: unallocated: " << name << endl;
    return false; // No data allocated
  }

  if(grow) {
    return file->write_rec(**(f->getData()), name, mesh->ngx, mesh->ngy, mesh->ngz);
  }else {
    return file->write(**(f->getData()), name, mesh->ngx, mesh->ngy, mesh->ngz);
  }
}

bool Datafile::varAdded(const string &name) {
  for(std::vector< VarStr<int> >::iterator it = int_arr.begin(); it != int_arr.end(); it++) {
    if(name == it->name)
      return true;
  }

  for(std::vector< VarStr<BoutReal> >::iterator it = BoutReal_arr.begin(); it != BoutReal_arr.end(); it++) {
    if(name == it->name)
      return true;
  }

  for(std::vector< VarStr<Field2D> >::iterator it = f2d_arr.begin(); it != f2d_arr.end(); it++) {
    if(name == it->name)
      return true;
  }
  
  for(std::vector< VarStr<Field3D> >::iterator it = f3d_arr.begin(); it != f3d_arr.end(); it++) {
    if(name == it->name)
      return true;
  }
  
  for(std::vector< VarStr<Vector2D> >::iterator it = v2d_arr.begin(); it != v2d_arr.end(); it++) {
    if(name == it->name)
      return true;
  }

  for(std::vector< VarStr<Vector3D> >::iterator it = v3d_arr.begin(); it != v3d_arr.end(); it++) {
    if(name == it->name)
      return true;
  }
  return false;
}
