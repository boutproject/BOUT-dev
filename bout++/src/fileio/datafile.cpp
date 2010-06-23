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
#include "mpi.h" // For MPI_Wtime()

#define DATAFILE_ORIGIN
#include "datafile.h"
#undef DATAFILE_ORIGIN

#include "globals.h"

#ifdef PDBF
#include "pdb_format.h"
#endif

#ifdef NCDF
#include "nc_format.h"
#endif

#include <string.h>

// Define a default file extension
#ifdef PDBF
char DEFAULT_FILE_EXT[] = "pdb";
#else
#ifdef NCDF
char DEFAULT_FILE_EXT[] = "nc";
#else

#error No file format available; aborting.

#endif // NCDF
#endif // PDBF

int match_string(const char *str, int n, const char **match)
{
  for(int i=0;i<n;i++)
    if(strcasecmp(str, match[i]) == 0)
      return i;
  return -1;
}

// Work out which data format to use for given filename
DataFormat *data_format(const char *filename)
{
  if(filename == NULL) {
    // Return default file format

#ifdef PDBF
    //output.write("\tUsing default format (PDB)\n");
    return new PdbFormat;
#else

#ifdef NCDF
    //output.write("\tUsing default format (NetCDF)\n");
    return new NcFormat;
#else

#error No file format available; aborting.

#endif // NCDF
#endif // PDBF
  }

  // Extract the file extension

  int len = strlen(filename);

  int ind = len-1;  
  while((ind != -1) && (filename[ind] != '.')) {
    ind--;
  }
  
  const char *s = filename + ind+1;

  // Match strings
  
#ifdef PDBF
  const char *pdb_match[] = {"pdb"};
  if(match_string(s, 1, pdb_match) != -1) {
    output.write("\tUsing PDB format for file '%s'\n", filename);
    return new PdbFormat;
  }
#endif

#ifdef NCDF
  const char *ncdf_match[] = {"cdl", "nc", "ncdf"};
  if(match_string(s, 3, ncdf_match) != -1) {
    output.write("\tUsing NetCDF format for file '%s'\n", filename);
    return new NcFormat;
  }
#endif

  output.write("\tFile extension not recognised for '%s'\n", filename);
  // Set to the default
  return data_format(NULL);
}

///////////////////////////////////////
// Global variables, shared between Datafile objects
bool Datafile::enabled = true;
real Datafile::wtime = 0.0;

Datafile::Datafile()
{
  low_prec = false;
  file = NULL;
  setFormat(data_format()); // Set default format
}

Datafile::Datafile(DataFormat *format)
{
  low_prec = false;
  file = NULL;
  setFormat(format);
}

Datafile::~Datafile()
{
  
}

void Datafile::setFormat(DataFormat *format)
{
  if(file != NULL)
    delete file;
  
  file = format;
  
  if(low_prec)
    file->setLowPrecision();
}

void Datafile::setLowPrecision()
{
  low_prec = true;
  file->setLowPrecision();
}

void Datafile::add(int &i, const char *name, int grow)
{
  VarStr<int> d;

  d.ptr = &i;
  d.name = string(name);
  d.grow = (grow > 0) ? true : false;
  
  int_arr.push_back(d);
}

void Datafile::add(real &r, const char *name, int grow)
{
  VarStr<real> d;

  d.ptr = &r;
  d.name = string(name);
  d.grow = (grow > 0) ? true : false;
  
  real_arr.push_back(d);
}

void Datafile::add(Field2D &f, const char *name, int grow)
{
  VarStr<Field2D> d;

  d.ptr = &f;
  d.name = string(name);
  d.grow = (grow > 0) ? true : false;
  
  f2d_arr.push_back(d);
}

void Datafile::add(Field3D &f, const char *name, int grow)
{
  VarStr<Field3D> d;

  d.ptr = &f;
  d.name = string(name);
  d.grow = (grow > 0) ? true : false;
  
  f3d_arr.push_back(d);
}

void Datafile::add(Vector2D &f, const char *name, int grow)
{
  VarStr<Vector2D> d;

  d.ptr = &f;
  d.name = string(name);
  d.grow = (grow > 0) ? true : false;
  d.covar = f.covariant;
  
  v2d_arr.push_back(d);
}

void Datafile::add(Vector3D &f, const char *name, int grow)
{
  VarStr<Vector3D> d;

  d.ptr = &f;
  d.name = string(name);
  d.grow = (grow > 0) ? true : false;
  d.covar = f.covariant;
  
  v3d_arr.push_back(d);
}

int Datafile::read(const char *format, ...)
{
  va_list ap;  // List of arguments
  
  if(format == (const char*) NULL)
    return 1;

  char filename[512];
  va_start(ap, format);
    vsprintf(filename, format, ap);
  va_end(ap);

  // Record starting time
  real tstart = MPI_Wtime();
  
  // Open the file
  
  if(!file->openr(filename))
    return 1;

  if(!file->is_valid())
    return 1;  

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

  // Read reals

  for(std::vector< VarStr<real> >::iterator it = real_arr.begin(); it != real_arr.end(); it++) {
    if(it->grow) {
      if(!file->read_rec(it->ptr, it->name)) {
	output.write("\tWARNING: Could not read real %s. Setting to zero\n", it->name.c_str());
	*(it->ptr) = 0;
	continue;
      }
    }else {
      if(!file->read(it->ptr, it->name)) {
	output.write("\tWARNING: Could not read real %s. Setting to zero\n", it->name.c_str());
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
  
  file->close();

  wtime += MPI_Wtime() - tstart;

  return 0;
}

int Datafile::write(const char *format, ...)
{
  va_list ap;  // List of arguments
  
  if(format == (const char*) NULL)
    return 1;

  char filename[512];
  va_start(ap, format);
    vsprintf(filename, format, ap);
  va_end(ap);
  
  if(write(string(filename), false))
    return 0;
  
  return 1;
}

int Datafile::append(const char *format, ...)
{
  va_list ap;  // List of arguments
  
  if(format == (const char*) NULL)
    return 1;
  
  char filename[512];
  va_start(ap, format);
    vsprintf(filename, format, ap);
  va_end(ap);
  
  if(write(string(filename), true))
    return 0;
  
  return 1;
}

bool Datafile::write(const string &filename, bool append)
{
  if(!enabled)
    return true; // Just pretend it worked
  
  // Record starting time
  real tstart = MPI_Wtime();

  if(!file->openw(filename, append))
    return false;

  if(!file->is_valid())
    return false;
  
  file->setRecord(-1); // Latest record

  // Write integers
  for(std::vector< VarStr<int> >::iterator it = int_arr.begin(); it != int_arr.end(); it++) {
    if(it->grow) {
      file->write_rec(it->ptr, it->name);
    }else {
      file->write(it->ptr, it->name);
    }
  }
  
  // Write reals
  for(std::vector< VarStr<real> >::iterator it = real_arr.begin(); it != real_arr.end(); it++) {
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

  file->close();

  wtime += MPI_Wtime() - tstart;

  return true;
}

/////////////////////////////////////////////////////////////

bool Datafile::read_f2d(const string &name, Field2D *f, bool grow)
{
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

bool Datafile::read_f3d(const string &name, Field3D *f, bool grow)
{
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

bool Datafile::write_f2d(const string &name, Field2D *f, bool grow)
{
  if(!f->isAllocated())
    return false; // No data allocated
  
  if(grow) {
    return file->write_rec(*(f->getData()), name, mesh->ngx, mesh->ngy);
  }else {
    return file->write(*(f->getData()), name, mesh->ngx, mesh->ngy);
  }
}

bool Datafile::write_f3d(const string &name, Field3D *f, bool grow)
{
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

