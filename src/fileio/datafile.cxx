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
#include <bout/mesh.hxx>
#include <bout/sys/timer.hxx>
#include "bout_types.hxx"
#include <datafile.hxx>
#include <boutexception.hxx>
#include <options.hxx>
#include <output.hxx>
#include <boutcomm.hxx>
#include <utils.hxx>
#include <msg_stack.hxx>
#include <cstring>
#include "formatfactory.hxx"

Datafile::Datafile(Options* opt, Mesh* mesh_in)
    : mesh(mesh_in == nullptr ? bout::globals::mesh : mesh_in), file(nullptr) {
  
  if (opt == nullptr) {
    return; // To allow static initialisation
  }
  // Read options
  
  OPTION(opt, parallel, false); // By default no parallel formats for now
  OPTION(opt, flush, true);     // Safer. Disable explicitly if required
  OPTION(opt, guards, true);    // Compatible with old behavior
  OPTION(opt, floats, false); // High precision by default
  OPTION(opt, openclose, true); // Open and close every write or read
  OPTION(opt, enabled, true);
  OPTION(opt, init_missing, false); // Initialise missing variables?
  OPTION(opt, shiftoutput, false); // Do we want to write 3D fields in shifted space?
  OPTION(opt, shiftinput, false); // Do we want to read 3D fields in shifted space?
  OPTION(opt, flushfrequency, 1); // How frequently do we flush the file
}

Datafile::Datafile(Datafile&& other) noexcept
    : mesh(other.mesh), parallel(other.parallel), flush(other.flush),
      guards(other.guards), floats(other.floats), openclose(other.openclose),
      Lx(other.Lx), Ly(other.Ly), Lz(other.Lz), enabled(other.enabled),
      init_missing(other.init_missing), shiftoutput(other.shiftoutput),
      shiftinput(other.shiftinput), flushFrequencyCounter(other.flushFrequencyCounter),
      flushfrequency(other.flushfrequency), file(std::move(other.file)),
      filename(std::move(other.filename)), writable(other.writable),
      appending(other.appending), first_time(other.first_time),
      int_arr(std::move(other.int_arr)), int_vec_arr(std::move(other.int_vec_arr)),
      string_arr(std::move(other.string_arr)),
      BoutReal_arr(std::move(other.BoutReal_arr)), bool_arr(std::move(other.bool_arr)),
      f2d_arr(std::move(other.f2d_arr)), f3d_arr(std::move(other.f3d_arr)),
      v2d_arr(std::move(other.v2d_arr)), v3d_arr(std::move(other.v3d_arr)) {
  other.file = nullptr;
}

Datafile::Datafile(const Datafile& other)
    : mesh(other.mesh), parallel(other.parallel), flush(other.flush),
      guards(other.guards), floats(other.floats), openclose(other.openclose),
      Lx(other.Lx), Ly(other.Ly), Lz(other.Lz), enabled(other.enabled),
      init_missing(other.init_missing), shiftoutput(other.shiftoutput),
      shiftinput(other.shiftinput), flushFrequencyCounter(other.flushFrequencyCounter),
      flushfrequency(other.flushfrequency), file(nullptr), filename(other.filename),
      writable(other.writable), appending(other.appending), first_time(other.first_time),
      int_arr(other.int_arr), int_vec_arr(other.int_vec_arr),
      string_arr(other.string_arr), BoutReal_arr(other.BoutReal_arr),
      bool_arr(other.bool_arr), f2d_arr(other.f2d_arr), f3d_arr(other.f3d_arr),
      v2d_arr(other.v2d_arr), v3d_arr(other.v3d_arr) {}

Datafile& Datafile::operator=(Datafile &&rhs) noexcept {
  mesh         = rhs.mesh;
  parallel     = rhs.parallel;
  flush        = rhs.flush;
  guards       = rhs.guards;
  floats       = rhs.floats;
  openclose    = rhs.openclose;
  enabled      = rhs.enabled;
  init_missing = rhs.init_missing;
  shiftoutput  = rhs.shiftoutput;
  shiftinput   = rhs.shiftinput;
  flushFrequencyCounter = 0;
  flushfrequency = rhs.flushfrequency;
  file         = std::move(rhs.file);
  filename     = std::move(rhs.filename);
  writable     = rhs.writable;
  appending    = rhs.appending;
  first_time   = rhs.first_time;
  int_arr      = std::move(rhs.int_arr);
  int_vec_arr  = std::move(rhs.int_vec_arr);
  string_arr   = std::move(rhs.string_arr);
  BoutReal_arr = std::move(rhs.BoutReal_arr);
  bool_arr     = std::move(rhs.bool_arr);
  f2d_arr      = std::move(rhs.f2d_arr);
  f3d_arr      = std::move(rhs.f3d_arr);
  v2d_arr      = std::move(rhs.v2d_arr);
  v3d_arr      = std::move(rhs.v3d_arr);
  return *this;
}

bool Datafile::openr(const std::string& filename_) {
  this->filename = filename_;

  if (filename.empty()) {
    throw BoutException("Datafile::open: No argument given for opening file!");
  }

  // Get the data format
  file = FormatFactory::getInstance()->createDataFormat(filename.c_str(), parallel);
  
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
    if(!file->openr(filename, BoutComm::rank()))
      throw BoutException("Datafile::open: Failed to open file {:s} for reading!",
                          filename);
  }

  writable = false;
  
  return true;
}

bool Datafile::openw(const std::string& filename_) {
  if(!enabled)
    return true;
  
  this->filename = filename_;

  if (filename.empty()) {
    throw BoutException("Datafile::open: No argument given for opening file!");
  }

  // Get the data format
  file = FormatFactory::getInstance()->createDataFormat(filename.c_str(), parallel, mesh);
  
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
  // Open the file
  if(!file->openw(filename, BoutComm::rank(), appending))
    throw BoutException("Datafile::open: Failed to open file {:s} for writing!",
                        filename);

  appending = true;
  first_time = true; // newly opened file, so write attributes when variables are first written

  // Add variables to file
  // Add integers
  for(const auto& var : int_arr) {
    if (!file->addVarInt(var.name, var.save_repeat)) {
      throw BoutException("Failed to add int variable {:s} to Datafile", var.name);
    }
  }

  // Add vectors of integers
  for(const auto& var : int_vec_arr) {
    if (!file->addVarIntVec(var.name, var.save_repeat, var.ptr->size())) {
      throw BoutException("Failed to add int vector variable {:s} to Datafile", var.name);
    }
  }
  
  // Add strings
  for(const auto& var : string_arr) {
    if (!file->addVarString(var.name, var.save_repeat, var.ptr->size())) {
      throw BoutException("Failed to add string variable {:s} to Datafile", var.name);
    }
  }

  // Add BoutReals
  for(const auto& var : BoutReal_arr) {
    if (!file->addVarBoutReal(var.name, var.save_repeat)) {
      throw BoutException("Failed to add BoutReal variable {:s} to Datafile", var.name);
    }
  }

  // Add bools
  for(const auto& var : bool_arr) {
    if (!file->addVarInt(var.name, var.save_repeat)) {
      throw BoutException("Failed to add bool variable {:s} to Datafile", var.name);
    }
  }

  // Add 2D fields
  for (const auto& var : f2d_arr) {
    if (!file->addVarField2D(var.name, var.save_repeat)) {
      throw BoutException("Failed to add Field2D variable {:s} to Datafile", var.name);
    }
  }

  // Add 3D fields
  for (const auto& var : f3d_arr) {
    if (!file->addVarField3D(var.name, var.save_repeat)) {
      throw BoutException("Failed to add Field3D variable {:s} to Datafile", var.name);
    }
  }

  // Add FieldPerps
  for (const auto& var : fperp_arr) {
    if (!file->addVarFieldPerp(var.name, var.save_repeat)) {
      throw BoutException("Failed to add FieldPerp variable {:s} to Datafile", var.name);
    }
  }

  // 2D vectors
  for(const auto& var : v2d_arr) {
    auto name = var.covar ? var.name + "_" : var.name;
    if (!file->addVarField2D(name + "x", var.save_repeat)) {
      throw BoutException("Failed to add Vector2D variable {:s} to Datafile", name);
    }
    if (!file->addVarField2D(name + "y", var.save_repeat)) {
      throw BoutException("Failed to add Vector2D variable {:s} to Datafile", name);
    }
    if (!file->addVarField2D(name + "z", var.save_repeat)) {
      throw BoutException("Failed to add Vector2D variable {:s} to Datafile", name);
    }
  }

  // 3D vectors
  for(const auto& var : v3d_arr) {
    auto name = var.covar ? var.name + "_" : var.name;
    if (!file->addVarField3D(name + "x", var.save_repeat)) {
      throw BoutException("Failed to add Vector3D variable {:s} to Datafile", name);
    }
    if (!file->addVarField3D(name + "y", var.save_repeat)) {
      throw BoutException("Failed to add Vector3D variable {:s} to Datafile", name);
    }
    if (!file->addVarField3D(name + "z", var.save_repeat)) {
      throw BoutException("Failed to add Vector3D variable {:s} to Datafile", name);
    }
  }
  if (openclose) {
    file->close();
  }

  writable = true;

  return true;
}

bool Datafile::opena(const std::string& filename_) {
  if(!enabled)
    return true;

  this->filename = filename_;

  if (filename.empty()) {
    throw BoutException("Datafile::open: No argument given for opening file!");
  }

  // Get the data format
  file = FormatFactory::getInstance()->createDataFormat(filename.c_str(), parallel);
  
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
  // Open the file
  if(!file->openw(filename, BoutComm::rank(), true))
    throw BoutException("Datafile::open: Failed to open file {:s} for appending!",
                        filename);

  first_time = true; // newly opened file, so write attributes when variables are first written

  // Add variables to file
  // Add integers
  for(const auto& var : int_arr) {
    if (!file->addVarInt(var.name, var.save_repeat)) {
      throw BoutException("Failed to add int variable {:s} to Datafile", var.name);
    }
  }

  // Add vectors of integers
  for(const auto& var : int_vec_arr) {
    if (!file->addVarIntVec(var.name, var.save_repeat, var.ptr->size())) {
      throw BoutException("Failed to add int vector variable {:s} to Datafile", var.name);
    }
  }

  // Add strings
  for(const auto& var : string_arr) {
    if (!file->addVarString(var.name, var.save_repeat, var.ptr->size())) {
      throw BoutException("Failed to add string variable {:s} to Datafile", var.name);
    }
  }

  // Add BoutReals
  for(const auto& var : BoutReal_arr) {
    if (!file->addVarBoutReal(var.name, var.save_repeat)) {
      throw BoutException("Failed to add BoutReal variable {:s} to Datafile", var.name);
    }
  }

  // Add bools
  for(const auto& var : bool_arr) {
    if (!file->addVarInt(var.name, var.save_repeat)) {
      throw BoutException("Failed to add bool variable {:s} to Datafile", var.name);
    }
  }

  // Add 2D fields
  for (const auto& var : f2d_arr) {
    if (!file->addVarField2D(var.name, var.save_repeat)) {
      throw BoutException("Failed to add Field2D variable {:s} to Datafile", var.name);
    }
  }

  // Add 3D fields
  for (const auto& var : f3d_arr) {
    if (!file->addVarField3D(var.name, var.save_repeat)) {
      throw BoutException("Failed to add Field3D variable {:s} to Datafile", var.name);
    }
  }

  // Add FieldPerps
  for (const auto& var : fperp_arr) {
    if (!file->addVarFieldPerp(var.name, var.save_repeat)) {
      throw BoutException("Failed to add FieldPerp variable {:s} to Datafile", var.name);
    }
  }

  // 2D vectors
  for(const auto& var : v2d_arr) {
    auto name = var.covar ? var.name + "_" : var.name;
    if (!file->addVarField2D(name + "x", var.save_repeat)) {
      throw BoutException("Failed to add Vector2D variable {:s} to Datafile", name);
    }
    if (!file->addVarField2D(name + "y", var.save_repeat)) {
      throw BoutException("Failed to add Vector2D variable {:s} to Datafile", name);
    }
    if (!file->addVarField2D(name + "z", var.save_repeat)) {
      throw BoutException("Failed to add Vector2D variable {:s} to Datafile", name);
    }
  }

  // 3D vectors
  for(const auto& var : v3d_arr) {
    auto name = var.covar ? var.name + "_" : var.name;
    if (!file->addVarField3D(name + "x", var.save_repeat)) {
      throw BoutException("Failed to add Vector3D variable {:s} to Datafile", name);
    }
    if (!file->addVarField3D(name + "y", var.save_repeat)) {
      throw BoutException("Failed to add Vector3D variable {:s} to Datafile", name);
    }
    if (!file->addVarField3D(name + "z", var.save_repeat)) {
      throw BoutException("Failed to add Vector3D variable {:s} to Datafile", name);
    }
  }
  if (openclose) {
    file->close();
  }

  writable = true;

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
  writable = false;
}

void Datafile::setLowPrecision() {
  if(!enabled)
    return;
  floats = true;
  file->setLowPrecision();
}

void Datafile::add(int &i, const char *name, bool save_repeat, const std::string &description) {
  TRACE("DataFile::add(int)");
  if (!enabled)
    return;
  if (varAdded(name)) {
    // Check if it's the same variable
    if (&i == varPtr(name)) {
      output_warn.write("WARNING: variable '{:s}' already added to Datafile, skipping...\n",
                        name);
      return;
    } else {
      throw BoutException("Variable with name '{:s}' already added to Datafile", name);
    }
  }

  VarStr<int> d;

  d.ptr = &i;
  d.name = name;
  d.save_repeat = save_repeat;
  d.covar = false;
  d.description = description;
  
  int_arr.push_back(d);

  if (writable) {
    // Otherwise will add variables when Datafile is opened for writing/appending
    if (openclose) {
      // Open the file
      // Check filename has been set
      if (filename.empty())
        throw BoutException("Datafile::add: Filename has not been set");
      if(!file->openw(filename, BoutComm::rank(), appending)) {
        if (appending) {
          throw BoutException("Datafile::add: Failed to open file {:s} for appending!",
                              filename);
        } else {
          throw BoutException("Datafile::add: Failed to open file {:s} for writing!",
                              filename);
        }
      }
      appending = true;
    }

    if(!file->is_valid())
      throw BoutException("Datafile::add: File is not valid!");

    // Add variable to file
    if (!file->addVarInt(name, save_repeat)) {
      throw BoutException("Failed to add int variable {:s} to Datafile", name);
    }

    if(openclose) {
      file->close();
    }
  }
}

void Datafile::add(std::vector<int> &i, const char *name, bool save_repeat, const std::string &description) {
  TRACE("DataFile::add(std::vector<int>)");
  if (!enabled)
    return;
  if (varAdded(name)) {
    // Check if it's the same variable
    if (&i == varPtr(name)) {
      output_warn.write("WARNING: variable '{:s}' already added to Datafile, skipping...\n",
                        name);
      return;
    } else {
      throw BoutException("Variable with name '{:s}' already added to Datafile", name);
    }
  }

  VarStr<std::vector<int>> d;

  d.ptr = &i;
  d.name = name;
  d.save_repeat = save_repeat;
  d.covar = false;
  d.size = i.size();
  d.description = description;

  int_vec_arr.push_back(d);

  if (writable) {
    // Otherwise will add variables when Datafile is opened for writing/appending
    if (openclose) {
      // Open the file
      // Check filename has been set
      if (filename.empty())
        throw BoutException("Datafile::add: Filename has not been set");
      if(!file->openw(filename, BoutComm::rank(), appending)) {
        if (appending) {
          throw BoutException("Datafile::add: Failed to open file {:s} for appending!",
                              filename);
        } else {
          throw BoutException("Datafile::add: Failed to open file {:s} for writing!",
                              filename);
        }
      }
      appending = true;
    }

    if(!file->is_valid())
      throw BoutException("Datafile::add: File is not valid!");

    // Add variable to file
    if (!file->addVarIntVec(name, save_repeat, i.size())) {
      throw BoutException("Failed to add int vector variable {:s} to Datafile", name);
    }

    if(openclose) {
      file->close();
    }
  }
}

void Datafile::add(std::string &s, const char *name, bool save_repeat, const std::string &description) {
  TRACE("DataFile::add(std::string)");
  if (!enabled) {
    return;
  }
  if (varAdded(name)) {
    // Check if it's the same variable
    if (&s == varPtr(name)) {
      output_warn.write("WARNING: variable '{:s}' already added to Datafile, skipping...\n",
                        name);
      return;
    } else {
      throw BoutException("Variable with name '{:s}' already added to Datafile", name);
    }
  }

  VarStr<std::string> d;

  d.ptr = &s;
  d.name = name;
  d.save_repeat = save_repeat;
  d.covar = false;
  d.size = s.size();
  d.description = description;

  string_arr.push_back(d);

  if (writable) {
    // Otherwise will add variables when Datafile is opened for writing/appending
    if (openclose) {
      // Open the file
      // Check filename has been set
      if (filename.empty()) {
        throw BoutException("Datafile::add: Filename has not been set");
      }
      if(!file->openw(filename, BoutComm::rank(), appending)) {
        if (appending) {
          throw BoutException("Datafile::add: Failed to open file {:s} for appending!",
                              filename);
        } else {
          throw BoutException("Datafile::add: Failed to open file {:s} for writing!",
                              filename);
        }
      }
      appending = true;
    }

    if (!file->is_valid()) {
      throw BoutException("Datafile::add: File is not valid!");
    }

    // Add variable to file
    if (!file->addVarString(name, save_repeat, s.size())) {
      throw BoutException("Failed to add string variable {:s} to Datafile", name);
    }

    if (openclose) {
      file->close();
    }
  }
}

void Datafile::add(BoutReal &r, const char *name, bool save_repeat, const std::string &description) {
  TRACE("DataFile::add(BoutReal)");
  if (!enabled)
    return;
  if (varAdded(name)) {
    // Check if it's the same variable
    if (&r == varPtr(name)) {
      output_warn.write("WARNING: variable '{:s}' already added to Datafile, skipping...\n",
                        name);
      return;
    } else {
      throw BoutException("Variable with name '{:s}' already added to Datafile", name);
    }
  }

  VarStr<BoutReal> d;

  d.ptr = &r;
  d.name = name;
  d.save_repeat = save_repeat;
  d.covar = false;
  d.description = description;
  
  BoutReal_arr.push_back(d);

  if (writable) {
    // Otherwise will add variables when Datafile is opened for writing/appending
    if (openclose) {
      // Open the file
      if (filename.empty())
        throw BoutException("Datafile::add: Filename has not been set");
      if(!file->openw(filename, BoutComm::rank(), appending)) {
        if (appending) {
          throw BoutException("Datafile::add: Failed to open file {:s} for appending!",
                              filename);
        } else {
          throw BoutException("Datafile::add: Failed to open file {:s} for writing!",
                              filename);
        }
      }
      appending = true;
    }

    if(!file->is_valid())
      throw BoutException("Datafile::add: File is not valid!");

    if(floats)
      file->setLowPrecision();

    // Add variable to file
    if (!file->addVarBoutReal(name, save_repeat)) {
      throw BoutException("Failed to add BoutReal variable {:s} to Datafile", name);
    }

    if(openclose) {
      file->close();
    }
  }
}

void Datafile::add(bool &b, const char *name, bool save_repeat, const std::string &description) {
  TRACE("DataFile::add(bool)");
  if (!enabled)
    return;
  if (varAdded(name)) {
    // Check if it's the same variable
    if (&b == varPtr(name)) {
      output_warn.write("WARNING: variable '{:s}' already added to Datafile, skipping...\n",
                        name);
      return;
    } else {
      throw BoutException("Variable with name '{:s}' already added to Datafile", name);
    }
  }

  VarStr<bool> d;

  d.ptr = &b;
  d.name = name;
  d.save_repeat = save_repeat;
  d.covar = false;
  d.description = description;

  bool_arr.push_back(d);

  if (writable) {
    // Otherwise will add variables when Datafile is opened for writing/appending
    if (openclose) {
      // Open the file
      // Check filename has been set
      if (filename.empty())
        throw BoutException("Datafile::add: Filename has not been set");
      if(!file->openw(filename, BoutComm::rank(), appending)) {
        if (appending) {
          throw BoutException("Datafile::add: Failed to open file {:s} for appending!",
                              filename);
        } else {
          throw BoutException("Datafile::add: Failed to open file {:s} for writing!",
                              filename);
        }
      }
      appending = true;
    }

    if(!file->is_valid())
      throw BoutException("Datafile::add: File is not valid!");

    // Add variable to file
    if (!file->addVarInt(name, save_repeat)) {
      throw BoutException("Failed to add bool variable {:s} to Datafile", name);
    }

    if(openclose) {
      file->close();
    }
  }
}

void Datafile::add(Field2D &f, const char *name, bool save_repeat, const std::string &description) {
  TRACE("DataFile::add(Field2D)");
  if (!enabled)
    return;
  if (varAdded(name)) {
    // Check if it's the same variable
    if (&f == varPtr(name)) {
      output_warn.write("WARNING: variable '{:s}' already added to Datafile, skipping...\n",
                        name);
      return;
    } else {
      throw BoutException("Variable with name '{:s}' already added to Datafile", name);
    }
  }

  VarStr<Field2D> d;

  d.ptr = &f;
  d.name = name;
  d.save_repeat = save_repeat;
  d.covar = false;
  d.description = description;
  
  f2d_arr.push_back(d);

  if (writable) {
    // Otherwise will add variables when Datafile is opened for writing/appending
    if (openclose) {
      // Open the file
      if (filename.empty())
        throw BoutException("Datafile::add: Filename has not been set");
      if(!file->openw(filename, BoutComm::rank(), appending)) {
        if (appending) {
          throw BoutException("Datafile::add: Failed to open file {:s} for appending!",
                              filename);
        } else {
          throw BoutException("Datafile::add: Failed to open file {:s} for writing!",
                              filename);
        }
      }
      appending = true;
    }

    if(!file->is_valid())
      throw BoutException("Datafile::add: File is not valid!");

    if(floats)
      file->setLowPrecision();

    // Add variable to file
    if (!file->addVarField2D(name, save_repeat)) {
      throw BoutException("Failed to add Field2D variable {:s} to Datafile", name);
    }

    if(openclose) {
      file->close();
    }
  }
}

void Datafile::add(Field3D &f, const char *name, bool save_repeat, const std::string &description) {
  TRACE("DataFile::add(Field3D)");
  if (!enabled)
    return;
  if (varAdded(name)) {
    // Check if it's the same variable
    if (&f == varPtr(name)) {
      output_warn.write("WARNING: variable '{:s}' already added to Datafile, skipping...\n",
                        name);
      return;
    } else {
      throw BoutException("Variable with name '{:s}' already added to Datafile", name);
    }
  }

  VarStr<Field3D> d;

  d.ptr = &f;
  d.name = name;
  d.save_repeat = save_repeat;
  d.covar = false;
  d.description = description;
  
  f3d_arr.push_back(d);

  if (writable) {
    // Otherwise will add variables when Datafile is opened for writing/appending
    if (openclose) {
      // Open the file
      if (filename.empty())
        throw BoutException("Datafile::add: Filename has not been set");
      if(!file->openw(filename, BoutComm::rank(), appending)) {
        if (appending) {
          throw BoutException("Datafile::add: Failed to open file {:s} for appending!",
                              filename);
        } else {
          throw BoutException("Datafile::add: Failed to open file {:s} for writing!",
                              filename);
        }
      }
      appending = true;
    }

    if(!file->is_valid())
      throw BoutException("Datafile::add: File is not valid!");

    if(floats)
      file->setLowPrecision();

    // Add variable to file
    if (!file->addVarField3D(name, save_repeat)) {
      throw BoutException("Failed to add Field3D variable {:s} to Datafile", name);
    }

    if(openclose) {
      file->close();
    }
  }
}

void Datafile::add(FieldPerp &f, const char *name, bool save_repeat, const std::string &description) {
  AUTO_TRACE();
  if (!enabled)
    return;
  if (varAdded(name)) {
    // Check if it's the same variable
    if (&f == varPtr(name)) {
      output_warn.write("WARNING: variable '{:s}' already added to Datafile, skipping...\n",
                        name);
      return;
    } else {
      throw BoutException("Variable with name '{:s}' already added to Datafile", name);
    }
  }

  VarStr<FieldPerp> d;

  d.ptr = &f;
  d.name = name;
  d.save_repeat = save_repeat;
  d.covar = false;
  d.description = description;

  fperp_arr.push_back(d);

  if (writable) {
    // Otherwise will add variables when Datafile is opened for writing/appending
    if (openclose) {
      // Open the file
      if (filename.empty())
        throw BoutException("Datafile::add: Filename has not been set");
      if(!file->openw(filename, BoutComm::rank(), appending)) {
        if (appending) {
          throw BoutException("Datafile::add: Failed to open file {:s} for appending!",
                              filename);
        } else {
          throw BoutException("Datafile::add: Failed to open file {:s} for writing!",
                              filename);
        }
      }
      appending = true;
    }

    if(!file->is_valid())
      throw BoutException("Datafile::add: File is not valid!");

    if(floats)
      file->setLowPrecision();

    // Add variable to file
    if (!file->addVarFieldPerp(name, save_repeat)) {
      throw BoutException("Failed to add FieldPerp variable {:s} to Datafile", name);
    }

    if(openclose) {
      file->close();
    }
  }
}

void Datafile::add(Vector2D &f, const char *name, bool save_repeat, const std::string &description) {
  TRACE("DataFile::add(Vector2D)");
  if (!enabled)
    return;
  if (varAdded(name)) {
    // Check if it's the same variable
    if (&f == varPtr(name)) {
      output_warn.write("WARNING: variable '{:s}' already added to Datafile, skipping...\n",
                        name);
      return;
    } else {
      throw BoutException("Variable with name '{:s}' already added to Datafile", name);
    }
  }

  VarStr<Vector2D> d;

  d.ptr = &f;
  d.name = name;
  d.save_repeat = save_repeat;
  d.covar = f.covariant;
  d.description = description;

  v2d_arr.push_back(d);

  if (writable) {
    // Otherwise will add variables when Datafile is opened for writing/appending
    if (openclose) {
      // Open the file
      if (filename.empty())
        throw BoutException("Datafile::add: Filename has not been set");
      if(!file->openw(filename, BoutComm::rank(), appending)) {
        if (appending) {
          throw BoutException("Datafile::add: Failed to open file {:s} for appending!",
                              filename);
        } else {
          throw BoutException("Datafile::add: Failed to open file {:s} for writing!",
                              filename);
        }
      }
      appending = true;
    }

    if(!file->is_valid())
      throw BoutException("Datafile::add: File is not valid!");

    if(floats)
      file->setLowPrecision();

    // Add variables to file
    auto dname = d.covar ? d.name + "_" : d.name;
#if BOUT_USE_METRIC_3D
    // Add variables to file
    if (!file->addVarField3D(dname + "x", save_repeat)) {
      throw BoutException("Failed to add Vector2D variable {:s} to Datafile", dname);
    }
    if (!file->addVarField3D(dname + "y", save_repeat)) {
      throw BoutException("Failed to add Vector2D variable {:s} to Datafile", dname);
    }
    if (!file->addVarField3D(dname + "z", save_repeat)) {
      throw BoutException("Failed to add Vector2D variable {:s} to Datafile", dname);
    }
#else
    if (!file->addVarField2D(dname + "x", save_repeat)) {
      throw BoutException("Failed to add Vector2D variable {:s} to Datafile", dname);
    }
    if (!file->addVarField2D(dname + "y", save_repeat)) {
      throw BoutException("Failed to add Vector2D variable {:s} to Datafile", dname);
    }
    if (!file->addVarField2D(dname + "z", save_repeat)) {
      throw BoutException("Failed to add Vector2D variable {:s} to Datafile", dname);
    }
#endif

    if(openclose) {
      file->close();
    }
  }
}

void Datafile::add(Vector3D &f, const char *name, bool save_repeat, const std::string &description) {
  TRACE("DataFile::add(Vector3D)");
  if (!enabled)
    return;
  if (varAdded(name)) {
    // Check if it's the same variable
    if (&f == varPtr(name)) {
      output_warn.write("WARNING: variable '{:s}' already added to Datafile, skipping...\n",
                        name);
      return;
    } else {
      throw BoutException("Variable with name '{:s}' already added to Datafile", name);
    }
  }

  VarStr<Vector3D> d;

  d.ptr = &f;
  d.name = name;
  d.save_repeat = save_repeat;
  d.covar = f.covariant;
  d.description = description;

  v3d_arr.push_back(d);

  if (writable) {
    // Otherwise will add variables when Datafile is opened for writing/appending
    if (openclose) {
      // Open the file
      if (filename.empty())
        throw BoutException("Datafile::add: Filename has not been set");
      if(!file->openw(filename, BoutComm::rank(), appending)) {
        if (appending) {
          throw BoutException("Datafile::add: Failed to open file {:s} for appending!",
                              filename);
        } else {
          throw BoutException("Datafile::add: Failed to open file {:s} for writing!",
                              filename);
        }
      }
      appending = true;
    }

    if(!file->is_valid())
      throw BoutException("Datafile::add: File is not valid!");

    if(floats)
      file->setLowPrecision();

    // Add variables to file
    auto dname = d.covar ? d.name + "_" : d.name;
    if (!file->addVarField3D(dname + "x", save_repeat)) {
      throw BoutException("Failed to add Vector3D variable {:s} to Datafile", dname);
    }
    if (!file->addVarField3D(dname + "y", save_repeat)) {
      throw BoutException("Failed to add Vector3D variable {:s} to Datafile", dname);
    }
    if (!file->addVarField3D(dname + "z", save_repeat)) {
      throw BoutException("Failed to add Vector3D variable {:s} to Datafile", dname);
    }

    if(openclose) {
      file->close();
    }
  }
}

namespace {
// Read a value from file and check it matches reference_value, throw if not
void checkGridValue(DataFormat* file, const std::string& name,
                    const std::string& filename, const int reference_value) {
  int file_value;
  if (!file->read(&file_value, name)) {
    throw BoutException("Could not read {} from file '{}'", name, filename);
  }

  if (file_value != reference_value) {
    throw BoutException("{} ({}) in file '{}' does not match value in mesh ({})", name,
                        file_value, filename, reference_value);
  }
}

// Check that the array sizes in \p file match those in existing \p mesh
void checkFileGrid(DataFormat* file, const std::string& filename, const Mesh* mesh) {
  checkGridValue(file, "MXG", filename, mesh->xstart);
  checkGridValue(file, "MYG", filename, mesh->ystart);
  checkGridValue(file, "MZG", filename, mesh->zstart);
  // nx includes boundaries
  checkGridValue(file, "nx", filename, mesh->GlobalNx);
  checkGridValue(file, "ny", filename, mesh->GlobalNyNoBoundaries);
  checkGridValue(file, "nz", filename, mesh->GlobalNzNoBoundaries);
}
} // namespace

bool Datafile::read() {
  Timer timer("io");  ///< Start timer. Stops when goes out of scope

  if(openclose) {
    // Open the file
    if(!file->openr(filename, BoutComm::rank())) {
      throw BoutException("Datafile::read: Failed to open file {:s} for reading!",
                          filename);
    }
  }
  
  if(!file->is_valid())
    throw BoutException("Datafile::read: File is not valid!");

  checkFileGrid(file.get(), filename, mesh);

  file->setRecord(-1); // Read the latest record

  // Read integers
  for(const auto& var : int_arr) {
    if(var.save_repeat) {
      if(!file->read_rec(var.ptr, var.name.c_str())) {
        if(!init_missing) {
          throw BoutException(
              "Missing data for {:s} in input. Set init_missing=true to set to zero.",
              var.name);
        }
        output_warn.write("\tWARNING: Could not read integer {:s}. Setting to zero\n", var.name);
        *(var.ptr) = 0;
        continue;
      }
    } else {
      if(!file->read(var.ptr, var.name.c_str())) {
        if(!init_missing) {
          throw BoutException(
              "Missing data for {:s} in input. Set init_missing=true to set to zero.",
              var.name);
        }
        output_warn.write("\tWARNING: Could not read integer {:s}. Setting to zero\n", var.name);
        *(var.ptr) = 0;
        continue;
      }
    }
  }

  // Read vectors of integers
  for(const auto& var : int_vec_arr) {
    if (var.ptr->size() != var.size) {
      throw BoutException("Size of std::vector<int> '{:s}' has changed since being added "
                          "to Datafile. Cannot read.", var.name);
    }
    if(var.save_repeat) {
      if(!file->read_rec(&(*var.ptr)[0], var.name.c_str(), var.ptr->size())) {
        if(!init_missing) {
          throw BoutException(
              "Missing data for {:s} in input. Set init_missing=true to create empty vector.",
              var.name);
        }
        output_warn.write("\tWARNING: Could not read integer vector {:s}. Creating empty vector\n", var.name);
        *(var.ptr) = {};
        continue;
      }
    } else {
      if(!file->read(&(*var.ptr)[0], var.name.c_str(), var.ptr->size())) {
        if(!init_missing) {
          throw BoutException(
              "Missing data for {:s} in input. Set init_missing=true to create empty vector.",
              var.name);
        }
        output_warn.write("\tWARNING: Could not read integer vector {:s}. Creating empty vector\n", var.name);
        *(var.ptr) = {};
        continue;
      }
    }
  }

  // Read strings
  for (const auto& var : string_arr) {
    if (var.ptr->size() != var.size) {
      throw BoutException("Size of std::string '{:s}' has changed since being "
                          "added to Datafile. Cannot read.", var.name);
    }
    var.ptr->resize(var.size);
    if (var.save_repeat) {
      if (!file->read_rec(&(*var.ptr)[0], var.name.c_str(), var.size)) {
        if (!init_missing) {
          throw BoutException(
              "Missing data for {:s} in input. Set init_missing=true to create empty string.",
              var.name);
        }
        output_warn.write("\tWARNING: Could not read string {:s}. Creating empty string\n", var.name);
      }
    } else {
      if (!file->read(&(*var.ptr)[0], var.name.c_str(), var.size)) {
        if (!init_missing) {
          throw BoutException(
              "Missing data for {:s} in input. Set init_missing=true to create empty string.",
              var.name);
        }
        output_warn.write("\tWARNING: Could not read string {:s}. Creating empty string\n", var.name);
      }
    }
  }

  // Read BoutReals
  for(const auto& var : BoutReal_arr) {
    if(var.save_repeat) {
      if(!file->read_rec(var.ptr, var.name)) {
        if(!init_missing) {
          throw BoutException(
              "Missing data for {:s} in input. Set init_missing=true to set to zero.",
              var.name);
        }
        output_warn.write("\tWARNING: Could not read BoutReal {:s}. Setting to zero\n", var.name);
        *(var.ptr) = 0;
        continue;
      }
    } else {
      if(!file->read(var.ptr, var.name)) {
        if(!init_missing) {
          throw BoutException(
              "Missing data for {:s} in input. Set init_missing=true to set to zero.",
              var.name);
        }
        output_warn.write("\tWARNING: Could not read BoutReal {:s}. Setting to zero\n", var.name);
        *(var.ptr) = 0;
        continue;
      }
    }
  }
  
  // Read bools
  for(const auto& var : bool_arr) {
    int var_as_int = 0;
    if(var.save_repeat) {
      if(!file->read_rec(&var_as_int, var.name.c_str())) {
        if(!init_missing) {
          throw BoutException(
              "Missing data for {:s} in input. Set init_missing=true to set to false.",
              var.name);
        }
        output_warn.write("\tWARNING: Could not read bool {:s}. Setting to false\n", var.name);
        *(var.ptr) = false;
        continue;
      }
    } else {
      if(!file->read(&var_as_int, var.name.c_str())) {
        if(!init_missing) {
          throw BoutException(
              "Missing data for {:s} in input. Set init_missing=true to set to false.",
              var.name);
        }
        output_warn.write("\tWARNING: Could not read bool {:s}. Setting to false\n", var.name);
        *(var.ptr) = false;
        continue;
      }
    }
    *var.ptr = bool(var_as_int);
  }

  // Read 2D fields
  for(const auto& var : f2d_arr) {
    read_f2d(var.name, var.ptr, var.save_repeat);
  }

  // Read 3D fields
  for(const auto& var : f3d_arr) {
    read_f3d(var.name, var.ptr, var.save_repeat);
  }

  // Read FieldPerps
  for(const auto& var : fperp_arr) {
    read_fperp(var.name, var.ptr, var.save_repeat);
  }

  // 2D vectors
#if BOUT_USE_METRIC_3D
  for (const auto& var : v2d_arr) {
    if (var.covar) {
      // Reading covariant vector
      read_f3d(var.name + "_x", &(var.ptr->x), var.save_repeat);
      read_f3d(var.name + "_y", &(var.ptr->y), var.save_repeat);
      read_f3d(var.name + "_z", &(var.ptr->z), var.save_repeat);
    } else {
      read_f3d(var.name + "x", &(var.ptr->x), var.save_repeat);
      read_f3d(var.name + "y", &(var.ptr->y), var.save_repeat);
      read_f3d(var.name + "z", &(var.ptr->z), var.save_repeat);
    }
#else
  for(const auto& var : v2d_arr) {
    if(var.covar) {
      // Reading covariant vector
      read_f2d(var.name + "_x", &(var.ptr->x), var.save_repeat);
      read_f2d(var.name + "_y", &(var.ptr->y), var.save_repeat);
      read_f2d(var.name + "_z", &(var.ptr->z), var.save_repeat);
    } else {
      read_f2d(var.name + "x", &(var.ptr->x), var.save_repeat);
      read_f2d(var.name + "y", &(var.ptr->y), var.save_repeat);
      read_f2d(var.name + "z", &(var.ptr->z), var.save_repeat);
    }
#endif
    var.ptr->covariant = var.covar;
  }

  // 3D vectors
  for(const auto& var : v3d_arr) {
    if(var.covar) {
      // Reading covariant vector
      read_f3d(var.name + "_x", &(var.ptr->x), var.save_repeat);
      read_f3d(var.name + "_y", &(var.ptr->y), var.save_repeat);
      read_f3d(var.name + "_z", &(var.ptr->z), var.save_repeat);
    } else {
      read_f3d(var.name + "x", &(var.ptr->x), var.save_repeat);
      read_f3d(var.name + "y", &(var.ptr->y), var.save_repeat);
      read_f3d(var.name + "z", &(var.ptr->z), var.save_repeat);
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

  if(openclose && (flushFrequencyCounter % flushfrequency == 0)) {
    // Open the file
    if(!file->openw(filename, BoutComm::rank(), appending)) {
      if (appending) {
        throw BoutException("Datafile::add: Failed to open file {:s} for appending!",
                            filename);
      } else {
        throw BoutException("Datafile::add: Failed to open file {:s} for writing!",
                            filename);
      }
    }
    appending = true;
    flushFrequencyCounter = 0;
  }
  
  if(!file->is_valid())
    throw BoutException("Datafile::open: File is not valid!");

  if(floats)
    file->setLowPrecision();
  
  Timer timer("io");
  
  file->setRecord(-1); // Latest record

  if (first_time) {
    first_time = false;

    // Set the field attributes from field meta-data.
    // Attributes must have been set for all fields before the first time
    // output is written, since this happens after the first rhs evaluation

    // Integer variables
    for(const auto& var : int_arr) {
      if (not var.description.empty()) {
        file->setAttribute(var.name, "description", var.description);
      }
    }

    // Vectors of integers
    for(const auto& var : int_vec_arr) {
      if (not var.description.empty()) {
        file->setAttribute(var.name, "description", var.description);
      }
    }

    // String variables
    for (const auto& var : string_arr) {
      if (not var.description.empty()) {
        file->setAttribute(var.name, "description", var.description);
      }
    }

    // BoutReal variables
    for(const auto& var : BoutReal_arr) {
      if (not var.description.empty()) {
        file->setAttribute(var.name, "description", var.description);
      }
    }

    // bool variables
    for(const auto& var : bool_arr) {
      if (not var.description.empty()) {
        file->setAttribute(var.name, "description", var.description);
      }
    }

    // 2D fields
    for (const auto& var : f2d_arr) {
      file->writeFieldAttributes(var.name, *var.ptr);
      if (not var.description.empty()) {
        file->setAttribute(var.name, "description", var.description);
      }
    }

    // 3D fields
    for (const auto& var : f3d_arr) {
      file->writeFieldAttributes(var.name, *var.ptr, shiftoutput);
      if (not var.description.empty()) {
        file->setAttribute(var.name, "description", var.description);
      }
    }

    // FieldPerps
    for (const auto& var : fperp_arr) {
      file->writeFieldAttributes(var.name, *var.ptr, shiftoutput);
      if (not var.description.empty()) {
        file->setAttribute(var.name, "description", var.description);
      }
    }

    // 2D vectors
    for(const auto& var : v2d_arr) {
      Vector2D v  = *(var.ptr);
      auto name = var.covar ? var.name + "_" : var.name;
      file->writeFieldAttributes(name+"x", v.x);
      file->writeFieldAttributes(name+"y", v.y);
      file->writeFieldAttributes(name+"z", v.z);
      if (not var.description.empty()) {
        file->setAttribute(var.name, "description", var.description);
      }
    }

    // 3D vectors
    for(const auto& var : v3d_arr) {
      Vector3D v  = *(var.ptr);
      auto name = var.covar ? var.name + "_" : var.name;
      file->writeFieldAttributes(name+"x", v.x, shiftoutput);
      file->writeFieldAttributes(name+"y", v.y, shiftoutput);
      file->writeFieldAttributes(name+"z", v.z, shiftoutput);
      if (not var.description.empty()) {
        file->setAttribute(var.name, "description", var.description);
      }
    }
  }

  // Write integers
  for(const auto& var : int_arr) {
    write_int(var.name, var.ptr, var.save_repeat);
  }
  
  // Write vectors of integers
  for(const auto& var : int_vec_arr) {
    if (var.ptr->size() != var.size) {
      throw BoutException("Size of std::vector<int> '{:s}' has changed since being added "
                          "to Datafile. Cannot write.", var.name);
    }
    write_int_vec(var.name, var.ptr, var.save_repeat);
  }

  // Write strings
  for (const auto& var : string_arr) {
    if (var.ptr->size() != var.size) {
      throw BoutException("Size of string '{:s}' has changed since being "
                          "added to Datafile. Cannot write.", var.name);
    }
    write_string(var.name, var.ptr, var.save_repeat);
  }

  // Write BoutReals
  for(const auto& var : BoutReal_arr) {
    write_real(var.name, var.ptr, var.save_repeat);
  }

  // Write bools
  for(const auto& var : bool_arr) {
    int var_as_int = int(*var.ptr);
    write_int(var.name, &var_as_int, var.save_repeat);
  }

  // Write 2D fields
  for (const auto& var : f2d_arr) {
    write_f2d(var.name, var.ptr, var.save_repeat);
  }

  // Write 3D fields
  for (const auto& var : f3d_arr) {
    write_f3d(var.name, var.ptr, var.save_repeat);
  }
  
  // Write FieldPerps
  for (const auto& var : fperp_arr) {
    write_fperp(var.name, var.ptr, var.save_repeat);
  }

  // 2D vectors
  for(const auto& var : v2d_arr) {
    Vector2D v  = *(var.ptr);
    auto name = var.name;

    if(var.covar) {
      // Writing covariant vector
      v.toCovariant();
      name += "_";
    } else {
      // Writing contravariant vector
      v.toContravariant();
    }

#if BOUT_USE_METRIC_3D
    write_f3d(name + "x", &(v.x), var.save_repeat);
    write_f3d(name + "y", &(v.y), var.save_repeat);
    write_f3d(name + "z", &(v.z), var.save_repeat);
#else
    write_f2d(name+"x", &(v.x), var.save_repeat);
    write_f2d(name+"y", &(v.y), var.save_repeat);
    write_f2d(name+"z", &(v.z), var.save_repeat);
#endif
  }

  // 3D vectors
  for(const auto& var : v3d_arr) {
    Vector3D v  = *(var.ptr);
    auto name = var.name;

    if(var.covar) {
      // Writing covariant vector
      v.toCovariant();
      name += "_";
    } else {
      // Writing contravariant vector
      v.toContravariant();
    }

    write_f3d(name+"x", &(v.x), var.save_repeat);
    write_f3d(name+"y", &(v.y), var.save_repeat);
    write_f3d(name+"z", &(v.z), var.save_repeat);
  }
  
  if(openclose  && (flushFrequencyCounter+1 % flushfrequency == 0)){
    file->close();
  }
  flushFrequencyCounter++;
  return true;
}

bool Datafile::write(const std::string& filename_) const {
  if(!enabled)
    return true;

  if (filename_.empty()) {
    throw BoutException("Datafile::write: No argument given!");
  }

  // Create a new datafile
  Datafile tmp(*this);

  tmp.openw(filename_);
  bool ret = tmp.write();
  tmp.close();

  return ret;
}

void Datafile::setAttribute(const std::string &varname, const std::string &attrname, const std::string &text) {

  TRACE("Datafile::setAttribute(string, string, string)");

  Timer timer("io");

  if (!enabled) {
    return;
  }

  if(!file)
    throw BoutException("Datafile::write: File is not valid!");

  if(openclose && (flushFrequencyCounter % flushfrequency == 0)) {
    // Open the file
    if(!file->openw(filename, BoutComm::rank(), appending)) {
      if (appending) {
        throw BoutException("Datafile::add: Failed to open file {:s} for appending!",
                            filename);
      } else {
        throw BoutException("Datafile::add: Failed to open file {:s} for writing!",
                            filename);
      }
    }
    appending = true;
    flushFrequencyCounter = 0;
  }

  if(!file->is_valid())
    throw BoutException("Datafile::setAttribute: File is not valid!");

  file->setAttribute(varname, attrname, text);

  if (openclose) {
    file->close();
  }
}

void Datafile::setAttribute(const std::string &varname, const std::string &attrname, int value) {

  TRACE("Datafile::setAttribute(string, string, int)");

  Timer timer("io");

  if (!enabled) {
    return;
  }

  if(!file)
    throw BoutException("Datafile::write: File is not valid!");

  if(openclose && (flushFrequencyCounter % flushfrequency == 0)) {
    // Open the file
    if(!file->openw(filename, BoutComm::rank(), appending)) {
      if (appending) {
        throw BoutException("Datafile::add: Failed to open file {:s} for appending!",
                            filename);
      } else {
        throw BoutException("Datafile::add: Failed to open file {:s} for writing!",
                            filename);
      }
    }
    appending = true;
    flushFrequencyCounter = 0;
  }

  if(!file->is_valid())
    throw BoutException("Datafile::setAttribute: File is not valid!");

  file->setAttribute(varname, attrname, value);

  if (openclose) {
    file->close();
  }
}

void Datafile::setAttribute(const std::string &varname, const std::string &attrname, BoutReal value) {

  TRACE("Datafile::setAttribute(string, string, BoutReal)");

  Timer timer("io");

  if (!enabled) {
    return;
  }

  if(!file)
    throw BoutException("Datafile::write: File is not valid!");

  if(openclose && (flushFrequencyCounter % flushfrequency == 0)) {
    // Open the file
    if(!file->openw(filename, BoutComm::rank(), appending)) {
      if (appending) {
        throw BoutException("Datafile::add: Failed to open file {:s} for appending!",
                            filename);
      } else {
        throw BoutException("Datafile::add: Failed to open file {:s} for writing!",
                            filename);
      }
    }
    appending = true;
    flushFrequencyCounter = 0;
  }

  if(!file->is_valid())
    throw BoutException("Datafile::setAttribute: File is not valid!");

  file->setAttribute(varname, attrname, value);

  if (openclose) {
    file->close();
  }
}

/////////////////////////////////////////////////////////////

bool Datafile::read_f2d(const std::string &name, Field2D *f, bool save_repeat) {
  try {
    file->readFieldAttributes(name, *f);
  } catch (const BoutException&) {
    if (init_missing) {
      output_warn.write("\tWARNING: Could not read 2D field {:s} attributes.\n", name);
    } else {
      throw;
    }
  }

  f->allocate();
  
  if(save_repeat) {
    if(!file->read_rec(&((*f)(0,0)), name, mesh->LocalNx, mesh->LocalNy)) {
      if(init_missing) {
        output_warn.write("\tWARNING: Could not read 2D field {:s}. Setting to zero\n", name);
        *f = 0.0;
      } else {
        throw BoutException("Missing 2D evolving field {:s} in input. Set "
                            "init_missing=true to set to zero.",
                            name);
      }
      return false;
    }
  }else {
    if(!file->read(&((*f)(0,0)), name, mesh->LocalNx, mesh->LocalNy)) {
      if(init_missing) {
        output_warn.write("\tWARNING: Could not read 2D field {:s}. Setting to zero\n", name);
        *f = 0.0;
      } else {
        throw BoutException(
            "Missing 2D field {:s} in input. Set init_missing=true to set to zero.",
            name);
      }
      return false;
    }
  }
  
  return true;
}

bool Datafile::read_f3d(const std::string &name, Field3D *f, bool save_repeat) {
  try {
    file->readFieldAttributes(name, *f);
  } catch (const BoutException&) {
    if (init_missing) {
      output_warn.write("\tWARNING: Could not read 3D field {:s} attributes.\n", name);
    } else {
      throw;
    }
  }

  f->allocate();
  
  if(save_repeat) {
    if(!file->read_rec(&((*f)(0,0,0)), name, mesh->LocalNx, mesh->LocalNy, mesh->LocalNz)) {
      if(init_missing) {
        output_warn.write("\tWARNING: Could not read 3D field {:s}. Setting to zero\n", name);
        *f = 0.0;
      }else {
        throw BoutException("Missing 3D evolving field {:s} in input. Set "
                            "init_missing=true to set to zero.",
                            name);
      }
      return false;
    }
  }else {
    if(!file->read(&((*f)(0,0,0)), name, mesh->LocalNx, mesh->LocalNy, mesh->LocalNz)) {
      if(init_missing) {
        output_warn.write("\tWARNING: Could not read 3D field {:s}. Setting to zero\n", name);
        *f = 0.0;
      }else {
        throw BoutException(
            "Missing 3D field {:s} in input. Set init_missing=true to set to zero.",
            name);
      }
      return false;
    }
  }

  
  if (shiftinput) {
    // Input file is in field-aligned coordinates e.g. BOUT++ 3.x restart file
    *f = fromFieldAligned(*f, "RGN_ALL");
  }
  
  return true;
}

bool Datafile::read_fperp(const std::string &name, FieldPerp *f, bool save_repeat) {
  try {
    file->readFieldAttributes(name, *f);
  } catch (const BoutException&) {
    if (init_missing) {
      output_warn.write("\tWARNING: Could not read FieldPerp {:s} attributes.\n", name);
    } else {
      throw;
    }
  }

  int yindex = f->getIndex();
  if (yindex >= 0 and yindex < mesh->LocalNy) {
    // yindex is in the range of this processor, so read FieldPerp

    f->allocate();

    if(save_repeat) {
      if(!file->read_rec_perp(&((*f)(0,0)), name, mesh->LocalNx, mesh->LocalNz)) {
        if(init_missing) {
          output_warn.write("\tWARNING: Could not read FieldPerp {:s}. Setting to zero\n", name);
          *f = 0.0;
        }else {
          throw BoutException("Missing evolving FieldPerp {:s} in input. Set "
                              "init_missing=true to set to zero.",
                              name);
        }
        return false;
      }
    }else {
      if(!file->read_perp(&((*f)(0,0)), name, mesh->LocalNx, mesh->LocalNz)) {
        if(init_missing) {
          output_warn.write("\tWARNING: Could not read FieldPerp {:s}. Setting to zero\n", name);
          *f = 0.0;
        }else {
          throw BoutException(
              "Missing FieldPerp {:s} in input. Set init_missing=true to set to zero.",
              name);
        }
        return false;
      }
    }

    if (shiftinput) {
      // Input file is in field-aligned coordinates e.g. BOUT++ 3.x restart file
      *f = fromFieldAligned(*f, "RGN_ALL");
    }
  }

  return true;
}

bool Datafile::write_int(const std::string &name, int *f, bool save_repeat) {
  if(save_repeat) {
    return file->write_rec(f, name);
  }else {
    return file->write(f, name);
  }
}

bool Datafile::write_int_vec(const std::string &name, std::vector<int> *f, bool save_repeat) {
  if(save_repeat) {
    return file->write_rec(&(*f)[0], name, f->size());
  }else {
    return file->write(&(*f)[0], name, f->size());
  }
}

bool Datafile::write_string(const std::string &name, std::string *f, bool save_repeat) {
  if (save_repeat) {
    return file->write_rec(&(*f)[0], name, f->size());
  } else {
    return file->write(&(*f)[0], name, f->size());
  }
}

bool Datafile::write_real(const std::string &name, BoutReal *f, bool save_repeat) {
  if(save_repeat) {
    return file->write_rec(f, name);
  }else {
    return file->write(f, name);
  }
}

bool Datafile::write_f2d(const std::string &name, Field2D *f, bool save_repeat) {
  if (!f->isAllocated()) {
    throw BoutException("Datafile::write_f2d: Field2D '{:s}' is not allocated!", name);
  }
  if (save_repeat) {
    if (!file->write_rec(&((*f)(0, 0)), name, mesh->LocalNx, mesh->LocalNy)) {
      throw BoutException("Datafile::write_f2d: Failed to write {:s}!", name);
    }
  } else {
    if (!file->write(&((*f)(0, 0)), name, mesh->LocalNx, mesh->LocalNy)) {
      throw BoutException("Datafile::write_f2d: Failed to write {:s}!", name);
    }
  }
  return true;
}

bool Datafile::write_f3d(const std::string &name, Field3D *f, bool save_repeat) {
  if (!f->isAllocated()) {
    throw BoutException("Datafile::write_f3d: Field3D '{:s}' is not allocated!", name);
  }

  //Deal with shifting the output
  Field3D f_out{emptyFrom(*f)};
  if(shiftoutput and not (f->getDirectionY() == YDirectionType::Aligned)) {
    f_out = toFieldAligned(*f);
  }else {
    f_out = *f;
  }

  if(save_repeat) {
    return file->write_rec(&(f_out(0,0,0)), name, mesh->LocalNx, mesh->LocalNy, mesh->LocalNz);
  }else {
    return file->write(&(f_out(0,0,0)), name, mesh->LocalNx, mesh->LocalNy, mesh->LocalNz);
  }
}

bool Datafile::write_fperp(const std::string &name, FieldPerp *f, bool save_repeat) {
  int yindex = f->getIndex();
  if (yindex >= 0 and yindex < mesh->LocalNy) {
    if (!f->isAllocated()) {
      throw BoutException("Datafile::write_fperp: FieldPerp '{:s}' is not allocated!",
                          name);
    }

    //Deal with shifting the output
    FieldPerp f_out{emptyFrom(*f)};
    if(shiftoutput and not (f->getDirectionY() == YDirectionType::Aligned)) {
      f_out = toFieldAligned(*f);
    }else {
      f_out = *f;
    }

    if(save_repeat) {
      return file->write_rec_perp(&(f_out(0,0)), name, mesh->LocalNx, mesh->LocalNz);
    }else {
      return file->write_perp(&(f_out(0,0)), name, mesh->LocalNx, mesh->LocalNz);
    }
  }

  // Don't need to write f as it's y-index is not on this processor. Return
  // without doing anything.
  return true;
}

bool Datafile::varAdded(const std::string &name) {
  for(const auto& var : int_arr ) {
    if(name == var.name)
      return true;
  }

  for(const auto& var : int_vec_arr ) {
    if(name == var.name)
      return true;
  }

  for (const auto& var : string_arr) {
    if(name == var.name) {
      return true;
    }
  }

  for(const auto& var : BoutReal_arr ) {
    if(name == var.name)
      return true;
  }

  for(const auto& var : bool_arr ) {
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
  
  for(const auto& var : fperp_arr ) {
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

void *Datafile::varPtr(const std::string &name) {
  for (const auto &var : int_arr) {
    if (name == var.name) {
      return static_cast<void *>(var.ptr);
    }
  }

  for (const auto &var : int_vec_arr) {
    if (name == var.name) {
      return static_cast<void *>(var.ptr);
    }
  }

  for (const auto &var : string_arr) {
    if (name == var.name) {
      return static_cast<void *>(var.ptr);
    }
  }

  for (const auto &var : BoutReal_arr) {
    if (name == var.name) {
      return static_cast<void *>(var.ptr);
    }
  }

  for (const auto &var : bool_arr) {
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

  for (const auto &var : fperp_arr) {
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

