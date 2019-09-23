/*!************************************************************************
 * Provides access to the PETSc library, handling initialisation and
 * finalisation. 
 *
 * Usage
 * -----
 *
 * #include <bout/petsclib.hxx>
 * 
 * class MyClass {
 *   public:
 *   
 *   private:
 *     PetscLib lib;
 * };
 * 
 *
 * This will then automatically initialise Petsc the first time an object
 * is created, and finalise it when the last object is destroyed.
 * 
 **************************************************************************
 * Copyright 2012 B.D.Dudson, S.Farley, M.V.Umansky, X.Q.Xu
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

class PetscLib;

#ifndef __PETSCLIB_H__
#define __PETSCLIB_H__

#ifdef BOUT_HAS_PETSC

#include "options.hxx"

#include <petsc.h>
#include <petscversion.h>

/*!
 * Handles initialisation and finalisation of PETSc library.
 * The first instance which is created initialises PETSc
 * Keeps a count of the number of how many instances exist
 * When the last instance is destroyed it finalises PETSc.
 */ 
class PetscLib {
public:
  /*!
   * Ensure that PETSc has been initialised
   */
  PetscLib(Options* opt = nullptr);
  
  /*!
   * Calls PetscFinalize when all PetscLib instances are destroyed
   */ 
  ~PetscLib();
  
  /*!
   * This is called once to set the command-line options.
   * Should be done early in the program, before any instances of
   * PetscLib are created.
   * The arguments will be passed to PetscInitialize()
   */ 
  static void setArgs(int &c, char** &v) { pargc = &c; pargv = &v;}

  /// Create a KSP linear solver that uses the options specific to this PetscLib
  void createKSPWithOptions(MPI_Comm& comm, KSP& ksp);

  /*!
   * Force cleanup. This will call PetscFinalize, printing a warning
   * if any instances of PetscLib still exist
   */ 
  static void cleanup(); 
private:
  static int count; ///< How many instances?
  static char help[]; ///< Help string
  
  // Command-line arguments
  static int* pargc;
  static char*** pargv;
  
  // Prefix for object-specific options
  std::string options_prefix;

  static PetscLogEvent USER_EVENT;

  void setPetscOptions(Options& options, std::string pass_options_prefix);
};

#ifndef PETSC_VERSION_GE
// Newer versions of PETSc define these symbols for testing library version
// This is a re-implementation of the PETSc BSD-licensed code

#define PETSC_VERSION_GE(MAJOR,MINOR,SUBMINOR)                          \
  ( (PETSC_VERSION_MAJOR > MAJOR) ||                                    \
    ( (PETSC_VERSION_MAJOR == MAJOR) && ( (PETSC_VERSION_MINOR > MINOR) || \
                                          ( (PETSC_VERSION_MINOR == MINOR) && (PETSC_VERSION_SUBMINOR >= SUBMINOR)))))

#endif // PETSC_VERSION_GE

#else // BOUT_HAS_PETSC

#include "unused.hxx"

class PetscLib {
public:
  PetscLib(Options* UNUSED(opt) = nullptr) {}
  ~PetscLib() {}
  
  static void setArgs(int &UNUSED(c), char** &UNUSED(v)) {}
  
  void createKSPWithOptions(MPI_Comm& comm, KSP& ksp);

  static void cleanup() {}
};

#endif // BOUT_HAS_PETSC


#endif //  __PETSCLIB_H__
