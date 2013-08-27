/**************************************************************************
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

// Fix error in PETSC_DEPRECATED("Use SNESGetLineSearch()") on Hopper (PETSc-3.4)
#define PETSC_DEPRECATED(a)

#include <petsc.h>

class PetscLib {
public:
  PetscLib();
  ~PetscLib();
  
  static void setArgs(int &c, char** &v) { pargc = &c; pargv = &v;}
  
  static void cleanup(); // Force cleanup
private:
  static int count; // How many instances?
  static char help[]; // Help string
  
  // Command-line arguments
  static int* pargc;
  static char*** pargv;
  
  static PetscLogEvent USER_EVENT;
};

#else // BOUT_HAS_PETSC

class PetscLib {
public:
  PetscLib() {}
  ~PetscLib() {}
  
  static void setArgs(int &c, char** &v) {}
  
  static void cleanup() {}
};

#endif // BOUT_HAS_PETSC


#endif //  __PETSCLIB_H__
