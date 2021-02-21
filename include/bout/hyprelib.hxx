/*!************************************************************************
 * Provides access to the Hypre library, handling initialisation and
 * finalisation. 
 *
 * Usage
 * -----
 *
 * #include <bout/hyprelib.hxx>
 * 
 * class MyClass {
 *   public:
 *   
 *   private:
 *     HypreLib lib;
 * };
 * 
 *
 * This will then automatically initialise Hypre the first time an object
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

class HypreLib;

#ifndef __HYPRELIB_H__
#define __HYPRELIB_H__

#include "bout/build_config.hxx"

class Options;

#if BOUT_HAS_HYPRE

#include <HYPRE.h>
#include "HYPRE_utilities.h"
#include "_hypre_utilities.h"


/*!
 * Handles initialisation and finalisation of Hypre library.
 * The first instance which is created initialises Hypre
 * Keeps a count of the number of how many instances exist
 * When the last instance is destroyed it finalises Hypre.
 */ 
class HypreLib {
public:
  /*!
   * Ensure that Hypre has been initialised
   */
  explicit HypreLib();
  
  /*!
   * Calls PetscFinalize when all HypreLib instances are destroyed
   */ 
  ~HypreLib();
  
  /*!
   * Force cleanup. This will call PetscFinalize, printing a warning
   * if any instances of HypreLib still exist
   */ 
  static void cleanup(); 
private:
  static int count; ///< How many instances?
};


#else // BOUT_HAS_HYPRE

#include "unused.hxx"


class HypreLib {
public:
  explicit HypreLib(Options* UNUSED(opt) = nullptr) {}
  ~HypreLib() {}
  
  static void cleanup() {}
};

#endif // BOUT_HAS_HYPRE


#endif //  __HYPRELIB_H__
