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
 * This header tries to workaround some annoying PETSc features, and
 * so it *must* be included before *any* PETSc header.
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

#ifndef BOUT_PETSCLIB_H
#define BOUT_PETSCLIB_H

#include "bout/build_defines.hxx"

class Options;

#if BOUT_HAS_PETSC

// PETSc "helpfully" defines macros for MPI functions that clobber the
// real names, and short of `#undef`-ing all of them in every file
// that includes any PETSc header, we can define the following macro
// which should disable them, which I'm sure will work forever. This
// means we _must_ `#include` this header _before_ any PETSc header!
#define PETSC_HAVE_BROKEN_RECURSIVE_MACRO

#include <petsc.h> // IWYU pragma: export
#include <petscversion.h>

#include "bout/boutexception.hxx"

// NOLINTNEXTLINE(cppcoreguidelines-macro-usage)
#define BOUT_DO_PETSC(cmd) PetscLib::assertIerr((cmd), #cmd)

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
  explicit PetscLib(Options* opt = nullptr);

  PetscLib(const PetscLib&) = default;
  PetscLib(PetscLib&&) = default;
  PetscLib& operator=(const PetscLib&) = default;
  PetscLib& operator=(PetscLib&&) = default;

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
  static void setArgs(int& argc, char**& argv) {
    pargc = &argc;
    pargv = &argv;
  }

  /// Set options for a KSP linear solver that uses the options specific to this PetscLib,
  /// by setting an options prefix for the KSP, and adding that prefix to all the options
  /// set in the [petsc] section, or [petsc] subsection of the options, if non-null 'opt'
  /// was passed to the constructor.
  void setOptionsFromInputFile(KSP& ksp);

  /// Set options for a SNES linear solver that uses the options specific to this PetscLib,
  /// by setting an options prefix for the SNES, and adding that prefix to all the options
  /// set in the [petsc] section, or [petsc] subsection of the options, if non-null 'opt'
  /// was passed to the constructor.
  void setOptionsFromInputFile(SNES& snes);

  /*!
   * Force cleanup. This will call PetscFinalize, printing a warning
   * if any instances of PetscLib still exist
   */
  static void cleanup();

  static inline void assertIerr(PetscErrorCode ierr,
                                const std::string& petsc_op = "PETSc operation") {
    if (ierr != 0) {
      throw BoutException("{:s} failed with {:d}", petsc_op, ierr);
    }
  }

  static BoutException SNESFailure(SNES& snes);

private:
  // NOLINTBEGIN(cppcoreguidelines-avoid-non-const-global-variables)
  static inline int count = 0; ///< How many instances?

  // Command-line arguments
  static inline int* pargc = nullptr;
  static inline char*** pargv = nullptr;

  // Prefix for object-specific options
  std::string options_prefix;

  static inline PetscLogEvent USER_EVENT;
  // NOLINTEND(cppcoreguidelines-avoid-non-const-global-variables)
};

#ifndef PETSC_VERSION_GE
// Newer versions of PETSc define these symbols for testing library version
// This is a re-implementation of the PETSc BSD-licensed code

#define PETSC_VERSION_GE(MAJOR, MINOR, SUBMINOR) \
  ((PETSC_VERSION_MAJOR > MAJOR)                 \
   || ((PETSC_VERSION_MAJOR == MAJOR)            \
       && ((PETSC_VERSION_MINOR > MINOR)         \
           || ((PETSC_VERSION_MINOR == MINOR)    \
               && (PETSC_VERSION_SUBMINOR >= SUBMINOR)))))

#endif // PETSC_VERSION_GE

#else // BOUT_HAS_PETSC

#include "bout/unused.hxx"

// PETSc not available, so KSP and SNES not already defined. KSP and SNES should never be
// called, so forward declaration OK here.
class KSP;
class SNES;

class PetscLib {
public:
  explicit PetscLib(Options* UNUSED(opt) = nullptr) {}
  ~PetscLib() {}

  static void setArgs(int& UNUSED(c), char**& UNUSED(v)) {}

  void setOptionsFromInputFile(KSP& UNUSED(ksp)) {}
  void setOptionsFromInputFile(SNES& UNUSED(snes)) {}

  static void cleanup() {}
};

#endif // BOUT_HAS_PETSC

#endif // BOUT_PETSCLIB_H
