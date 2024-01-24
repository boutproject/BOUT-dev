/*! \file ******************************************************************
 * Provides access to the Hypre library, handling initialisation and
 * finalisation.
 *
 * Usage
 * -----
 *
 *     #include <bout/hyprelib.hxx>
 *
 *     class MyClass {
 *       HypreLib lib;
 *     public:
 *       // ...
 *     };
 *
 *
 * This will then automatically ensure Hypre is initialised the first
 * time an instance is created, and finalise it when the last instances is
 * destroyed (across all objects using `HypreLib`)
 */
/**************************************************************************
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

#ifndef BOUT_HYPRELIB_H
#define BOUT_HYPRELIB_H

#include "bout/build_defines.hxx"

namespace bout {
/*!
 * Handles initialisation and finalisation of Hypre library.
 *
 * The first instance which is created initialises Hypre. Keeps a
 * count of the number of how many instances exist, and when the last
 * instance is destroyed it finalises Hypre.
 */
#if BOUT_HAS_HYPRE
class HypreLib {
public:
  explicit HypreLib();

  HypreLib(const HypreLib& other) noexcept;
  HypreLib(HypreLib&& other) noexcept;
  HypreLib& operator=(const HypreLib& other) = default;
  HypreLib& operator=(HypreLib&& other) = default;

  ~HypreLib();

  /// Immediately finalise Hypre library
  static void cleanup();

private:
  static int count; ///< Current number of instances
};

#else // BOUT_HAS_HYPRE

class HypreLib {
public:
  static void cleanup() {}
};

#endif // BOUT_HAS_HYPRE
} // namespace bout
#endif //  BOUT_HYPRELIB_H
