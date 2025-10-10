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

#ifndef BOUT_HYPRELIB_H
#define BOUT_HYPRELIB_H

#include "bout/build_defines.hxx"

namespace bout {
#if BOUT_HAS_HYPRE
/*!
 * Handles initialisation and finalisation of Hypre library.
 * The first instance which is created initialises Hypre.
 * Keeps a count of the number of how many instances exist
 * when the last instance is destroyed it finalises Hypre.
 */
class HypreLib {
public:
  explicit HypreLib();

  HypreLib(const HypreLib& other) noexcept;
  HypreLib(HypreLib&& other) noexcept;
  HypreLib& operator=(const HypreLib& other) = default;
  HypreLib& operator=(HypreLib&& other) = default;

  ~HypreLib();

  static void cleanup();

private:
  static int count; ///< How many instances?
};

#else // BOUT_HAS_HYPRE

class HypreLib {
public:
  static void cleanup() {}
};

#endif // BOUT_HAS_HYPRE
} // namespace bout
#endif //  BOUT_HYPRELIB_H
