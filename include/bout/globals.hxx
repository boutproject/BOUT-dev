/************************************************************************/ /**
 * \brief Global variables for BOUT++
 * 
 * 
 **************************************************************************
 * Copyright 2010 - 2025 BOUT++ contributors
 *
 * Contact: Ben Dudson, dudson2@llnl.gov
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

#ifndef BOUT_GLOBALS_H
#define BOUT_GLOBALS_H

#include "bout/macro_for_each.hxx"

class MpiWrapper;

namespace bout {
namespace globals {
#ifndef GLOBALORIGIN
#define GLOBAL extern
#define SETTING(name, val) extern name
#else
#define GLOBAL
#define SETTING(name, val) name = val
#endif

SETTING(MpiWrapper* mpi, nullptr); ///< The MPI wrapper object

///////////////////////////////////////////////////////////////

#undef GLOBAL
#undef SETTING
} // namespace globals
} // namespace bout

#endif // BOUT_GLOBALS_H
