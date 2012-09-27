/************************************************************************//**
 * \brief Global variables for BOUT++
 * 
 * 
 **************************************************************************
 * Copyright 2010 B.D.Dudson, S.Farley, M.V.Umansky, X.Q.Xu
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

#ifndef __GLOBALS_H__
#define __GLOBALS_H__

#include "datafile.hxx"
#include "bout/mesh.hxx"

#ifndef GLOBALORIGIN
#define GLOBAL extern
#define SETTING(name, val) extern name;
#else
#define GLOBAL
#define SETTING(name, val) name = val;
#endif

SETTING(Mesh *mesh, NULL); ///< The mesh object

///////////////////////////////////////////////////////////////

/// Dump file object
GLOBAL Datafile dump;

// Error handling (bout++.cpp)
void bout_error();
void bout_error(const char *str);

#undef GLOBAL
#undef SETTING

#endif // __GLOBALS_H__
