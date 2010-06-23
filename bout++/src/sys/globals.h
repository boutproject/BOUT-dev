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

#include "mpi.h"

#include "bout_types.h"
#include "field2d.h"
#include "options.h"
#include "output.h"
#include "msg_stack.h" 

#include "datafile.h"
#include "grid.h"
#include "mesh.h"

#ifndef GLOBALORIGIN
#define GLOBAL extern
#define SETTING(name, val) extern name;
#else
#define GLOBAL
#define SETTING(name, val) name = val;
#endif

GLOBAL Mesh *mesh; ///< The mesh object

const real PI = 3.141592653589793;
const real TWOPI = 6.2831853071795;

GLOBAL int MYPE_IN_CORE; // 1 if processor in core

///////////////////////////////////////////////////////////////

/// Options file object
GLOBAL OptionFile options;

/// Define for reading options which passes the variable name
#define OPTION(var, def) options.get(#var, var, def)

/// Output object
GLOBAL Output output;

/// Dump file object
GLOBAL Datafile dump;

/// Status message stack. Used for debugging messages
GLOBAL MsgStack msg_stack;

/// Define for reading a variable from the grid
#define GRID_LOAD(var) mesh->get(var, #var)

// Settings

// Timing information
GLOBAL real wtime_invert; //< Time spent performing inversions

GLOBAL bool non_uniform;  // Use corrections for non-uniform meshes

// Error handling (bout++.cpp)
void bout_error();
void bout_error(const char *str);

#undef GLOBAL
#undef SETTING


/// Concise way to write time-derivatives
#define ddt(f) (*((f).timeDeriv()))

#endif // __GLOBALS_H__
