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
#include "bout/macro_for_each.hxx"

class Mesh;
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

SETTING(Mesh *mesh, nullptr); ///< The mesh object
SETTING(MpiWrapper *mpi, nullptr); ///< The MPI wrapper object
  
/// Define for reading a variable from the grid
#define GRID_LOAD1(var) mesh->get(var, #var)
#define GRID_LOAD2(var1, var2) {\
    mesh->get(var1, #var1); \
    mesh->get(var2, #var2);}
#define GRID_LOAD3(var1, var2, var3) {\
    mesh->get(var1, #var1); \
    mesh->get(var2, #var2); \
    mesh->get(var3, #var3);}
#define GRID_LOAD4(var1, var2, var3, var4) { \
    mesh->get(var1, #var1); \
    mesh->get(var2, #var2); \
    mesh->get(var3, #var3); \
    mesh->get(var4, #var4); }
#define GRID_LOAD5(var1, var2, var3, var4, var5) {\
    mesh->get(var1, #var1); \
    mesh->get(var2, #var2); \
    mesh->get(var3, #var3); \
    mesh->get(var4, #var4); \
    mesh->get(var5, #var5);}
#define GRID_LOAD6(var1, var2, var3, var4, var5, var6) {\
    mesh->get(var1, #var1); \
    mesh->get(var2, #var2); \
    mesh->get(var3, #var3); \
    mesh->get(var4, #var4); \
    mesh->get(var5, #var5); \
    mesh->get(var6, #var6);}

/// Read fields from the global mesh
/// The name of the variable will be used as the name
/// in the input.
/// This should accept up to 10 arguments
#define GRID_LOAD(...) \
  { MACRO_FOR_EACH_FN(GRID_LOAD1, __VA_ARGS__) }

///////////////////////////////////////////////////////////////

/// Dump file object
GLOBAL Datafile dump;

#undef GLOBAL
#undef SETTING
} // namespace globals
} // namespace bout

#endif // __GLOBALS_H__
