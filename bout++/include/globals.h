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

const BoutReal PI = 3.141592653589793;
const BoutReal TWOPI = 6.2831853071795;

GLOBAL int MYPE_IN_CORE; // 1 if processor in core

///////////////////////////////////////////////////////////////

/// Options file object
GLOBAL OptionFile options;

/// Define for reading options which passes the variable name
#define OPTION(var, def) options.get(#var, var, def)

#define OPTION2(var1, var2, def){ \
    options.get(#var1, var1, def);  \
    options.get(#var2, var2, def);}

#define OPTION3(var1, var2, var3, def){          \
    options.get(#var1, var1, def);               \
    options.get(#var2, var2, def);               \
    options.get(#var3, var3, def);}

#define OPTION4(var1, var2, var3, var4, def){    \
    options.get(#var1, var1, def);               \
    options.get(#var2, var2, def);               \
    options.get(#var3, var3, def);               \
    options.get(#var4, var4, def);}

#define OPTION5(var1, var2, var3, var4, var5, def){     \
    options.get(#var1, var1, def);                      \
    options.get(#var2, var2, def);                      \
    options.get(#var3, var3, def);                      \
    options.get(#var4, var4, def);                      \
    options.get(#var5, var5, def);}

#define OPTION6(var1, var2, var3, var4, var5, var6, def){        \
    options.get(#var1, var1, def);                               \
    options.get(#var2, var2, def);                               \
    options.get(#var3, var3, def);                               \
    options.get(#var4, var4, def);                               \
    options.get(#var5, var5, def);                               \
    options.get(#var6, var6, def);}

/// Macro to replace bout_solve, passing variable name
#define SOLVE_FOR(var) bout_solve(var, #var)
#define SOLVE_FOR2(var1, var2) { \
  bout_solve(var1, #var1);        \
  bout_solve(var2, #var2);}
#define SOLVE_FOR3(var1, var2, var3) {  \
  bout_solve(var1, #var1);             \
  bout_solve(var2, #var2);             \
  bout_solve(var3, #var3);}
#define SOLVE_FOR4(var1, var2, var3, var4) { \
  bout_solve(var1, #var1);             \
  bout_solve(var2, #var2);             \
  bout_solve(var3, #var3);             \
  bout_solve(var4, #var4);}
#define SOLVE_FOR5(var1, var2, var3, var4, var5) {    \
  bout_solve(var1, #var1);             \
  bout_solve(var2, #var2);             \
  bout_solve(var3, #var3);             \
  bout_solve(var4, #var4);             \
  bout_solve(var5, #var5);}
#define SOLVE_FOR6(var1, var2, var3, var4, var5, var6) {      \
  bout_solve(var1, #var1);             \
  bout_solve(var2, #var2);             \
  bout_solve(var3, #var3);             \
  bout_solve(var4, #var4);             \
  bout_solve(var5, #var5);             \
  bout_solve(var6, #var6);}

/// Output object
GLOBAL Output output;

/// Dump file object
GLOBAL Datafile dump;

/// Write this variable once to the grid file
#define SAVE_ONCE(var) dump.add(var, #var, 0)
#define SAVE_ONCE2(var1, var2) { \
    dump.add(var1, #var1, 0); \
    dump.add(var2, #var2, 0);}
#define SAVE_ONCE3(var1, var2, var3) {\
    dump.add(var1, #var1, 0); \
    dump.add(var2, #var2, 0); \
    dump.add(var3, #var3, 0);}
#define SAVE_ONCE4(var1, var2, var3, var4) { \
    dump.add(var1, #var1, 0); \
    dump.add(var2, #var2, 0); \
    dump.add(var3, #var3, 0); \
    dump.add(var4, #var4, 0);}
#define SAVE_ONCE5(var1, var2, var3, var4, var5) {\
    dump.add(var1, #var1, 0); \
    dump.add(var2, #var2, 0); \
    dump.add(var3, #var3, 0); \
    dump.add(var4, #var4, 0); \
    dump.add(var5, #var5, 0);}
#define SAVE_ONCE6(var1, var2, var3, var4, var5, var6) {\
    dump.add(var1, #var1, 0); \
    dump.add(var2, #var2, 0); \
    dump.add(var3, #var3, 0); \
    dump.add(var4, #var4, 0); \
    dump.add(var5, #var5, 0); \
    dump.add(var6, #var6, 0);}

/// Status message stack. Used for debugging messages
GLOBAL MsgStack msg_stack;

/// Define for reading a variable from the grid
#define GRID_LOAD(var) mesh->get(var, #var)
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
    mesh->get(var4, #var4);}
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

// Settings

// Timing information
GLOBAL BoutReal wtime_invert; //< Time spent performing inversions

GLOBAL bool non_uniform;  // Use corrections for non-uniform meshes

// Error handling (bout++.cpp)
void bout_error();
void bout_error(const char *str);

#undef GLOBAL
#undef SETTING


/// Concise way to write time-derivatives
#define ddt(f) (*((f).timeDeriv()))

#endif // __GLOBALS_H__
