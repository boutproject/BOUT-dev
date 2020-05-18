/************************************************************************//**
 * \brief Openmp utilitiy wrappers
 * 
 * 
 **************************************************************************
 * Copyright 2017
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

#ifndef __OPENMPWRAP_H__
#define __OPENMPWRAP_H__

//Some helpers for indirection -- required so that the _Pragma gets "omp <x>"
//where <x> is any number of valid omp options/environments (e.g. atomic, critical etc.)
#define INDIRECT0(a) #a
#define INDIRECT1(a) INDIRECT0(omp a)
#define INDIRECT2(b) INDIRECT1(b)

//Define a macro wrapper to the use of `#pragma omp` to avoid unknown pragma
//warnings when compiling without openmp support.
//#ifdef _OPENMP
//#define BOUT_OMP(...) _Pragma(INDIRECT2(__VA_ARGS__))
//#else
#define BOUT_OMP(...)
//#endif

//Perhaps want to cleanup local helpers with below, but DON'T!
//This would cause uses of BOUT_OMP to break
// #undef INDIRECT0
// #undef INDIRECT1
// #undef INDIRECT2
#endif

