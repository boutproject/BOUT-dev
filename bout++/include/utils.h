/**************************************************************************
 * Memory allocation routines
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

#ifndef __UTILS_H__
#define __UTILS_H__

#include "bout_types.h"
#include "dcomplex.h"

#include <string>
#include <list>

BoutReal *rvector(int size);
BoutReal *rvresize(BoutReal *v, int newsize);
int *ivector(int size);
int *ivresize(int *v, int newsize);
BoutReal **rmatrix(int xsize, int ysize);
int **imatrix(int xsize, int ysize);
void free_rmatrix(BoutReal **m);
void free_imatrix(int **m);
BoutReal ***r3tensor(int nrow, int ncol, int ndep);
void free_r3tensor(BoutReal ***m);

dcomplex **cmatrix(int nrow, int ncol);
void free_cmatrix(dcomplex** cm);

BoutReal SQ(BoutReal x);
int ROUND(BoutReal x);
void SWAP(BoutReal &a, BoutReal &b);
void SWAP(BoutReal* &a, BoutReal* &b);
void SWAP(dcomplex &a, dcomplex &b);
void SWAP(int &a, int &b);
int MAX(int a, int b);
bool is_pow2(int x); // Check if a number is a power of 2

/*
BoutReal operator^(BoutReal lhs, int rhs);
BoutReal operator^(BoutReal lhs, const BoutReal &rhs);
*/

char* copy_string(const char* s);
char *strconcat(const char* left, const char *right);

/// Convert a string to lower case
const string lowercase(const string &str);
const string lowercasequote(const string &str);

/// Split a string on a given delimiter
std::list<std::string> &strsplit(const std::string &s, char delim, std::list<std::string> &elems);
std::list<std::string> strsplit(const std::string &s, char delim);

#endif // __UTILS_H__
