/*************************************************************************
 *
 *                   BOUNDARY CONDITIONS
 * 
 * Similar to the derivs.cpp system, here names and codes are mapped
 * to functions using lookup tables.
 *
 * Changelog:
 *
 * 2010-06 Ben Dudson <bd512@york.ac.uk>
 *      * Adaping to use Mesh interface
 *
 *
 **************************************************************************
 * Copyright 2010 B.D.Dudson, S.Farley, M.V.Umansky, X.Q.Xu
 *
 * Contact Ben Dudson, bd512@york.ac.uk
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

#include "globals.h"
#include "boundary.h"
#include "utils.h"
#include "fft.h"
#include "dcomplex.h"
#include "invert_laplace.h" // For Laplacian coefficients

#include <math.h>

// Define types for boundary functions
typedef void (*bndry_func3d)(Field3D &);
typedef void (*bndry_func2d)(Field2D &);

/*******************************************************************************
 * Boundary functions. These just calculate values in the boundary region
 * based on values in the domain. More complicated boundaries (e.g. relaxation)
 * are then built on top of these.
 *******************************************************************************/

///////////////////////// INNER X ////////////////////////////

void bndry_inner_zero(Field2D &var)
{
  if(!mesh->first_x())
    return;

  for(int jx=0;jx<mesh->xstart;jx++)
    for(int jy=0;jy<mesh->ngy;jy++)
      var[jx][jy] = 0.0;
}

void bndry_inner_zero(Field3D &var)
{
  if(!mesh->first_x())
    return;
  
  for(int jx=0;jx<mesh->xstart;jx++)
    for(int jy=0;jy<mesh->ngy;jy++)
      for(int jz=0;jz<mesh->ngz;jz++)
	var[jx][jy][jz] = 0.0;
}

void bndry_inner_flat(Field2D &var)
{
  if(!mesh->first_x())
    return;
  
  for(int jx=mesh->xstart-1;jx>=0;jx--)
    for(int jy=0;jy<mesh->ngy;jy++)
      var[jx][jy] = var[jx+1][jy]; // Setting this way for relaxed boundaries
  
}

void bndry_inner_flat(Field3D &var)
{
  if(!mesh->first_x())
    return;
  
  if(mesh->ShiftXderivs) // Shift into real space
    var = var.ShiftZ(true);
  
  for(int jx=mesh->xstart-1;jx>=0;jx--)
    for(int jy=0;jy<mesh->ngy;jy++)
      for(int jz=0;jz<mesh->ngz;jz++)
	var[jx][jy][jz] = var[mesh->xstart][jy][jz];
    
  if(mesh->ShiftXderivs) // Shift back
    var = var.ShiftZ(false);
}

///////////////////////// OUTER X ////////////////////////////

void bndry_outer_zero(Field2D &var)
{
  if(!mesh->last_x())
    return;

  for(int jx=mesh->ngx-1;jx>mesh->xend;jx--)
    for(int jy=0;jy<mesh->ngy;jy++)
      var[jx][jy] = 0.0;
}

void bndry_outer_zero(Field3D &var)
{
  if(!mesh->last_x())
    return;

  for(int jx=mesh->ngx-1;jx>mesh->xend;jx--) {
    for(int jy=0;jy<mesh->ngy;jy++) {
      for(int jz=0;jz<mesh->ngz;jz++) {
	var[jx][jy][jz] = 0.0;
      }
    }
  }
}

///////////////////////// UPPER Y ////////////////////////////

void bndry_yup_flat(Field2D &var)
{
  RangeIter* xr = mesh->iterateBndryUpperY();
  for(xr->first(); !xr->isDone(); xr->next())
    for(int jy=mesh->yend+1;jy<mesh->ngy;jy++) {
      
      var[xr->ind][jy] = var[xr->ind][mesh->yend];
      
    }
  delete xr;
}

void bndry_yup_zero(Field2D &var)
{
  RangeIter* xr = mesh->iterateBndryUpperY();
  for(xr->first(); !xr->isDone(); xr->next())
    for(int jy=mesh->yend+1;jy<mesh->ngy;jy++) {
      var[xr->ind][jy] = 0.0;
    }
  delete xr;
}

void bndry_yup_flat(Field3D &var)
{
  RangeIter* xr = mesh->iterateBndryUpperY();
  for(xr->first(); !xr->isDone(); xr->next())
    for(int jy=mesh->yend+1;jy<mesh->ngy;jy++)
      for(int jz=0;jz<mesh->ngz;jz++) {
	var[xr->ind][jy][jz] = var[xr->ind][mesh->yend][jz];
      }
  delete xr;
}

void bndry_yup_zero(Field3D &var)
{
  RangeIter* xr = mesh->iterateBndryUpperY();
  for(xr->first(); !xr->isDone(); xr->next())
    for(int jy=mesh->yend+1;jy<mesh->ngy;jy++)
      for(int jz=0;jz<mesh->ngz;jz++) {
	var[xr->ind][jy][jz] = 0.0;
      }
  delete xr;
}

///////////////////////// LOWER Y ////////////////////////////

void bndry_ydown_flat(Field2D &var)
{
  RangeIter* xr = mesh->iterateBndryLowerY();
  for(xr->first(); !xr->isDone(); xr->next())
    for(int jy=0;jy<mesh->ystart;jy++)
      var[xr->ind][jy] = var[xr->ind][mesh->ystart];
  delete xr;
}

void bndry_ydown_zero(Field2D &var)
{
  RangeIter* xr = mesh->iterateBndryLowerY();
  for(xr->first(); !xr->isDone(); xr->next())
    for(int jy=0;jy<mesh->ystart;jy++)
      var[xr->ind][jy] = 0.0;
  delete xr;
}

void bndry_ydown_flat(Field3D &var)
{
  RangeIter* xr = mesh->iterateBndryLowerY();
  for(xr->first(); !xr->isDone(); xr->next())
    for(int jy=0;jy<mesh->ystart;jy++)
      for(int jz=0;jz<mesh->ngz;jz++) {
	var[xr->ind][jy][jz] = var[xr->ind][mesh->ystart][jz];
      }
  delete xr;
}

void bndry_ydown_zero(Field3D &var)
{
  RangeIter* xr = mesh->iterateBndryLowerY();
  for(xr->first(); !xr->isDone(); xr->next())
    for(int jy=0;jy<mesh->ystart;jy++)
      for(int jz=0;jz<mesh->ngz;jz++) {
	var[xr->ind][jy][jz] = 0.0;
      }
  delete xr;
}

/*******************************************************************************
 * Lookup tables of functions. Map between names, codes and functions
 *******************************************************************************/

/// Enum for boundary type
enum BNDRY_TYPE {BNDRY_NULL,  // Null boundary
		 BNDRY_Z,  // Zero value
		 BNDRY_ZERO_GRAD   // Zero gradient
};

/// Enum for boundary location. Declared as a bit-field so can be or'd together
enum BNDRY_LOC {INNER_X = 1, OUTER_X = 2, UPPER_Y = 4, LOWER_Y = 8};

/// Translate between short names, long names, and BNDRY_TYPE
struct BndryNameLookup {
  BNDRY_TYPE type;
  const char *label; ///< Short name used in BOUT.inp file
  const char *name;  ///< Long descriptive name
};

/// Find a function for a particular boundary type
struct BndryLookup {
  BNDRY_TYPE type;
  bndry_func3d func3d; ///< 3D function (can be NULL)
  bndry_func2d func2d; ///< 2D function (can be NULL)
};

/////////////// LOOKUP TABLES ////////////////

/// Translate between short names, long names and BNDRY_TYPE codes
static BndryNameLookup BndryNameTable[] = { {BNDRY_Z,         "ZERO",      "Zero value"},
					    {BNDRY_ZERO_GRAD, "ZERO_GRAD", "Zero gradient"},
					    {BNDRY_NULL} }; // Terminate list

static BndryLookup InnerBndryTable[] = { {BNDRY_Z,         bndry_inner_zero, bndry_inner_zero},
					 {BNDRY_ZERO_GRAD, bndry_inner_flat, bndry_inner_flat},
					 {BNDRY_NULL} };

static BndryLookup OuterBndryTable[] = { {BNDRY_Z, bndry_inner_zero, bndry_inner_zero},
					 {BNDRY_NULL} };

static BndryLookup LowerBndryTable[] = { {BNDRY_Z,          bndry_ydown_zero, bndry_ydown_zero},
					 {BNDRY_ZERO_GRAD,  bndry_ydown_flat, bndry_ydown_flat},
					 {BNDRY_NULL} };

static BndryLookup UpperBndryTable[] = { {BNDRY_Z,          bndry_yup_zero, bndry_yup_zero},
					 {BNDRY_ZERO_GRAD,  bndry_yup_flat, bndry_yup_flat},
					 {BNDRY_NULL} };
/////////////// BOUNDARY NAMES ////////////////

/// Describes a boundary
struct BndryDesc {
  BndryLookup *table; ///< The lookup table
  BNDRY_LOC loc;      ///< Boundary location code
  const char *general; ///< The general name used in BOUT.inp (doesn't change)
  char *code; ///< Short code for BOUT.inp file
  char *desc; ///< Longer description
};

/// One description for each boundary
static BndryDesc InnerBndry = {InnerBndryTable, INNER_X, "xinner", NULL, NULL}, 
  OuterBndry = {OuterBndryTable, OUTER_X, "xouter", NULL, NULL}, 
  UpperBndry = {UpperBndryTable, UPPER_Y, "yupper", NULL, NULL}, 
  LowerBndry = {LowerBndryTable, LOWER_Y, "ylower", NULL, NULL};

/// Set boundary codes. Allows a simulation to customise
/// different regions of the domain. E.g. for tokamaks, specify 'PF' or 'CORE'
void BndrySetName(BNDRY_LOC loc, const char *code, const char *desc = NULL)
{
  BndryDesc *b;
  
  switch(loc) {
  case INNER_X:
    b = &InnerBndry;
    break;
  case OUTER_X:
    b = &OuterBndry;
    break;
  case UPPER_Y:
    b = &UpperBndry;
    break;
  case LOWER_Y:
    b = &LowerBndry;
    break;
  default:
    bout_error("Unrecognised boundary to BndrySet");
    return;
  }
  
  b->code = copy_string(code);
  b->desc = copy_string(desc);
}

/*******************************************************************************
 * Routines using the above tables to map between names, codes, and functions
 *******************************************************************************/

/// Lookup a function operating on Field3D
bndry_func3d BndryLookup3D(BndryLookup *table, BNDRY_TYPE type)
{
  int i = 0;
  
  if(type == BNDRY_NULL)
    return table[0].func3d; // Return the first entry
  
  do {
    if(table[i].type == type) {
      // Matched type. Check if function non-NULL
      if(table[i].func3d != NULL)
	return table[i].func3d;
      // otherwise, return the first entry
      return table[0].func3d;
    }
    i++;
  }while(table[i].type != BNDRY_NULL);
  // Not found in list. Return the first entry
  
  return table[0].func3d;
}

bndry_func3d BndryLookup3D(BNDRY_LOC loc, BNDRY_TYPE type)
{
  BndryLookup *table;
  
  switch(loc) {
  case INNER_X: {
    table = InnerBndryTable;
    break;
  }
  case OUTER_X: {
    table = OuterBndryTable;
    break;
  }
  case UPPER_Y: {
    break;
  }
  case LOWER_Y: {
    break;
  }
  default: {
    // Don't know what to do
    bout_error("BndryLookup3D(BNDRY_LOC, BNDRY_TYPE) passed multiple locations");
    exit(1);
  }
  }

  return BndryLookup3D(table, type);
}

/// Lookup a function operating on Field2D
bndry_func2d BndryLookup2D(BndryLookup *table, BNDRY_TYPE type)
{
  int i = 0;
  
  if(type == BNDRY_NULL)
    return table[0].func2d; // Return the first entry
  
  do {
    if(table[i].type == type) {
      // Matched type. Check if function non-NULL
      if(table[i].func2d != NULL)
	return table[i].func2d;
      // otherwise, return the first entry
      return table[0].func2d;
    }
    i++;
  }while(table[i].type != BNDRY_NULL);
  // Not found in list. Return the first entry
  
  return table[0].func2d;
}

bndry_func2d BndryLookup2D(BNDRY_LOC loc, BNDRY_TYPE type)
{
  BndryLookup *table;
  
  switch(loc) {
  case INNER_X: {
    table = InnerBndryTable;
    break;
  }
  case OUTER_X: {
    table = OuterBndryTable;
    break;
  }
  case UPPER_Y: {
    break;
  }
  case LOWER_Y: {
    break;
  }
  default: {
    // Don't know what to do
    bout_error("BndryLookup2D(BNDRY_LOC, BNDRY_TYPE) passed multiple locations");
    exit(1);
  }
  }

  return BndryLookup2D(table, type);
}

/// Test if a boundary method is implemented for a given boundary
bool isImplemented(BndryLookup *table, BNDRY_TYPE type)
{
  int i = 0;
  do {
    if(table[i].type == type)
      return true;
    i++;
  }while(table[i].type != BNDRY_NULL);
  
  return false;
}

/// Look up a boundary code in the tables
BNDRY_TYPE BndryLookup3D(BndryLookup *table, const char *label, bool dummy = false)
{
  if(label == NULL)
    return table[0].type;
  if(label[0] == '\0')
    return table[0].type;

  // Loop through the name lookup table
  int i = 0;
  do {
    
  }while(BndryNameTable[i].type != BNDRY_NULL);

  return BNDRY_NULL; // Not implemented yet
}

/// Look up a boundary code in the tables
BNDRY_TYPE BndryLookup2D(BndryLookup *table, const char *label, bool dummy = false)
{
  return BNDRY_NULL; // Not implemented yet
}

/*******************************************************************************
 * Relaxation boundary conditions. Convert a non-relaxing boundary into a
 * relaxing one. Probably not the most efficient system, but general.
 *******************************************************************************/

/// Turn an inner boundary condition into a relaxing boundary condition
/// NOTE: CURRENTLY THIS WON'T WORK AS ADVERTISED
void BndryRelax(Field3D &var, Field3D &F_var, real tconst, bndry_func3d func, BNDRY_LOC loc)
{
  Field3D tmpvar;

  if(tconst <= 0.0) {
    // Not relaxing. It's unlikely that the user wants a negative time-constant.
    // Just apply the boundary to the time-derivative
    func(F_var);
    return;
  }
  // Make a copy of the variable
  tmpvar = var;
  
  // Apply the boundary condition to the variable copy
  func(tmpvar);
  
  // Now make the variable relax to this value
  
  int jx,jy,jz;

  switch(loc) {
  case INNER_X: {
    for(jx=0;jx<mesh->xstart;jx++)
      for(jy=0;jy<mesh->ngy;jy++)
	for(jz=0;jz<mesh->ngz;jz++)
	  F_var[jx][jy][jz] = (tmpvar[jx][jy][jz] - var[jx][jy][jz]) / tconst;
    break;
  }
  case OUTER_X: {
    for(jx=mesh->ngx-1;jx>mesh->xend;jx--)
      for(jy=0;jy<mesh->ngy;jy++)
	for(jz=0;jz<mesh->ngz;jz++)
	  F_var[jx][jy][jz] = (tmpvar[jx][jy][jz] - var[jx][jy][jz]) / tconst;
  }
  case UPPER_Y: {
    RangeIter *xi = mesh->iterateBndryUpperY();
    /// Iterate over boundary
    for(xi->first(); !xi->isDone(); xi->next())
      for(jy=mesh->ngy-1;jy > mesh->yend;jy--)
	for(jz=0;jz<mesh->ngz;jz++)
	  F_var[xi->ind][jy][jz] = (tmpvar[xi->ind][jy][jz] - var[xi->ind][jy][jz]) / tconst;
    delete xi;
    break;
  }
  case LOWER_Y: {
    RangeIter *xi = mesh->iterateBndryLowerY();
    /// Iterate over boundary
    for(xi->first(); !xi->isDone(); xi->next())
      for(jy=0;jy<mesh->ystart;jy++)
	for(jz=0;jz<mesh->ngz;jz++)
	  F_var[xi->ind][jy][jz] = (tmpvar[xi->ind][jy][jz] - var[xi->ind][jy][jz]) / tconst;
    delete xi;
    break;
  }
  }
}

void BndryRelax(Field2D &var, Field2D &F_var, real tconst, bndry_func2d func, BNDRY_LOC loc)
{
  Field2D tmpvar;

  if(tconst <= 0.0) {
    // Not relaxing. It's unlikely that the user wants a negative time-constant.
    // Just apply the boundary to the variable
    func(var);
  }
  // Make a copy of the variable
  tmpvar = var;
  
  // Apply the boundary condition to the copy
  func(tmpvar);
  
  // Now make the variable relax to this value
  
  int jx,jy;

  switch(loc) {
  case INNER_X: {
    for(jx=0;jx<mesh->xstart;jx++)
      for(jy=0;jy<mesh->ngy;jy++)
	F_var[jx][jy] = (tmpvar[jx][jy] - var[jx][jy]) / tconst;
    break;
  }
  case OUTER_X: {
    for(jx=mesh->ngx-1;jx>mesh->xend;jx--)
      for(jy=0;jy<mesh->ngy;jy++)
	F_var[jx][jy] = (tmpvar[jx][jy] - var[jx][jy]) / tconst;
  }
  case UPPER_Y: {
    RangeIter *xi = mesh->iterateBndryUpperY();
    /// Iterate over boundary
    for(xi->first(); !xi->isDone(); xi->next())
      for(jy=mesh->ngy-1;jy > mesh->xend;jy--)
	F_var[xi->ind][jy] = (tmpvar[xi->ind][jy] - var[xi->ind][jy]) / tconst;
    delete xi;
    break;
  }
  case LOWER_Y: {
    RangeIter *xi = mesh->iterateBndryLowerY();
    /// Iterate over boundary
    for(xi->first(); !xi->isDone(); xi->next())
      for(jy=0;jy<mesh->ystart;jy++)
	F_var[xi->ind][jy] = (tmpvar[xi->ind][jy] - var[xi->ind][jy]) / tconst;
    delete xi;
    break;
  }
  }
}

/*******************************************************************************
 * Intermediate-level interface functions
 * 
 * These take a boundary location and type code, and apply to 2D and 3D fields.
 *******************************************************************************/

void apply_boundary(Field3D &var, Field3D &F_var, BNDRY_LOC loc, BNDRY_TYPE type, real tconst)
{
  if(type == BNDRY_NULL)
    return;

  BndryLookup *table; // Table of functions to use

  switch(loc) {
  case INNER_X: {
    table = InnerBndryTable;
    break;
  }
  case OUTER_X: {
    table = OuterBndryTable;
    break;
  }
  case UPPER_Y: {
    table = UpperBndryTable;
    break;
  }
  case LOWER_Y: {
    table = LowerBndryTable;
    break;
  }
  default: {
    // More than one boundary or'd together
    int mask = 1;
    for(int i=0;i<3;i++) {
      if(loc & mask) {
	// Apply this boundary
	apply_boundary(var, F_var, (BNDRY_LOC) mask, type, tconst);
      }
      mask *= 2;
    }
    return;
  }
  }
  
  // Get the function to apply
  bndry_func3d func = BndryLookup3D(table, type);

  // Apply the boundary condition
  BndryRelax(var, F_var, tconst, func, loc);
}

void apply_boundary(Field2D &var, Field2D &F_var, BNDRY_LOC loc, BNDRY_TYPE type, real tconst)
{
  if(type == BNDRY_NULL)
    return;

  BndryLookup *table; // Table of functions to use

  switch(loc) {
  case INNER_X: {
    table = InnerBndryTable;
    break;
  }
  case OUTER_X: {
    table = OuterBndryTable;
    break;
  }
  case UPPER_Y: {
    table = UpperBndryTable;
    break;
  }
  case LOWER_Y: {
    table = LowerBndryTable;
    break;
  }
  default: {
    // More than one boundary or'd together
    int mask = 1;
    for(int i=0;i<3;i++) {
      if(loc & mask) {
	// Apply this boundary
	apply_boundary(var, F_var, (BNDRY_LOC) mask, type, tconst);
      }
      mask *= 2;
    }
    return;
  }
  }
  
  // Get the function to apply
  bndry_func2d func = BndryLookup2D(table, type);

  // Apply the boundary condition
  BndryRelax(var, F_var, tconst, func, loc);
}

/*******************************************************************************
 * High-level interface (ones actually used by rest of code)
 * 
 * Each boundary can be given a code
 *
 * Inner boundary
 *   [name] / <inner code>  (e.g. [P] / pf    )
 *   [name] / inner         (e.g. [P] / inner )
 *   [All]  / <inner code>  (e.g. [All] / pf  )
 *   [All]  / inner
 *  None found  -  No boundary condition applied
 *
 * Relaxation time-constants are also read from the input file, 
 * with the following precedance
 * 
 *   [name] / <inner code>_tconst   (e.g. [P] / pf_tconst    )
 *   [name] / inner_tconst          (e.g. [P] / inner_tconst )
 *   [name] / bndry_tconst          (e.g. [P] / bndry_tconst )
 *   [All]  / <inner code>_tconst   (e.g. [All] / pf_tconst  )
 *   [All]  / inner_tconst          (e.g. [All] / inner_tconst )
 *   [All]  / bndry_tconst          (e.g. [All] / bndry_tconst )
 *  None found - constant set to -1
 *
 * Setting a negative relaxation time-constant corresponds to a non-relaxing boundary condition.
 *
 * There are also (optional) alternative names. These are for vectors, 
 * so for V_x (for example), the order would be
 *  - Check options under [V_x]
 *  - Check under [V]
 *  - Check under [All]
 * 
 * Similarly for the other boundaries
 *******************************************************************************/

/// Get a boundary string from the options file. Starts looking at the very specific, then the general.
char* getBndryString(const char *name, const char *altname, BndryDesc *bndry, bool verbose = false)
{
  char *str;
  if((str = options.getString(name, bndry->code)) != NULL) {
    if(verbose) output.write("\tOption %s / %s = %s\n", name, bndry->code, str);
  }else if((str = options.getString(name, bndry->general)) != NULL) {
    if(verbose) output.write("\tOption %s / %s = %s\n", name, bndry->general, str);

  }else if((str = options.getString(altname, bndry->code)) != NULL) {
    if(verbose) output.write("\tOption %s / %s = %s\n", altname, bndry->code, str);
  }else if((str = options.getString(altname, bndry->general)) != NULL) {
    if(verbose) output.write("\tOption %s / %s = %s\n", altname, bndry->general, str);

  }else if((str = options.getString("All", bndry->code)) != NULL) {
    if(verbose) output.write("\tOption All / %s = %s\n", bndry->code, str);
  }else if((str = options.getString("All", bndry->general)) != NULL) {
    if(verbose) output.write("\tOption All / %s = %s\n", bndry->general, str);
  }
  
  return str;
}

/// Get relaxation constant for the boundary
real getBndryRelax(const char *name, const char *altname, BndryDesc *bndry, bool verbose = false)
{
  real tconst = -1.0;
  char *s;
  
  if(options.getReal(name, s = strconcat(bndry->code, "_tconst"), tconst) == 0) {
    if(verbose) output.write("\t  Option %s / %s = %e\n", name, s, tconst);
  }else if(options.getReal(name, s = strconcat(bndry->general, "_tconst"), tconst) == 0) {
    if(verbose) output.write("\t  Option %s / %s = %e\n", name, s, tconst);
  }else if(options.getReal(name, "bndry_tconst", tconst) == 0) {
    if(verbose) output.write("\t  Option %s / bndry_tconst = %e\n", name, tconst); 

  }else if(options.getReal(altname, s = strconcat(bndry->code, "_tconst"), tconst) == 0) {
    if(verbose) output.write("\t  Option %s / %s = %e\n", altname, s, tconst);
  }else if(options.getReal(altname, s = strconcat(bndry->general, "_tconst"), tconst) == 0) {
    if(verbose) output.write("\t  Option %s / %s = %e\n", altname, s, tconst);
  }else if(options.getReal(altname, "bndry_tconst", tconst) == 0) {
    if(verbose) output.write("\t  Option %s / bndry_tconst = %e\n", altname, tconst); 

  }else if(options.getReal("All", s = strconcat(bndry->code, "_tconst"), tconst) == 0) {
    if(verbose) output.write("\t  Option All / %s = %e\n", s, tconst);
  }else if(options.getReal("All", s = strconcat(bndry->general, "_tconst"), tconst) == 0) {
    if(verbose) output.write("\t  Option All / %s = %e\n", s, tconst);
  }else if(options.getReal("All", "bndry_tconst", tconst) == 0) {
    if(verbose) output.write("\t  Option All / bndry_tconst = %e\n", tconst); 

  }
  
  return tconst;
}

//////// 3D fields /////////

void apply_boundary(Field3D &var, Field3D &F_var, const char* name, const char* altname,
		    BndryDesc *bndry, bool dummy = false)
{
  const char *str;
  real tconst;
  BNDRY_TYPE type;
  
  if(dummy) {
    // Find a string to refer to this boundary
    if((str = bndry->desc) == NULL)
      if((str = bndry->code) == NULL)
	str = bndry->general;
      
    output.write("\t %s boundary\n", str);
  }
  
  // Get the boundary condition string from options file
  if((str = getBndryString(name, altname, bndry, dummy)) == NULL) {
    if(dummy) output.write("\t**WARNING: No boundary condition applied\n");
  }else {
    
    // Look up boundary type code
    type = BndryLookup3D(bndry->table, str, dummy);
    
    // Get time constant
    if( ((tconst = getBndryRelax(name, altname, bndry, dummy)) > 0.) && dummy )
      output.write("\t->Relaxing boundary\n");
    
    if(!dummy) {
      // Not a dummy run - apply the boundary condition
      
      apply_boundary(var, F_var, bndry->loc, type, tconst);
    }
  }
}

/// Apply a boundary condition to a variable, depending on input options.
/// Setting dummy=true prints out the boundary, rather than applying it.
void apply_boundary(Field3D &var, Field3D &F_var, const char* name, const char* altname, bool dummy = false)
{
  if(dummy) output.write("\tBoundary condition for '%s'\n", name);

  /// Inner boundary
  apply_boundary(var, F_var, name, altname, &InnerBndry, dummy);
  
  /// Outer boundary
  apply_boundary(var, F_var, name, altname, &OuterBndry, dummy);
  
  /// Upper boundary
  apply_boundary(var, F_var, name, altname, &UpperBndry, dummy);
  
  /// Lower boundary
  apply_boundary(var, F_var, name, altname, &LowerBndry, dummy);
}

void apply_boundary(Field3D &var, Field3D &F_var, const char* name, bool dummy = false)
{
  apply_boundary(var, F_var, name, name, dummy);
}

//////// 2D fields /////////

void apply_boundary(Field2D &var, Field2D &F_var, const char* name, const char* altname,
		    BndryDesc *bndry, bool dummy = false)
{
  const char *str;
  real tconst;
  BNDRY_TYPE type;
  
  if(dummy) {
    // Find a string to refer to this boundary
    if((str = bndry->desc) == NULL)
      if((str = bndry->code) == NULL)
	str = bndry->general;
      
    output.write("\t %s boundary\n", str);
  }
  
  // Get the boundary condition string from options file
  if((str = getBndryString(name, altname, bndry, dummy)) == NULL) {
    if(dummy) output.write("\t**WARNING: No boundary condition applied\n");
  }else {
    
    // Look up boundary type code
    type = BndryLookup2D(bndry->table, str, dummy);
    
    // Get time constant
    if( ((tconst = getBndryRelax(name, altname, bndry, dummy)) > 0.) && dummy )
      output.write("\t->Relaxing boundary\n");
    
    if(!dummy) {
      // Not a dummy run - apply the boundary condition
      
      apply_boundary(var, F_var, bndry->loc, type, tconst);
    }
  }
}

/// Apply a boundary condition to a variable, depending on input options.
/// Setting dummy=true prints out the boundary, rather than applying it.
void apply_boundary(Field2D &var, Field2D &F_var, const char* name, const char* altname, bool dummy = false)
{
  if(dummy) output.write("\tBoundary condition for '%s'\n", name);

  /// Inner boundary
  apply_boundary(var, F_var, name, altname, &InnerBndry, dummy);
  
  /// Outer boundary
  apply_boundary(var, F_var, name, altname, &OuterBndry, dummy);
  
  /// Upper boundary
  apply_boundary(var, F_var, name, altname, &UpperBndry, dummy);
  
  /// Lower boundary
  apply_boundary(var, F_var, name, altname, &LowerBndry, dummy);
}

void apply_boundary(Field2D &var, Field2D &F_var, const char* name, bool dummy = false)
{
  apply_boundary(var, F_var, name, name, dummy);
}

//////// 3D Vectors /////////

void apply_boundary(Vector3D &var, Vector3D &F_var, const char* name, bool dummy = false)
{
  
}

/*******************************************************************************
 * OBSOLETE CODE. KEPT FOR BACKWARDS-COMPATIBILITY (FOR NOW)
 *******************************************************************************/

// Dummy run - just prints out the boundary conditions
void print_boundary(const char* fullname, const char* shortname)
{
  int opt;

  if(options.getInt(fullname, "xinner", opt)) // Look in the section with the variable name
    if(options.getInt(shortname, "xinner", opt))
      if(options.getInt("All", "xinner", opt)) // Otherwise look for global settings
	opt = 0; // Do nothing

  switch(opt) {
  case BNDRY_NONE:
    output.write("%s inner x boundary: NONE\n", fullname);
    break;
  case BNDRY_ZERO: {
    output.write("%s inner x boundary: Zero value\n", fullname);
    break;
  }case BNDRY_GRADIENT: {
    output.write("%s inner x boundary: Zero gradient\n", fullname);
    break;
  }
  case BNDRY_LAPLACE: {
    output.write("%s inner x boundary: Zero Laplacian\n", fullname);
    break;
  }
  case BNDRY_LAPLACE_GRAD: {
    output.write("%s inner x boundary: Zero Laplacian + zero gradient\n", fullname);
    break;
  }
  case BNDRY_DIVCURL: {
    output.write("%s inner x boundary: Div = 0, Curl = 0\n", fullname);
    break;
  }
  case BNDRY_LAPLACE_ZERO: {
    output.write("%s inner x boundary: Zero laplacian + zero value\n", fullname);
    break;
  }
  case BNDRY_LAPLACE_DECAY: {
    output.write("%s inner x boundary: Zero Laplacian decaying solution\n", fullname);
    break;
  }
  case BNDRY_C_LAPLACE_DECAY: {
    output.write("%s inner x boundary: Constant Laplacian decaying solution\n", fullname);
    break;
  }
  default: {
    output.write("Error: Invalid option %d for inner x boundary\n", opt);
    bout_error("Aborting\n");
  }
  };

  // Get outer x option
  if(options.getInt(fullname, "xouter", opt)) // Look in the section with the variable name
    if(options.getInt(shortname, "xouter", opt))
      if(options.getInt("All", "xouter", opt)) // Otherwise look for global settings
	opt = 0; // Do nothing

  switch(opt) {
  case BNDRY_NONE:
    output.write("%s outer x boundary: NONE\n", fullname);
    break;
  case BNDRY_ZERO: { 
    output.write("%s outer x boundary: Zero value\n", fullname);
    break;
  }case BNDRY_GRADIENT: {
    output.write("%s outer x boundary: Zero gradient\n", fullname);
    break;
  }
  case BNDRY_LAPLACE: {
    output.write("%s outer x boundary: Zero Laplacian\n", fullname);
    break;
  }
  case BNDRY_DIVCURL: {
    output.write("%s outer x boundary: Div = 0, Curl = 0\n", fullname);
    break;
  }
  case BNDRY_LAPLACE_DECAY: {
    output.write("%s outer x boundary: Zero Laplacian decaying solution\n", fullname);
    break;
  }
  case BNDRY_C_LAPLACE_DECAY: {
    output.write("%s outer x boundary: Constant Laplacian decaying solution\n", fullname);
    break;
  }
  default: {
    output.write("Error: Invalid option %d for outer x boundary\n", opt);
    exit(1);
  }
  };

  // Get lower y option
  if(options.getInt(fullname, "ylower", opt)) // Look in the section with the variable name
    if(options.getInt(shortname, "ylower", opt))
      if(options.getInt("All", "ylower", opt)) // Otherwise look for global settings
	opt = 0; // Do nothing

  switch(opt) {
  case BNDRY_NONE:
    output.write("%s lower y boundary: NONE\n", fullname);
    break;
  case BNDRY_ZERO: { 
    output.write("%s lower y boundary: Zero value\n", fullname);
    break;
  }case BNDRY_GRADIENT: {
    output.write("%s lower y boundary: Zero Gradient\n", fullname);
    break;
  }
  case BNDRY_ROTATE: {
    output.write("%s lower y boundary: Rotate 180 degrees\n", fullname);
    break;
  }
  case BNDRY_ZAVERAGE: {
    output.write("%s lower y boundary: Z Average\n", fullname);
    break;
  }
  case BNDRY_ROTATE_NEG: {
    output.write("%s lower y boundary: Rotate 180 degrees and reverse sign\n", fullname);
    break;
  }
  default: {
    output.write("Error: Invalid option %d for lower y boundary\n", opt);
    exit(1);
  }
  };

  // Get upper y option
  if(options.getInt(fullname, "yupper", opt)) // Look in the section with the variable name
    if(options.getInt(shortname, "yupper", opt))
      if(options.getInt("All", "yupper", opt)) // Otherwise look for global settings
	opt = 0; // Do nothing

  switch(opt) {
  case BNDRY_NONE:
    output.write("%s upper y boundary: NONE\n", fullname);
    break;
  case BNDRY_ZERO: { 
    output.write("%s upper y boundary: Zero value\n", fullname);
    break;
  }case BNDRY_GRADIENT: {
    output.write("%s upper y boundary: Zero gradient\n", fullname);
    break;
  }
  default: {
    output.write("Error: Invalid option %d for upper y boundary\n", opt);
    exit(1);
  }
  };
}

// For fields
void print_boundary(const char* name)
{
  print_boundary(name, name);
}

// For vectors
void print_boundary(const char* name, bool covariant)
{
  static char buffer[128];
  // X component

  if(covariant) {
    sprintf(buffer, "%s_x", name);
  }else
    sprintf(buffer, "%sx", name);

  print_boundary(buffer, name);

  // Y component
  
  if(covariant) {
    sprintf(buffer, "%s_y", name);
  }else
    sprintf(buffer, "%sy", name);

  print_boundary(buffer, name);

  // Z component
  
  if(covariant) {
    sprintf(buffer, "%s_z", name);
  }else
    sprintf(buffer, "%sz", name);

  print_boundary(buffer, name);
}

/// Applies a boundary condition, depending on setting in BOUT.inp
void apply_boundary(Field3D &var, const char* fullname, const char* shortname)
{
  int opt;

  //output.write("bndry...");

  // Get inner x option
  if(options.getInt(fullname, "xinner", opt)) // Look in the section with the variable name
    if(options.getInt(shortname, "xinner", opt))
      if(options.getInt("All", "xinner", opt)) // Otherwise look for global settings
	opt = BNDRY_NONE; // Do nothing

  switch(opt) {
  case BNDRY_NONE:
    break;
  case BNDRY_ZERO: { 
    bndry_inner_zero(var);
    break;
  }case BNDRY_GRADIENT: {
    bndry_inner_flat(var);
    break;
  }
  case BNDRY_LAPLACE: {
    bndry_core_laplace2(var);
    bndry_pf_laplace(var);
    break;
  }
  case BNDRY_LAPLACE_GRAD: {
    bndry_inner_laplace(var);
    break;
  }
  case BNDRY_DIVCURL: {
    break;
  }
  case BNDRY_LAPLACE_ZERO: {
    bndry_inner_zero_laplace(var);
    break;
  }
  case BNDRY_LAPLACE_DECAY: {
    bndry_inner_laplace_decay(var);
    break;
  }
  case BNDRY_C_LAPLACE_DECAY: {
    bndry_inner_const_laplace_decay(var);
    break;
  }
  default: {
    output.write("Error: Invalid option %d for inner x boundary\n", opt);
    exit(1);
  }
  };

  // Get outer x option
  if(options.getInt(fullname, "xouter", opt)) // Look in the section with the variable name
    if(options.getInt(shortname, "xouter", opt))
      if(options.getInt("All", "xouter", opt)) // Otherwise look for global settings
	opt = BNDRY_NONE; // Do nothing

  switch(opt) {
  case BNDRY_NONE:
    break;
  case BNDRY_ZERO: { 
    bndry_sol_zero(var);
    break;
  }case BNDRY_GRADIENT: {
    bndry_sol_flat(var);
    break;
  }
  case BNDRY_LAPLACE: {
    bndry_sol_laplace(var);
    break;
  }
  case BNDRY_LAPLACE_DECAY: {
    bndry_outer_laplace_decay(var);
    break;
  }
  case BNDRY_C_LAPLACE_DECAY: {
    bndry_outer_laplace_decay(var);
    break;
  }
  case BNDRY_DIVCURL: {
    break;
  }
  default: {
    output.write("Error: Invalid option %d for outer x boundary\n", opt);
    bout_error("Aborting");
  }
  };

  // Get lower y option
  if(options.getInt(fullname, "ylower", opt)) // Look in the section with the variable name
    if(options.getInt(shortname, "ylower", opt))
      if(options.getInt("All", "ylower", opt)) // Otherwise look for global settings
	opt = BNDRY_NONE; // Do nothing

  switch(opt) {
  case BNDRY_NONE:
    break;
  case BNDRY_ZERO: { 
    bndry_ydown_zero(var);
    break;
  }case BNDRY_GRADIENT: {
    bndry_ydown_flat(var);
    break;
  }
  case BNDRY_ROTATE: {
    bndry_ydown_rotate(var);
    break;
  }
  case BNDRY_ZAVERAGE: {
    bndry_ydown_zaverage(var);
    break;
  }
  case BNDRY_ROTATE_NEG: {
    bndry_ydown_rotate(var, true);
    break;
  }
  default: {
    output.write("Error: Invalid option %d for lower y boundary\n", opt);
    exit(1);
  }
  };

  //output.write("ly ");

  // Get upper y option
  if(options.getInt(fullname, "yupper", opt)) // Look in the section with the variable name
    if(options.getInt(shortname, "yupper", opt))
      if(options.getInt("All", "yupper", opt)) // Otherwise look for global settings
	opt = BNDRY_NONE; // Do nothing

  switch(opt) {
  case BNDRY_NONE:
    break;
  case BNDRY_ZERO: { 
    bndry_yup_zero(var);
    break;
  }case BNDRY_GRADIENT: {
    bndry_yup_flat(var);
    break;
  }
  default: {
    output.write("Error: Invalid option %d for upper y boundary\n", opt);
    exit(1);
  }
  };

  //output.write("uy ");

  bndry_toroidal(var);

  //output.write("done\n");
}

void apply_boundary(Field2D &var, const char* fullname, const char* shortname)
{
  int opt;

  // Get inner x option
  if(options.getInt(fullname, "xinner", opt)) // Look in the section with the variable name
    if(options.getInt(shortname, "xinner", opt))
      if(options.getInt("All", "xinner", opt)) // Otherwise look for global settings
	opt = BNDRY_NONE; // Do nothing

  switch(opt) {
  case BNDRY_NONE:
    break;
  case BNDRY_ZERO: { 
    bndry_inner_zero(var);
    break;
  }case BNDRY_GRADIENT: {
    bndry_inner_flat(var);
    break;
  }
  default: {
    output.write("Error: Invalid option %d for inner x boundary of %s\n", opt, fullname);
    exit(1);
  }
  };

  // Get outer x option
  if(options.getInt(fullname, "xouter", opt)) // Look in the section with the variable name
    if(options.getInt(shortname, "xouter", opt))
      if(options.getInt("All", "xouter", opt)) // Otherwise look for global settings
	opt = BNDRY_NONE; // Do nothing

  switch(opt) {
  case BNDRY_NONE:
    break;
  case BNDRY_ZERO: { 
    bndry_sol_zero(var);
    break;
  }case BNDRY_GRADIENT: {
    bndry_sol_flat(var);
    break;
  }
  default: {
    output.write("Error: Invalid option %d for outer x boundary of %s\n", opt, fullname);
    exit(1);
  }
  };
  
  // Get lower y option
  if(options.getInt(fullname, "ylower", opt)) // Look in the section with the variable name
    if(options.getInt(shortname, "ylower", opt))
      if(options.getInt("All", "ylower", opt)) // Otherwise look for global settings
	opt = BNDRY_NONE; // Do nothing

  switch(opt) {
  case BNDRY_NONE:
    break;
  case BNDRY_ZERO: { 
    bndry_ydown_zero(var);
    break;
  }case BNDRY_GRADIENT: {
    bndry_ydown_flat(var);
    break;
  }
  default: {
    output.write("Error: Invalid option %d for lower y boundary of %s\n", opt, fullname);
    exit(1);
  }
  };

  // Get upper y option
  if(options.getInt(fullname, "yupper", opt)) // Look in the section with the variable name
    if(options.getInt(shortname, "yupper", opt))
      if(options.getInt("All", "yupper", opt)) // Otherwise look for global settings
	opt = BNDRY_NONE; // Do nothing

  switch(opt) {
  case BNDRY_NONE:
    break;
  case BNDRY_ZERO: { 
    bndry_yup_zero(var);
    break;
  }case BNDRY_GRADIENT: {
    bndry_yup_flat(var);
    break;
  }
  default: {
    output.write("Error: Invalid option %d for upper y boundary of %s\n", opt, fullname);
    exit(1);
  }
  };
}

void apply_boundary(Field2D &var, const char* name)
{
  apply_boundary(var, name, name);
}

void apply_boundary(Field3D &var, const char* name)
{
  apply_boundary(var, name, name);
}

void apply_boundary(Vector3D &var, const char* name)
{
  static char buffer[128];
  int opt;
  // X component
  
  if(var.covariant) {
    sprintf(buffer, "%s_x", name);
  }else
    sprintf(buffer, "%sx", name);

  apply_boundary(var.x, buffer, name);

  // Y component
  
  if(var.covariant) {
    sprintf(buffer, "%s_y", name);
  }else
    sprintf(buffer, "%sy", name);

  apply_boundary(var.y, buffer, name);

  // Z component
  
  if(var.covariant) {
    sprintf(buffer, "%s_z", name);
  }else
    sprintf(buffer, "%sz", name);

  apply_boundary(var.z, buffer, name);

  // Check for vector boundary

  if(options.getInt(name, "xinner", opt))
    if(options.getInt("All", "xinner", opt))
      opt = BNDRY_NONE; // Do nothing

  if(opt == BNDRY_DIVCURL) {
    bndry_inner_divcurl(var);
  }

  if(options.getInt(name, "xouter", opt))
    if(options.getInt("All", "xouter", opt))
      opt = BNDRY_NONE; // Do nothing

  if(opt == BNDRY_DIVCURL) {
    bndry_sol_divcurl(var);
  }

}

/////////////////// X BOUNDARIES /////////////////////

void bndry_core_zero(Field2D &var)
{
  int jx, jy;

  if((MYPE_IN_CORE == 1) && mesh->first_x())
    for(jx=0;jx<mesh->xstart;jx++)
      for(jy=0;jy<mesh->ngy;jy++)
	var[jx][jy] = 0.0;
  
}

// Inner x core boundary
void bndry_core_zero(Field3D &var)
{
  // Set boundary to zero
  int jx, jy, jz;
  
  if((MYPE_IN_CORE == 1) && mesh->first_x()) {
    for(jx=0;jx<mesh->xstart;jx++) {
      for(jy=0;jy<mesh->ngy;jy++) {
	for(jz=0;jz<mesh->ngz;jz++) {
	  var[jx][jy][jz] = 0.0;
	}
      }
    }
  }
}

void bndry_core_zero(Vector3D &var)
{
  bndry_core_zero(var.x);
  bndry_core_zero(var.y);
  bndry_core_zero(var.z);
}

void bndry_core_flat(Field2D &var)
{
  int jx, jy;

  if((MYPE_IN_CORE == 1) && mesh->first_x())
    for(jx=0;jx<mesh->xstart;jx++)
      for(jy=0;jy<mesh->ngy;jy++)
	var[jx][jy] = var[mesh->xstart][jy];
  
}

void bndry_core_flat(Field3D &var)
{
  int jx, jy, jz;
  
  if((MYPE_IN_CORE == 1) && mesh->first_x()) {
    if(mesh->ShiftXderivs)
      var = var.ShiftZ(true);

    for(jx=0;jx<mesh->xstart;jx++) {
      for(jy=0;jy<mesh->ngy;jy++) {
	for(jz=0;jz<mesh->ngz;jz++) {
	  var[jx][jy][jz] = var[mesh->xstart][jy][jz];
	}
      }
    }
    
    if(mesh->ShiftXderivs)
      var = var.ShiftZ(false);
  }
}

void bndry_core_flat(Vector3D &var)
{
  bndry_core_flat(var.x);
  bndry_core_flat(var.y);
  bndry_core_flat(var.z);
}


void bndry_core_laplace(Field3D &var)
{
  static dcomplex *c0 = (dcomplex*) NULL, *c1, *c2;
  int jy, jz;
  real kwave;
  real coef1, coef2, coef3;
  dcomplex a, b, c;

  if((MYPE_IN_CORE != 1) || mesh->first_x())
    return;

  int ncz = mesh->ngz-1;
  
  if(c0 == (dcomplex*) NULL) {
    // Allocate memory
    c0 = new dcomplex[ncz/2 + 1];
    c1 = new dcomplex[ncz/2 + 1];
    c2 = new dcomplex[ncz/2 + 1];
  }

  for(jy=0;jy<mesh->ngy;jy++) {

    // Take FFT
    ZFFT(var[2][jy], mesh->zShift[2][jy], c2);

    coef1=mesh->g11[1][jy]/(SQ(mesh->dx[1][jy]));
    coef2=mesh->g33[1][jy];
    coef3=mesh->g13[1][jy]/(2.0*mesh->dx[1][jy]);

    for(jz=0;jz<=ncz/2;jz++) {
      kwave=jz*2.0*PI/mesh->zlength; // wave number is 1/[rad]
      
      a = dcomplex(coef1,-kwave*coef3);
      b = dcomplex(-2.0*coef1 - SQ(kwave)*coef2,0.0);
      c = dcomplex(coef1,kwave*coef3);

      // a*c0 + b*c1 + c*c2 = 0, c0 = c1
      c0[jz] = c1[jz] = -c * c2[jz] / (a+b);
    }

    // Reverse FFT
    ZFFT_rev(c0, mesh->zShift[0][jy], var[0][jy]);
    ZFFT_rev(c1, mesh->zShift[1][jy], var[1][jy]);
    
    var[0][jy][ncz] = var[0][jy][0];
    var[1][jy][ncz] = var[0][jy][0];
  }
}

void bndry_core_laplace2(Field3D &var)
{
  static dcomplex *c0 = (dcomplex*) NULL, *c1, *c2;
  int jx, jy, jz;
  real kwave;
  real coef1, coef2, coef3;
  dcomplex a, b, c;

  if((MYPE_IN_CORE != 1) || mesh->first_x())
    return;

  int ncz = mesh->ngz-1;

  if(c0 == (dcomplex*) NULL) {
    // Allocate memory
    c0 = new dcomplex[ncz/2 + 1];
    c1 = new dcomplex[ncz/2 + 1];
    c2 = new dcomplex[ncz/2 + 1];
  }

  for(jy=0;jy<mesh->ngy;jy++) {
    for(jx=mesh->xstart-1;jx>=0;jx--) {

      // Take FFT
      ZFFT(var[jx+1][jy], mesh->zShift[jx+1][jy], c1);
      ZFFT(var[jx+2][jy], mesh->zShift[jx+2][jy], c2);

      coef1=mesh->g11[jx+1][jy]/(SQ(mesh->dx[jx+1][jy]));
      coef2=mesh->g33[jx+1][jy];
      coef3=mesh->g13[jx+1][jy]/(2.0*mesh->dx[jx+1][jy]);

      for(jz=0;jz<=ncz/2;jz++) {
	kwave=jz*2.0*PI/mesh->zlength; // wave number is 1/[rad]
	
	a = dcomplex(coef1,-kwave*coef3);
	b = dcomplex(-2.0*coef1 - SQ(kwave)*coef2,0.0);
	c = dcomplex(coef1,kwave*coef3);
	
	// a*c0 + b*c1 + c*c2 = 0
	c0[jz] = -(b*c1[jz] + c*c2[jz])/a;
      }

      // Reverse FFT
      ZFFT_rev(c0, mesh->zShift[jx][jy], var[jx][jy]);
      
      var[jx][jy][ncz] = var[jx][jy][0];
    }
  }
}

// Inner x PF boundary
void bndry_pf_zero(Field2D &var)
{
  // Set boundary to zero
  int jx, jy;
  
  if((MYPE_IN_CORE == 0) && mesh->first_x())
    for(jx=0;jx<mesh->xstart;jx++)
      for(jy=0;jy<mesh->ngy;jy++)
	var[jx][jy] = 0.0;
  
}

void bndry_pf_zero(Field3D &var)
{
  // Set boundary to zero
  int jx, jy, jz;
  
  if((MYPE_IN_CORE == 0) && mesh->first_x()) {
    for(jx=0;jx<mesh->xstart;jx++) {
      for(jy=0;jy<mesh->ngy;jy++) {
	for(jz=0;jz<mesh->ngz;jz++) {
	  var[jx][jy][jz] = 0.0;
	}
      }
    }
  }
}

void bndry_pf_zero(Vector3D &var)
{
  bndry_pf_zero(var.x);
  bndry_pf_zero(var.y);
  bndry_pf_zero(var.z);
}

void bndry_pf_flat(Field2D &var)
{
  // Set boundary to zero
  int jx, jy;
  
  if((MYPE_IN_CORE == 0) && mesh->first_x())
    for(jx=0;jx<mesh->xstart;jx++)
      for(jy=0;jy<mesh->ngy;jy++)
	var[jx][jy] = var[mesh->xstart][jy];
  
}

void bndry_pf_flat(Field3D &var)
{
  int jx, jy, jz;

  if((MYPE_IN_CORE == 0) && mesh->first_x()) {
    if(mesh->ShiftXderivs)
      var = var.ShiftZ(true);

    for(jx=0;jx<mesh->xstart;jx++) {
      for(jy=0;jy<mesh->ngy;jy++) {
	for(jz=0;jz<mesh->ngz;jz++) {
	  var[jx][jy][jz] = var[mesh->xstart][jy][jz];
	}
      }
    }
    
    if(mesh->ShiftXderivs)
      var = var.ShiftZ(false);
  }
}

void bndry_pf_flat(Vector3D &var)
{
  bndry_pf_flat(var.x);
  bndry_pf_flat(var.y);
  bndry_pf_flat(var.z);
}

void bndry_pf_laplace(Field3D &var)
{
  static dcomplex *c0 = (dcomplex*) NULL, *c1, *c2;
  int jy, jz;
  real kwave;
  real coef1, coef2, coef3;
  dcomplex a, b, c;

  int ncz = mesh->ngz-1;

  if((MYPE_IN_CORE != 0) || (!mesh->first_x()))
    return;

  if(c0 == (dcomplex*) NULL) {
    // Allocate memory
    c0 = new dcomplex[ncz/2 + 1];
    c1 = new dcomplex[ncz/2 + 1];
    c2 = new dcomplex[ncz/2 + 1];
  }

  for(jy=0;jy<mesh->ngy;jy++) {

    // Take FFT
    ZFFT(var[2][jy], mesh->zShift[2][jy], c2);

    coef1=mesh->g11[1][jy]/(SQ(mesh->dx[1][jy]));
    coef2=mesh->g33[1][jy];
    coef3=mesh->g13[1][jy]/(2.0*mesh->dx[1][jy]);

    for(jz=0;jz<=ncz/2;jz++) {
      kwave=jz*2.0*PI/mesh->zlength; // wave number is 1/[rad]
      
      a = dcomplex(coef1,-kwave*coef3);
      b = dcomplex(-2.0*coef1 - SQ(kwave)*coef2,0.0);
      c = dcomplex(coef1,kwave*coef3);

      // a*c0 + b*c1 + c*c2 = 0, c0 = c1
      c0[jz] = c1[jz] = -c * c2[jz] / (a+b);
    }

    // Reverse FFT
    ZFFT_rev(c0, mesh->zShift[0][jy], var[0][jy]);
    ZFFT_rev(c1, mesh->zShift[1][jy], var[1][jy]);
    
    var[0][jy][ncz] = var[0][jy][0];
    var[1][jy][ncz] = var[0][jy][0];
  }
}

// Inner x (core and PF) boundary

void bndry_inner_zero(Vector3D &var)
{
  bndry_inner_zero(var.x);
  bndry_inner_zero(var.y);
  bndry_inner_zero(var.z);
}

void bndry_inner_flat(Vector3D &var)
{
  bndry_inner_flat(var.x);
  bndry_inner_flat(var.y);
  bndry_inner_flat(var.z);
}

void bndry_inner_laplace(Field3D &var)
{
  bndry_core_laplace(var);
  bndry_pf_laplace(var);
}


// Outer x boundary
void bndry_sol_zero(Field2D &var)
{
  int jx, jy;
  
  if(!mesh->last_x())
    return;

  for(jx=mesh->ngx;jx>mesh->xend;jx--)
    for(jy=0;jy<mesh->ngy;jy++)
      var[jx][jy] = 0.0;
}

void bndry_sol_zero(Field3D &var)
{
  int jx, jy, jz;
  
  if(!mesh->last_x())
    return;

  for(jx=mesh->ngx-1;jx>mesh->xend;jx--) {
    for(jy=0;jy<mesh->ngy;jy++) {
      for(jz=0;jz<mesh->ngz;jz++) {
	var[jx][jy][jz] = 0.0;
      }
    }
  }
}

void bndry_sol_zero(Vector3D &var)
{
  bndry_sol_zero(var.x);
  bndry_sol_zero(var.y);
  bndry_sol_zero(var.z);
}

void bndry_sol_flat(Field2D &var)
{
  int jx, jy;
  
  if(!mesh->last_x())
    return;

  for(jx=mesh->ngx-1;jx>mesh->xend;jx--)
    for(jy=0;jy<mesh->ngy;jy++)
      var[jx][jy] = var[mesh->xend][jy];
}

void bndry_sol_flat(Field3D &var)
{
  int jx, jy, jz;

  if(!mesh->last_x())
    return;
  
  if(mesh->ShiftXderivs)
    var = var.ShiftZ(true);

  for(jx=mesh->ngx-1;jx>mesh->xend;jx--) {
    for(jy=0;jy<mesh->ngy;jy++) {
      for(jz=0;jz<mesh->ngz;jz++) {
	var[jx][jy][jz] = var[mesh->xend][jy][jz];
      }
    }
  }

  if(mesh->ShiftXderivs)
    var = var.ShiftZ(false);
}

void bndry_sol_flat(Vector3D &var)
{
  bndry_sol_flat(var.x);
  bndry_sol_flat(var.y);
  bndry_sol_flat(var.z);
}

void bndry_sol_laplace(Field3D &var)
{
  static dcomplex *c0 = (dcomplex*) NULL, *c1, *c2;
  int jx, jy, jz;
  real kwave;
  real coef1, coef2, coef3;
  dcomplex a, b, c;

  if(!mesh->last_x())
    return;

  int ncz = mesh->ngz-1;

  if(c0 == (dcomplex*) NULL) {
    // Allocate memory
    c0 = new dcomplex[ncz/2 + 1];
    c1 = new dcomplex[ncz/2 + 1];
    c2 = new dcomplex[ncz/2 + 1];
  }

  for(jy=0;jy<mesh->ngy;jy++) {
    for(jx=mesh->xend+1;jx<mesh->ngx;jx++) {

      // Take FFT
      ZFFT(var[jx-2][jy], mesh->zShift[jx-2][jy], c0);
      ZFFT(var[jx-1][jy], mesh->zShift[jx-1][jy], c1);

      coef1=mesh->g11[jx-1][jy]/(SQ(mesh->dx[jx-1][jy]));
      coef2=mesh->g33[jx-1][jy];
      coef3=mesh->g13[jx-1][jy]/(2.0*mesh->dx[jx-1][jy]);

      for(jz=0;jz<=ncz/2;jz++) {
	kwave=jz*2.0*PI/mesh->zlength; // wave number is 1/[rad]
	
	a = dcomplex(coef1,-kwave*coef3);
	b = dcomplex(-2.0*coef1 - SQ(kwave)*coef2,0.0);
	c = dcomplex(coef1,kwave*coef3);
	
	// a*c0 + b*c1 + c*c2 = 0
	c2[jz] = -(a*c0[jz] + b*c1[jz])/c;
      }

      // Reverse FFT
      ZFFT_rev(c2, mesh->zShift[jx][jy], var[jx][jy]);
      
      var[jx][jy][ncz] = var[jx][jy][0];
    }
  }
}

//////////// BOUNDARY CONDITION FOR B //////////////
// Boundary condition for B fields: Curl B = Div B = 0

void bndry_inner_divcurl(Vector3D &var)
{
  
}

void bndry_sol_divcurl(Vector3D &var)
{
  int jx, jy, jz, jzp, jzm;
  real tmp;
  
  int ncz = mesh->ngz-1;

  if(!mesh->last_x())
    return;

  var.to_covariant();

  if(mesh->xstart > 2) {
    output.write("Error: Div = Curl = 0 boundary condition doesn't work for MXG > 2. Sorry\n");
    exit(1);
  }

  jx = mesh->xend+1;
  for(jy=1;jy<mesh->ngy-1;jy++) {
    for(jz=0;jz<ncz;jz++) {
      jzp = (jz+1) % ncz;
      jzm = (jz - 1 + ncz) % ncz;

      // dB_y / dx = dB_x / dy
      
      // dB_x / dy
      tmp = (var.x[jx-1][jy+1][jz] - var.x[jx-1][jy-1][jz]) / (mesh->dy[jx-1][jy-1] + mesh->dy[jx-1][jy]);
      
      var.y[jx][jy][jz] = var.y[jx-2][jy][jz] + (mesh->dx[jx-2][jy] + mesh->dx[jx-1][jy]) * tmp;
      if(mesh->xstart == 2)
	// 4th order to get last point
	var.y[jx+1][jy][jz] = var.y[jx-3][jy][jz] + 4.*mesh->dx[jx][jy]*tmp;

      // dB_z / dx = dB_x / dz

      tmp = (var.x[jx-1][jy][jzp] - var.x[jx-1][jy][jzm]) / (2.*mesh->dz);
      
      var.z[jx][jy][jz] = var.z[jx-2][jy][jz] + (mesh->dx[jx-2][jy] + mesh->dx[jx-1][jy]) * tmp;
      if(mesh->xstart == 2)
	var.z[jx+1][jy][jz] = var.z[jx-3][jy][jz] + 4.*mesh->dx[jx][jy]*tmp;

      // d/dx( Jmesh->g11 B_x ) = - d/dx( Jmesh->g12 B_y + Jmesh->g13 B_z) 
      //                    - d/dy( JB^y ) - d/dz( JB^z )
	
      tmp = -( mesh->J[jx][jy]*mesh->g12[jx][jy]*var.y[jx][jy][jz] + mesh->J[jx][jy]*mesh->g13[jx][jy]*var.z[jx][jy][jz]
	       - mesh->J[jx-2][jy]*mesh->g12[jx-2][jy]*var.y[jx-2][jy][jz] + mesh->J[jx-2][jy]*mesh->g13[jx-2][jy]*var.z[jx-2][jy][jz] )
	/ (mesh->dx[jx-2][jy] + mesh->dx[jx-1][jy]); // First term (d/dx) using vals calculated above
      tmp -= (mesh->J[jx-1][jy+1]*mesh->g12[jx-1][jy+1]*var.x[jx-1][jy+1][jz] - mesh->J[jx-1][jy-1]*mesh->g12[jx-1][jy-1]*var.x[jx-1][jy-1][jz]
	      + mesh->J[jx-1][jy+1]*mesh->g22[jx-1][jy+1]*var.y[jx-1][jy+1][jz] - mesh->J[jx-1][jy-1]*mesh->g22[jx-1][jy-1]*var.y[jx-1][jy-1][jz]
	      + mesh->J[jx-1][jy+1]*mesh->g23[jx-1][jy+1]*var.z[jx-1][jy+1][jz] - mesh->J[jx-1][jy-1]*mesh->g23[jx-1][jy-1]*var.z[jx-1][jy-1][jz])
	/ (mesh->dy[jx-1][jy-1] + mesh->dy[jx-1][jy]); // second (d/dy)
      tmp -= (mesh->J[jx-1][jy]*mesh->g13[jx-1][jy]*(var.x[jx-1][jy][jzp] - var.x[jx-1][jy][jzm]) +
	      mesh->J[jx-1][jy]*mesh->g23[jx-1][jy]*(var.y[jx-1][jy][jzp] - var.y[jx-1][jy][jzm]) +
	      mesh->J[jx-1][jy]*mesh->g33[jx-1][jy]*(var.z[jx-1][jy][jzp] - var.z[jx-1][jy][jzm])) / (2.*mesh->dz);
      
      var.x[jx][jy][jz] = ( mesh->J[jx-2][jy]*mesh->g11[jx-2][jy]*var.x[jx-2][jy][jz] + 
			    (mesh->dx[jx-2][jy] + mesh->dx[jx-1][jy]) * tmp ) / mesh->J[jx][jy]*mesh->g11[jx][jy];
      if(mesh->xstart == 2)
	var.x[jx+1][jy][jz] = ( mesh->J[jx-3][jy]*mesh->g11[jx-3][jy]*var.x[jx-3][jy][jz] + 
				4.*mesh->dx[jx][jy]*tmp ) / mesh->J[jx+1][jy]*mesh->g11[jx+1][jy];
    }
  }
}

// Relaxing X boundary conditions

void bndry_core_relax_val(Field3D &F_var, const Field3D &var, real value, real rate)
{
  if(MYPE_IN_CORE)
    bndry_inner_relax_val(F_var, var, value, rate);
}

void bndry_pf_relax_val(Field3D &F_var, const Field3D &var, real value, real rate)
{
  if(!MYPE_IN_CORE)
    bndry_inner_relax_val(F_var, var, value, rate);
}

void bndry_inner_relax_val(Field3D &F_var, const Field3D &var, real value, real rate)
{
  int jx, jy, jz;

  if(!mesh->first_x())
    return;
  
  rate = fabs(rate);

  for(jx=0;jx<mesh->xstart;jx++)
    for(jy=0;jy<mesh->ngy;jy++)
      for(jz=0;jz<mesh->ngz-1;jz++)
	F_var[jx][jy][jz] = rate*(value - var[jx][jy][jz]);
}

// Combination of relaxing zero-gradient and zero-value
void bndry_inner_relax_val2(Field3D &F_var, const Field3D &var, real value, real rate)
{
  int jx, jy, jz;

  if(!mesh->first_x())
    return;
  
  rate = fabs(rate);

  for(jx=mesh->xstart-1;jx>=0;jx--)
    for(jy=0;jy<mesh->ngy;jy++)
      for(jz=0;jz<mesh->ngz-1;jz++) {
	F_var[jx][jy][jz] = F_var[jx+1][jy][jz] + rate*(value - var[jx][jy][jz]);
      }
}

void bndry_sol_relax_val(Field3D &F_var, const Field3D &var, real value, real rate)
{
  int jx, jy, jz;
  
  if(!mesh->last_x())
    return;

  rate = fabs(rate);

  for(jx=mesh->xend+1;jx<mesh->ngx;jx++)
    for(jy=0;jy<mesh->ngy;jy++)
      for(jz=0;jz<mesh->ngz-1;jz++)
	F_var[jx][jy][jz] = rate*(value - var[jx][jy][jz]);
}

void bndry_sol_relax_val2(Field3D &F_var, const Field3D &var, real value, real rate)
{
  int jx, jy, jz;
  
  if(!mesh->last_x())
    return;

  rate = fabs(rate);

  for(jx=mesh->xend+1;jx<mesh->ngx;jx++)
    for(jy=0;jy<mesh->ngy;jy++)
      for(jz=0;jz<mesh->ngz-1;jz++) {
	F_var[jx][jy][jz] = F_var[jx-1][jy][jz] + rate*(value - var[jx][jy][jz]);
      }
}

/// RELAX TO ZERO-GRADIENT

void bndry_core_relax_flat(Field3D &F_var, const Field3D &var, real rate)
{
  if(MYPE_IN_CORE)
    bndry_inner_relax_flat(F_var, var, rate);
}

void bndry_pf_relax_flat(Field3D &F_var, const Field3D &var, real rate)
{
  if(!MYPE_IN_CORE)
    bndry_inner_relax_flat(F_var, var, rate);
}

void bndry_inner_relax_flat(Field3D &F_var, const Field3D &var, real rate)
{
  int jx, jy, jz;
  
  if(!mesh->first_x())
    return;
  
  rate = fabs(rate);

  for(jx=0;jx<mesh->xstart;jx++)
    for(jy=0;jy<mesh->ngy;jy++)
      for(jz=0;jz<mesh->ngz-1;jz++)
	F_var[jx][jy][jz] = rate*(var[jx+1][jy][jz] - var[jx][jy][jz]);
}

void bndry_sol_relax_flat(Field3D &F_var, const Field3D &var, real rate)
{
  int jx, jy, jz;
  
  if(!mesh->last_x())
    return;

  rate = fabs(rate);

  for(jx=mesh->xend+1;jx<mesh->ngx;jx++)
    for(jy=0;jy<mesh->ngy;jy++)
      for(jz=0;jz<mesh->ngz-1;jz++)
	F_var[jx][jy][jz] = rate*(var[jx-1][jy][jz] - var[jx][jy][jz]);
}

// Symmetric boundary

void bndry_core_sym(Field3D &var)
{
  if(MYPE_IN_CORE)
    bndry_inner_sym(var);
}

void bndry_pf_sym(Field3D &var)
{
  if(!MYPE_IN_CORE)
    bndry_inner_sym(var);
}

void bndry_inner_sym(Field3D &var)
{
  int jx, jy, jz, xb;

  if(!mesh->first_x())
    return;

  if(mesh->ShiftXderivs)
    var = var.ShiftZ(true);

  xb = 2*mesh->xstart;
  if(!mesh->BoundaryOnCell)
    xb--;  // Move boundary to between cells

  for(jx=0;jx<mesh->xstart;jx++)
    for(jy=0;jy<mesh->ngy;jy++)
      for(jz=0;jz<mesh->ngz;jz++)
	var[jx][jy][jz] = var[xb-jx][jy][jz];
  
  if(mesh->ShiftXderivs)
    var = var.ShiftZ(false);
}

void bndry_sol_sym(Field3D &var)
{
  int jx, jy, jz, xb;

  if(!mesh->last_x())
    return;
  
  if(mesh->ShiftXderivs)
    var = var.ShiftZ(true);

  xb = mesh->xend;
  if(mesh->BoundaryOnCell)
    xb--;

  for(jx=0;jx<mesh->xstart;jx++)
    for(jy=0;jy<mesh->ngy;jy++)
      for(jz=0;jz<mesh->ngz;jz++)
	var[mesh->xend+1+jx][jy][jz] = var[xb-jx][jy][jz];

  if(mesh->ShiftXderivs)
    var = var.ShiftZ(false);
}

// Relax to symmetric boundary

void bndry_core_relax_sym(Field3D &F_var, const Field3D &var, real rate)
{
  if(MYPE_IN_CORE)
    bndry_inner_relax_sym(F_var, var, rate);
}

void bndry_pf_relax_sym(Field3D &F_var, const Field3D &var, real rate)
{
  if(!MYPE_IN_CORE)
    bndry_inner_relax_sym(F_var, var, rate);
}

void bndry_inner_relax_sym(Field3D &F_var, const Field3D &var1, real rate)
{
  int jx, jy, jz, xb;

  if(!mesh->first_x())
    return;

  rate = fabs(rate);

  Field3D var;
  var = var1;
  if(mesh->ShiftXderivs)
    var = var.ShiftZ(true);

  xb = 2*mesh->xstart;
  if(!mesh->BoundaryOnCell)
    xb--;  // Move boundary to between cells

  for(jx=0;jx<mesh->xstart;jx++)
    for(jy=0;jy<mesh->ngy;jy++)
      for(jz=0;jz<mesh->ngz;jz++)
	F_var[jx][jy][jz] = rate*(var[xb-jx][jy][jz] - var[jx][jy][jz]);
  
  if(mesh->ShiftXderivs)
    F_var = F_var.ShiftZ(false);
}

void bndry_sol_relax_sym(Field3D &F_var, const Field3D &var1, real rate)
{
  int jx, jy, jz, xb;

  if(!mesh->last_x())
    return;

  rate = fabs(rate);
  
  Field3D var;
  var = var1; 
  if(mesh->ShiftXderivs)
    var = var.ShiftZ(true);

  xb = mesh->xend;
  if(mesh->BoundaryOnCell)
    xb--;

  for(jx=0;jx<mesh->xstart;jx++)
    for(jy=0;jy<mesh->ngy;jy++)
      for(jz=0;jz<mesh->ngz;jz++)
	F_var[mesh->xend+1+jx][jy][jz] = rate*(var[xb-jx][jy][jz] - var[mesh->xend+1+jx][jy][jz]);

  if(mesh->ShiftXderivs)
    var = var.ShiftZ(false);
}

/////////////////// Y BOUNDARIES ////////////////




void bndry_ydown_zero(Vector3D &var)
{
  bndry_ydown_zero(var.x);
  bndry_ydown_zero(var.y);
  bndry_ydown_zero(var.z);
}

void bndry_yup_zero(Vector3D &var)
{
  bndry_yup_zero(var.x);
  bndry_yup_zero(var.y);
  bndry_yup_zero(var.z);
}

// Rotates 180 degrees
void bndry_ydown_rotate(Field3D &var, bool reverse)
{
  int jx, jy, jy2, jz;
  real kwave;
  
  static dcomplex *cv = (dcomplex*) NULL;

  int ncz = mesh->ngz-1;

  if(cv == (dcomplex*) NULL)
    cv = new dcomplex[ncz/2 + 1];

  for(jx=0;jx<mesh->ngx;jx++) {
    for(jy=(mesh->ystart-1); jy >= 0;jy--) {
      jy2 = 2*mesh->ystart - jy - 1;
      
      // jy2 is rotated 180 degrees and put into jy
      
      rfft(var[jx][jy2], ncz, cv);

      for(jz=1;jz<=ncz/2;jz++) {
	// Rotate by 180 degrees
	kwave=jz*2.0*PI/mesh->zlength; // wave number is 1/[rad]
	cv[jz] *= dcomplex(cos(kwave*PI) , sin(kwave*PI));
      }
      
      irfft(cv, ncz, var[jx][jy]);
      
      if(reverse) {
	for(jz=0;jz<ncz;jz++)
	  var[jx][jy][jz] *= -1.0;
      }
      var[jx][jy][ncz] = var[jx][jy][0];
    }
  }
}

// For the centre of a circle
void bndry_ydown_zaverage(Field3D &var)
{
  int jx, jy, jz;
  real avg, w;

  // centre point is average of last point
  // others are just a linear interpolation

  for(jx=0;jx<mesh->ngx;jx++) {
    avg = 0.0;
    for(jz=0;jz<mesh->ngz-1;jz++)
      avg += var[jx][mesh->ystart][jz];
    avg /= (real) mesh->ngz-1;
    for(jz=0;jz<mesh->ngz-1;jz++) {
      var[jx][0][jz] = avg;
      for(jy=1;jy<mesh->ystart;jy++) {
        w = ((real) jy) / ((real) mesh->ystart);
        var[jx][jy][jz] = w*var[jx][mesh->ystart][jz] + (1.0 - w)*avg;
      }
    }
  }
  
}

// Relax towards a given value

void bndry_ydown_relax_val(Field3D &F_var, const Field3D &var, real value, real rate)
{
  rate = fabs(rate); // RATE IS ALWAYS POSITIVE

  RangeIter *xi = mesh->iterateBndryLowerY();
  /// Iterate over boundary
  for(xi->first(); !xi->isDone(); xi->next())
    for(int jy=mesh->ystart-1;jy>=0;jy--) 
      for(int jz=0;jz<mesh->ngz;jz++)
	F_var[xi->ind][jy][jz] = F_var[xi->ind][jy+1][jz] + rate*(value - var[xi->ind][jy][jz]);
  delete xi;
}

void bndry_yup_relax_val(Field3D &F_var, const Field3D &var, real value, real rate)
{
  rate = fabs(rate); // RATE IS ALWAYS POSITIVE

  // If not communicating with another processor
  RangeIter *xi = mesh->iterateBndryUpperY();
  /// Iterate over boundary
  for(xi->first(); !xi->isDone(); xi->next())
    for(int jy=mesh->xend+1;jy<mesh->ngy;jy++) 
      for(int jz=0;jz<mesh->ngz;jz++)
	F_var[xi->ind][jy][jz] = F_var[xi->ind][jy-1][jz] + rate*(value - var[xi->ind][jy][jz]);
  delete xi;
}

/// Relax to zero gradient

void bndry_ydown_relax_flat(Field3D &F_var, const Field3D &var, real rate)
{
  rate = fabs(rate); // RATE IS ALWAYS POSITIVE

  RangeIter *xi = mesh->iterateBndryLowerY();
  /// Iterate over boundary
  for(xi->first(); !xi->isDone(); xi->next())
    for(int jy=mesh->xstart-1;jy>=0;jy--) 
      for(int jz=0;jz<mesh->ngz;jz++)
	F_var[xi->ind][jy][jz] = F_var[xi->ind][jy+1][jz]
	  + rate*(var[xi->ind][jy+1][jz] - var[xi->ind][jy][jz]);
  delete xi;
}

void bndry_yup_relax_flat(Field3D &F_var, const Field3D &var, real rate)
{
  rate = fabs(rate); // RATE IS ALWAYS POSITIVE

  RangeIter *xi = mesh->iterateBndryUpperY();
  /// Iterate over boundary
  for(xi->first(); !xi->isDone(); xi->next())
    for(int jy=mesh->xend+1;jy<mesh->ngy;jy++) 
      for(int jz=0;jz<mesh->ngz;jz++)
	F_var[xi->ind][jy][jz] = F_var[xi->ind][jy-1][jz]
	  + rate*(var[xi->ind][jy-1][jz] - var[xi->ind][jy][jz]);
  delete xi;
}

// Symmetric boundary condition

void bndry_ydown_sym(Field3D &var)
{
  int yb = 2*mesh->ystart;
  if(!mesh->BoundaryOnCell)
    yb--;
  
  RangeIter *xi = mesh->iterateBndryLowerY();
  /// Iterate over boundary
  for(xi->first(); !xi->isDone(); xi->next())
    for(int jy=0;jy<mesh->ystart;jy++)
      for(int jz=0;jz<mesh->ngz;jz++)
	var[xi->ind][jy][jz] = var[xi->ind][yb-jy][jz];
}

void bndry_yup_sym(Field3D &var)
{
  int yb = mesh->yend;
  if(mesh->BoundaryOnCell)
    yb--;
  
  RangeIter *xi = mesh->iterateBndryUpperY();
  /// Iterate over boundary
  for(xi->first(); !xi->isDone(); xi->next())
    for(int jy=mesh->yend+1;jy<mesh->ngy;jy++)
      for(int jz=0;jz<mesh->ngz;jz++)
	var[xi->ind][jy][jz] = var[xi->ind][yb-jy+mesh->yend+1][jz];
  delete xi;
}

// Relax to symmetric boundaries

void bndry_ydown_relax_sym(Field3D &F_var, const Field3D &var, real rate)
{
  rate = fabs(rate);

  int yb = 2*mesh->ystart;
  if(!mesh->BoundaryOnCell)
    yb--;
  
  RangeIter *xi = mesh->iterateBndryLowerY();
  /// Iterate over boundary
  for(xi->first(); !xi->isDone(); xi->next())
    for(int jy=mesh->ystart-1;jy>=0;jy--)
      for(int jz=0;jz<mesh->ngz;jz++)
	F_var[xi->ind][jy][jz] = F_var[xi->ind][jy+1][jz]
	  + rate*(var[xi->ind][yb-jy][jz] - var[xi->ind][jy][jz]);
  delete xi;
}

void bndry_yup_relax_sym(Field3D &F_var, const Field3D &var, real rate)
{
  int yb = mesh->yend;
  if(mesh->BoundaryOnCell)
    yb--;

  rate = fabs(rate);
  
  RangeIter *xi = mesh->iterateBndryUpperY();
  /// Iterate over boundary
  for(xi->first(); !xi->isDone(); xi->next())
    for(int jy=0;jy<mesh->ystart;jy++)
      for(int jz=0;jz<mesh->ngz;jz++)
	F_var[xi->ind][mesh->yend+1+jy][jz] = F_var[xi->ind][mesh->yend+jy][jz]
	  + rate*(var[xi->ind][yb-jy][jz] - var[xi->ind][mesh->yend+1+jy][jz]);
  delete xi;
}

////////////////// Z BOUNDARIES ///////////////

void bndry_toroidal(Field3D &var)
{
  int jx, jy;

  for(jx=0; jx<mesh->ngx; jx++) {
    for(jy=0; jy<mesh->ngy; jy++) {
      var[jx][jy][mesh->ngz-1] = var[jx][jy][0];
    }
  }
}

void bndry_toroidal(Vector3D &var)
{
  bndry_toroidal(var.x);
  bndry_toroidal(var.y);
  bndry_toroidal(var.z);
}

/////////////////// ZERO LAPLACE ///////////////

/// Sets innermost point to zero, and sets laplacian = 0 on last "real" point
void bndry_inner_zero_laplace(Field3D &var)
{
  static dcomplex *c1 = (dcomplex*) NULL, *c2, *c3;
  int jy, jz;
  dcomplex d, e, f;

  if(!mesh->first_x())
    return;

  int ncz = mesh->ngz-1;

  if(c1 == (dcomplex*) NULL) {
    // Allocate memory
    c1 = new dcomplex[ncz/2 + 1];
    c2 = new dcomplex[ncz/2 + 1];
    c3 = new dcomplex[ncz/2 + 1];
  }

  for(jy=0;jy<mesh->ngy;jy++) {

    // Take FFT
    ZFFT(var[2][jy], mesh->zShift[2][jy], c2);
    ZFFT(var[3][jy], mesh->zShift[3][jy], c3);

    for(jz=0;jz<=ncz/2;jz++) {
      
      laplace_tridag_coefs(2, jy, jz, d, e, f);

      // d*c1 + e*c2 + f*c3 = 0
      
      c1[jz] = - (e*c2[jz] + f*c3[jz]) / d;
    }

    // Reverse FFT 
    ZFFT_rev(c1, mesh->zShift[1][jy], var[1][jy]);

    for(jz=0;jz<ncz;jz++)
      var[0][jy][jz] = 0.; // Innermost point set to zero
    
    var[0][jy][ncz] = var[0][jy][0];
    var[1][jy][ncz] = var[0][jy][0];
  }
}

void bndry_core_zero_laplace(Field3D &var)
{
  if(MYPE_IN_CORE == 1)
    bndry_inner_zero_laplace(var);
}

void bndry_pf_zero_laplace(Field3D &var)
{
  if(MYPE_IN_CORE == 0)
    bndry_inner_zero_laplace(var);
}



/////////////

void bndry_outer_zero_laplace(Field3D &var)
{
  static dcomplex *c1 = (dcomplex*) NULL, *c2, *c3;
  int jy, jz;
  dcomplex d, e, f;

  if(!mesh->last_x())
    return;
  
  int ncz = mesh->ngz-1;
  
  if(c1 == (dcomplex*) NULL) {
    // Allocate memory
    c1 = new dcomplex[ncz/2 + 1];
    c2 = new dcomplex[ncz/2 + 1];
    c3 = new dcomplex[ncz/2 + 1];
  }

  for(jy=0;jy<mesh->ngy;jy++) {

    // Take FFT
    ZFFT(var[mesh->ngx-4][jy], mesh->zShift[mesh->ngx-4][jy], c1);
    ZFFT(var[mesh->ngx-3][jy], mesh->zShift[mesh->ngx-3][jy], c2);
    
    for(jz=0;jz<=ncz/2;jz++) {
      
      laplace_tridag_coefs(mesh->ngx-3, jy, jz, d, e, f);
      
      // d*c1 + e*c2 + f*c3 = 0
      
      c3[jz] = - (d*c1[jz] + e*c2[jz]) / f;
    }

    // Reverse FFT to get next-to-last point
    ZFFT_rev(c3, mesh->zShift[mesh->ngx-2][jy], var[mesh->ngx-2][jy]);
    
    for(jz=0;jz<ncz;jz++)
      var[mesh->ngx-1][jy][jz] = 0.; // Outermost point set to zero

    var[mesh->ngx-1][jy][ncz] = var[mesh->ngx-1][jy][0];
    var[mesh->ngx-2][jy][ncz] = var[mesh->ngx-2][jy][0];
  }
}

/////////////////////////////////////////////////////////////

void bndry_inner_laplace_decay(Field3D &var)
{
  static dcomplex *c0 = (dcomplex*) NULL, *c1;

  if(!mesh->first_x())
    return;

  int ncz = mesh->ngz-1;

  if(c0 == (dcomplex*) NULL) {
    // Allocate memory
    c0 = new dcomplex[ncz/2 + 1];
    c1 = new dcomplex[ncz/2 + 1];
  }
  
  for(int jy=0;jy<mesh->ngy;jy++) {
    // Take FFT of first "real" points
    ZFFT(var[mesh->xstart][jy], mesh->zShift[mesh->xstart][jy], c0);
    ZFFT(var[mesh->xstart+1][jy], mesh->zShift[mesh->xstart+1][jy], c1);
    c1[0] -= c0[0]; // Only need DC component
    
    // Solve  mesh->g11*d2f/dx2 - mesh->g33*kz^2f = 0
    // Assume mesh->g11, mesh->g33 constant -> exponential growth or decay

    for(int jx=mesh->xstart-1;jx>=0;jx--) {
      // kz = 0 solution
      c0[0] -= c1[0];  // Straight line
      
      // kz != 0 solution
      real coef = -1.0*sqrt(mesh->g33[jx][jy] / mesh->g11[jx][jy])*mesh->dx[jx][jy];
      for(int jz=1;jz<=ncz/2;jz++) {
	real kwave=jz*2.0*PI/mesh->zlength; // wavenumber in [rad^-1]
	c0[jz] *= exp(coef*kwave); // The decaying solution only
      }
      // Reverse FFT
      ZFFT_rev(c0, mesh->zShift[jx][jy], var[jx][jy]);
    }
  }
}

void bndry_outer_laplace_decay(Field3D &var)
{
  static dcomplex *c0 = (dcomplex*) NULL, *c1;

  if(!mesh->last_x())
    return;

  int ncz = mesh->ngz-1;

  if(c0 == (dcomplex*) NULL) {
    // Allocate memory
    c0 = new dcomplex[ncz/2 + 1];
    c1 = new dcomplex[ncz/2 + 1];
  }
  
  for(int jy=0;jy<mesh->ngy;jy++) {
    // Take FFT of final "real" points
    ZFFT(var[mesh->xend][jy], mesh->zShift[mesh->xend][jy], c0);
    ZFFT(var[mesh->xend-1][jy], mesh->zShift[mesh->xend-1][jy], c1);
    c1[0] = c0[0] - c1[0]; // Only used for DC component
    
    // Solve  mesh->g11*d2f/dx2 - mesh->g33*kz^2f = 0
    // Assume mesh->g11, mesh->g33 constant -> exponential growth or decay

    for(int jx=mesh->xend+1;jx<mesh->ngx;jx++) {
      // kz = 0 solution
      c0[0] += c1[0]; // Just a straight line
      
      // kz != 0 solution
      real coef = -1.0*sqrt(mesh->g33[jx-1][jy] / mesh->g11[jx-1][jy])*mesh->dx[jx-1][jy];
      for(int jz=1;jz<=ncz/2;jz++) {
	real kwave=jz*2.0*PI/mesh->zlength; // wavenumber in [rad^-1]
	c0[jz] *= exp(coef*kwave); // The decaying solution only
      }
      // Reverse FFT
      ZFFT_rev(c0, mesh->zShift[jx][jy], var[jx][jy]);
    }
  }
}

/////////////////////////////////////////////////////////////

void bndry_inner_const_laplace_decay(Field3D &var)
{
  static dcomplex *c0 = (dcomplex*) NULL, *c1, *c2;

  if(!mesh->first_x())
    return;

  int ncz = mesh->ngz-1;

  if(c0 == (dcomplex*) NULL) {
    // Allocate memory
    c0 = new dcomplex[ncz/2 + 1];
    c1 = new dcomplex[ncz/2 + 1];
    c2 = new dcomplex[ncz/2 + 1];
  }
  
  for(int jy=0;jy<mesh->ngy;jy++) {
    // Take FFT of first three "real" points
    ZFFT(var[mesh->xstart][jy],   mesh->zShift[mesh->xstart][jy],   c0);
    ZFFT(var[mesh->xstart+1][jy], mesh->zShift[mesh->xstart+1][jy], c1);
    ZFFT(var[mesh->xstart+2][jy], mesh->zShift[mesh->xstart+2][jy], c2);
    
    dcomplex k0lin = (c1[0] - c0[0])/mesh->dx[mesh->xstart][jy]; // for kz=0 solution
    
    // Calculate Delp2 on point MXG+1 (and put into c1)
    for(int jz=0;jz<=ncz/2;jz++) {
      dcomplex d,e,f;
      laplace_tridag_coefs(mesh->xstart+1, jy, jz, d, e, f);
      c1[jz] = d*c0[jz] + e*c1[jz] + f*c2[jz];
    }
    // Solve  mesh->g11*d2f/dx2 - mesh->g33*kz^2f = 0
    // Assume mesh->g11, mesh->g33 constant -> exponential growth or decay

    real xpos = 0.0;
    for(int jx=mesh->xstart-1;jx>=0;jx--) {
      // kz = 0 solution
      xpos -= mesh->dx[jx][jy];
      c2[0] = c0[0] + k0lin*xpos + 0.5*c1[0]*xpos*xpos/mesh->g11[jx+1][jy];
      // kz != 0 solution
      real coef = -1.0*sqrt(mesh->g33[jx+1][jy] / mesh->g11[jx+1][jy])*mesh->dx[jx+1][jy];
      for(int jz=1;jz<=ncz/2;jz++) {
	real kwave=jz*2.0*PI/mesh->zlength; // wavenumber in [rad^-1]
	c0[jz] *= exp(coef*kwave); // The decaying solution only
	// Add the particular solution
	c2[jz] = c0[jz] - c1[jz]/(mesh->g33[jx+1][jy]*kwave*kwave); 
      }
      // Reverse FFT
      ZFFT_rev(c2, mesh->zShift[jx][jy], var[jx][jy]);
    }
  }
}

void bndry_outer_const_laplace_decay(Field3D &var)
{
  static dcomplex *c0 = (dcomplex*) NULL, *c1, *c2;

  if(!mesh->first_x())
    return;

  int ncz = mesh->ngz-1;

  if(c0 == (dcomplex*) NULL) {
    // Allocate memory
    c0 = new dcomplex[ncz/2 + 1];
    c1 = new dcomplex[ncz/2 + 1];
    c2 = new dcomplex[ncz/2 + 1];
  }
  
  for(int jy=0;jy<mesh->ngy;jy++) {
    // Take FFT of last three "real" points
    
    ZFFT(var[mesh->xend][jy], mesh->zShift[mesh->xend][jy], c0);
    ZFFT(var[mesh->xend-1][jy], mesh->zShift[mesh->xend-1][jy], c1);
    ZFFT(var[mesh->xend-2][jy], mesh->zShift[mesh->xend-2][jy], c1);
    
    dcomplex k0lin = (c0[0] - c1[0])/mesh->dx[mesh->xend][jy]; // for kz=0 solution
    
    // Calculate Delp2 on point MXG+1 (and put into c1)
    for(int jz=0;jz<=ncz/2;jz++) {
      dcomplex d,e,f;
      laplace_tridag_coefs(mesh->xend-1, jy, jz, d, e, f);
      c1[jz] = d*c2[jz] + e*c1[jz] + f*c0[jz];
    }
    // Solve  mesh->g11*d2f/dx2 - mesh->g33*kz^2f = 0
    // Assume mesh->g11, mesh->g33 constant -> exponential growth or decay

    real xpos = 0.0;
    for(int jx=mesh->xend+1;jx<mesh->ngx;jx++) {
      // kz = 0 solution: a + bx + 0.5*jpar*x^2/mesh->g11
      xpos += mesh->dx[jx][jy];
      c2[0] = c0[0] + k0lin*xpos + 0.5*c1[0]*xpos*xpos/mesh->g11[jx-1][jy]; 
      // kz != 0 solution
      real coef = -1.0*sqrt(mesh->g33[jx-1][jy] / mesh->g11[jx-1][jy])*mesh->dx[jx-1][jy];
      for(int jz=1;jz<=ncz/2;jz++) {
	real kwave=jz*2.0*PI/mesh->zlength; // wavenumber in [rad^-1]
	c0[jz] *= exp(coef*kwave); // The decaying solution only
	// Add the particular solution
	c2[jz] = c0[jz] - c1[jz]/(mesh->g33[jx][jy]*kwave*kwave); 
      }
      // Reverse FFT
      ZFFT_rev(c2, mesh->zShift[jx][jy], var[jx][jy]);
    }
  }
}
