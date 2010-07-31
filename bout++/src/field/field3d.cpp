/**************************************************************************
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

#include "globals.h"

#include "mpi.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "field3d.h"
#include "utils.h"
#include "fft.h"
#include "dcomplex.h"
#include "interpolation.h"
#include "boundary_op.h"
#include "boutexception.h"

/// Constructor
Field3D::Field3D() : background(NULL)
{
#ifdef MEMDEBUG
  output.write("Field3D %u: constructor\n", (unsigned int) this);
#endif
#ifdef TRACK
  name = "<F3D>";
#endif

  /// Mark data as unallocated
  block = NULL;

  location = CELL_CENTRE; // Cell centred variable by default
  
  ddt = NULL;
}

/// Doesn't copy any data, just create a new reference to the same data (copy on change later)
Field3D::Field3D(const Field3D& f) : background(NULL)
{
#ifdef MEMDEBUG
  output.write("Field3D %u: Copy constructor from %u\n", (unsigned int) this, (unsigned int) &f);
#endif

#ifdef CHECK
  msg_stack.push("Field3D: Copy constructor");
  f.checkData();
#endif

  /// Copy a reference to the block
  block = f.block;
  /// Increase reference count
  block->refs++;

  location = f.location;

  ddt = NULL;
  
#ifdef CHECK
  msg_stack.pop();
#endif
}

Field3D::Field3D(const Field2D& f) : background(NULL)
{
#ifdef CHECK
  msg_stack.push("Field3D: Copy constructor from Field2D");
#endif
  
  /// Mark data as unallocated
  block = NULL;

  location = CELL_CENTRE; // Cell centred variable by default
  
  ddt = NULL;
  
  *this = f;
  
#ifdef CHECK
  msg_stack.pop();
#endif
}

Field3D::Field3D(const BoutReal val) : background(NULL)
{
#ifdef CHECK
  msg_stack.push("Field3D: Copy constructor from value");
#endif
  
  /// Mark data as unallocated
  block = NULL;

  location = CELL_CENTRE; // Cell centred variable by default
  
  ddt = NULL;
  
  *this = val;
  
#ifdef CHECK
  msg_stack.pop();
#endif
}

Field3D::~Field3D()
{
  /// free the block of data if allocated
  freeData();
  
  /// Delete the time derivative variable if allocated
  if(ddt != NULL)
    delete ddt;
}

Field3D* Field3D::clone() const
{
  return new Field3D(*this);
}

/// Make sure data is allocated and unique to this object
/// NOTE: Logically constant, but internals mutable
void Field3D::allocate() const
{
  allocData();
}

BoutReal*** Field3D::getData() const
{
#ifdef CHECK
  if(block ==  NULL)
    error("Field3D: getData() returning null pointer\n");

  // NOTE: Can't check data here since may not be set yet.
  //       Using this function circumvents the checking
#endif
  
  // User might alter data, so need to make unique
  allocate();

  return(block->data);
}

Field3D* Field3D::timeDeriv()
{
  if(ddt == NULL)
    ddt = new Field3D();
  
  return ddt;
}

const Field2D Field3D::DC()
{
  Field2D result;
  int jx, jy, jz;
  BoutReal **d;

#ifdef CHECK
  msg_stack.push("Field3D: DC");
  checkData();
#endif

#ifdef TRACK
  result.name = "DC(" + name + ")";
#endif

  result = 0.0;
  d = result.getData();

  BoutReal inv_n = 1. / (BoutReal) (mesh->ngz-1);

  for(jx=0;jx<mesh->ngx;jx++)
    for(jy=0;jy<mesh->ngy;jy++) {
      for(jz=0;jz<(mesh->ngz-1);jz++)
	d[jx][jy] += block->data[jx][jy][jz];
      d[jx][jy] *= inv_n;
    }
  
#ifdef CHECK
  msg_stack.pop();
#endif

  return(result);
}

void Field3D::setLocation(CELL_LOC loc)
{
  if(loc == CELL_VSHIFT) {
    error("Field3D: CELL_VSHIFT cell location only makes sense for vectors");
  }
  
  if(loc == CELL_DEFAULT)
    loc = CELL_CENTRE;
  
  location = loc;
}

CELL_LOC Field3D::getLocation() const
{
  return location;
}

/***************************************************************
 *                         OPERATORS 
 ***************************************************************/

BoutReal** Field3D::operator[](int jx) const
{
#ifdef CHECK
  if(block == NULL) {
    error("Field3D: [] operator on empty data");
    exit(1);
  }
  
  if((jx < 0) || (jx >= mesh->ngx)) {
    error("Field3D: [%d] operator out of bounds", jx);
    exit(1);
  }
#endif
 
  // User might alter data, so need to make unique
 allocate();
 
  return(block->data[jx]);
}

BoutReal& Field3D::operator[](bindex &bx) const
{
#ifdef CHECK
  if(block == NULL) {
    error("Field3D: [bindex] operator on empty data");
    exit(1);
  }
  if((bx.jx < 0) || (bx.jx >= mesh->ngx)) {
    error("Field3D: [bindex.jx = %d] out of range", bx.jx);
  }
  if((bx.jy < 0) || (bx.jy >= mesh->ngy)) {
    error("Field3D: [bindex.jy = %d] out of range", bx.jy);
  }
  if((bx.jz < 0) || (bx.jz >= mesh->ngz)) {
    error("Field3D: [bindex.jz = %d] out of range", bx.jz);
  }
#endif

  return block->data[bx.jx][bx.jy][bx.jz];
}

/////////////////// ASSIGNMENT ////////////////////

Field3D & Field3D::operator=(const Field3D &rhs)
{
  /// Check for self-assignment
  if(this == &rhs)
    return(*this); // skip this assignment

#ifdef CHECK

  msg_stack.push("Field3D: Assignment from Field3D");
  rhs.checkData(true);

#endif

#ifdef TRACK
  name = rhs.name;
#endif

  if(block != NULL)
    freeData(); // get rid of the current data

  /// Copy reference, don't copy data
  block = rhs.block;
  block->refs++;
  
  location = rhs.location;
  
#ifdef CHECK
  msg_stack.pop();
#endif

  return(*this);
}

Field3D & Field3D::operator=(const Field2D &rhs)
{
  int jx, jy, jz;
  BoutReal **d;

#ifdef CHECK
  msg_stack.push("Field3D: Assignment from Field2D");
  rhs.checkData(true);
#endif

  d = rhs.getData();

#ifdef TRACK
  name = "F3D("+rhs.name+")";
#endif

  allocate();

  /// Copy data

  for(jx=0;jx<mesh->ngx;jx++)
    for(jy=0;jy<mesh->ngy;jy++)
      for(jz=0;jz<mesh->ngz;jz++)
	block->data[jx][jy][jz] = d[jx][jy];

  /// Only 3D fields have locations
  //location = CELL_CENTRE;

#ifdef CHECK
  msg_stack.pop();
#endif

  return(*this);
}

Field3D & Field3D::operator=(const FieldPerp &rhs)
{
  int jx, jy, jz; 
  BoutReal **d;
  
  jy = rhs.getIndex();
  
  d = rhs.getData();

#ifdef CHECK
  if(d == (BoutReal**) NULL) {
    // No data
    error("Field3D: No data in assignment from FieldPerp");
  }
  
  /// Test rhs values
  for(jx=mesh->xstart;jx<=mesh->xend;jx++)
    for(jz=0;jz<mesh->ngz-1;jz++)
      if(!finite(d[jx][jz])) {
	error("Field3D: Assignment from non-finite FieldPerp data at (%d,%d,%d)\n", jx,jy,jz);
      }
#endif

#ifdef TRACK
  name = "F3D("+rhs.name+")";
#endif

 allocate();

  /// Copy data

  for(jx=0;jx<mesh->ngx;jx++)
    for(jz=0;jz<mesh->ngz;jz++)
      block->data[jx][jy][jz] = d[jx][jz];

  return(*this);
}

const bvalue & Field3D::operator=(const bvalue &bv)
{
 allocate();

#ifdef CHECK
  if(!finite(bv.val)) {
    output.write("Field3D: assignment from non-finite value at (%d,%d,%d)\n", 
	   bv.jx, bv.jy,bv.jz);
    exit(1);
  }
#endif

#ifdef TRACK
  name = "<bv3D>";
#endif

  block->data[bv.jx][bv.jy][bv.jz] = bv.val;
  
  return(bv);
}

BoutReal Field3D::operator=(const BoutReal val)
{
  int jx, jy, jz;
  
 allocate();

#ifdef CHECK
  if(!finite(val)) {
    error("Field3D: Assignment from non-finite BoutReal\n");
  }
#endif

#ifdef TRACK
  name = "<r3D>";
#endif

  for(jx=0;jx<mesh->ngx;jx++)
    for(jy=0;jy<mesh->ngy;jy++)
      for(jz=0;jz<mesh->ngz;jz++)
	block->data[jx][jy][jz] = val;

  // Only 3D fields have locations
  //location = CELL_CENTRE;
  // DON'T RE-SET LOCATION

  return(val);
}

////////////////// ADDITION //////////////////////

Field3D & Field3D::operator+=(const Field3D &rhs)
{
#ifdef CHECK
  msg_stack.push("Field3D: += Field3D");
  
  rhs.checkData();
  checkData();
#endif

  if(mesh->StaggerGrids && (rhs.location != location)) {
    // Interpolate and call again
    return (*this) += interp_to(rhs, location);
  }

#ifdef TRACK
  name = "(" + name + "+" + rhs.name + ")";
#endif

  if(block->refs == 1) {
    // This is the only reference to this data
    for(int jx=0;jx<mesh->ngx;jx++)
      for(int jy=0;jy<mesh->ngy;jy++)
	for(int jz=0;jz<mesh->ngz;jz++)
	  block->data[jx][jy][jz] += rhs.block->data[jx][jy][jz];
  }else {
    // Need to put result in a new block

    memblock3d *nb = newBlock();

    for(int jx=0;jx<mesh->ngx;jx++)
      for(int jy=0;jy<mesh->ngy;jy++)
	for(int jz=0;jz<mesh->ngz;jz++)
	  nb->data[jx][jy][jz] = block->data[jx][jy][jz] + rhs.block->data[jx][jy][jz];

    block->refs--;
    block = nb;
  }

#ifdef CHECK
  msg_stack.pop();
#endif

  return(*this);
}

Field3D & Field3D::operator+=(const Field2D &rhs)
{
  int jx, jy, jz;
  BoutReal **d;

#ifdef CHECK
  msg_stack.push("Field3D: += ( Field2D )");

  checkData();
  rhs.checkData();
#endif

  d = rhs.getData();

#ifdef TRACK
  name = "(" + name + "+" + rhs.name + ")";
#endif

  if(block->refs == 1) {
    for(jx=0;jx<mesh->ngx;jx++)
      for(jy=0;jy<mesh->ngy;jy++)
	for(jz=0;jz<mesh->ngz;jz++)
	  block->data[jx][jy][jz] += d[jx][jy];
  }else {
    memblock3d *nb = newBlock();
    
    for(jx=0;jx<mesh->ngx;jx++)
      for(jy=0;jy<mesh->ngy;jy++)
	for(jz=0;jz<mesh->ngz;jz++)
	  nb->data[jx][jy][jz] = block->data[jx][jy][jz] + d[jx][jy];

    block->refs--;
    block = nb;
  }

#ifdef CHECK
  msg_stack.pop();
#endif

  return(*this);
}

Field3D & Field3D::operator+=(const BoutReal &rhs)
{
  int jx, jy, jz;
#ifdef CHECK
  msg_stack.push("Field3D: += ( BoutReal )");

  checkData();

  if(!finite(rhs)) {
    error("Field3D: += operator passed non-finite BoutReal number");
  }
#endif
  
#ifdef TRACK
  name = "(" + name + "+BoutReal)";
#endif

  if(block->refs == 1) {
    for(jx=0;jx<mesh->ngx;jx++)
      for(jy=0;jy<mesh->ngy;jy++)
	for(jz=0;jz<mesh->ngz;jz++)
	  block->data[jx][jy][jz] += rhs;
  }else {
    memblock3d *nb = newBlock();
    
    for(jx=0;jx<mesh->ngx;jx++)
      for(jy=0;jy<mesh->ngy;jy++)
	for(jz=0;jz<mesh->ngz;jz++)
	  nb->data[jx][jy][jz] = block->data[jx][jy][jz] + rhs;

    block->refs--;
    block = nb;
  }

#ifdef CHECK
  msg_stack.pop();
#endif
  
  return *this;
}

/////////////////// SUBTRACTION /////////////////////////

Field3D & Field3D::operator-=(const Field3D &rhs)
{
  int jx, jy, jz;

#ifdef CHECK
  msg_stack.push("Field3D: -= ( Field3D )");
  rhs.checkData();
  checkData();
#endif

  if(mesh->StaggerGrids && (rhs.location != location)) {
    // Interpolate and call again
    return (*this) -= interp_to(rhs, location);
  }

#ifdef TRACK
  name = "(" + name + "-" + rhs.name + ")";
#endif

  if(block->refs == 1) {
    for(jx=0;jx<mesh->ngx;jx++)
      for(jy=0;jy<mesh->ngy;jy++)
	for(jz=0;jz<mesh->ngz;jz++)
	  block->data[jx][jy][jz] -= rhs.block->data[jx][jy][jz];
  }else {
    memblock3d *nb = newBlock();
    
    for(jx=0;jx<mesh->ngx;jx++)
      for(jy=0;jy<mesh->ngy;jy++)
	for(jz=0;jz<mesh->ngz;jz++)
	  nb->data[jx][jy][jz] = block->data[jx][jy][jz] - rhs.block->data[jx][jy][jz];

    block->refs--;
    block = nb;
  }
  
#ifdef CHECK
  msg_stack.pop();
#endif
  
  return(*this);
}

Field3D & Field3D::operator-=(const Field2D &rhs)
{
  int jx, jy, jz;
  BoutReal **d;

#ifdef CHECK
  msg_stack.push("Field3D: -= ( Field2D )");
  rhs.checkData();
  checkData();
#endif

  d = rhs.getData();

#ifdef TRACK
  name = "(" + name + "-" + rhs.name + ")";
#endif

  if(block->refs == 1) {
    for(jx=0;jx<mesh->ngx;jx++)
      for(jy=0;jy<mesh->ngy;jy++)
	for(jz=0;jz<mesh->ngz;jz++)
	  block->data[jx][jy][jz] -= d[jx][jy];

  }else {
    memblock3d *nb = newBlock();

    for(jx=0;jx<mesh->ngx;jx++)
      for(jy=0;jy<mesh->ngy;jy++)
	for(jz=0;jz<mesh->ngz;jz++)
	  nb->data[jx][jy][jz] = block->data[jx][jy][jz] - d[jx][jy];

    block->refs--;
    block = nb;
  }

#ifdef CHECK
  msg_stack.pop();
#endif

  return(*this);
}

Field3D & Field3D::operator-=(const BoutReal &rhs)
{
  int jx, jy, jz;

#ifdef CHECK
  msg_stack.push("Field3D: -= ( BoutReal )");
  checkData();

  if(!finite(rhs)) {
    error("Field3D: -= operator passed non-finite BoutReal number");
  }
#endif

#ifdef TRACK
  name = "(" + name + "-BoutReal)";
#endif
  
  if(block->refs == 1) {
    for(jx=0;jx<mesh->ngx;jx++)
      for(jy=0;jy<mesh->ngy;jy++)
	for(jz=0;jz<mesh->ngz;jz++)
	  block->data[jx][jy][jz] -= rhs;
  }else {
    memblock3d *nb = newBlock();
    
    for(jx=0;jx<mesh->ngx;jx++)
      for(jy=0;jy<mesh->ngy;jy++)
	for(jz=0;jz<mesh->ngz;jz++)
	  nb->data[jx][jy][jz] = block->data[jx][jy][jz] - rhs;

    block->refs--;
    block = nb;
  }
  
#ifdef CHECK
  msg_stack.pop();
#endif

  return *this;
}

/////////////////// MULTIPLICATION ///////////////////////

Field3D & Field3D::operator*=(const Field3D &rhs)
{
  int jx, jy, jz;

#ifdef CHECK
  msg_stack.push("Field3D: *= ( Field3D )");

  rhs.checkData();
  checkData();
#endif

  if(mesh->StaggerGrids && (rhs.location != location)) {
    // Interpolate and call again
    return (*this) *= interp_to(rhs, location);
  }

#ifdef TRACK
  name = "(" + name + "*" + rhs.name + ")";
#endif

  if(block->refs == 1) {
    for(jx=0;jx<mesh->ngx;jx++)
      for(jy=0;jy<mesh->ngy;jy++)
	for(jz=0;jz<mesh->ngz;jz++)
	  block->data[jx][jy][jz] *= rhs.block->data[jx][jy][jz];
  }else {
    memblock3d *nb = newBlock();
    
    for(jx=0;jx<mesh->ngx;jx++)
      for(jy=0;jy<mesh->ngy;jy++)
	for(jz=0;jz<mesh->ngz;jz++)
	  nb->data[jx][jy][jz] = block->data[jx][jy][jz] * rhs.block->data[jx][jy][jz];

    block->refs--;
    block = nb;
  }

#ifdef CHECK
  msg_stack.pop();
#endif

  return(*this);
}

Field3D & Field3D::operator*=(const Field2D &rhs)
{
  int jx, jy, jz;
  BoutReal **d;

#ifdef CHECK
  msg_stack.push("Field3D: *= ( Field2D )");
  rhs.checkData();
  checkData();
#endif

  d = rhs.getData();

#ifdef TRACK
  name = "(" + name + "*"+rhs.name+")";
#endif

  if(block->refs == 1) {
    for(jx=0;jx<mesh->ngx;jx++)
      for(jy=0;jy<mesh->ngy;jy++)
	for(jz=0;jz<mesh->ngz;jz++)
	  block->data[jx][jy][jz] *= d[jx][jy];
  }else {
    memblock3d *nb = newBlock();

    for(jx=0;jx<mesh->ngx;jx++)
      for(jy=0;jy<mesh->ngy;jy++)
	for(jz=0;jz<mesh->ngz;jz++)
	  nb->data[jx][jy][jz] = block->data[jx][jy][jz] * d[jx][jy];

    block->refs--;
    block = nb;
  }

#ifdef CHECK
  msg_stack.pop();
#endif
  
  return(*this);
}

Field3D & Field3D::operator*=(const BoutReal rhs)
{
  int jx, jy, jz;
  
#ifdef CHECK
  msg_stack.push("Field3D: *= ( BoutReal )");
  checkData();

  if(!finite(rhs)) {
    error("Field3D: *= operator passed non-finite BoutReal number");
  }
#endif

#ifdef TRACK
  name = "(" + name + "*BoutReal)";
#endif

  if(block->refs == 1) {
    for(jx=0;jx<mesh->ngx;jx++)
      for(jy=0;jy<mesh->ngy;jy++)
	for(jz=0;jz<mesh->ngz;jz++)
	  block->data[jx][jy][jz] *= rhs;

  }else {
    memblock3d *nb = newBlock();

    for(jx=0;jx<mesh->ngx;jx++)
      for(jy=0;jy<mesh->ngy;jy++)
	for(jz=0;jz<mesh->ngz;jz++)
	  nb->data[jx][jy][jz] = block->data[jx][jy][jz] * rhs;

    block->refs--;
    block = nb;
  }

#ifdef CHECK
  msg_stack.pop();
#endif
  
  return(*this);
}

//////////////////// DIVISION /////////////////////

Field3D & Field3D::operator/=(const Field3D &rhs)
{
  int jx, jy, jz;

  if(mesh->StaggerGrids && (rhs.location != location)) {
    // Interpolate and call again
    
#ifdef CHECK
    msg_stack.push("Field3D /= Interpolating");
#endif
    (*this) /= interp_to(rhs, location);
#ifdef CHECK
    msg_stack.pop();
#endif
    return *this;
  }

#ifdef CHECK
  msg_stack.push("Field3D: /= ( Field3D )");
  rhs.checkData();
  checkData();
#endif

#ifdef TRACK
  name = "(" + name + "/" + rhs.name+")";
#endif

  if(block->refs == 1) {
    for(jx=0;jx<mesh->ngx;jx++)
      for(jy=0;jy<mesh->ngy;jy++)
	for(jz=0;jz<mesh->ngz;jz++)
	  block->data[jx][jy][jz] /= rhs.block->data[jx][jy][jz];
    
  }else {
    memblock3d *nb = newBlock();

    for(jx=0;jx<mesh->ngx;jx++)
      for(jy=0;jy<mesh->ngy;jy++)
	for(jz=0;jz<mesh->ngz;jz++)
	  nb->data[jx][jy][jz] = block->data[jx][jy][jz] / rhs.block->data[jx][jy][jz];

    block->refs--;
    block = nb;
  }

#ifdef CHECK
  msg_stack.pop();
#endif

  return(*this);
}

Field3D & Field3D::operator/=(const Field2D &rhs)
{
  int jx, jy, jz;
  BoutReal **d;

#ifdef CHECK
  msg_stack.push("Field3D: /= ( Field2D )");
  rhs.checkData();
  checkData();
#endif

  d = rhs.getData();

#ifdef TRACK
  name = "(" + name + "/" + rhs.name+")";
#endif

  /// NOTE: The faster version of this uses multiplications rather than division
  /// This seems to be enough to make a difference to the number of CVODE steps
  /// Hence for now straight division is used

  if(block->refs == 1) {
    for(jx=0;jx<mesh->ngx;jx++)
      for(jy=0;jy<mesh->ngy;jy++) {
	BoutReal val = 1.0L / d[jx][jy]; // Because multiplications are faster than divisions
	for(jz=0;jz<mesh->ngz;jz++)
	  block->data[jx][jy][jz] *= val;
	  //block->data[jx][jy][jz] /= d[jx][jy];
      }
  }else {
    memblock3d *nb = newBlock();

    for(jx=0;jx<mesh->ngx;jx++)
      for(jy=0;jy<mesh->ngy;jy++) {
	BoutReal val = 1.0L / d[jx][jy];
	for(jz=0;jz<mesh->ngz;jz++)
	  nb->data[jx][jy][jz] = block->data[jx][jy][jz] * val;
	//nb->data[jx][jy][jz] = block->data[jx][jy][jz] / d[jx][jy];
      }

    block->refs--;
    block = nb;
  }
  
#ifdef CHECK
  msg_stack.pop();
#endif

  return(*this);
}

Field3D & Field3D::operator/=(const BoutReal rhs)
{
  int jx, jy, jz;
  
#ifdef CHECK
  msg_stack.push("Field3D: /= ( BoutReal )");
  checkData();

  if(!finite(rhs)) {
    error("Field3D: /= operator passed non-finite BoutReal number");
  }
#endif

#ifdef TRACK
  name = "(" + name + "/BoutReal)";
#endif

  BoutReal val = 1.0 / rhs; // Because multiplication faster than division

  if(block->refs == 1) {
    for(jx=0;jx<mesh->ngx;jx++)
      for(jy=0;jy<mesh->ngy;jy++)
	for(jz=0;jz<mesh->ngz;jz++)
	  block->data[jx][jy][jz] *= val;
  }else {
    memblock3d *nb = newBlock();
    
    for(jx=0;jx<mesh->ngx;jx++)
      for(jy=0;jy<mesh->ngy;jy++)
	for(jz=0;jz<mesh->ngz;jz++)
	  nb->data[jx][jy][jz] = block->data[jx][jy][jz] * val;

    block->refs--;
    block = nb;
  }
  
#ifdef CHECK
  msg_stack.pop();
#endif

  return(*this);
}

////////////////// EXPONENTIATION ///////////////////////

Field3D & Field3D::operator^=(const Field3D &rhs)
{
  int jx, jy, jz;

  if(mesh->StaggerGrids && (rhs.location != location)) {
    // Interpolate and call again
    
#ifdef CHECK
    msg_stack.push("Field3D ^= Interpolating");
#endif
    (*this) ^= interp_to(rhs, location);
#ifdef CHECK
    msg_stack.pop();
#endif
    return *this;
  }

#ifdef CHECK
  msg_stack.push("Field3D: ^= ( Field3D )");
  rhs.checkData();
  checkData();
#endif

#ifdef TRACK
  name = "(" + name + "^"+rhs.name+")";
#endif

  if(block->refs == 1) {
    for(jx=0;jx<mesh->ngx;jx++)
      for(jy=0;jy<mesh->ngy;jy++)
	for(jz=0;jz<mesh->ngz;jz++)
	  block->data[jx][jy][jz] = pow(block->data[jx][jy][jz], rhs.block->data[jx][jy][jz]);

  }else {
    memblock3d *nb = newBlock();
    
    for(jx=0;jx<mesh->ngx;jx++)
      for(jy=0;jy<mesh->ngy;jy++)
	for(jz=0;jz<mesh->ngz;jz++)
	  nb->data[jx][jy][jz] = pow(block->data[jx][jy][jz], rhs.block->data[jx][jy][jz]);
    
    block->refs--;
    block = nb;
  }

#ifdef CHECK
  msg_stack.pop();
#endif

  return(*this);
}

Field3D & Field3D::operator^=(const Field2D &rhs)
{
  int jx, jy, jz;
  BoutReal **d;

#ifdef CHECK
  msg_stack.push("Field3D: ^= ( Field2D )");
  rhs.checkData();
  checkData();
#endif

  d = rhs.getData();

#ifdef TRACK
  name = "(" + name + "^"+rhs.name+")";
#endif

  if(block->refs == 1) {
    for(jx=0;jx<mesh->ngx;jx++)
      for(jy=0;jy<mesh->ngy;jy++)
	for(jz=0;jz<mesh->ngz;jz++)
	  block->data[jx][jy][jz] = pow(block->data[jx][jy][jz], d[jx][jy]);

  }else {
    memblock3d *nb = newBlock();

    for(jx=0;jx<mesh->ngx;jx++)
      for(jy=0;jy<mesh->ngy;jy++)
	for(jz=0;jz<mesh->ngz;jz++)
	  nb->data[jx][jy][jz] = pow(block->data[jx][jy][jz], d[jx][jy]);

    block->refs--;
    block = nb;
  }
  
#ifdef CHECK
  msg_stack.pop();
#endif

  return(*this);
}

Field3D & Field3D::operator^=(const BoutReal rhs)
{
  int jx, jy, jz;

#ifdef CHECK
  msg_stack.push("Field3D: ^= ( BoutReal )");
  checkData();

  if(!finite(rhs)) {
    error("Field3D: ^= operator passed non-finite BoutReal number");
  }
#endif

#ifdef TRACK
  name = "(" + name + "^BoutReal)";
#endif

  if(block->refs == 1) {
    for(jx=0;jx<mesh->ngx;jx++)
      for(jy=0;jy<mesh->ngy;jy++)
	for(jz=0;jz<mesh->ngz;jz++)
	  block->data[jx][jy][jz] = pow(block->data[jx][jy][jz], rhs);

  }else {
    memblock3d *nb = newBlock();

    for(jx=0;jx<mesh->ngx;jx++)
      for(jy=0;jy<mesh->ngy;jy++)
	for(jz=0;jz<mesh->ngz;jz++)
	  nb->data[jx][jy][jz] = pow(block->data[jx][jy][jz], rhs);

    block->refs--;
    block = nb;
  }
  
#ifdef CHECK
  msg_stack.pop();
#endif

  return(*this);
}


/***************************************************************
 *                      BINARY OPERATORS 
 ***************************************************************/

/////////////////// ADDITION ///////////////////

const Field3D Field3D::operator+() const
{
  Field3D result = *this;

#ifdef TRACK
  result.name = name;
#endif
  
  return result;
}


const Field3D Field3D::operator+(const Field3D &other) const
{
  Field3D result = *this;
  result += other;
  return(result);
}

const Field3D Field3D::operator+(const Field2D &other) const
{
  Field3D result = *this;
  result += other;
  return(result);
}

const FieldPerp Field3D::operator+(const FieldPerp &other) const
{
  FieldPerp result = other;
  result += (*this);
  return(result);
}

const Field3D Field3D::operator+(const BoutReal &rhs) const
{
  Field3D result = *this;
  result += rhs;
  return(result);
}

/////////////////// SUBTRACTION ////////////////

const Field3D Field3D::operator-() const
{
  Field3D result = *this;

  result *= -1.0;

#ifdef TRACK
  result.name = "(-"+name+")";
#endif
  
  return result;
}

const Field3D Field3D::operator-(const Field3D &other) const
{
  Field3D result = *this;
  result -= other;
  return(result);
}

const Field3D Field3D::operator-(const Field2D &other) const
{
  Field3D result = *this;
  result -= other;
  return(result);
}

const FieldPerp Field3D::operator-(const FieldPerp &other) const
{
  BoutReal **d;
  int jx, jz;
  int jy = other.getIndex();
  FieldPerp result = other;

#ifdef CHECK
  if(block == NULL)
    error("Field3D: - FieldPerp operates on empty data");
#endif

  d = result.getData();
  for(jx=0;jx<mesh->ngx;jx++)
    for(jz=0;jz<mesh->ngz;jz++)
      d[jx][jz] = block->data[jx][jy][jz] - d[jx][jz];
  
  return(result);
}

const Field3D Field3D::operator-(const BoutReal &rhs) const
{
  Field3D result = *this;
  result -= rhs;
  return(result);
}

///////////////// MULTIPLICATION ///////////////

const Field3D Field3D::operator*(const Field3D &other) const
{
  Field3D result = *this;
  result *= other;
  return(result);
}

const Field3D Field3D::operator*(const Field2D &other) const
{
  Field3D result = *this;
  result *= other;
  return(result);
}

const FieldPerp Field3D::operator*(const FieldPerp &other) const
{
  FieldPerp result = other;
  result *= (*this);
  return(result);
}

const Field3D Field3D::operator*(const BoutReal rhs) const
{
  Field3D result = *this;
  result *= rhs;
  return(result);
}

//////////////////// DIVISION ////////////////////

const Field3D Field3D::operator/(const Field3D &other) const
{
  Field3D result = *this;
  result /= other;
  return(result);
}

const Field3D Field3D::operator/(const Field2D &other) const
{
  Field3D result = *this;
  result /= other;
  return(result);
}

const FieldPerp Field3D::operator/(const FieldPerp &other) const
{
  BoutReal **d;
  int jx, jz;
  int jy = other.getIndex();
  FieldPerp result = other;
  
#ifdef CHECK
  if(block == NULL)
    error("Field3D: / FieldPerp operates on empty data");
#endif

  d = result.getData();
  for(jx=0;jx<mesh->ngx;jx++)
    for(jz=0;jz<mesh->ngz;jz++)
      d[jx][jz] = block->data[jx][jy][jz] / d[jx][jz];
  
#ifdef TRACK
  result.name = "(" + name + "/" + other.name + ")";
#endif  
  
  return(result);
}

const Field3D Field3D::operator/(const BoutReal rhs) const
{
  Field3D result = *this;
  result /= rhs;
  return(result);
}

////////////// EXPONENTIATION /////////////////

const Field3D Field3D::operator^(const Field3D &other) const
{
  Field3D result = *this;
  result ^= other;
  return(result);
}

const Field3D Field3D::operator^(const Field2D &other) const
{
  Field3D result = *this;
  result ^= other;
  return(result);
}

const FieldPerp Field3D::operator^(const FieldPerp &other) const
{
  BoutReal **d;
  int jx, jz;
  int jy = other.getIndex();
  FieldPerp result = other;
  
#ifdef CHECK
  if(block == NULL)
    error("Field3D: ^ FieldPerp operates on empty data");
#endif

  d = result.getData();
  for(jx=0;jx<mesh->ngx;jx++)
    for(jz=0;jz<mesh->ngz;jz++)
      d[jx][jz] = pow(block->data[jx][jy][jz], d[jx][jz]);
  
#ifdef TRACK
  result.name = "("+name+"^"+other.name + ")";
#endif

  return(result);
}

const Field3D Field3D::operator^(const BoutReal rhs) const
{
  Field3D result = *this;
  result ^= rhs;
  return(result);
}

/***************************************************************
 *                         STENCILS
 ***************************************************************/

/*!
  18 Aug 2008: Added need_x argument to disable Z interpolation when not needed
  since interpZ is taking ~25% of run time (on single processor)
*/
void Field3D::setStencil(bstencil *fval, bindex *bx, bool need_x) const
{
  fval->jx = bx->jx;
  fval->jy = bx->jy;
  fval->jz = bx->jz;
  
#ifdef CHECK
  // Check data set
  if(block == NULL) {
    error("Field3D: Setting stencil for empty data\n");
    exit(1);
  }
#endif

  fval->cc = block->data[bx->jx][bx->jy][bx->jz];

  if(need_x) {
    if(mesh->ShiftXderivs) {
      fval->xp = interpZ(bx->jxp, bx->jy, bx->jz, bx->xp_offset, mesh->ShiftOrder);
      fval->xm = interpZ(bx->jxm, bx->jy, bx->jz, bx->xm_offset, mesh->ShiftOrder);
      fval->x2p = interpZ(bx->jx2p, bx->jy, bx->jz, bx->x2p_offset, mesh->ShiftOrder);
      fval->x2m = interpZ(bx->jx2m, bx->jy, bx->jz, bx->x2m_offset, mesh->ShiftOrder);
    }else {
      // No shift in the z direction
      fval->xp = block->data[bx->jxp][bx->jy][bx->jz];
      fval->xm = block->data[bx->jxm][bx->jy][bx->jz];
      fval->x2p = block->data[bx->jx2p][bx->jy][bx->jz];
      fval->x2m = block->data[bx->jx2m][bx->jy][bx->jz];
    }
  }

  // TWIST-SHIFT CONDITION
  if(bx->yp_shift) {
    fval->yp = interpZ(bx->jx, bx->jyp, bx->jz, bx->yp_offset, mesh->TwistOrder);
  }else
    fval->yp = block->data[bx->jx][bx->jyp][bx->jz];
  
  if(bx->ym_shift) {
    fval->ym = interpZ(bx->jx, bx->jym, bx->jz, bx->ym_offset, mesh->TwistOrder);
  }else
    fval->ym = block->data[bx->jx][bx->jym][bx->jz];

  if(bx->y2p_shift) {
    fval->y2p = interpZ(bx->jx, bx->jy2p, bx->jz, bx->yp_offset, mesh->TwistOrder);
  }else
    fval->y2p = block->data[bx->jx][bx->jy2p][bx->jz];

  if(bx->y2m_shift) {
    fval->y2m = interpZ(bx->jx, bx->jy2m, bx->jz, bx->ym_offset, mesh->TwistOrder);
  }else
    fval->y2m = block->data[bx->jx][bx->jy2m][bx->jz];

  // Z neighbours
  
  fval->zp = block->data[bx->jx][bx->jy][bx->jzp];
  fval->zm = block->data[bx->jx][bx->jy][bx->jzm];
  fval->z2p = block->data[bx->jx][bx->jy][bx->jz2p];
  fval->z2m = block->data[bx->jx][bx->jy][bx->jz2m];
}

void Field3D::setXStencil(stencil &fval, const bindex &bx, CELL_LOC loc) const
{
  fval.jx = bx.jx;
  fval.jy = bx.jy;
  fval.jz = bx.jz;
  
#ifdef CHECK
  // Check data set
  if(block == NULL) {
    error("Field3D: Setting X stencil for empty data\n");
    exit(1);
  }
#endif

  fval.c = block->data[bx.jx][bx.jy][bx.jz];

  if(mesh->ShiftXderivs && (mesh->ShiftOrder != 0)) {
    fval.p = interpZ(bx.jxp, bx.jy, bx.jz, bx.xp_offset, mesh->ShiftOrder);
    fval.m = interpZ(bx.jxm, bx.jy, bx.jz, bx.xm_offset, mesh->ShiftOrder);
    fval.pp = interpZ(bx.jxp, bx.jy, bx.jz, bx.x2p_offset, mesh->ShiftOrder);
    fval.mm = interpZ(bx.jxm, bx.jy, bx.jz, bx.x2m_offset, mesh->ShiftOrder);
  }else {
    // No shift in the z direction
    fval.p = block->data[bx.jxp][bx.jy][bx.jz];
    fval.m = block->data[bx.jxm][bx.jy][bx.jz];
    fval.pp = block->data[bx.jx2p][bx.jy][bx.jz];
    fval.mm = block->data[bx.jx2m][bx.jy][bx.jz];
  }

  if(mesh->StaggerGrids && (loc != CELL_DEFAULT) && (loc != location)) {
    // Non-centred stencil

    if((location == CELL_CENTRE) && (loc == CELL_XLOW)) {
      // Producing a stencil centred around a lower X value
      fval.pp = fval.p;
      fval.p  = fval.c;
      
    }else if(location == CELL_XLOW) {
      // Stencil centred around a cell centre
      
      fval.mm = fval.m;
      fval.m  = fval.c;
    }
    // Shifted in one direction -> shift in another
    // Could produce warning
  }
}

void Field3D::setYStencil(stencil &fval, const bindex &bx, CELL_LOC loc) const
{
  fval.jx = bx.jx;
  fval.jy = bx.jy;
  fval.jz = bx.jz;
  
#ifdef CHECK
  // Check data set
  if(block == NULL) {
    error("Field3D: Setting Y stencil for empty data\n");
    exit(1);
  }
#endif

  fval.c = block->data[bx.jx][bx.jy][bx.jz];

  
  //if((!TwistShift) || (mesh->TwistOrder == 0)) {
    // Either no twist-shift, or already done in communicator
    
    fval.p = block->data[bx.jx][bx.jyp][bx.jz];
    fval.m = block->data[bx.jx][bx.jym][bx.jz];
    fval.pp = block->data[bx.jx][bx.jy2p][bx.jz];
    fval.mm = block->data[bx.jx][bx.jy2m][bx.jz];
    /*
  }else {
    // TWIST-SHIFT CONDITION
    if(bx.yp_shift) {
      fval.p = interpZ(bx.jx, bx.jyp, bx.jz, bx.yp_offset, mesh->TwistOrder);
    }else
      fval.p = block->data[bx.jx][bx.jyp][bx.jz];
    
    if(bx.ym_shift) {
      fval.m = interpZ(bx.jx, bx.jym, bx.jz, bx.ym_offset, mesh->TwistOrder);
    }else
      fval.m = block->data[bx.jx][bx.jym][bx.jz];
    
    if(bx.y2p_shift) {
      fval.pp = interpZ(bx.jx, bx.jy2p, bx.jz, bx.yp_offset, mesh->TwistOrder);
    }else
      fval.pp = block->data[bx.jx][bx.jy2p][bx.jz];
    
    if(bx.y2m_shift) {
      fval.mm = interpZ(bx.jx, bx.jy2m, bx.jz, bx.ym_offset, mesh->TwistOrder);
    }else
      fval.mm = block->data[bx.jx][bx.jy2m][bx.jz];
  }
    */

  if(mesh->StaggerGrids && (loc != CELL_DEFAULT) && (loc != location)) {
    // Non-centred stencil

    if((location == CELL_CENTRE) && (loc == CELL_YLOW)) {
      // Producing a stencil centred around a lower Y value
      fval.pp = fval.p;
      fval.p  = fval.c;
    }else if(location == CELL_YLOW) {
      // Stencil centred around a cell centre
      
      fval.mm = fval.m;
      fval.m  = fval.c;
    }
    // Shifted in one direction -> shift in another
    // Could produce warning
  }
}

void Field3D::setZStencil(stencil &fval, const bindex &bx, CELL_LOC loc) const
{
  fval.jx = bx.jx;
  fval.jy = bx.jy;
  fval.jz = bx.jz;
  
#ifdef CHECK
  // Check data set
  if(block == NULL) {
    error("Field3D: Setting stencil for empty data\n");
    exit(1);
  }
#endif

  fval.c = block->data[bx.jx][bx.jy][bx.jz];

  fval.p = block->data[bx.jx][bx.jy][bx.jzp];
  fval.m = block->data[bx.jx][bx.jy][bx.jzm];
  fval.pp = block->data[bx.jx][bx.jy][bx.jz2p];
  fval.mm = block->data[bx.jx][bx.jy][bx.jz2m];

  if(mesh->StaggerGrids && (loc != CELL_DEFAULT) && (loc != location)) {
    // Non-centred stencil

    if((location == CELL_CENTRE) && (loc == CELL_ZLOW)) {
      // Producing a stencil centred around a lower Z value
      fval.pp = fval.p;
      fval.p  = fval.c;
      
    }else if(location == CELL_ZLOW) {
      // Stencil centred around a cell centre
      
      fval.mm = fval.m;
      fval.m  = fval.c;
    }
    // Shifted in one direction -> shift in another
    // Could produce warning
  }
}

BoutReal Field3D::interpZ(int jx, int jy, int jz0, BoutReal zoffset, int order) const
{
  int zi;
  BoutReal result;
  int jzp, jzm, jz2p;

  zi = ROUND(zoffset);  // Find the nearest integer
  zoffset -= (BoutReal) zi; // Difference (-0.5 to +0.5)

  if((zoffset < 0.0) && (order > 1)) { // If order = 0 or 1, want closest
    // For higher-order interpolation, expect zoffset > 0
    zi--;
    zoffset += 1.0;
  }
  
  int ncz = mesh->ngz - 1;
  jz0 = (((jz0 + zi)%ncz) + ncz) % ncz;
  jzp = (jz0 + 1) % ncz;
  jz2p = (jz0 + 2) % ncz;
  jzm = (jz0 - 1 + ncz) % ncz;

  switch(order) {
  case 2: {
    // 2-point linear interpolation

    result = (1.0 - zoffset)*block->data[jx][jy][jz0] + zoffset*block->data[jx][jy][jzp];

    break;
  }
  case 3: {
    // 3-point Lagrange interpolation

    result = 0.5*zoffset*(zoffset-1.0)*block->data[jx][jy][jzm]
      + (1.0 - zoffset*zoffset)*block->data[jx][jy][jz0]
      + 0.5*zoffset*(zoffset + 1.0)*block->data[jx][jy][jzp];
    break;
  }
  case 4: {
    // 4-point Lagrange interpolation
    result = -zoffset*(zoffset-1.0)*(zoffset-2.0)*block->data[jx][jy][jzm]/6.0
      + 0.5*(zoffset*zoffset - 1.0)*(zoffset-2.0)*block->data[jx][jy][jz0]
      - 0.5*zoffset*(zoffset+1.0)*(zoffset-2.0)*block->data[jx][jy][jzp]
      + zoffset*(zoffset*zoffset - 1.0)*block->data[jx][jy][jz2p]/6.0;
    break;
  }
  default: {
    // Nearest neighbour
    result = block->data[jx][jy][jz0];
  }
  };
  return result;
}

void Field3D::shiftZ(int jx, int jy, double zangle)
{
  static dcomplex *v = (dcomplex*) NULL;
  int jz;
  BoutReal kwave;
  
#ifdef CHECK
  // Check data set
  if(block == NULL) {
    error("Field3D: Shifting in Z an empty data set\n");
    exit(1);
  }
#endif

  int ncz = mesh->ngz-1;

  if(ncz == 1)
    return;

 allocate();

  if(v == (dcomplex*) NULL) {
    //allocate memory
    v = new dcomplex[ncz/2 + 1];
  }
  
  rfft(block->data[jx][jy], ncz, v); // Forward FFT

  // Apply phase shift
  for(jz=1;jz<=ncz/2;jz++) {
    kwave=jz*2.0*PI/mesh->zlength; // wave number is 1/[rad]
    v[jz] *= dcomplex(cos(kwave*zangle) , -sin(kwave*zangle));
  }

  irfft(v, ncz, block->data[jx][jy]); // Reverse FFT

  block->data[jx][jy][ncz] = block->data[jx][jy][0];
}

const Field3D Field3D::shiftZ(const Field2D zangle) const
{
  Field3D result;
  int jx, jy;

#ifdef CHECK
  msg_stack.push("Field3D: shiftZ ( Field2D )");
  checkData();
#endif

  result = *this;

  for(jx=0;jx<mesh->ngx;jx++) {
    for(jy=0;jy<mesh->ngy;jy++) {
      result.shiftZ(jx, jy, zangle[jx][jy]);
    }
  }

#ifdef CHECK
  msg_stack.pop();
#endif

  return result;
}

const Field3D Field3D::shiftZ(const BoutReal zangle) const
{
  Field3D result;
  int jx, jy;

#ifdef CHECK
  msg_stack.push("Field3D: shiftZ ( BoutReal )");
  checkData();
#endif

  result = *this;

  for(jx=0;jx<mesh->ngx;jx++) {
    for(jy=0;jy<mesh->ngy;jy++) {
      result.shiftZ(jx, jy, zangle);
    }
  }

#ifdef CHECK
  msg_stack.pop();
#endif

  return result;
}

const Field3D Field3D::shiftZ(bool toBoutReal) const
{
  if(toBoutReal) {
    return shiftZ(mesh->zShift);
  }
  return shiftZ(-mesh->zShift);
}


/***************************************************************
 *                         SLICING
 ***************************************************************/

void Field3D::getXArray(int y, int z, rvec &xv) const
{
#ifdef CHECK
  // Check data set
  if(block == NULL) {
    error("Field3D: getXArray on an empty data set\n");
    exit(1);
  }
#endif

  xv.resize(mesh->ngx);
  
  for(int x=0;x<mesh->ngx;x++)
    xv[x] = block->data[x][y][z];
}

void Field3D::getYArray(int x, int z, rvec &yv) const
{
#ifdef CHECK
  // Check data set
  if(block == NULL) {
    error("Field3D: getYArray on an empty data set\n");
    exit(1);
  }
#endif

  yv.resize(mesh->ngy);
  
  for(int y=0;y<mesh->ngy;y++)
    yv[y] = block->data[x][y][z];
}

void Field3D::getZArray(int x, int y, rvec &zv) const
{
#ifdef CHECK
  // Check data set
  if(block == NULL) {
    error("Field3D: getZArray on an empty data set\n");
    exit(1);
  }
#endif

  zv.resize(mesh->ngz-1);
  
  for(int z=0;z<mesh->ngz-1;z++)
    zv[z] = block->data[x][y][z];
}

void Field3D::setXArray(int y, int z, const rvec &xv)
{
 allocate();

#ifdef CHECK
  // Check that vector is correct size
  if(xv.capacity() != (unsigned int) mesh->ngx) {
    error("Field3D: setXArray has incorrect size\n");
    exit(1);
  }
#endif

  for(int x=0;x<mesh->ngx;x++)
    block->data[x][y][z] = xv[x];
}

void Field3D::setYArray(int x, int z, const rvec &yv)
{
  allocate();

#ifdef CHECK
  // Check that vector is correct size
  if(yv.capacity() != (unsigned int) mesh->ngy) {
    error("Field3D: setYArray has incorrect size\n");
    exit(1);
  }
#endif

  for(int y=0;y<mesh->ngy;y++)
    block->data[x][y][z] = yv[y];
}

void Field3D::setZArray(int x, int y, const rvec &zv)
{
  allocate();

#ifdef CHECK
  // Check that vector is correct size
  if(zv.capacity() != (unsigned int) (mesh->ngz-1)) {
    error("Field3D: setZArray has incorrect size\n");
    exit(1);
  }
#endif

  for(int z=0;z<(mesh->ngz-1);z++)
    block->data[x][y][z] = zv[z];
}

const FieldPerp Field3D::slice(int y) const
{
  FieldPerp result;

  result.set(*this, y);

#ifdef TRACK
  result.name = "Slice("+name+")";
#endif

  return(result);
}

/***************************************************************
 *                      MATH FUNCTIONS
 ***************************************************************/

const Field3D Field3D::sqrt() const
{
  int jx, jy, jz;
  Field3D result;

#ifdef CHECK
  msg_stack.push("Field3D: Sqrt()");

  // Check data set
  if(block == NULL) {
    error("Field3D: Taking sqrt of empty data\n");
    exit(1);
  }
    
  // Test values
  for(jx=mesh->xstart;jx<=mesh->xend;jx++)
    for(jy=mesh->ystart;jy<=mesh->yend;jy++) 
      for(jz=0;jz<mesh->ngz-1;jz++) {
	if(block->data[jx][jy][jz] < 0.0) {
	  error("Field3D: Sqrt operates on negative value at [%d,%d,%d]\n", jx, jy, jz);
	}
      }
#endif

#ifdef TRACK
  result.name = "Sqrt("+name+")";
#endif

  result.allocate();

  for(jx=0;jx<mesh->ngx;jx++)
    for(jy=0;jy<mesh->ngy;jy++) 
      for(jz=0;jz<mesh->ngz;jz++)
	result.block->data[jx][jy][jz] = ::sqrt(block->data[jx][jy][jz]);

#ifdef CHECK
  msg_stack.pop();
#endif

  result.location = location;
  
  return result;
}

const Field3D Field3D::abs() const
{
  int jx, jy, jz;
  Field3D result;

#ifdef CHECK
  // Check data set
  if(block == NULL) {
    error("Field3D: Taking abs of empty data\n");
    exit(1);
  }
#endif

#ifdef TRACK
  result.name = "Abs("+name+")";
#endif

  result.allocate();

  for(jx=0;jx<mesh->ngx;jx++)
    for(jy=0;jy<mesh->ngy;jy++) 
      for(jz=0;jz<mesh->ngz;jz++)
	result.block->data[jx][jy][jz] = fabs(block->data[jx][jy][jz]);

  result.location = location;

  return result;
}

BoutReal Field3D::min(bool allpe) const
{
  int jx, jy, jz;
  BoutReal result;

#ifdef CHECK
  if(block == NULL) {
    error("Field3D: min() method on empty data");
    exit(1);
  }

  if(allpe) {
    msg_stack.push("Field3D::Min() over all PEs");
  }else
    msg_stack.push("Field3D::Min()");
#endif

  result = block->data[0][0][0];

  for(jx=0;jx<mesh->ngx;jx++)
    for(jy=0;jy<mesh->ngy;jy++)
      for(jz=0;jz<mesh->ngz;jz++)
	if(block->data[jx][jy][jz] < result)
	  result = block->data[jx][jy][jz];

  if(allpe) {
    // MPI reduce
    BoutReal localresult = result;
    MPI_Allreduce(&localresult, &result, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
  }

#ifdef CHECK
  msg_stack.pop();
#endif

  return result;
}

BoutReal Field3D::max(bool allpe) const
{
  int jx, jy, jz;
  BoutReal result;

#ifdef CHECK
  if(block == NULL) {
    error("Field3D: max() method on empty data");
    exit(1);
  }
  if(allpe) {
    msg_stack.push("Field3D::Max() over all PEs");
  }else
    msg_stack.push("Field3D::Max()");
#endif
  
  result = block->data[0][0][0];

  for(jx=0;jx<mesh->ngx;jx++)
    for(jy=0;jy<mesh->ngy;jy++)
      for(jz=0;jz<mesh->ngz;jz++)
	if(block->data[jx][jy][jz] > result)
	  result = block->data[jx][jy][jz];
  
  if(allpe) {
    // MPI reduce
    BoutReal localresult = result;
    MPI_Allreduce(&localresult, &result, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  }
  
#ifdef CHECK
  msg_stack.pop();
#endif

  return result;
}

///////////////////// FieldData VIRTUAL FUNCTIONS //////////

int Field3D::getData(int x, int y, int z, void *vptr) const
{
#ifdef CHECK
  // Check data set
  if(block ==  NULL) {
    error("Field3D: getData on empty data\n");
    exit(1);
  }
  
  // check ranges
  if((x < 0) || (x >= mesh->ngx) || (y < 0) || (y >= mesh->ngy) || (z < 0) || (z >= mesh->ngz)) {
    error("Field3D: getData (%d,%d,%d) out of bounds\n", x, y, z);
    exit(1);
  }
#endif
  BoutReal *ptr = (BoutReal*) vptr;
  *ptr = block->data[x][y][z];
  
  return sizeof(BoutReal);
}

int Field3D::getData(int x, int y, int z, BoutReal *rptr) const
{
#ifdef CHECK
  // Check data set
  if(block == NULL) {
    error("Field3D: getData on empty data\n");
    exit(1);
  }
  
  // check ranges
  if((x < 0) || (x >= mesh->ngx) || (y < 0) || (y >= mesh->ngy) || (z < 0) || (z >= mesh->ngz)) {
    error("Field3D: getData (%d,%d,%d) out of bounds\n", x, y, z);
    exit(1);
  }
#endif

  *rptr = block->data[x][y][z];
  return 1;
}

int Field3D::setData(int x, int y, int z, void *vptr)
{
  allocate();
#ifdef CHECK
  // check ranges
  if((x < 0) || (x >= mesh->ngx) || (y < 0) || (y >= mesh->ngy) || (z < 0) || (z >= mesh->ngz)) {
    error("Field3D: fillArray (%d,%d,%d) out of bounds\n", x, y, z);
    exit(1);
  }
#endif
  BoutReal *ptr = (BoutReal*) vptr;
  block->data[x][y][z] = *ptr;
  
  return sizeof(BoutReal);
}

int Field3D::setData(int x, int y, int z, BoutReal *rptr)
{
  allocate();
#ifdef CHECK
  // check ranges
  if((x < 0) || (x >= mesh->ngx) || (y < 0) || (y >= mesh->ngy) || (z < 0) || (z >= mesh->ngz)) {
    error("Field3D: setData (%d,%d,%d) out of bounds\n", x, y, z);
    exit(1);
  }
#endif

  block->data[x][y][z] = *rptr;
  return 1;
}

#ifdef CHECK
/// Check if the data is valid
bool Field3D::checkData(bool vital) const
{
  if(block ==  NULL) {
    error("Field3D: Operation on empty data\n");
    return true;
  }

  if( vital || ( CHECK > 2 ) ) { 
    // Do full checks
    // Normally this is done just for some operations (equalities)
    // If CHECKS > 2, all operations perform checks
    
    int jx, jy, jz;
    
    for(jx=mesh->xstart;jx<=mesh->xend;jx++)
      for(jy=mesh->ystart;jy<=mesh->yend;jy++)
	for(jz=0;jz<mesh->ngz-1;jz++)
	  if(!finite(block->data[jx][jy][jz])) {
	    error("Field3D: Operation on non-finite data at [%d][%d][%d]\n", jx, jy, jz);
	    return true;
	  }
  }
    
  return false;
}
#endif

void Field3D::cleanup()
{
  while(blocklist != NULL) {
    memblock3d *nb = blocklist->all_next;
    
    // Free the 3D data
    free_r3tensor(blocklist->data);
    // Delete the structure
    delete blocklist;
    // Move to the next one
    blocklist = nb;
    nblocks--;
  }
  
  // Reset to starting
  nblocks = 0;
  free_block = NULL;
}

///////////////////// BOUNDARY CONDITIONS //////////////////

void Field3D::setBackground(const Field2D &f2d)
{
  background = &f2d;
}

void Field3D::applyBoundary()
{
#ifdef CHECK
  msg_stack.push("Field3D::applyBoundary()");
  
  if(block == NULL)
    output << "WARNING: Applying boundary condition to empty data" << endl;
#endif
  
  if(block == NULL)
    return;
  
  if(background != NULL) {
    // Apply boundary to the total of this and background
    
    Field3D tot = *this + (*background);
    tot.applyBoundary();
    *this = tot - (*background);
  }else {
    // Apply boundary to this field
    for(vector<BoundaryOp*>::iterator it = bndry_op.begin(); it != bndry_op.end(); it++)
      (*it)->apply(*this);
  }
  
  // Set the corners to zero
  for(int jx=0;jx<mesh->xstart;jx++) {
    for(int jy=0;jy<mesh->ystart;jy++) {
      for(int jz=0;jz<mesh->ngz;jz++)
        block->data[jx][jy][jz] = 0.;
    }
    for(int jy=mesh->yend+1;jy<mesh->ngy;jy++) {
      for(int jz=0;jz<mesh->ngz;jz++)
        block->data[jx][jy][jz] = 0.;
    }
  }
  for(int jx=mesh->xend+1;jx<mesh->ngx;jx++) {
    for(int jy=0;jy<mesh->ystart;jy++) {
      for(int jz=0;jz<mesh->ngz;jz++)
        block->data[jx][jy][jz] = 0.;
    }
    for(int jy=mesh->yend+1;jy<mesh->ngy;jy++) {
      for(int jz=0;jz<mesh->ngz;jz++)
        block->data[jx][jy][jz] = 0.;
    }
  }

#ifdef CHECK
  msg_stack.pop();
#endif
}

void Field3D::applyTDerivBoundary()
{
#ifdef CHECK
  msg_stack.push("Field3D::applyTDerivBoundary()");

  if(block == NULL)
    output << "WARNING: Applying boundary condition to empty data" << endl;
#endif
  
  if(block == NULL)
    return;
  
  if(background != NULL)
    *this += *background;
    
  for(vector<BoundaryOp*>::iterator it = bndry_op.begin(); it != bndry_op.end(); it++)
    (*it)->apply_ddt(*this);
  
  if(background != NULL)
    *this -= *background;

  // Set the corners to zero
  for(int jx=0;jx<mesh->xstart;jx++) {
    for(int jy=0;jy<mesh->ystart;jy++) {
      for(int jz=0;jz<mesh->ngz;jz++)
        ddt->block->data[jx][jy][jz] = 0.;
    }
    for(int jy=mesh->yend+1;jy<mesh->ngy;jy++) {
      for(int jz=0;jz<mesh->ngz;jz++)
        ddt->block->data[jx][jy][jz] = 0.;
    }
  }
  for(int jx=mesh->xend+1;jx<mesh->ngx;jx++) {
    for(int jy=0;jy<mesh->ystart;jy++) {
      for(int jz=0;jz<mesh->ngz;jz++)
        ddt->block->data[jx][jy][jz] = 0.;
    }
    for(int jy=mesh->yend+1;jy<mesh->ngy;jy++) {
      for(int jz=0;jz<mesh->ngz;jz++)
        ddt->block->data[jx][jy][jz] = 0.;
    }
  }

#ifdef CHECK
  msg_stack.pop();
#endif
}

/***************************************************************
 *                     PRIVATE FUNCTIONS
 ***************************************************************/

// GLOBAL VARS

int Field3D::nblocks = 0;
memblock3d* Field3D::blocklist = NULL;
memblock3d* Field3D::free_block = NULL;

/// Get a new block of data, either from free list or allocate
memblock3d *Field3D::newBlock() const
{
  memblock3d *nb;

  if(free_block != NULL) {
    // just pop off the top of the stack
    nb = free_block;
    free_block = nb->next;
    nb->next = NULL;
    nb->refs = 1;
  }else {
    // No more blocks left - allocate a new block
    nb = new memblock3d;

    nb->data = r3tensor(mesh->ngx, mesh->ngy, mesh->ngz);
    nb->refs = 1;
    nb->next = NULL;
    
    // add to the global list
    nb->all_next = blocklist;
    blocklist = nb;
    
    nblocks++;
  }

  return nb;
}


/// Makes sure data is allocated and only referenced by this object
void Field3D::allocData() const
{
  /// Check if any data associated with this object
  if(block != (memblock3d*) NULL) {
    // Already has a block of data
    
    /// Check if data shared with other objects
    if(block->refs > 1) {
      // Need to get a new block and copy across

      memblock3d* nb = newBlock();

      for(int jx=0;jx<mesh->ngx;jx++)
	for(int jy=0;jy<mesh->ngy;jy++)
	  for(int jz=0;jz<mesh->ngz;jz++)
	    nb->data[jx][jy][jz] = block->data[jx][jy][jz];

      block->refs--;
      block = nb;
    }
  }else {
    // No data - get a new block

    block = newBlock();
    
  }
}

void Field3D::freeData()
{
  // put data block onto stack

  // Need to check for either no data, or all data has been cleared
  if((block == NULL) || (nblocks == 0))
    return;

  block->refs--;

  if(block->refs == 0) {
    // No more references to this data - put on free list

#ifdef DISABLE_FREELIST
    // For debugging, free memory
    free_r3tensor(block->data);
    delete block;
#else
    block->next = free_block;
    free_block = block;
#endif
  }

  block = NULL;
}

/***************************************************************
 *               NON-MEMBER OVERLOADED OPERATORS
 ***************************************************************/

const Field3D operator-(const BoutReal &lhs, const Field3D &rhs)
{
  Field3D result;
  int jx, jy, jz;

#ifdef TRACK
  result.name = "(BoutReal-"+rhs.name+")";
#endif

  result.allocate();
  
  for(jx=0;jx<mesh->ngx;jx++)
    for(jy=0;jy<mesh->ngy;jy++)
      for(jz=0;jz<mesh->ngz;jz++)
	result.block->data[jx][jy][jz] = lhs - rhs.block->data[jx][jy][jz];

  result.location = rhs.location;

  return result;
}

const Field3D operator+(const BoutReal &lhs, const Field3D &rhs)
{
  return rhs+lhs;
}

const Field3D operator*(const BoutReal lhs, const Field3D &rhs)
{
  return(rhs * lhs);
}

const Field3D operator/(const BoutReal lhs, const Field3D &rhs)
{
  Field3D result = rhs;
  int jx, jy, jz;
  BoutReal ***d;

  d = result.getData();
#ifdef CHECK
  if(d == (BoutReal***) NULL) {
    bout_error("Field3D: left / operator has invalid Field3D argument");
  }
#endif

#ifdef TRACK
  result.name = "(BoutReal/"+rhs.name+")";
#endif
  
  for(jx=0;jx<mesh->ngx;jx++)
    for(jy=0;jy<mesh->ngy;jy++)
      for(jz=0;jz<mesh->ngz;jz++)
	d[jx][jy][jz] = lhs / d[jx][jy][jz];

  result.setLocation( rhs.getLocation() );

  return(result);
}

const Field3D operator^(const BoutReal lhs, const Field3D &rhs)
{
  Field3D result = rhs;
  int jx, jy, jz;
  BoutReal ***d;

  d = result.getData();

#ifdef CHECK
  if(d == (BoutReal***) NULL) {
    output.write("Field3D: left ^ operator has invalid Field3D argument");
    exit(1);
  }
#endif

#ifdef TRACK
  result.name = "(BoutReal^"+rhs.name+")";
#endif
  
  for(jx=0;jx<mesh->ngx;jx++)
    for(jy=0;jy<mesh->ngy;jy++)
      for(jz=0;jz<mesh->ngz;jz++)
	d[jx][jy][jz] = pow(lhs, d[jx][jy][jz]);

  result.setLocation( rhs.getLocation() );

  return(result);
}

//////////////// NON-MEMBER FUNCTIONS //////////////////

const Field3D sqrt(const Field3D &f)
{
  return f.sqrt();
}

const Field3D abs(const Field3D &f)
{
  return f.abs();
}

BoutReal min(const Field3D &f, bool allpe)
{
  return f.min(allpe);
}

BoutReal max(const Field3D &f, bool allpe)
{
  return f.max(allpe);
}

/////////////////////////////////////////////////////////////////////
// Friend functions

const Field3D exp(const Field3D &f)
{
#ifdef CHECK
  msg_stack.push("exp(Field3D)");

  if(f.block == NULL)
    throw BoutException("Argument to exp(Field3D) is empty\n");
#endif

  Field3D result;
  result.allocate();
  
  for(int jx=0;jx<mesh->ngx;jx++)
    for(int jy=0;jy<mesh->ngy;jy++)
      for(int jz=0;jz<mesh->ngz;jz++)
        result.block->data[jx][jy][jz] = exp(f.block->data[jx][jy][jz]);
  
#ifdef CHECK
  msg_stack.pop();
#endif
  return result;
}

const Field3D log(const Field3D &f)
{
#ifdef CHECK
  msg_stack.push("log(Field3D)");

  if(f.block == NULL)
    throw BoutException("Argument to log(Field3D) is empty\n");
#endif

  Field3D result;
  result.allocate();
  
  for(int jx=0;jx<mesh->ngx;jx++)
    for(int jy=0;jy<mesh->ngy;jy++)
      for(int jz=0;jz<mesh->ngz;jz++) {
#ifdef CHECK
        if(f.block->data[jx][jy][jz] < 0.)
          throw BoutException("log(Field3D) has negative argument at [%d][%d][%d]\n", jx, jy, jz);
#endif
        result.block->data[jx][jy][jz] = log(f.block->data[jx][jy][jz]);
      }
  
#ifdef CHECK
  msg_stack.pop();
#endif
  return result;
}

const Field3D sin(const Field3D &f)
{
  Field3D result;
  int jx, jy, jz;
  
  result.allocate();
  
  for(jx=0;jx<mesh->ngx;jx++)
    for(jy=0;jy<mesh->ngy;jy++)
      for(jz=0;jz<mesh->ngz;jz++)
	result.block->data[jx][jy][jz] = sin(f.block->data[jx][jy][jz]);

#ifdef TRACK
  result.name = "sin("+f.name+")";
#endif

  result.location = f.location;

  return result;
}

const Field3D cos(const Field3D &f)
{
  Field3D result;
  int jx, jy, jz;
  
  result.allocate();
  
  for(jx=0;jx<mesh->ngx;jx++)
    for(jy=0;jy<mesh->ngy;jy++)
      for(jz=0;jz<mesh->ngz;jz++)
	result.block->data[jx][jy][jz] = cos(f.block->data[jx][jy][jz]);

#ifdef TRACK
  result.name = "cos("+f.name+")";
#endif

  result.location = f.location;

  return result;
}

const Field3D tan(const Field3D &f)
{
  Field3D result;
  int jx, jy, jz;
  
  result.allocate();
  
  for(jx=0;jx<mesh->ngx;jx++)
    for(jy=0;jy<mesh->ngy;jy++)
      for(jz=0;jz<mesh->ngz;jz++)
	result.block->data[jx][jy][jz] = tan(f.block->data[jx][jy][jz]);

#ifdef TRACK
  result.name = "tan("+f.name+")";
#endif

  result.location = f.location;

  return result;
}

const Field3D sinh(const Field3D &f)
{
  Field3D result;
  int jx, jy, jz;
  
  result.allocate();
  
  for(jx=0;jx<mesh->ngx;jx++)
    for(jy=0;jy<mesh->ngy;jy++)
      for(jz=0;jz<mesh->ngz;jz++)
	result.block->data[jx][jy][jz] = sinh(f.block->data[jx][jy][jz]);

#ifdef TRACK
  result.name = "sinh("+f.name+")";
#endif

  result.location = f.location;

  return result;
}

const Field3D cosh(const Field3D &f)
{
  Field3D result;
  int jx, jy, jz;
  
  result.allocate();
  
  for(jx=0;jx<mesh->ngx;jx++)
    for(jy=0;jy<mesh->ngy;jy++)
      for(jz=0;jz<mesh->ngz;jz++)
	result.block->data[jx][jy][jz] = cosh(f.block->data[jx][jy][jz]);

#ifdef TRACK
  result.name = "cosh("+f.name+")";
#endif

  result.location = f.location;

  return result;
}

const Field3D tanh(const Field3D &f)
{
  Field3D result;
  int jx, jy, jz;
  
  result.allocate();
  
  for(jx=0;jx<mesh->ngx;jx++)
    for(jy=0;jy<mesh->ngy;jy++)
      for(jz=0;jz<mesh->ngz;jz++)
	result.block->data[jx][jy][jz] = tanh(f.block->data[jx][jy][jz]);

#ifdef TRACK
  result.name = "tanh("+f.name+")";
#endif

  result.location = f.location;

  return result;
}

const Field3D filter(const Field3D &var, int N0)
{
  Field3D result;
  static dcomplex *f = (dcomplex*) NULL;
  int jx, jy, jz;
  
  int ncz = mesh->ngz-1;

  if(f == (dcomplex*) NULL) {
    //.allocate memory
    f = new dcomplex[ncz/2 + 1];
  }

  result.allocate();

  for(jx=0;jx<mesh->ngx;jx++) {
    for(jy=0;jy<mesh->ngy;jy++) {

      rfft(var.block->data[jx][jy], ncz, f); // Forward FFT

      for(jz=0;jz<=ncz/2;jz++) {
	
	if(jz != N0) {
	  // Zero this component
	  f[jz] = 0.0;
	}
      }

      irfft(f, ncz, result.block->data[jx][jy]); // Reverse FFT

      result.block->data[jx][jy][ncz] = result.block->data[jx][jy][0];
    }
  }
  
#ifdef TRACK
  result.name = "filter("+var.name+")";
#endif
  
  result.location = var.location;

  return result;
}

// Smooths a field in Fourier space
// DOESN'T WORK VERY WELL
/*
const Field3D smooth(const Field3D &var, BoutReal zmax, BoutReal xmax)
{
  Field3D result;
  static dcomplex **f = NULL, *fx;
  int jx, jy, jz, zmi, xmi;

  if(f == NULL) {
    f = cmatrix(mesh->ngx, ncz/2 + 1); 
    fx = new dcomplex[2*mesh->ngx];
  }
  
  if((zmax > 1.0) || (xmax > 1.0)) {
    // Removed everyting
    result = 0.0;
    return result;
  }
  
  result.allocate();

  zmi = ncz/2;
  xmi = mesh->ngx;

  if(zmax > 0.0)
    zmi = (int) ((1.0 - zmax)*((BoutReal) (ncz/2)));

  if(xmax > 0.0)
    xmi = (int) ((1.0 - xmax)*((BoutReal) mesh->ngx));

  //output.write("filter: %d, %d\n", xmi, zmi);

  for(jy=0;jy<mesh->ngy;jy++) {

    for(jx=0;jx<mesh->ngx;jx++) {
      // Take FFT in the Z direction, shifting into BoutReal space
      ZFFT(var.block->data[jx][jy], mesh->zShift[jx][jy], f[jx]);
    }

    if(zmax > 0.0) {
      // filter in z
      for(jx=0;jx<mesh->ngx;jx++) {
	for(jz=zmi+1;jz<=ncz/2;jz++) {
	  f[jx][jz] = 0.0;
	}
      }
    }

    if(is_pow2(mesh->ngx) && (xmax > 0.0)) {
      // mesh->ngx is a power of 2 - filter in x too
      for(jz=0;jz<=zmi;jz++) { // Go through non-zero toroidal modes
	for(jx=0;jx<mesh->ngx;jx++) {
	  fx[jx] = f[jx][jz];
	  fx[2*mesh->ngx - 1 - jx] = f[jx][jz]; // NOTE:SYMMETRIC
	}
	
	// FFT in X direction
	
	cfft(fx, 2*mesh->ngx, -1); // FFT
	
	for(jx=xmi+1; jx<=mesh->ngx; jx++) {
	  fx[jx] = 0.0;
	  fx[2*mesh->ngx-jx] = 0.0;
	}
	
	// Reverse X FFT
	cfft(fx, 2*mesh->ngx, 1);

	for(jx=0;jx<mesh->ngx;jx++)
	  f[jx][jz] = fx[jx];
	
      }
    }

    // Reverse Z FFT
    for(jx=0;jx<mesh->ngx;jx++) {
      ZFFT_rev(f[jx], mesh->zShift[jx][jy], result.block->data[jx][jy]);
      result.block->data[jx][jy][ncz] = result.block->data[jx][jy][0];
    }
  }
  
  return result;
}
*/

// Fourier filter in z
const Field3D lowPass(const Field3D &var, int zmax)
{
  Field3D result;
  static dcomplex *f = NULL;
  int jx, jy, jz;

#ifdef CHECK
  msg_stack.push("low_pass(Field3D, %d)", zmax);
#endif

  int ncz = mesh->ngz-1;
  
  if(!var.isAllocated())
    return var;

  if(f == NULL)
    f = new dcomplex[ncz/2 + 1];
 
  if((zmax >= ncz/2) || (zmax < 0)) {
    // Removing nothing
    return var;
  }
  
  result.allocate();

  for(jx=0;jx<mesh->ngx;jx++) {
    for(jy=0;jy<mesh->ngy;jy++) {
      // Take FFT in the Z direction
      rfft(var.block->data[jx][jy], ncz, f);
      
      // Filter in z
      for(jz=zmax+1;jz<=ncz/2;jz++)
	f[jz] = 0.0;

      irfft(f, ncz, result.block->data[jx][jy]); // Reverse FFT
      result.block->data[jx][jy][ncz] = result.block->data[jx][jy][0];
    }
  }
  
  result.location = var.location;

#ifdef CHECK
  msg_stack.pop();
#endif
  
  return result;
}
// Fourier filter in z with zmin
const Field3D lowPass(const Field3D &var, int zmax, int zmin)
{
  Field3D result;
  static dcomplex *f = NULL;
  int jx, jy, jz;

#ifdef CHECK
  msg_stack.push("low_pass(Field3D, %d, %d)", zmax, zmin);
#endif

  if(!var.isAllocated())
    return var;

  int ncz = mesh->ngz-1;

  if(f == NULL)
    f = new dcomplex[ncz/2 + 1];
 
  if((zmax >= ncz/2) || (zmax < 0)) {
    // Removing nothing
    return result;
  }
  
  result.allocate();

  for(jx=0;jx<mesh->ngx;jx++) {
    for(jy=0;jy<mesh->ngy;jy++) {
      // Take FFT in the Z direction
      rfft(var.block->data[jx][jy], ncz, f);
      
      // Filter in z
      for(jz=zmax+1;jz<=ncz/2;jz++)
	f[jz] = 0.0;

      // Filter zonal mode
      if(zmin==0) {
	f[0] = 0.0;
      }
      irfft(f, ncz, result.block->data[jx][jy]); // Reverse FFT
      result.block->data[jx][jy][ncz] = result.block->data[jx][jy][0];
    }
  }
  
  result.location = var.location;

#ifdef CHECK
  msg_stack.pop();
#endif
  
  return result;
}

bool finite(const Field3D &f)
{
#ifdef CHECK
  msg_stack.push("finite( Field3D )");
#endif

  if(!f.isAllocated()) {
#ifdef CHECK
    msg_stack.pop();
#endif
    return false;
  }
  
  for(int jx=0;jx<mesh->ngx;jx++)
    for(int jy=0;jy<mesh->ngy;jy++)
      for(int jz=0;jz<mesh->ngz-1;jz++)
	if(!finite(f.block->data[jx][jy][jz])) {
#ifdef CHECK
	  msg_stack.pop();
#endif
	  return false;
	}

#ifdef CHECK
  msg_stack.pop();
#endif

  return true;
}
