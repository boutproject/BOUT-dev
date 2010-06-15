/**************************************************************************
 * Class for 2D X-Y profiles
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

#include "mpi.h"

#include "globals.h"

#include "field2d.h"

#include "utils.h"

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

Field2D::Field2D()
{
  data = (real**) NULL;
  is_const = false;

#ifdef TRACK
  name = "<F2D>";
#endif
}

Field2D::Field2D(const Field2D& f)
{
  data = (real**) NULL;
  is_const = false;
  *this = f;
}

Field2D::Field2D(real val)
{
  data = (real**) NULL;
  *this = val;
}

Field2D::~Field2D()
{
  free_data();
}

Field2D* Field2D::clone() const
{
  return new Field2D(*this);
}

void Field2D::Allocate()
{
  alloc_data();
}

real **Field2D::getData() const
{
#ifdef CHECK
  if(data == NULL)
    error("Field2D::getData returning null pointer\n");
#endif  
  return data;
}

///////////// OPERATORS ////////////////

Field2D & Field2D::operator=(const Field2D &rhs)
{
  int jx, jy;

  // Check for self-assignment
  if(this == &rhs)
    return(*this); // skip this assignment

#ifdef CHECK
  msg_stack.push("Field2D: Assignment from Field2D");
  
  rhs.check_data(true);
#endif
  
#ifdef TRACK
  name = rhs.name;
#endif

  alloc_data(); // Make sure data is allocated

  // Copy data across

  for(jx=0;jx<mesh->ngx;jx++)
    for(jy=0;jy<mesh->ngy;jy++)
      data[jx][jy] = rhs.data[jx][jy];

#ifdef CHECK
  msg_stack.pop();
#endif

  return(*this);
}

Field2D & Field2D::operator=(const real rhs)
{
  int jx, jy;
  
#ifdef TRACK
  name = "<r2D>";
#endif

  alloc_data(); // Make sure data is allocated

  for(jx=0;jx<mesh->ngx;jx++)
    for(jy=0;jy<mesh->ngy;jy++)
      data[jx][jy] = rhs;
  return(*this);
}

real* Field2D::operator[](int jx) const
{

#ifdef CHECK

  if(data == (real**) NULL) {
    error("Field2D: [] operator on empty data");
    exit(1);
  }
  
  if((jx < 0) || (jx >= mesh->ngx)) {
    error("Field2D: [] operator out of bounds (%d , %d)\n", jx, mesh->ngx);
    exit(1);
  }
#endif
  
  return(data[jx]);
}

Field2D & Field2D::operator+=(const Field2D &rhs)
{
  int jx, jy;

#ifdef CHECK
  msg_stack.push("Field2D: += ( Field2D )");
  rhs.check_data();
  check_data();
#endif
  
#ifdef TRACK
  name = "("+name + "+" + rhs.name + ")";
#endif

  for(jx=0;jx<mesh->ngx;jx++)
    for(jy=0;jy<mesh->ngy;jy++)
      data[jx][jy] += rhs.data[jx][jy];

#ifdef CHECK
  msg_stack.pop();
#endif

  return(*this);
}

Field2D & Field2D::operator+=(const real rhs)
{
  int jx, jy;

#ifdef CHECK
  if(data == (real**) NULL) {
    error("Field2D: += operates on empty data");
  }
#endif

#ifdef TRACK
  name = "("+name + "+real)";
#endif

  for(jx=0;jx<mesh->ngx;jx++)
    for(jy=0;jy<mesh->ngy;jy++)
      data[jx][jy] += rhs;

  return(*this);
}

Field2D & Field2D::operator-=(const Field2D &rhs)
{
  int jx, jy;

#ifdef CHECK
  if(rhs.data == (real**) NULL) {
    // Invalid data
    error("Field2D: - operator has invalid Field2D argument");
  }

  if(data == (real**) NULL) {
    error("Field2D: -= operates on empty data");
  }
#endif

#ifdef TRACK
  name = "("+name + "-" + rhs.name + ")";
#endif

  for(jx=0;jx<mesh->ngx;jx++)
    for(jy=0;jy<mesh->ngy;jy++)
      data[jx][jy] -= rhs.data[jx][jy];

  return(*this);
}

Field2D & Field2D::operator-=(const real rhs)
{
  int jx, jy;

#ifdef CHECK
  if(data == (real**) NULL) {
    error("Field2D: -= operates on empty data");
  }
#endif

#ifdef TRACK
  name = "("+name + "-real)";
#endif

  for(jx=0;jx<mesh->ngx;jx++)
    for(jy=0;jy<mesh->ngy;jy++)
      data[jx][jy] -= rhs;

  return(*this);
}

Field2D & Field2D::operator*=(const Field2D &rhs)
{
  int jx, jy;

#ifdef CHECK
  if(rhs.data == (real**) NULL) {
    // Invalid data
    error("Field2D: * operator has invalid Field2D argument");
  }

  if(data == (real**) NULL) {
    error("Field2D: *= operates on empty data");
  }
#endif

#ifdef TRACK
  name = "("+name + "*" + rhs.name + ")";
#endif
  
  for(jx=0;jx<mesh->ngx;jx++)
    for(jy=0;jy<mesh->ngy;jy++)
      data[jx][jy] *= rhs.data[jx][jy];

  return(*this);
}

Field2D & Field2D::operator*=(const real rhs)
{
  int jx, jy;

#ifdef CHECK
  if(data == (real**) NULL) {
    error("Field2D: *= operates on empty data");
  }
#endif

#ifdef TRACK
  name = "("+name + "*real)";
#endif

  for(jx=0;jx<mesh->ngx;jx++) {
    for(jy=0;jy<mesh->ngy;jy++) {
      data[jx][jy] *= rhs;
    }
  }
  
  return(*this);
}

Field2D & Field2D::operator/=(const Field2D &rhs)
{
  int jx, jy;

#ifdef CHECK
  if(rhs.data == (real**) NULL) {
    // Invalid data
    error("Field2D: / operator has invalid Field2D argument");
  }

  if(data == (real**) NULL) {
    error("Field2D: /= operates on empty data");
  }
#endif

#ifdef TRACK
  name = "("+name + "/" + rhs.name + ")";
#endif
  
  for(jx=0;jx<mesh->ngx;jx++)
    for(jy=0;jy<mesh->ngy;jy++)
      data[jx][jy] /= rhs.data[jx][jy];

  return(*this);
}

Field2D & Field2D::operator/=(const real rhs)
{
  int jx, jy;

#ifdef CHECK
  if(data == (real**) NULL) {
    error("Field2D: /= operates on empty data");
  }
#endif

#ifdef TRACK
  name = "("+name + "/real)";
#endif

  real inv_rhs = 1. / rhs; // Multiplication faster than division

  for(jx=0;jx<mesh->ngx;jx++) {
    for(jy=0;jy<mesh->ngy;jy++) {
      data[jx][jy] *= inv_rhs;
    }
  }
  
  return(*this);
}

Field2D & Field2D::operator^=(const Field2D &rhs)
{
  int jx, jy;

#ifdef CHECK
  if(rhs.data == (real**) NULL) {
    // Invalid data
    error("Field2D: ^ operator has invalid Field2D argument");
  }

  if(data == (real**) NULL) {
    error("Field2D: *= operates on empty data");
  }
#endif

#ifdef TRACK
  name = "("+name + "^" + rhs.name + ")";
#endif
  
  for(jx=0;jx<mesh->ngx;jx++)
    for(jy=0;jy<mesh->ngy;jy++)
      data[jx][jy] = pow(data[jx][jy], rhs.data[jx][jy]);

  return(*this);
}

Field2D & Field2D::operator^=(const real rhs)
{
  int jx, jy;

#ifdef CHECK
  if(data == (real**) NULL) {
    error("Field2D: *= operates on empty data");
  }
#endif

#ifdef TRACK
  name = "("+name + "^real)";
#endif

  for(jx=0;jx<mesh->ngx;jx++) {
    for(jy=0;jy<mesh->ngy;jy++) {
      data[jx][jy] = pow(data[jx][jy], rhs);
    }
  }
  
  return(*this);
}


///////////////// BINARY OPERATORS ////////////////////


const Field2D Field2D::operator+(const Field2D &other) const
{
  Field2D result = *this;
  result += other;
  return(result);
}

const Field2D Field2D::operator+(const real rhs) const
{
  Field2D result = *this;
  result += rhs;
  return result;
}

const Field2D Field2D::operator-() const
{
  Field2D result = *this;

  result *= -1.0;

#ifdef TRACK
  result.name = "(-" + name + ")";
#endif

  return result;
}

const Field2D Field2D::operator-(const Field2D &other) const
{
  Field2D result = *this;
  result -= other;
  return(result);
}

const Field2D Field2D::operator-(const real rhs) const
{
  Field2D result = *this;
  result -= rhs;
  return result;
}

const Field2D Field2D::operator*(const Field2D &other) const
{
  Field2D result = *this;
  result *= other;
  return(result);
}

const Field2D Field2D::operator*(const real rhs) const
{
  Field2D result = *this;
  result *= rhs;
  return(result);
}

const Field2D Field2D::operator/(const Field2D &other) const
{
  Field2D result = *this;
  result /= other;
  return(result);
}

const Field2D Field2D::operator/(const real rhs) const
{
  Field2D result = *this;
  result /= rhs;
  return(result);
}

const Field2D Field2D::operator^(const Field2D &other) const
{
  Field2D result = *this;
  result ^= other;
  return(result);
}

const Field2D Field2D::operator^(const real rhs) const
{
  Field2D result = *this;
  result ^= rhs;
  return(result);
}

///////////// Left binary operators ////////////////

const Field3D Field2D::operator+(const Field3D &other) const
{
  // just turn operator around
  return(other + (*this));
}

const Field3D Field2D::operator-(const Field3D &other) const
{
  Field3D result = other;
  real ***d;
  int jx, jy, jz;

  d = result.getData();

#ifdef CHECK
  if(d == (real***) NULL)
    error("Field2D: - operator has invalid Fiel3D argument");
  if(data == (real**) NULL)
    error("Field2D: - operates on empty data");
#endif

#ifdef TRACK
  result.name = "(" + name + "-" + other.name + ")";
#endif  

  for(jx=0;jx<mesh->ngx;jx++)
    for(jy=0;jy<mesh->ngy;jy++)
      for(jz=0;jz<mesh->ngz;jz++)
	d[jx][jy][jz] = data[jx][jy] - d[jx][jy][jz];

  result.setLocation( other.getLocation() );

  return(result);
}

const Field3D Field2D::operator*(const Field3D &other) const
{
  // turn operator around
  return(other * (*this));
}

const Field3D Field2D::operator/(const Field3D &other) const
{
  Field3D result = other;
  real ***d;
  int jx, jy, jz;

  d = result.getData();

#ifdef CHECK
  if(d == (real***) NULL)
    error("Field2D: / operator has invalid Fiel3D argument");
  if(data == (real**) NULL)
    error("Field2D: / operates on empty data");
#endif  

#ifdef TRACK
  result.name = "(" + name + "/" + other.name + ")";
#endif

  for(jx=0;jx<mesh->ngx;jx++)
    for(jy=0;jy<mesh->ngy;jy++)
      for(jz=0;jz<mesh->ngz;jz++)
	d[jx][jy][jz] = data[jx][jy] / d[jx][jy][jz];

  result.setLocation( other.getLocation() );

  return(result);
}

const Field3D Field2D::operator^(const Field3D &other) const
{
  Field3D result = other;
  real ***d;
  int jx, jy, jz;

  d = result.getData();

#ifdef CHECK
  if(d == (real***) NULL)
    error("Field2D: ^ operator has invalid Fiel3D argument");
  if(data == (real**) NULL)
    error("Field2D: ^ operates on empty data");
#endif

#ifdef TRACK
  result.name = "(" + name + "^" + other.name + ")";
#endif  

  for(jx=0;jx<mesh->ngx;jx++)
    for(jy=0;jy<mesh->ngy;jy++)
      for(jz=0;jz<mesh->ngz;jz++)
	d[jx][jy][jz] = pow(data[jx][jy], d[jx][jy][jz]);

  result.setLocation( other.getLocation() );

  return(result);
}

const FieldPerp Field2D::operator+(const FieldPerp &other) const
{
  FieldPerp result = other;
  result += (*this);
  return(result);
}

const FieldPerp Field2D::operator-(const FieldPerp &other) const
{
  FieldPerp result = other;
  real **d = result.getData();
  int jx, jy, jz;

  jy = result.getIndex();

  for(jx=0;jx<mesh->ngx;jx++)
    for(jz=0;jz<mesh->ngz;jz++)
      d[jx][jz] = data[jx][jy] - d[jx][jz];
  
#ifdef TRACK
  result.name = "(" + name + "-" + other.name + ")";
#endif

  return(result);
}

const FieldPerp Field2D::operator*(const FieldPerp &other) const
{
  FieldPerp result = other;
  result *= (*this);
  return(result);
}

const FieldPerp Field2D::operator/(const FieldPerp &other) const
{
  FieldPerp result = other;
  real **d = result.getData();
  int jx, jy, jz;

  jy = result.getIndex();

  for(jx=0;jx<mesh->ngx;jx++)
    for(jz=0;jz<mesh->ngz;jz++)
      d[jx][jz] = data[jx][jy] / d[jx][jz];

#ifdef TRACK
  result.name = "(" + name + "/" + other.name + ")";
#endif

  return(result);
}

const FieldPerp Field2D::operator^(const FieldPerp &other) const
{
  FieldPerp result = other;
  real **d = result.getData();
  int jx, jy, jz;

  jy = result.getIndex();

  for(jx=0;jx<mesh->ngx;jx++)
    for(jz=0;jz<mesh->ngz;jz++)
      d[jx][jz] = pow(data[jx][jy], d[jx][jz]);

#ifdef TRACK
  result.name = "(" + name + "^" + other.name + ")";
#endif

  return(result);
}


////////////////////// STENCILS //////////////////////////

void Field2D::getXarray(int y, int z, rvec &xv) const
{
#ifdef CHECK
  // Check data set
  if(data == (real**) NULL) {
    error("Field2D: getXarray on an empty data set\n");
    exit(1);
  }
#endif

  xv.resize(mesh->ngx);
  
  for(int x=0;x<mesh->ngx;x++)
    xv[x] = data[x][y];
}

void Field2D::getYarray(int x, int z, rvec &yv) const
{
#ifdef CHECK
  // Check data set
  if(data == (real**) NULL) {
    error("Field2D: getYarray on an empty data set\n");
    exit(1);
  }
#endif

  yv.resize(mesh->ngy);
  
  for(int y=0;y<mesh->ngy;y++)
    yv[y] = data[x][y];
}

void Field2D::getZarray(int x, int y, rvec &zv) const
{
#ifdef CHECK
  // Check data set
  if(data == (real**) NULL) {
    error("Field2D: getZarray on an empty data set\n");
    exit(1);
  }
#endif

  zv.resize(mesh->ngz-1);
  
  for(int z=0;z<mesh->ngz-1;z++)
    zv[z] = data[x][y];
}

void Field2D::setXarray(int y, int z, const rvec &xv)
{
  alloc_data();

#ifdef CHECK
  // Check that vector is correct size
  if(xv.capacity() != (unsigned int) mesh->ngx) {
    error("Field2D: setXarray has incorrect size\n");
    exit(1);
  }
#endif

  for(int x=0;x<mesh->ngx;x++)
    data[x][y] = xv[x];
}

void Field2D::setYarray(int x, int z, const rvec &yv)
{
  alloc_data();

#ifdef CHECK
  // Check that vector is correct size
  if(yv.capacity() != (unsigned int) mesh->ngy) {
    error("Field2D: setYarray has incorrect size\n");
    exit(1);
  }
#endif

  for(int y=0;y<mesh->ngy;y++)
    data[x][y] = yv[y];
}

void Field2D::SetStencil(bstencil *fval, bindex *bx) const
{

  // Check data set
  if(data == (real**) NULL) {
    error("Field2D: Setting stencil for empty data\n");
    exit(1);
  }

  fval->jx = bx->jx;
  fval->jy = bx->jy;
  fval->jz = bx->jz;
  
  fval->cc = data[bx->jx][bx->jy];
  fval->xp = data[bx->jxp][bx->jy];
  fval->xm = data[bx->jxm][bx->jy];
  fval->yp = data[bx->jx][bx->jyp];
  fval->ym = data[bx->jx][bx->jym];
  fval->zp = data[bx->jx][bx->jy];
  fval->zm = data[bx->jx][bx->jy];

  fval->x2p = data[bx->jx2p][bx->jy];
  fval->x2m = data[bx->jx2m][bx->jy];
  fval->y2p = data[bx->jx][bx->jy2p];
  fval->y2m = data[bx->jx][bx->jy2m];
  fval->z2p = data[bx->jx][bx->jy];
  fval->z2m = data[bx->jx][bx->jy];
}

void Field2D::SetXStencil(stencil &fval, const bindex &bx, CELL_LOC loc) const
{
  fval.mm = data[bx.jx2m][bx.jy];
  fval.m  = data[bx.jxm][bx.jy];
  fval.c  = data[bx.jx][bx.jy];
  fval.p  = data[bx.jxp][bx.jy];
  fval.pp = data[bx.jx2p][bx.jy];
}

void Field2D::SetYStencil(stencil &fval, const bindex &bx, CELL_LOC loc) const
{
  fval.mm = data[bx.jx][bx.jy2m];
  fval.m  = data[bx.jx][bx.jym];
  fval.c  = data[bx.jx][bx.jy];
  fval.p  = data[bx.jx][bx.jyp];
  fval.pp = data[bx.jx][bx.jy2p];
}

void Field2D::SetZStencil(stencil &fval, const bindex &bx, CELL_LOC loc) const
{
  fval = data[bx.jx][bx.jy];
}

///////////////////// MATH FUNCTIONS ////////////////////


const Field2D Field2D::Sqrt() const
{
  int jx, jy;
  Field2D result;

  // Check data set
  if(data == (real**) NULL) {
    error("Field2D: Taking sqrt of empty data\n");
    exit(1);
  }

#ifdef CHECK
  // Test values
  for(jx=0;jx<mesh->ngx;jx++)
    for(jy=0;jy<mesh->ngy;jy++)
      if(data[jx][jy] < 0.0) {
	output.write("Field2D: Sqrt operates on negative value at [%d,%d]\n", jx, jy);
	exit(1);
      }
#endif

#ifdef TRACK
  result.name = "Sqrt("+name+")";
#endif

  result.Allocate();

  for(jx=0;jx<mesh->ngx;jx++)
    for(jy=0;jy<mesh->ngy;jy++)
      result.data[jx][jy] = sqrt(data[jx][jy]);

  return result;
}

const Field2D Field2D::Abs() const
{
  int jx, jy;
  Field2D result;

#ifdef CHECK
  // Check data set
  if(data == (real**) NULL) {
    error("Field2D: Taking abs of empty data\n");
    exit(1);
  }
#endif

#ifdef TRACK
  result.name = "Abs("+name+")";
#endif

  result.Allocate();

  for(jx=0;jx<mesh->ngx;jx++)
    for(jy=0;jy<mesh->ngy;jy++)
      result.data[jx][jy] = fabs(data[jx][jy]);

  return result;
}

real Field2D::Min(bool allpe) const
{
  int jx, jy;
  real result;

#ifdef CHECK
  // Check data set
  if(data == (real**) NULL) {
    error("Field2D: Taking min of empty data\n");
    exit(1);
  }
  if(allpe) {
    msg_stack.push("Field2D::Min() over all PEs");
  }else
    msg_stack.push("Field2D::Min()");
#endif

  result = data[0][0];

  for(jx=0;jx<mesh->ngx;jx++)
    for(jy=0;jy<mesh->ngy;jy++)
      if(data[jx][jy] < result)
	result = data[jx][jy];
  
  if(allpe) {
    // MPI reduce
    real localresult = result;
    MPI_Allreduce(&localresult, &result, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
  }
  
#ifdef CHECK
  msg_stack.pop();
#endif
  
  return result;
}

real Field2D::Max(bool allpe) const
{
  int jx, jy;
  real result;

#ifdef CHECK
  // Check data set
  if(data == (real**) NULL) {
    error("Field2D: Taking max of empty data\n");
    exit(1);
  }
  if(allpe) {
    msg_stack.push("Field2D::Max() over all PEs");
  }else
    msg_stack.push("Field2D::Max()");
#endif

  result = data[0][0];

  for(jx=0;jx<mesh->ngx;jx++)
    for(jy=0;jy<mesh->ngy;jy++)
      if(data[jx][jy] > result)
	result = data[jx][jy];
  
  if(allpe) {
    // MPI reduce
    real localresult = result;
    MPI_Allreduce(&localresult, &result, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  }
  
#ifdef CHECK
  msg_stack.pop();
#endif
  
  return result;
}

bool Field2D::Finite() const
{
  int jx, jy;

#ifdef CHECK
  // Check data set
  if(data == (real**) NULL) {
    error("Field2D: Taking finite of empty data\n");
    exit(1);
  }
#endif

  for(jx=0;jx<mesh->ngx;jx++)
    for(jy=0;jy<mesh->ngy;jy++)
      if(!finite(data[jx][jy]))
	return false;

  return true;
}

///////////////////// FieldData VIRTUAL FUNCTIONS //////////

int Field2D::getData(int x, int y, int z, void *vptr) const
{
#ifdef CHECK
  // Check data set
  if(data == (real**) NULL) {
    error("Field2D: getData on empty data\n");
    exit(1);
  }
  
  // check ranges
  if((x < 0) || (x > ncx) || (y < 0) || (y > ncy) || (z < 0) || (z >= ncz)) {
    error("Field2D: getData (%d,%d,%d) out of bounds\n", x, y, z);
    exit(1);
  }
#endif
  real *ptr = (real*) vptr;
  *ptr = data[x][y];
  
  return sizeof(real);
}

int Field2D::getData(int x, int y, int z, real *rptr) const
{
#ifdef CHECK
  // Check data set
  if(data == (real**) NULL) {
    error("Field2D: getData on empty data\n");
    exit(1);
  }
  
  // check ranges
  if((x < 0) || (x > ncx) || (y < 0) || (y > ncy) || (z < 0) || (z >= ncz)) {
    error("Field2D: getData (%d,%d,%d) out of bounds\n", x, y, z);
    exit(1);
  }
#endif

  *rptr = data[x][y];
  return 1;
}

int Field2D::setData(int x, int y, int z, void *vptr)
{
  Allocate();
#ifdef CHECK
  // check ranges
  if((x < 0) || (x > ncx) || (y < 0) || (y > ncy) || (z < 0) || (z >= ncz)) {
    error("Field2D: setData (%d,%d,%d) out of bounds\n", x, y, z);
    exit(1);
  }
#endif
  real *ptr = (real*) vptr;
  data[x][y] = *ptr;
  
  return sizeof(real);
}

int Field2D::setData(int x, int y, int z, real *rptr)
{
  Allocate();
#ifdef CHECK
  // check ranges
  if((x < 0) || (x > ncx) || (y < 0) || (y > ncy) || (z < 0) || (z >= ncz)) {
    error("Field2D: setData (%d,%d,%d) out of bounds\n", x, y, z);
    exit(1);
  }
#endif

  data[x][y] = *rptr;
  return 1;
}

#ifdef CHECK
/// Check if the data is valid
bool Field2D::check_data(bool vital) const
{
  if(data == (real**) NULL) {
    error("Field2D: Operation on empty data\n");
  }

  if( vital || ( CHECK > 2 ) ) { 
    // Do full checks
    int jx, jy;

    for(jx=MXG;jx<mesh->ngx-MXG;jx++)
      for(jy=MYG;jy<mesh->ngy-MYG;jy++)
	if(!finite(data[jx][jy])) {
	  error("Field2D: Operation on non-finite data at [%d][%d]\n", jx, jy);
	}
  }
  return false;
}
#endif

///////////////////// PRIVATE FUNCTIONS ////////////////////

// GLOBAL VARS

int Field2D::nblocks = 0;
int Field2D::max_blocks = 0;
real*** Field2D::block = (real***) NULL;

void Field2D::alloc_data()
{
  if(data != (real**) NULL)
    return; // already allocated

  if(nblocks > 0) {
    // Some free blocks

    nblocks--;
    data = block[nblocks];

  }else {
    // Need to create another block

    data = rmatrix(mesh->ngx, mesh->ngy);

  }
}

void Field2D::free_data()
{
  // put data block onto stack

  if(data == (real**) NULL)
    return; // No data

  if(nblocks == max_blocks) {
    // need to increase size of stack
    if(max_blocks == 0) {
      block = (real***) malloc(sizeof(real**));
    }else {
      block = (real***) realloc(block, sizeof(real**)*(max_blocks+1));
    }
    max_blocks++;
  }

  block[nblocks] = data;
  nblocks++;

  data = (real**) NULL;
}

////////////// NON-MEMBER OVERLOADED OPERATORS //////////////

const Field2D operator+(const real lhs, const Field2D &rhs)
{
  return rhs+lhs;
}

const Field2D operator-(const real lhs, const Field2D &rhs)
{
  return -1.0*(rhs - lhs);
}

const Field2D operator*(const real lhs, const Field2D &rhs)
{
  // can just turn this operator around
  return(rhs * lhs);
}

const Field2D operator/(const real lhs, const Field2D &rhs)
{
  Field2D result = rhs;
  int jx, jy;
  real **d;

  d = result.data;

#ifdef CHECK
  if(d == (real**) NULL) {
    output.write("Field2D: left / operator has invalid Field2D argument");
    exit(1);
  }
#endif
  
#ifdef TRACK
  result.name = "(real/"+rhs.name+")";
#endif

  for(jx=0;jx<mesh->ngx;jx++)
    for(jy=0;jy<mesh->ngy;jy++)
      d[jx][jy] = lhs / d[jx][jy];

  return(result);
}

const Field2D operator^(const real lhs, const Field2D &rhs)
{
  Field2D result = rhs;
  int jx, jy;
  real **d;

  d = result.data;

#ifdef CHECK
  if(d == (real**) NULL) {
    output.write("Field2D: left ^ operator has invalid Field2D argument");
    exit(1);
  }
#endif

#ifdef TRACK
  result.name = "(real^"+rhs.name+")";
#endif
  
  for(jx=0;jx<mesh->ngx;jx++)
    for(jy=0;jy<mesh->ngy;jy++)
      d[jx][jy] = pow(lhs, d[jx][jy]);

  return(result);
}

//////////////// NON-MEMBER FUNCTIONS //////////////////

const Field2D sqrt(const Field2D &f)
{
  return f.Sqrt();
}

const Field2D abs(const Field2D &f)
{
  return f.Abs();
}

real min(const Field2D &f, bool allpe)
{
  return f.Min(allpe);
}

real max(const Field2D &f, bool allpe)
{
  return f.Max(allpe);
}

bool finite(const Field2D &f)
{
  return f.Finite();
}

// Friend functions

const Field2D sin(const Field2D &f)
{
  Field2D result;
  int jx, jy;
  
#ifdef TRACK
  result.name = "sin("+f.name+")";
#endif

  result.Allocate();
  
  for(jx=0;jx<mesh->ngx;jx++)
    for(jy=0;jy<mesh->ngy;jy++)
      result.data[jx][jy] = sin(f.data[jx][jy]);

  return result;
}

const Field2D cos(const Field2D &f)
{
  Field2D result;
  int jx, jy;
  
#ifdef TRACK
  result.name = "cos("+f.name+")";
#endif

  result.Allocate();
  
  for(jx=0;jx<mesh->ngx;jx++)
    for(jy=0;jy<mesh->ngy;jy++)
      result.data[jx][jy] = cos(f.data[jx][jy]);

  return result;
}

const Field2D tan(const Field2D &f)
{
  Field2D result;
  int jx, jy;
  
#ifdef TRACK
  result.name = "tan("+f.name+")";
#endif

  result.Allocate();
  
  for(jx=0;jx<mesh->ngx;jx++)
    for(jy=0;jy<mesh->ngy;jy++)
      result.data[jx][jy] = tan(f.data[jx][jy]);

  return result;
}

const Field2D sinh(const Field2D &f)
{
  Field2D result;
  int jx, jy;
  
#ifdef TRACK
  result.name = "sinh("+f.name+")";
#endif

  result.Allocate();
  
  for(jx=0;jx<mesh->ngx;jx++)
    for(jy=0;jy<mesh->ngy;jy++)
      result.data[jx][jy] = sinh(f.data[jx][jy]);

  return result;
}

const Field2D cosh(const Field2D &f)
{
  Field2D result;
  int jx, jy;
  
#ifdef TRACK
  result.name = "cosh("+f.name+")";
#endif

  result.Allocate();
  
  for(jx=0;jx<mesh->ngx;jx++)
    for(jy=0;jy<mesh->ngy;jy++)
      result.data[jx][jy] = cosh(f.data[jx][jy]);

  return result;
}

const Field2D tanh(const Field2D &f)
{
  Field2D result;
  int jx, jy;
  
#ifdef TRACK
  result.name = "tanh("+f.name+")";
#endif

  result.Allocate();
  
  for(jx=0;jx<mesh->ngx;jx++)
    for(jy=0;jy<mesh->ngy;jy++)
      result.data[jx][jy] = tanh(f.data[jx][jy]);

  return result;
}
