/**************************************************************************
 * Class for 2D X-Z slices
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

#include <globals.hxx>

#include <stdlib.h>
#include <math.h>

#include <fieldperp.hxx>
#include <utils.hxx>
#include <boutexception.hxx>

extern BoutReal** rmatrix(int nx, int ny);

FieldPerp::FieldPerp() {
  data = (BoutReal**) NULL;
  yindex = -1;
}

FieldPerp::FieldPerp(const FieldPerp &f) {
  data = (BoutReal**) NULL;
  *this = f;
}

FieldPerp::~FieldPerp() {
  freeData();
}

FieldPerp* FieldPerp::clone() const {
  return new FieldPerp(*this);
}

void FieldPerp::set(const Field3D &f, int y) {
  int jx, jz;
  BoutReal ***d = f.getData();

  if(d == (BoutReal***) NULL) {
    error("FieldPerp: Setting from empty Field3D");
    return;
  }

  yindex = y;

  allocData();

  for(jx=0;jx<mesh->ngx;jx++)
    for(jz=0;jz<mesh->ngz;jz++)
      data[jx][jz] = d[jx][y][jz];
}

void FieldPerp::allocate() {
  allocData();
}

void FieldPerp::setIndex(int y) {
  if((y < 0) || (y >= mesh->ngy) )
    throw BoutException("FieldPerp setIndex to invalid value: %d", y);
  yindex = y;
}

int FieldPerp::getIndex() const {
  if((yindex < 0) || (yindex >= mesh->ngy) )
    throw BoutException("FieldPerp has invalid yindex: %d", yindex);
  return yindex;
}

/***************************************************************
 *                         OPERATORS 
 ***************************************************************/

BoutReal* FieldPerp::operator[](int jx) const {
#if CHECK > 2
  if(data == (BoutReal**) NULL) {
    throw BoutException("FieldPerp: [] operator on empty data\n");
  }
  
  if((jx < 0) || (jx >= mesh->ngx)) {
    throw BoutException("FieldPerp: [] operator out of bounds\n");
  }
#endif
  
  return data[jx];
}

//////////////// ASSIGNMENT //////////////////

FieldPerp& FieldPerp::operator=(const FieldPerp &rhs) {
  int jx, jz;
  
  // Check for self-assignment
  if(this == &rhs)
    return(*this); // skip this assignment

  if(rhs.data == (BoutReal**) NULL) {
    // No data
    throw BoutException("FieldPerp: No data in assignment from FieldPerp");
  }

  allocData();

  yindex = rhs.getIndex();

  for(jx=0;jx<mesh->ngx;jx++)
    for(jz=0;jz<mesh->ngz;jz++)
      data[jx][jz] = rhs.data[jx][jz];

  return(*this);
}

FieldPerp & FieldPerp::operator=(const BoutReal rhs) {
  int jx, jz;

  allocData();

  for(jx=0;jx<mesh->ngx;jx++)
    for(jz=0;jz<mesh->ngz;jz++)
      data[jx][jz] = rhs;

  return(*this);
}

////////////////// ADDITION //////////////////////

FieldPerp & FieldPerp::operator+=(const FieldPerp &rhs) {
  int jx, jz;
  
  if(rhs.data == (BoutReal**) NULL) {
    // No data
    error("FieldPerp: += operates on empty FieldPerp");
    return(*this);
  }
  if(data == (BoutReal**) NULL) {
    error("FieldPerp: += operates on empty data");
    return(*this);
  }

  for(jx=0;jx<mesh->ngx;jx++)
    for(jz=0;jz<mesh->ngz;jz++)
      data[jx][jz] += rhs.data[jx][jz];

  return(*this);
}

FieldPerp & FieldPerp::operator+=(const Field3D &rhs)
{
  int jx, jz;
  BoutReal ***d;

  d = rhs.getData();
  
  if(d == (BoutReal***) NULL) {
    // No data
    error("FieldPerp: += operates on empty Field3D");
    return(*this);
  }
  
  if(data == (BoutReal**) NULL) {
    error("FieldPerp: += operates on empty data");
    return(*this);
  }

  for(jx=0;jx<mesh->ngx;jx++)
    for(jz=0;jz<mesh->ngz;jz++)
      data[jx][jz] += d[jx][yindex][jz];

  return(*this);
}

FieldPerp & FieldPerp::operator+=(const Field2D &rhs)
{
  int jx, jz;
  BoutReal **d;

  d = rhs.getData();
  
  if(d == (BoutReal**) NULL) {
    // No data
    error("FieldPerp: += operates on empty Field2D");
    return(*this);
  }
  
  if(data == (BoutReal**) NULL) {
    error("FieldPerp: += operates on empty data");
    return(*this);
  }

  for(jx=0;jx<mesh->ngx;jx++)
    for(jz=0;jz<mesh->ngz;jz++)
      data[jx][jz] += d[jx][yindex];

  return(*this);
}

FieldPerp & FieldPerp::operator+=(const BoutReal rhs)
{
  int jx, jz;

  if(data == (BoutReal**) NULL) {
    error("FieldPerp: += operates on empty data");
    return(*this);
  }
  
  for(jx=0;jx<mesh->ngx;jx++)
    for(jz=0;jz<mesh->ngz;jz++)
      data[jx][jz] += rhs;

  return(*this);
}

/////////////////// SUBTRACTION /////////////////////////

FieldPerp & FieldPerp::operator-=(const FieldPerp &rhs)
{
  int jx, jz;
  
  if(rhs.data == (BoutReal**) NULL) {
    // No data
    error("FieldPerp: -= operates on empty FieldPerp");
    return(*this);
  }
  
  if(data == (BoutReal**) NULL) {
    error("FieldPerp: -= operates on empty data");
    return(*this);
  }

  for(jx=0;jx<mesh->ngx;jx++)
    for(jz=0;jz<mesh->ngz;jz++)
      data[jx][jz] -= rhs.data[jx][jz];

  return(*this);
}

FieldPerp & FieldPerp::operator-=(const Field3D &rhs)
{
  int jx, jz;
  BoutReal ***d;

  d = rhs.getData();
  
  if(d == (BoutReal***) NULL) {
    // No data
    error("FieldPerp: -= operates on empty Field3D");
    return(*this);
  }
  
  if(data == (BoutReal**) NULL) {
    error("FieldPerp: -= operates on empty data");
    return(*this);
  }

  for(jx=0;jx<mesh->ngx;jx++)
    for(jz=0;jz<mesh->ngz;jz++)
      data[jx][jz] -= d[jx][yindex][jz];

  return(*this);
}

FieldPerp & FieldPerp::operator-=(const Field2D &rhs)
{
  int jx, jz;
  BoutReal **d;

  d = rhs.getData();
  
  if(d == (BoutReal**) NULL) {
    // No data
    error("FieldPerp: -= operates on empty Field2D");
    return(*this);
  }
  
  if(data == (BoutReal**) NULL) {
    error("FieldPerp: -= operates on empty data");
    return(*this);
  }

  for(jx=0;jx<mesh->ngx;jx++)
    for(jz=0;jz<mesh->ngz;jz++)
      data[jx][jz] -= d[jx][yindex];

  return(*this);
}

FieldPerp & FieldPerp::operator-=(const BoutReal rhs)
{
  int jx, jz;

  if(data == (BoutReal**) NULL) {
    error("FieldPerp: -= operates on empty data");
    return(*this);
  }
  
  for(jx=0;jx<mesh->ngx;jx++)
    for(jz=0;jz<mesh->ngz;jz++)
      data[jx][jz] -= rhs;

  return(*this);
}

/////////////////// MULTIPLICATION ///////////////////////

FieldPerp & FieldPerp::operator*=(const FieldPerp &rhs)
{
  int jx, jz;
  
  if(rhs.data == (BoutReal**) NULL) {
    // No data
    error("FieldPerp: *= operates on empty FieldPerp");
    return(*this);
  }
  
  if(data == (BoutReal**) NULL) {
    error("FieldPerp: *= operates on empty data");
    return(*this);
  }

  for(jx=0;jx<mesh->ngx;jx++)
    for(jz=0;jz<mesh->ngz;jz++)
      data[jx][jz] *= rhs.data[jx][jz];

  return(*this);
}

FieldPerp & FieldPerp::operator*=(const Field3D &rhs)
{
  int jx, jz;
  BoutReal ***d;

  d = rhs.getData();
  
  if(d == (BoutReal***) NULL) {
    // No data
    error("FieldPerp: *= operates on empty Field3D");
    return(*this);
  }
  
  if(data == (BoutReal**) NULL) {
    error("FieldPerp: *= operates on empty data");
    return(*this);
  }

  for(jx=0;jx<mesh->ngx;jx++)
    for(jz=0;jz<mesh->ngz;jz++)
      data[jx][jz] *= d[jx][yindex][jz];

  return(*this);
}

FieldPerp & FieldPerp::operator*=(const Field2D &rhs)
{
  int jx, jz;
  BoutReal **d;

  d = rhs.getData();
  
  if(d == (BoutReal**) NULL) {
    // No data
    error("FieldPerp: *= operates on empty Field2D");
    return(*this);
  }
  
  if(data == (BoutReal**) NULL) {
    error("FieldPerp: *= operates on empty data");
    return(*this);
  }

  for(jx=0;jx<mesh->ngx;jx++)
    for(jz=0;jz<mesh->ngz;jz++)
      data[jx][jz] *= d[jx][yindex];

  return(*this);
}

FieldPerp & FieldPerp::operator*=(const BoutReal rhs)
{
  int jx, jz;

  if(data == (BoutReal**) NULL) {
    error("FieldPerp: *= operates on empty data");
    return(*this);
  }
  
  for(jx=0;jx<mesh->ngx;jx++)
    for(jz=0;jz<mesh->ngz;jz++)
      data[jx][jz] *= rhs;

  return(*this);
}

//////////////////// DIVISION /////////////////////

FieldPerp & FieldPerp::operator/=(const FieldPerp &rhs)
{
  int jx, jz;
  
  if(rhs.data == (BoutReal**) NULL) {
    // No data
    error("FieldPerp: /= operates on empty FieldPerp");
    return(*this);
  }
  
  if(data == (BoutReal**) NULL) {
    error("FieldPerp: /= operates on empty data");
    return(*this);
  }

  for(jx=0;jx<mesh->ngx;jx++)
    for(jz=0;jz<mesh->ngz;jz++)
      data[jx][jz] /= rhs.data[jx][jz];

  return(*this);
}

FieldPerp & FieldPerp::operator/=(const Field3D &rhs)
{
  int jx, jz;
  BoutReal ***d;

  d = rhs.getData();
  
  if(d == (BoutReal***) NULL) {
    // No data
    error("FieldPerp: /= operates on empty Field3D");
    return(*this);
  }

  if(data == (BoutReal**) NULL) {
    error("FieldPerp: /= operates on empty data");
    return(*this);
  }

  for(jx=0;jx<mesh->ngx;jx++)
    for(jz=0;jz<mesh->ngz;jz++)
      data[jx][jz] /= d[jx][yindex][jz];

  return(*this);
}

FieldPerp & FieldPerp::operator/=(const Field2D &rhs)
{
  int jx, jz;
  BoutReal **d;

  d = rhs.getData();
  
  if(d == (BoutReal**) NULL) {
    // No data
    error("FieldPerp: /= operates on empty Field2D");
    return(*this);
  }

  if(data == (BoutReal**) NULL) {
    error("FieldPerp: /= operates on empty data");
    return(*this);
  }
  
  for(jx=0;jx<mesh->ngx;jx++)
    for(jz=0;jz<mesh->ngz;jz++)
      data[jx][jz] /= d[jx][yindex];

  return(*this);
}

FieldPerp & FieldPerp::operator/=(const BoutReal rhs)
{
  int jx, jz;
  
  if(data == (BoutReal**) NULL) {
    error("FieldPerp: /= operates on empty data");
    return(*this);
  }

  for(jx=0;jx<mesh->ngx;jx++)
    for(jz=0;jz<mesh->ngz;jz++)
      data[jx][jz] /= rhs;

  return(*this);
}

///////////////// EXPONENTIATION //////////////////

FieldPerp & FieldPerp::operator^=(const FieldPerp &rhs)
{
  int jx, jz;
  
  if(rhs.data == (BoutReal**) NULL) {
    // No data
    error("FieldPerp: ^= operates on empty FieldPerp");
    return(*this);
  }

  if(data == (BoutReal**) NULL) {
    error("FieldPerp: ^= operates on empty data");
    return(*this);
  }

  for(jx=0;jx<mesh->ngx;jx++)
    for(jz=0;jz<mesh->ngz;jz++)
      data[jx][jz] = pow(data[jx][jz], rhs.data[jx][jz]);

  return(*this);
}

FieldPerp & FieldPerp::operator^=(const Field3D &rhs)
{
  int jx, jz;
  BoutReal ***d;

  d = rhs.getData();
  
  if(d == (BoutReal***) NULL) {
    // No data
    error("FieldPerp: ^= operates on empty Field3D");
    return(*this);
  }

  if(data == (BoutReal**) NULL) {
    error("FieldPerp: ^= operates on empty data");
    return(*this);
  }

  for(jx=0;jx<mesh->ngx;jx++)
    for(jz=0;jz<mesh->ngz;jz++)
      data[jx][jz] = pow(data[jx][jz], d[jx][yindex][jz]);

  return(*this);
}

FieldPerp & FieldPerp::operator^=(const Field2D &rhs)
{
  int jx, jz;
  BoutReal **d;

  d = rhs.getData();
  
  if(d == (BoutReal**) NULL) {
    // No data
    error("FieldPerp: ^= operates on empty Field2D");
    return(*this);
  }

  if(data == (BoutReal**) NULL) {
    error("FieldPerp: ^= operates on empty data");
    return(*this);
  }

  for(jx=0;jx<mesh->ngx;jx++)
    for(jz=0;jz<mesh->ngz;jz++)
      data[jx][jz] = pow(data[jx][jz], d[jx][yindex]);

  return(*this);
}

FieldPerp & FieldPerp::operator^=(const BoutReal rhs)
{
  int jx, jz;
  
  if(data == (BoutReal**) NULL) {
    error("FieldPerp: ^= operates on empty data");
    return(*this);
  }

  for(jx=0;jx<mesh->ngx;jx++)
    for(jz=0;jz<mesh->ngz;jz++)
      data[jx][jz] = pow(data[jx][jz], rhs);

  return(*this);
}

/***************************************************************
 *                      BINARY OPERATORS 
 ***************************************************************/

/////////////////// ADDITION ///////////////////

const FieldPerp FieldPerp::operator+(const FieldPerp &other) const
{
  FieldPerp result = *this;
  result += other;
  return(result);
}

const FieldPerp FieldPerp::operator+(const Field3D &other) const
{
  FieldPerp result = *this;
  result += other;
  return(result);
}

const FieldPerp FieldPerp::operator+(const Field2D &other) const
{
  FieldPerp result = *this;
  result += other;
  return(result);
}

/////////////////// SUBTRACTION ////////////////

const FieldPerp FieldPerp::operator-(const FieldPerp &other) const
{
  FieldPerp result = *this;
  result -= other;
  return(result);
}

const FieldPerp FieldPerp::operator-(const Field3D &other) const
{
  FieldPerp result = *this;
  result -= other;
  return(result);
}

const FieldPerp FieldPerp::operator-(const Field2D &other) const
{
  FieldPerp result = *this;
  result -= other;
  return(result);
}

///////////////// MULTIPLICATION ///////////////

const FieldPerp FieldPerp::operator*(const FieldPerp &other) const
{
  FieldPerp result = *this;
  result *= other;
  return(result);
}

const FieldPerp FieldPerp::operator*(const Field3D &other) const
{
  FieldPerp result = *this;
  result *= other;
  return(result);
}

const FieldPerp FieldPerp::operator*(const Field2D &other) const
{
  FieldPerp result = *this;
  result *= other;
  return(result);
}

const FieldPerp FieldPerp::operator*(const BoutReal other) const
{
  FieldPerp result = *this;
  result *= other;
  return(result);
}

//////////////////// DIVISION //////////////////

const FieldPerp FieldPerp::operator/(const FieldPerp &other) const
{
  FieldPerp result = *this;
  result /= other;
  return(result);
}

const FieldPerp FieldPerp::operator/(const Field3D &other) const
{
  FieldPerp result = *this;
  result /= other;
  return(result);
}

const FieldPerp FieldPerp::operator/(const Field2D &other) const
{
  FieldPerp result = *this;
  result /= other;
  return(result);
}

const FieldPerp FieldPerp::operator/(const BoutReal other) const
{
  FieldPerp result = *this;
  result /= other;
  return(result);
}

///////////////// EXPONENTIATION ///////////////

const FieldPerp FieldPerp::operator^(const FieldPerp &other) const
{
  FieldPerp result = *this;
  result ^= other;
  return(result);
}

const FieldPerp FieldPerp::operator^(const Field3D &other) const
{
  FieldPerp result = *this;
  result ^= other;
  return(result);
}

const FieldPerp FieldPerp::operator^(const Field2D &other) const
{
  FieldPerp result = *this;
  result ^= other;
  return(result);
}

const FieldPerp FieldPerp::operator^(const BoutReal other) const
{
  FieldPerp result = *this;
  result ^= other;
  return(result);
}

////////////////////// STENCILS //////////////////////////

void FieldPerp::setStencil(bstencil *fval, bindex *bx) const
{
  fval->cc = data[bx->jx][bx->jz];

  if(mesh->ShiftXderivs && (mesh->ShiftOrder != 0)) {
    fval->xp = interpZ(bx->jxp, bx->jz, bx->xp_offset, mesh->ShiftOrder);
    fval->xm = interpZ(bx->jxm, bx->jz, bx->xm_offset, mesh->ShiftOrder);
    fval->x2p = interpZ(bx->jxp, bx->jz, bx->x2p_offset, mesh->ShiftOrder);
    fval->x2m = interpZ(bx->jxm, bx->jz, bx->x2m_offset, mesh->ShiftOrder);
  }else {
    fval->xp = data[bx->jxp][bx->jz];
    fval->xm = data[bx->jxm][bx->jz];
    fval->x2p = data[bx->jx2p][bx->jz];
    fval->x2m = data[bx->jx2m][bx->jz];
  }

  fval->yp = data[bx->jx][bx->jz];
  fval->ym = data[bx->jx][bx->jz];
  fval->zp = data[bx->jx][bx->jzp];
  fval->zm = data[bx->jx][bx->jzm];

  fval->y2p = data[bx->jx][bx->jy];
  fval->y2m = data[bx->jx][bx->jy];
  fval->z2p = data[bx->jx][bx->jz2p];
  fval->z2m = data[bx->jx][bx->jz2m];
}

void FieldPerp::setXStencil(stencil &fval, const bindex &bx, CELL_LOC loc) const
{
  if(mesh->ShiftXderivs && (mesh->ShiftOrder != 0)) {
    fval.p = interpZ(bx.jxp, bx.jz, bx.xp_offset, mesh->ShiftOrder);
    fval.m = interpZ(bx.jxm, bx.jz, bx.xm_offset, mesh->ShiftOrder);
    fval.pp = interpZ(bx.jxp, bx.jz, bx.x2p_offset, mesh->ShiftOrder);
    fval.mm = interpZ(bx.jxm, bx.jz, bx.x2m_offset, mesh->ShiftOrder);
  }else {
    fval.p = data[bx.jxp][bx.jz];
    fval.m = data[bx.jxm][bx.jz];
    fval.pp = data[bx.jx2p][bx.jz];
    fval.mm = data[bx.jx2m][bx.jz];
  }
}

void FieldPerp::setYStencil(stencil &fval, const bindex &bx, CELL_LOC loc) const
{
  fval = data[bx.jx][bx.jz];
}

void FieldPerp::setZStencil(stencil &fval, const bindex &bx, CELL_LOC loc) const
{
  fval.p = data[bx.jx][bx.jzp];
  fval.m = data[bx.jx][bx.jzm];
  fval.pp = data[bx.jx][bx.jz2p];
  fval.mm = data[bx.jx][bx.jz2m];
}

BoutReal FieldPerp::interpZ(int jx, int jz0, BoutReal zoffset, int order) const
{
  int zi;
  BoutReal result;
  int jzp, jzm, jz2p;

  zi = ROUND(zoffset);  // Find the nearest integer
  zoffset -= (BoutReal) zi; // Difference (-0.5 to +0.5)

  if(zoffset < 0.0) {
    zi--;
    zoffset += 1.0;
  }
  
  int ncz = mesh->ngz-1;
  
  jz0 = (((jz0 + zi)%ncz) + ncz) % ncz;
  jzp = (jz0 + 1) % ncz;
  jz2p = (jz0 + 2) % ncz;
  jzm = (jz0 - 1 + ncz) % ncz;

  switch(order) {
  case 2: {
    // 2-point linear interpolation

    result = (1.0 - zoffset)*data[jx][jz0] + zoffset*data[jx][jzp];

    break;
  }
  case 3: {
    // 3-point Lagrange interpolation

    result = 0.5*zoffset*(zoffset-1.0)*data[jx][jzm]
      + (1.0 - zoffset*zoffset)*data[jx][jz0]
      + 0.5*zoffset*(zoffset + 1.0)*data[jx][jzp];
    break;
  }
  case 4: {
    // 4-point Lagrange interpolation
    result = -zoffset*(zoffset-1.0)*(zoffset-2.0)*data[jx][jzm]/6.0
      + 0.5*(zoffset*zoffset - 1.0)*(zoffset-2.0)*data[jx][jz0]
      - 0.5*zoffset*(zoffset+1.0)*(zoffset-2.0)*data[jx][jzp]
      + zoffset*(zoffset*zoffset - 1.0)*data[jx][jz2p]/6.0;
    break;
  }
  default: {
    // Nearest neighbour
    result = data[jx][jz0];
  }
  };
  return result;
}

///////////////////// PRIVATE FUNCTIONS ////////////////////

// GLOBAL VARS

int FieldPerp::nblocks = 0;
int FieldPerp::max_blocks = 0;
BoutReal*** FieldPerp::block = (BoutReal***) NULL;

void FieldPerp::allocData()
{
  if(data != (BoutReal**) NULL)
    return; // already allocated

  if(nblocks > 0) {
    // Some free blocks

    nblocks--;
    data = block[nblocks];

  }else {
    // Need to create another block

    data = rmatrix(mesh->ngx, mesh->ngz);

  }
}

void FieldPerp::freeData()
{
  // put data block onto stack

  if(data == (BoutReal**) NULL)
    return; // No data

  if(nblocks == max_blocks) {
    // need to increase size of stack
    if(max_blocks == 0) {
      block = (BoutReal***) malloc(sizeof(BoutReal**));
    }else {
      block = (BoutReal***) realloc(block, sizeof(BoutReal**)*(max_blocks+1));
    }
    max_blocks++;
  }

  block[nblocks] = data;
  nblocks++;

  data = (BoutReal**) NULL;
}

////////////// NON-MEMBER OVERLOADED OPERATORS //////////////

const FieldPerp operator*(const BoutReal lhs, const FieldPerp &rhs)
{
  return(rhs*lhs);
}

const FieldPerp operator/(const BoutReal lhs, const FieldPerp &rhs)
{
  int jx,jz;
  FieldPerp result = rhs;
  BoutReal **d = result.getData();

  for(jx=0;jx<mesh->ngx;jx++)
    for(jz=0;jz<mesh->ngz;jz++)
      d[jx][jz] = lhs/d[jx][jz];

  return(result);
}

const FieldPerp operator^(const BoutReal lhs, const FieldPerp &rhs)
{
  int jx,jz;
  FieldPerp result = rhs;
  BoutReal **d = result.getData();

  for(jx=0;jx<mesh->ngx;jx++)
    for(jz=0;jz<mesh->ngz;jz++)
      d[jx][jz] = pow(lhs, d[jx][jz]);

  return(result);
}

const FieldPerp copy(const FieldPerp &f) {
  FieldPerp fcopy = f;
  fcopy.allocate();
  return fcopy;
}
