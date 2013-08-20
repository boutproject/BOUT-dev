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

#include <boutcomm.hxx>
#include <bout/rvec.hxx>

#include <globals.hxx> // for mesh

#include <field2d.hxx>

#include <utils.hxx>

#include <boundary_op.hxx>
#include <boundary_factory.hxx>

#include <boutexception.hxx>
#include <msg_stack.hxx>

#include <cmath>
#include <output.hxx>

Field2D::Field2D() : data(NULL), deriv(NULL) { 
  is_const = false;

  boundaryIsSet = false;

#ifdef TRACK
  name = "<F2D>";
#endif
}

Field2D::Field2D(const Field2D& f) : data(NULL), deriv(NULL) {
  is_const = false;
  boundaryIsSet = false;
  *this = f;
}

Field2D::Field2D(BoutReal val) : data(NULL), deriv(NULL) {
  boundaryIsSet = false;
  *this = val;
}

Field2D::~Field2D() {
  freeData();
  
  if(deriv != NULL)
    delete deriv;
}

Field2D* Field2D::clone() const {
  return new Field2D(*this);
}

void Field2D::allocate()
{
  allocData();
}

Field2D* Field2D::timeDeriv()
{
  if(deriv == NULL)
    deriv = new Field2D();
  return deriv;
}

BoutReal **Field2D::getData() const
{
#ifdef CHECK
  if(data == NULL)
    throw BoutException("Field2D::getData returning null pointer\n");
#endif  
  return data;
}

///////////// OPERATORS ////////////////

Field2D & Field2D::operator=(const Field2D &rhs) {
  // Check for self-assignment
  if(this == &rhs)
    return(*this); // skip this assignment

#ifdef CHECK
  msg_stack.push("Field2D: Assignment from Field2D");
  
  rhs.checkData(true);
#endif
  
#ifdef TRACK
  name = rhs.name;
#endif

  allocData(); // Make sure data is allocated

  // Copy data across
  
  #pragma omp parallel for
  for(int j=0;j<mesh->ngx*mesh->ngy;j++)
      data[0][j] = rhs.data[0][j];

#ifdef CHECK
  msg_stack.pop();
#endif

  return(*this);
}

Field2D & Field2D::operator=(const BoutReal rhs) {
#ifdef TRACK
  name = "<r2D>";
#endif

  allocData(); // Make sure data is allocated

  #pragma omp parallel for
  for(int j=0;j<mesh->ngx*mesh->ngy;j++)
    data[0][j] = rhs;
  return(*this);
}

////////////// Indexing ///////////////////

BoutReal* Field2D::operator[](int jx) const {
#ifdef CHECK

  if(data == (BoutReal**) NULL)
    throw BoutException("Field2D: [] operator on empty data");
  
  if((jx < 0) || (jx >= mesh->ngx))
    throw BoutException("Field2D: [] operator out of bounds (%d , %d)\n", jx, mesh->ngx);
#endif
  
  return(data[jx]);
}

BoutReal& Field2D::operator()(int jx, int jy) {
#if CHECK > 2
  if(data == (BoutReal**) NULL)
    throw BoutException("Field2D: [] operator on empty data");
  if((jx < 0) || (jx >= mesh->ngx) || (jy < 0) || (jy >= mesh->ngy) )
    throw BoutException("Field2D: (%d, %d) index out of bounds (%d , %d)\n", 
                        jx, jy, mesh->ngx, mesh->ngy);
#endif
  
  return data[jx][jy];
}

const BoutReal& Field2D::operator()(int jx, int jy) const {
#if CHECK > 2
  if(data == (BoutReal**) NULL)
    throw BoutException("Field2D: [] operator on empty data");
  if((jx < 0) || (jx >= mesh->ngx) || (jy < 0) || (jy >= mesh->ngy) )
    throw BoutException("Field2D: (%d, %d) index out of bounds (%d , %d)\n", 
                        jx, jy, mesh->ngx, mesh->ngy);
#endif
  
  return data[jx][jy];
}

///////// Operators

Field2D & Field2D::operator+=(const Field2D &rhs) {
#ifdef CHECK
  msg_stack.push("Field2D: += ( Field2D )");
  rhs.checkData();
  checkData();
#endif
  
#ifdef TRACK
  name = "("+name + "+" + rhs.name + ")";
#endif
  
  #pragma omp parallel for
  for(int j=0;j<mesh->ngx*mesh->ngy;j++)
    data[0][j] += rhs.data[0][j];

#ifdef CHECK
  msg_stack.pop();
#endif

  return(*this);
}

Field2D & Field2D::operator+=(const BoutReal rhs) {
#ifdef CHECK
  if(data == (BoutReal**) NULL)
    throw BoutException("Field2D: += (%e) operates on empty data", rhs);
  
#endif

#ifdef TRACK
  name = "("+name + "+BoutReal)";
#endif
  
  #pragma omp parallel for
  for(int j=0;j<mesh->ngx*mesh->ngy;j++)
    data[0][j] += rhs;

  return(*this);
}

Field2D & Field2D::operator-=(const Field2D &rhs) {
#ifdef CHECK
  if(rhs.data == (BoutReal**) NULL) {
    // Invalid data
    throw BoutException("Field2D: - operator has invalid Field2D argument");
  }

  if(data == (BoutReal**) NULL)
    throw BoutException("Field2D: -= operates on empty data");
#endif

#ifdef TRACK
  name = "("+name + "-" + rhs.name + ")";
#endif

  #pragma omp parallel for
  for(int j=0;j<mesh->ngx*mesh->ngy;j++)
    data[0][j] -= rhs.data[0][j];

  return(*this);
}

Field2D & Field2D::operator-=(const BoutReal rhs) {
#ifdef CHECK
  if(data == (BoutReal**) NULL)
    throw BoutException("Field2D: -= operates on empty data");
#endif

#ifdef TRACK
  name = "("+name + "-BoutReal)";
#endif

  #pragma omp parallel for
  for(int j=0;j<mesh->ngx*mesh->ngy;j++)
    data[0][j] -= rhs;

  return(*this);
}

Field2D & Field2D::operator*=(const Field2D &rhs) {
#ifdef CHECK
  if(rhs.data == (BoutReal**) NULL) {
    // Invalid data
    throw BoutException("Field2D: * operator has invalid Field2D argument");
  }

  if(data == (BoutReal**) NULL)
    throw BoutException("Field2D: *= operates on empty data");
#endif

#ifdef TRACK
  name = "("+name + "*" + rhs.name + ")";
#endif
  
  #pragma omp parallel for
  for(int j=0;j<mesh->ngx*mesh->ngy;j++)
    data[0][j] *= rhs.data[0][j];

  return(*this);
}

Field2D & Field2D::operator*=(const BoutReal rhs) {
#ifdef CHECK
  if(data == (BoutReal**) NULL)
    throw BoutException("Field2D: *= operates on empty data");
#endif

#ifdef TRACK
  name = "("+name + "*BoutReal)";
#endif

  #pragma omp parallel for
  for(int j=0;j<mesh->ngx*mesh->ngy;j++) {
    data[0][j] *= rhs;
  }
  
  return(*this);
}

Field2D & Field2D::operator/=(const Field2D &rhs) {
#ifdef CHECK
  if(rhs.data == (BoutReal**) NULL) {
    // Invalid data
    throw BoutException("Field2D: / operator has invalid Field2D argument");
  }

  if(data == (BoutReal**) NULL) {
    throw BoutException("Field2D: /= operates on empty data");
  }
#endif

#ifdef TRACK
  name = "("+name + "/" + rhs.name + ")";
#endif
  
  #pragma omp parallel for
  for(int j=0;j<mesh->ngx*mesh->ngy;j++)
    data[0][j] /= rhs.data[0][j];

  return(*this);
}

Field2D & Field2D::operator/=(const BoutReal rhs) {
#ifdef CHECK
  if(data == (BoutReal**) NULL) {
    throw BoutException("Field2D: /= operates on empty data");
  }
#endif

#ifdef TRACK
  name = "("+name + "/BoutReal)";
#endif

  BoutReal inv_rhs = 1. / rhs; // Multiplication faster than division
  
  #pragma omp parallel for
  for(int j=0;j<mesh->ngx*mesh->ngy;j++) {
    data[0][j] *= inv_rhs;
  }
  
  return(*this);
}

Field2D & Field2D::operator^=(const Field2D &rhs) {
#ifdef CHECK
  if(rhs.data == (BoutReal**) NULL) {
    // Invalid data
    throw BoutException("Field2D: ^ operator has invalid Field2D argument");
  }

  if(data == (BoutReal**) NULL) {
    throw BoutException("Field2D: *= operates on empty data");
  }
#endif

#ifdef TRACK
  name = "("+name + "^" + rhs.name + ")";
#endif
  
  #pragma omp parallel for
  for(int j=0;j<mesh->ngx*mesh->ngy;j++)
    data[0][j] = pow(data[0][j], rhs.data[0][j]);

  return(*this);
}

Field2D & Field2D::operator^=(const BoutReal rhs) {
#ifdef CHECK
  if(data == (BoutReal**) NULL) {
    throw BoutException("Field2D: *= operates on empty data");
  }
#endif

#ifdef TRACK
  name = "("+name + "^BoutReal)";
#endif
  
  #pragma omp parallel for
  for(int j=0;j<mesh->ngx*mesh->ngy;j++) {
    data[0][j] = pow(data[0][j], rhs);
  }
  
  return(*this);
}


///////////////// BINARY OPERATORS ////////////////////


const Field2D Field2D::operator+(const Field2D &other) const {
  Field2D result = *this;
  result += other;
  return(result);
}

const Field2D Field2D::operator+(const BoutReal rhs) const {
  Field2D result = *this;
  result += rhs;
  return result;
}

const Field2D Field2D::operator-() const {
  Field2D result = *this;

  result *= -1.0;

#ifdef TRACK
  result.name = "(-" + name + ")";
#endif

  return result;
}

const Field2D Field2D::operator-(const Field2D &other) const {
  Field2D result = *this;
  result -= other;
  return(result);
}

const Field2D Field2D::operator-(const BoutReal rhs) const {
  Field2D result = *this;
  result -= rhs;
  return result;
}

const Field2D Field2D::operator*(const Field2D &other) const {
  Field2D result = *this;
  result *= other;
  return(result);
}

const Field2D Field2D::operator*(const BoutReal rhs) const {
  Field2D result = *this;
  result *= rhs;
  return(result);
}

const Field2D Field2D::operator/(const Field2D &other) const {
  Field2D result = *this;
  result /= other;
  return(result);
}

const Field2D Field2D::operator/(const BoutReal rhs) const {
  Field2D result = *this;
  result /= rhs;
  return(result);
}

const Field2D Field2D::operator^(const Field2D &other) const {
  Field2D result = *this;
  result ^= other;
  return(result);
}

const Field2D Field2D::operator^(const BoutReal rhs) const {
  Field2D result = *this;
  result ^= rhs;
  return(result);
}

///////////// Left binary operators ////////////////

const Field3D Field2D::operator+(const Field3D &other) const {
  // just turn operator around
  return(other + (*this));
}

const Field3D Field2D::operator-(const Field3D &other) const {
  Field3D result = other;
  BoutReal ***d;

  d = result.getData();

#ifdef CHECK
  if(d == (BoutReal***) NULL)
    throw BoutException("Field2D: - operator has invalid Fiel3D argument");
  if(data == (BoutReal**) NULL)
    throw BoutException("Field2D: - operates on empty data");
#endif

#ifdef TRACK
  result.name = "(" + name + "-" + other.name + ")";
#endif  

  #pragma omp parallel for
  for(int j=0;j<mesh->ngx*mesh->ngy;j++)
    for(int jz=0;jz<mesh->ngz;jz++)
      d[0][j][jz] = data[0][j] - d[0][j][jz];

  result.setLocation( other.getLocation() );

  return(result);
}

const Field3D Field2D::operator*(const Field3D &other) const {
  // turn operator around
  return(other * (*this));
}

const Field3D Field2D::operator/(const Field3D &other) const {
  Field3D result = other;
  BoutReal ***d;

  d = result.getData();

#ifdef CHECK
  if(d == (BoutReal***) NULL)
    throw BoutException("Field2D: / operator has invalid Fiel3D argument");
  if(data == (BoutReal**) NULL)
    throw BoutException("Field2D: / operates on empty data");
#endif  

#ifdef TRACK
  result.name = "(" + name + "/" + other.name + ")";
#endif

  #pragma omp parallel for
  for(int jx=0;jx<mesh->ngx;jx++)
    for(int jy=0;jy<mesh->ngy;jy++)
      for(int jz=0;jz<mesh->ngz;jz++)
	d[jx][jy][jz] = data[jx][jy] / d[jx][jy][jz];

  result.setLocation( other.getLocation() );

  return(result);
}

const Field3D Field2D::operator^(const Field3D &other) const {
  Field3D result = other;
  BoutReal ***d;

  d = result.getData();

#ifdef CHECK
  if(d == (BoutReal***) NULL)
    throw BoutException("Field2D: ^ operator has invalid Fiel3D argument");
  if(data == (BoutReal**) NULL)
    throw BoutException("Field2D: ^ operates on empty data");
#endif

#ifdef TRACK
  result.name = "(" + name + "^" + other.name + ")";
#endif  
  
  #pragma omp parallel for
  for(int jx=0;jx<mesh->ngx;jx++)
    for(int jy=0;jy<mesh->ngy;jy++)
      for(int jz=0;jz<mesh->ngz;jz++)
	d[jx][jy][jz] = pow(data[jx][jy], d[jx][jy][jz]);

  result.setLocation( other.getLocation() );

  return(result);
}

const FieldPerp Field2D::operator+(const FieldPerp &other) const {
  FieldPerp result = other;
  result += (*this);
  return(result);
}

const FieldPerp Field2D::operator-(const FieldPerp &other) const {
  FieldPerp result = other;
  BoutReal **d = result.getData();

  int jy = result.getIndex();

  #pragma omp parallel for
  for(int jx=0;jx<mesh->ngx;jx++)
    for(int jz=0;jz<mesh->ngz;jz++)
      d[jx][jz] = data[jx][jy] - d[jx][jz];
  
#ifdef TRACK
  result.name = "(" + name + "-" + other.name + ")";
#endif

  return(result);
}

const FieldPerp Field2D::operator*(const FieldPerp &other) const {
  FieldPerp result = other;
  result *= (*this);
  return(result);
}

const FieldPerp Field2D::operator/(const FieldPerp &other) const {
  FieldPerp result = other;
  BoutReal **d = result.getData();

  int jy = result.getIndex();

  #pragma omp parallel for
  for(int jx=0;jx<mesh->ngx;jx++)
    for(int jz=0;jz<mesh->ngz;jz++)
      d[jx][jz] = data[jx][jy] / d[jx][jz];

#ifdef TRACK
  result.name = "(" + name + "/" + other.name + ")";
#endif

  return(result);
}

const FieldPerp Field2D::operator^(const FieldPerp &other) const {
  FieldPerp result = other;
  BoutReal **d = result.getData();

  int jy = result.getIndex();

  #pragma omp parallel for
  for(int jx=0;jx<mesh->ngx;jx++)
    for(int jz=0;jz<mesh->ngz;jz++)
      d[jx][jz] = pow(data[jx][jy], d[jx][jz]);

#ifdef TRACK
  result.name = "(" + name + "^" + other.name + ")";
#endif

  return(result);
}


////////////////////// STENCILS //////////////////////////

void Field2D::getXArray(int y, int z, rvec &xv) const {
#ifdef CHECK
  // Check data set
  if(data == (BoutReal**) NULL)
    throw BoutException("Field2D: getXArray on an empty data set\n");
#endif

  xv.resize(mesh->ngx);
  
  for(int x=0;x<mesh->ngx;x++)
    xv[x] = data[x][y];
}

void Field2D::getYArray(int x, int z, rvec &yv) const {
#ifdef CHECK
  // Check data set
  if(data == (BoutReal**) NULL)
    throw BoutException("Field2D: getYArray on an empty data set\n");
#endif

  yv.resize(mesh->ngy);
  
  for(int y=0;y<mesh->ngy;y++)
    yv[y] = data[x][y];
}

void Field2D::getZArray(int x, int y, rvec &zv) const {
#ifdef CHECK
  // Check data set
  if(data == (BoutReal**) NULL)
    throw BoutException("Field2D: getZArray on an empty data set\n");
#endif

  zv.resize(mesh->ngz-1);
  
  for(int z=0;z<mesh->ngz-1;z++)
    zv[z] = data[x][y];
}

void Field2D::setXArray(int y, int z, const rvec &xv) {
  allocData();

#ifdef CHECK
  // Check that vector is correct size
  if(xv.capacity() != (unsigned int) mesh->ngx)
    throw BoutException("Field2D: setXArray has incorrect size\n");
#endif

  for(int x=0;x<mesh->ngx;x++)
    data[x][y] = xv[x];
}

void Field2D::setYArray(int x, int z, const rvec &yv) {
  allocData();

#ifdef CHECK
  // Check that vector is correct size
  if(yv.capacity() != (unsigned int) mesh->ngy)
    throw BoutException("Field2D: setYArray has incorrect size\n");
#endif

  for(int y=0;y<mesh->ngy;y++)
    data[x][y] = yv[y];
}

void Field2D::setStencil(bstencil *fval, bindex *bx) const {

  // Check data set
  if(data == (BoutReal**) NULL)
    throw BoutException("Field2D: Setting stencil for empty data\n");

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

void Field2D::setXStencil(stencil &fval, const bindex &bx, CELL_LOC loc) const {
  fval.mm = data[bx.jx2m][bx.jy];
  fval.m  = data[bx.jxm][bx.jy];
  fval.c  = data[bx.jx][bx.jy];
  fval.p  = data[bx.jxp][bx.jy];
  fval.pp = data[bx.jx2p][bx.jy];
}

void Field2D::setYStencil(stencil &fval, const bindex &bx, CELL_LOC loc) const {
  fval.mm = data[bx.jx][bx.jy2m];
  fval.m  = data[bx.jx][bx.jym];
  fval.c  = data[bx.jx][bx.jy];
  fval.p  = data[bx.jx][bx.jyp];
  fval.pp = data[bx.jx][bx.jy2p];
}

void Field2D::setZStencil(stencil &fval, const bindex &bx, CELL_LOC loc) const {
  fval = data[bx.jx][bx.jy];
}

///////////////////// MATH FUNCTIONS ////////////////////


const Field2D Field2D::sqrt() const {
  // Check data set
  if(data == (BoutReal**) NULL)
    throw BoutException("Field2D: Taking sqrt of empty data\n");

#ifdef CHECK
  // Test values
  for(int jx=0;jx<mesh->ngx;jx++)
    for(int jy=0;jy<mesh->ngy;jy++)
      if(data[jx][jy] < 0.0) {
	throw BoutException("Field2D: Sqrt operates on negative value at [%d,%d]\n", jx, jy);
      }
#endif
  
  Field2D result;
  result.allocate();

  #pragma omp parallel for
  for(int j=0;j<mesh->ngx*mesh->ngy;j++)
    result.data[0][j] = ::sqrt(data[0][j]);

#ifdef TRACK
  result.name = "sqrt("+name+")";
#endif

  return result;
}

const Field2D Field2D::abs() const
{
#ifdef CHECK
  // Check data set
  if(data == (BoutReal**) NULL)
    throw BoutException("Field2D: Taking abs of empty data\n");
#endif

#ifdef TRACK
  result.name = "abs("+name+")";
#endif

  Field2D result;
  result.allocate();

  #pragma omp parallel for
  for(int j=0;j<mesh->ngx*mesh->ngy;j++)
    result.data[0][j] = fabs(data[0][j]);

  return result;
}

BoutReal Field2D::min(bool allpe) const {
#ifdef CHECK
  // Check data set
  if(data == (BoutReal**) NULL)
    throw BoutException("Field2D: Taking min of empty data\n");
  if(allpe) {
    msg_stack.push("Field2D::Min() over all PEs");
  }else
    msg_stack.push("Field2D::Min()");
#endif

  BoutReal result = data[0][0];
  
  for(int jx=0;jx<mesh->ngx;jx++)
    for(int jy=0;jy<mesh->ngy;jy++)
      if(data[jx][jy] < result)
	result = data[jx][jy];
  
  if(allpe) {
    // MPI reduce
    BoutReal localresult = result;
    MPI_Allreduce(&localresult, &result, 1, MPI_DOUBLE, MPI_MIN, BoutComm::get());
  }
  
#ifdef CHECK
  msg_stack.pop();
#endif
  
  return result;
}

BoutReal Field2D::max(bool allpe) const {
#ifdef CHECK
  // Check data set
  if(data == (BoutReal**) NULL)
    throw BoutException("Field2D: Taking max of empty data\n");
  if(allpe) {
    msg_stack.push("Field2D::Max() over all PEs");
  }else
    msg_stack.push("Field2D::Max()");
#endif

  BoutReal result = data[mesh->xstart][mesh->ystart];

  for(int jx=mesh->xstart;jx<mesh->xend;jx++)
    for(int jy=mesh->ystart;jy<mesh->yend;jy++) {
      //if(!isfinite(data[jx][jy]))
      //  output.write("Non-finite number at %d,%d\n", jx, jy);
      if(data[jx][jy] > result)
	result = data[jx][jy];
    }
  
  if(allpe) {
    // MPI reduce
    BoutReal localresult = result;
    MPI_Allreduce(&localresult, &result, 1, MPI_DOUBLE, MPI_MAX, BoutComm::get());
  }
  
#ifdef CHECK
  msg_stack.pop();
#endif
  
  return result;
}

bool Field2D::finite() const {
#ifdef CHECK
  // Check data set
  if(data == (BoutReal**) NULL)
    throw BoutException("Field2D: Taking finite of empty data\n");
#endif
  
  for(int jx=0;jx<mesh->ngx;jx++)
    for(int jy=0;jy<mesh->ngy;jy++)
      if(!::finite(data[jx][jy]))
	return false;

  return true;
}

///////////////////// FieldData VIRTUAL FUNCTIONS //////////

int Field2D::getData(int x, int y, int z, void *vptr) const {
#ifdef CHECK
  // Check data set
  if(data == (BoutReal**) NULL)
    throw BoutException("Field2D: getData on empty data\n");
  
  // check ranges
  if((x < 0) || (x >= mesh->ngx) || (y < 0) || (y >= mesh->ngy) || (z < 0) || (z >= mesh->ngx)) {
    throw BoutException("Field2D: getData (%d,%d,%d) out of bounds\n", x, y, z);
  }
#endif
  BoutReal *ptr = (BoutReal*) vptr;
  *ptr = data[x][y];
  
  return sizeof(BoutReal);
}

int Field2D::getData(int x, int y, int z, BoutReal *rptr) const {
#ifdef CHECK
  // Check data set
  if(data == (BoutReal**) NULL) {
    throw BoutException("Field2D: getData on empty data\n");
  }
  
  // check ranges
  if((x < 0) || (x >= mesh->ngx) || (y < 0) || (y >= mesh->ngy) || (z < 0) || (z >= mesh->ngz)) {
    throw BoutException("Field2D: getData (%d,%d,%d) out of bounds\n", x, y, z);
  }
#endif

  *rptr = data[x][y];
  return 1;
}

int Field2D::setData(int x, int y, int z, void *vptr) {
  allocate();
#ifdef CHECK
  // check ranges
  if((x < 0) || (x >= mesh->ngx) || (y < 0) || (y >= mesh->ngy) || (z < 0) || (z >= mesh->ngz)) {
    throw BoutException("Field2D: setData (%d,%d,%d) out of bounds\n", x, y, z);
  }
#endif
  BoutReal *ptr = (BoutReal*) vptr;
  data[x][y] = *ptr;
  
  return sizeof(BoutReal);
}

int Field2D::setData(int x, int y, int z, BoutReal *rptr) {
  allocate();
#ifdef CHECK
  // check ranges
  if((x < 0) || (x >= mesh->ngx) || (y < 0) || (y >= mesh->ngy) || (z < 0) || (z >= mesh->ngz)) {
    throw BoutException("Field2D: setData (%d,%d,%d) out of bounds\n", x, y, z);
  }
#endif

  data[x][y] = *rptr;
  return 1;
}

#ifdef CHECK
/// Check if the data is valid
bool Field2D::checkData(bool vital) const {
  if(data == (BoutReal**) NULL) {
    throw BoutException("Field2D: Operation on empty data\n");
  }

  if( vital || ( CHECK > 2 ) ) { 
    // Do full checks
    int jx, jy;

    for(jx=mesh->xstart;jx<=mesh->xend;jx++)
      for(jy=mesh->ystart;jy<=mesh->yend;jy++)
	if(!::finite(data[jx][jy])) {
	  throw BoutException("Field2D: Operation on non-finite data at [%d][%d]\n", jx, jy);
	}
  }
  return false;
}
#endif

///////////////////// BOUNDARY CONDITIONS //////////////////

void Field2D::applyBoundary() {
#ifdef CHECK
  msg_stack.push("Field2D::applyBoundary()");
  if(!boundaryIsSet)
    output << "WARNING: Call to Field2D::applyBoundary(), but no boundary set" << endl;
#endif
  for(vector<BoundaryOp*>::iterator it = bndry_op.begin(); it != bndry_op.end(); it++)
    (*it)->apply(*this);
}

void Field2D::applyBoundary(const string &condition) {
#ifdef CHECK
  msg_stack.push("Field2D::applyBoundary(condition)");
  
  if(data == NULL)
    output << "WARNING: Empty data in Field2D::applyBoundary(condition)" << endl;
#endif
  
  if(data == NULL)
    return;

  /// Get the boundary factory (singleton)
  BoundaryFactory *bfact = BoundaryFactory::getInstance();
  
  /// Get the mesh boundary regions
  vector<BoundaryRegion*> reg = mesh->getBoundaries();
  
  /// Loop over the mesh boundary regions
  for(vector<BoundaryRegion*>::iterator it=reg.begin(); it != reg.end(); it++) {
    BoundaryOp* op = bfact->create(condition, (*it));
    op->apply(*this);
    delete op;
  }
  
  // Set the corners to zero
  for(int jx=0;jx<mesh->xstart;jx++) {
    for(int jy=0;jy<mesh->ystart;jy++) {
      data[jx][jy] = 0.;
    }
    for(int jy=mesh->yend+1;jy<mesh->ngy;jy++) {
      data[jx][jy] = 0.;
    }
  }
  for(int jx=mesh->xend+1;jx<mesh->ngx;jx++) {
    for(int jy=0;jy<mesh->ystart;jy++) {
      data[jx][jy] = 0.;
    }
    for(int jy=mesh->yend+1;jy<mesh->ngy;jy++) {
      data[jx][jy] = 0.;
    }
  }
#ifdef CHECK
  msg_stack.pop();
#endif
}

void Field2D::applyBoundary(const string &region, const string &condition) {
  if(data == NULL)
    return;

  /// Get the boundary factory (singleton)
  BoundaryFactory *bfact = BoundaryFactory::getInstance();
  
  /// Get the mesh boundary regions
  vector<BoundaryRegion*> reg = mesh->getBoundaries();
  
  /// Loop over the mesh boundary regions
  for(vector<BoundaryRegion*>::iterator it=reg.begin(); it != reg.end(); it++) {
    if((*it)->label.compare(region) == 0) {
      BoundaryOp* op = bfact->create(condition, (*it));
      op->apply(*this);
      delete op;
      break;
    }
  }
  
  // Set the corners to zero
  for(int jx=0;jx<mesh->xstart;jx++) {
    for(int jy=0;jy<mesh->ystart;jy++) {
      data[jx][jy] = 0.;
    }
    for(int jy=mesh->yend+1;jy<mesh->ngy;jy++) {
      data[jx][jy] = 0.;
    }
  }
  for(int jx=mesh->xend+1;jx<mesh->ngx;jx++) {
    for(int jy=0;jy<mesh->ystart;jy++) {
      data[jx][jy] = 0.;
    }
    for(int jy=mesh->yend+1;jy<mesh->ngy;jy++) {
      data[jx][jy] = 0.;
    }
  }
}

void Field2D::applyTDerivBoundary() {
  for(vector<BoundaryOp*>::iterator it = bndry_op.begin(); it != bndry_op.end(); it++)
    (*it)->apply_ddt(*this);
}

void Field2D::setBoundaryTo(const Field2D &f2d) {
  allocate(); // Make sure data allocated
#ifdef CHECK
  msg_stack.push("Field2D::setBoundary(const Field2D&)");
  
  if(f2d.data == NULL)
    throw BoutException("Setting boundary condition to empty data\n");
#endif

  /// Get the mesh boundary regions
  vector<BoundaryRegion*> reg = mesh->getBoundaries();
  
  /// Loop over boundary regions
  for(vector<BoundaryRegion*>::iterator it = reg.begin(); it != reg.end(); it++) {
    BoundaryRegion* bndry= *it;
    /// Loop within each region
    for(bndry->first(); !bndry->isDone(); bndry->next())
      data[bndry->x][bndry->y] = f2d.data[bndry->x][bndry->y];
  }
#ifdef CHECK
  msg_stack.pop();
#endif
}

void Field2D::cleanup() {
  // Delete free blocks
  while(!block.empty()) {
    free_rmatrix(block.top());
    block.pop();
  }
  recycle = false; // Free each remaining block as objects are deleted
}

///////////////////// PRIVATE FUNCTIONS ////////////////////

// GLOBAL VARS
stack<BoutReal**> Field2D::block;
bool Field2D::recycle = true;

void Field2D::allocData() {
  if(data != (BoutReal**) NULL)
    return; // already allocated
  
  if(!block.empty()) {
    // Some free blocks
    
    data = block.top();
    block.pop();

  }else {
    // Need to create another block
    if(mesh == NULL)
      throw BoutException("Assignment to Field2D before mesh is created");
    data = rmatrix(mesh->ngx, mesh->ngy);
  }
}

void Field2D::freeData() {
  // put data block onto stack

  if(data == (BoutReal**) NULL)
    return; // No data

  if(recycle) {
    // Put block on stack
    block.push(data);
  }else {
    // Free the memory
    free_rmatrix(data); 
  }

  data = (BoutReal**) NULL;
}

////////////// NON-MEMBER OVERLOADED OPERATORS //////////////

const Field2D operator+(const BoutReal lhs, const Field2D &rhs) {
  return rhs+lhs;
}

const Field2D operator-(const BoutReal lhs, const Field2D &rhs) {
  return -1.0*(rhs - lhs);
}

const Field2D operator*(const BoutReal lhs, const Field2D &rhs) {
  // can just turn this operator around
  return(rhs * lhs);
}

const Field2D operator/(const BoutReal lhs, const Field2D &rhs) {
  Field2D result = rhs;
  BoutReal **d;

  d = result.data;

#ifdef CHECK
  if(d == (BoutReal**) NULL)
    throw BoutException("Field2D: left / operator has invalid Field2D argument");
#endif
  
#ifdef TRACK
  result.name = "(BoutReal/"+rhs.name+")";
#endif

  #pragma omp parallel for
  for(int jx=0;jx<mesh->ngx;jx++)
    for(int jy=0;jy<mesh->ngy;jy++)
      d[jx][jy] = lhs / d[jx][jy];

  return(result);
}

const Field2D operator^(const BoutReal lhs, const Field2D &rhs) {
  Field2D result = rhs;
  BoutReal **d;

  d = result.data;

#ifdef CHECK
  if(d == (BoutReal**) NULL)
    throw BoutException("Field2D: left ^ operator has invalid Field2D argument");
#endif

#ifdef TRACK
  result.name = "(BoutReal^"+rhs.name+")";
#endif
  
  #pragma omp parallel for
  for(int jx=0;jx<mesh->ngx;jx++)
    for(int jy=0;jy<mesh->ngy;jy++)
      d[jx][jy] = pow(lhs, d[jx][jy]);

  return(result);
}

//////////////// NON-MEMBER FUNCTIONS //////////////////

const Field2D SQ(const Field2D &f) {
  return f*f;
}

const Field2D sqrt(const Field2D &f) {
  return f.sqrt();
}

const Field2D abs(const Field2D &f) {
  return f.abs();
}

BoutReal min(const Field2D &f, bool allpe) {
  return f.min(allpe);
}

BoutReal max(const Field2D &f, bool allpe) {
  return f.max(allpe);
}

bool finite(const Field2D &f) {
  return f.finite();
}

/////////////////////////////////////////////////
// Friend functions

const Field2D exp(const Field2D &f) {
#ifdef CHECK
  msg_stack.push("exp(Field2D)");
#endif

  Field2D result;
  result.allocate();
  
  #pragma omp parallel for
  for(int jx=0;jx<mesh->ngx;jx++)
    for(int jy=0;jy<mesh->ngy;jy++)
      result.data[jx][jy] = ::exp(f.data[jx][jy]);
  
#ifdef CHECK
  msg_stack.pop();
#endif
  return result;
}

const Field2D log(const Field2D &f) {
#ifdef CHECK
  msg_stack.push("log(Field2D)");
#endif

  Field2D result;
  result.allocate();
  
  for(int jx=0;jx<mesh->ngx;jx++)
    for(int jy=0;jy<mesh->ngy;jy++) {
#ifdef CHECK
      if(f.data[jx][jy] < 0.)
        throw BoutException("log(Field2D) has negative argument at [%d][%d]\n", jx, jy);
#endif
      result.data[jx][jy] = ::log(f.data[jx][jy]);
    }
  
#ifdef CHECK
  msg_stack.pop();
#endif
  return result;
}

const Field2D sin(const Field2D &f) {
  Field2D result;
  
#ifdef TRACK
  result.name = "sin("+f.name+")";
#endif

  result.allocate();
  
  #pragma omp parallel for
  for(int jx=0;jx<mesh->ngx;jx++)
    for(int jy=0;jy<mesh->ngy;jy++)
      result.data[jx][jy] = ::sin(f.data[jx][jy]);

  return result;
}

const Field2D cos(const Field2D &f) {
  Field2D result;
  
#ifdef TRACK
  result.name = "cos("+f.name+")";
#endif

  result.allocate();
  
  #pragma omp parallel for
  for(int jx=0;jx<mesh->ngx;jx++)
    for(int jy=0;jy<mesh->ngy;jy++)
      result.data[jx][jy] = ::cos(f.data[jx][jy]);

  return result;
}

const Field2D tan(const Field2D &f) {
  Field2D result;
  
#ifdef TRACK
  result.name = "tan("+f.name+")";
#endif

  result.allocate();
  
  #pragma omp parallel for
  for(int jx=0;jx<mesh->ngx;jx++)
    for(int jy=0;jy<mesh->ngy;jy++)
      result.data[jx][jy] = ::tan(f.data[jx][jy]);

  return result;
}

const Field2D sinh(const Field2D &f) {
  Field2D result;
  
#ifdef TRACK
  result.name = "sinh("+f.name+")";
#endif

  result.allocate();
  
  #pragma omp parallel for
  for(int jx=0;jx<mesh->ngx;jx++)
    for(int jy=0;jy<mesh->ngy;jy++)
      result.data[jx][jy] = ::sinh(f.data[jx][jy]);

  return result;
}

const Field2D cosh(const Field2D &f) {
  Field2D result;
  
#ifdef TRACK
  result.name = "cosh("+f.name+")";
#endif

  result.allocate();
  
  #pragma omp parallel for
  for(int jx=0;jx<mesh->ngx;jx++)
    for(int jy=0;jy<mesh->ngy;jy++)
      result.data[jx][jy] = ::cosh(f.data[jx][jy]);

  return result;
}

const Field2D tanh(const Field2D &f) {
  Field2D result;
  
#ifdef TRACK
  result.name = "tanh("+f.name+")";
#endif

  result.allocate();
  
  #pragma omp parallel for
  for(int jx=0;jx<mesh->ngx;jx++)
    for(int jy=0;jy<mesh->ngy;jy++)
      result.data[jx][jy] = ::tanh(f.data[jx][jy]);

  return result;
}

const Field2D copy(const Field2D &f) {
  Field2D result = f;
  result.allocate();
  return result;
}

const Field2D floor(const Field2D &var, BoutReal f) {
  Field2D result = copy(var);

  for(int jx=0;jx<mesh->ngx;jx++)
    for(int jy=0;jy<mesh->ngy;jy++) {
      if(result(jx, jy) < f)
        result(jx, jy) = f;
    }
  return result;
}
