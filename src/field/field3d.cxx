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

#include <boutcomm.hxx>
#include <globals.hxx>

#include <cmath>

#include <field3d.hxx>
#include <utils.hxx>
#include <fft.hxx>
#include <dcomplex.hxx>
#include <interpolation.hxx>
#include <boundary_op.hxx>
#include <boundary_factory.hxx>
#include <boutexception.hxx>
#include <output.hxx>
#include <msg_stack.hxx>
#include <bout/constants.hxx>
#include <bout/assert.hxx>

/// Constructor
Field3D::Field3D(Mesh *msh) : background(NULL), fieldmesh(msh), deriv(NULL), yup_field(this), ydown_field(this) {
#ifdef TRACK
  name = "<F3D>";
#endif

  if(fieldmesh) {
    nx = fieldmesh->ngx;
    ny = fieldmesh->ngy;
    nz = fieldmesh->ngz;
  }
  
  location = CELL_CENTRE; // Cell centred variable by default

  boundaryIsSet = false;
}

/// Doesn't copy any data, just create a new reference to the same data (copy on change later)
Field3D::Field3D(const Field3D& f) : background(NULL),
				     fieldmesh(f.fieldmesh), // The mesh containing array sizes
				     data(f.data),   // This handles references to the data array
				     deriv(NULL),
				     yup_field(this), ydown_field(this) {

  MsgStackItem trace("Field3D(Field3D&)");
  
#if CHECK > 2
  checkData(f);
#endif

  if(fieldmesh) {
    nx = fieldmesh->ngx;
    ny = fieldmesh->ngy;
    nz = fieldmesh->ngz;
  }

  location = f.location;
 
  boundaryIsSet = false;
}

Field3D::Field3D(const Field2D& f) : background(NULL), fieldmesh(nullptr), deriv(NULL), yup_field(this), ydown_field(this) {
  
  MsgStackItem trace("Field3D: Copy constructor from Field2D");
  
  location = CELL_CENTRE; // Cell centred variable by default
  
  boundaryIsSet = false;

  fieldmesh = mesh;
  nx = fieldmesh->ngx;
  ny = fieldmesh->ngy;
  nz = fieldmesh->ngz;
  
  *this = f;
}

Field3D::Field3D(const BoutReal val) : background(NULL), fieldmesh(nullptr), deriv(NULL), yup_field(this), ydown_field(this) {
  
  MsgStackItem trace("Field3D: Copy constructor from value");

  location = CELL_CENTRE; // Cell centred variable by default
  
  boundaryIsSet = false;

  fieldmesh = mesh;
  nx = fieldmesh->ngx;
  ny = fieldmesh->ngy;
  nz = fieldmesh->ngz;
  
  *this = val;
}

Field3D::~Field3D() {
  /// Delete the time derivative variable if allocated
  if(deriv != NULL)
    delete deriv;
  
  if(yup_field != this)
    delete yup_field;
  
  if(ydown_field != this)
    delete ydown_field;
}

void Field3D::allocate() {
  if(data.empty()) {
    if(!fieldmesh) {
      /// If no mesh, use the global
      fieldmesh = mesh;
      nx = fieldmesh->ngx;
      ny = fieldmesh->ngy;
      nz = fieldmesh->ngz;
    }
    data = Array<BoutReal>(nx*ny*nz);
  }else
    data.ensureUnique();
}

Field3D* Field3D::timeDeriv() {
  if(deriv == NULL)
    deriv = new Field3D(fieldmesh);

  return deriv;
}

void Field3D::splitYupYdown() {
  MsgStackItem trace("Field3D::splitYupYdown");
  
  if(yup_field != this)
    return;

  // yup_field and ydown_field null
  yup_field = new Field3D(fieldmesh);
  ydown_field = new Field3D(fieldmesh);
}

void Field3D::mergeYupYdown() {
  MsgStackItem trace("Field3D::mergeYupYdown");
  
  if(yup_field == this)
    return;

  delete yup_field;
  delete ydown_field;

  yup_field = this;
  ydown_field = this;
}

void Field3D::setLocation(CELL_LOC loc) {
  if(loc == CELL_VSHIFT)
    throw BoutException("Field3D: CELL_VSHIFT cell location only makes sense for vectors");
  
  if(loc == CELL_DEFAULT)
    loc = CELL_CENTRE;
  
  location = loc;
}

CELL_LOC Field3D::getLocation() const {
  return location;
}

/***************************************************************
 *                         OPERATORS 
 ***************************************************************/

const DataIterator Field3D::iterator() const {
  return DataIterator(0, nx-1, 
                      0, ny-1,
                      0, nz-1);
}

const DataIterator Field3D::begin() const {
  return DataIterator(0, 0, nx-1, 
                      0, 0, ny-1,
                      0, 0, nz-1);
}

const DataIterator Field3D::end() const {
  // end() iterator should be one past the last element
  return ++DataIterator(nx-1, 0, nx-1,
			ny-1, 0, ny-1,
			nz-1, 0, nz-1);
}

/////////////////// ASSIGNMENT ////////////////////

Field3D & Field3D::operator=(const Field3D &rhs) {
  /// Check for self-assignment
  if(this == &rhs)
    return(*this); // skip this assignment

  MsgStackItem trace("Field3D: Assignment from Field3D");
  
  /// Check that the data is valid
  checkData(rhs);
  
  // Copy the data and data sizes
  fieldmesh = rhs.fieldmesh;
  nx = rhs.nx; ny = rhs.ny; nz = rhs.nz; 
  
  data = rhs.data;
  
  location = rhs.location;
  
  return *this;
}

Field3D & Field3D::operator=(const Field2D &rhs) {
  MsgStackItem trace("Field3D = Field2D");
  
  ASSERT1(rhs.isAllocated());
  
  /// Check that the data is valid
  rhs.checkData();
 
  /// Make sure there's a unique array to copy data into
  allocate();

  /// Copy data
  for(auto i : (*this))
    (*this)[i] = rhs[i];
  
  /// Only 3D fields have locations for now
  //location = CELL_CENTRE;
  
  return *this;
}

Field3D & Field3D::operator=(const FieldPerp &rhs) {
  ASSERT1(rhs.isAllocated());
  
  /// Make sure there's a unique array to copy data into
  allocate();

  /// Copy data
  for(auto i : rhs) {
    (*this)[i] = rhs[i];
  }

  return *this;
}

const bvalue & Field3D::operator=(const bvalue &bv) {
  MsgStackItem trace("Field3D = bvalue");
  
  allocate();

#ifdef CHECK
  if(!finite(bv.val))
    throw BoutException("Field3D: assignment from non-finite value at (%d,%d,%d)\n", 
			bv.jx, bv.jy,bv.jz);
#endif

  operator()(bv.jx, bv.jy,bv.jz) = bv.val;
  
  return bv;
}

BoutReal Field3D::operator=(const BoutReal val) {
  MsgStackItem trace("Field3D = BoutReal");
  allocate();

#ifdef CHECK
  if(!finite(val))
    throw BoutException("Field3D: Assignment from non-finite BoutReal\n");
#endif
  for(auto i : (*this))
    (*this)[i] = val;

  // Only 3D fields have locations
  //location = CELL_CENTRE;
  // DON'T RE-SET LOCATION

  return val;
}

/////////////////////////////////////////////////////////////////////

#define F3D_UPDATE_FIELD(op,bop,ftype)                       \
  Field3D & Field3D::operator op(const ftype &rhs) {         \
    msg_stack.push("Field3D: %s %s", #op, #ftype);           \
    checkData(rhs) ;                                         \
    checkData(*this);                                        \
    if(data.unique()) {                                      \
      /* This is the only reference to this data */          \
      for(auto i : (*this))                                  \
        (*this)[i] op rhs[i];                                \
    }else {                                                  \
      /* Shared data */                                      \
      (*this) = (*this) bop rhs;                             \
    }                                                        \
    msg_stack.pop();                                         \
    return *this;                                            \
  }

F3D_UPDATE_FIELD(+=, +, Field3D);    // operator+= Field3D
F3D_UPDATE_FIELD(-=, -, Field3D);    // operator-= Field3D
F3D_UPDATE_FIELD(*=, *, Field3D);    // operator*= Field3D
F3D_UPDATE_FIELD(/=, /, Field3D);    // operator/= Field3D

F3D_UPDATE_FIELD(+=, +, Field2D);    // operator+= Field2D
F3D_UPDATE_FIELD(-=, -, Field2D);    // operator-= Field2D
F3D_UPDATE_FIELD(*=, *, Field2D);    // operator*= Field2D
F3D_UPDATE_FIELD(/=, /, Field2D);    // operator/= Field2D

#define F3D_UPDATE_REAL(op,bop)                              \
  Field3D & Field3D::operator op(const BoutReal &rhs) {      \
    msg_stack.push("Field3D: %s Field3D", #op);              \
    if(!finite(rhs))                                         \
      throw BoutException("Field3D: %s operator passed non-finite BoutReal number", #op); \
    checkData(*this);                                        \
                                                             \
    if(data.unique()) {                                      \
      /* This is the only reference to this data */          \
      for(auto i : (*this))                                  \
        (*this)[i] op rhs;                                   \
    }else {                                                  \
      /* Need to put result in a new block */                \
      (*this) = (*this) bop rhs;                             \
    }                                                        \
    msg_stack.pop();                                         \
    return *this;                                            \
  }

F3D_UPDATE_REAL(+=,+);    // operator+= Field2D
F3D_UPDATE_REAL(-=,-);    // operator-= Field2D
F3D_UPDATE_REAL(*=,*);    // operator*= Field2D
F3D_UPDATE_REAL(/=,/);    // operator/= Field2D

/***************************************************************
 *                         STENCILS
 ***************************************************************/

void Field3D::setXStencil(stencil &fval, const bindex &bx, CELL_LOC loc) const {
  fval.jx = bx.jx;
  fval.jy = bx.jy;
  fval.jz = bx.jz;
  
#ifdef CHECK
  // Check data set
  if(data.empty())
    throw BoutException("Field3D: Setting X stencil for empty data\n");
#endif
  
  fval.c  = operator()(bx.jx,  bx.jy, bx.jz);
  fval.p  = operator()(bx.jxp, bx.jy, bx.jz);
  fval.m  = operator()(bx.jxm, bx.jy, bx.jz);
  fval.pp = operator()(bx.jx2p, bx.jy, bx.jz);
  fval.mm = operator()(bx.jx2m, bx.jy, bx.jz);

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

void Field3D::setXStencil(forward_stencil &fval, const bindex &bx, CELL_LOC loc) const
{
  fval.jx = bx.jx;
  fval.jy = bx.jy;
  fval.jz = bx.jz;
  
#ifdef CHECK
  // Check data set
  if(data.empty())
    throw BoutException("Field3D: Setting X stencil for empty data\n");
#endif
  
  if(mesh->StaggerGrids && (loc != CELL_DEFAULT) && (loc != location)) {
    // Non-centred stencil

    if((location == CELL_CENTRE) && (loc == CELL_XLOW)) {
      // Producing a stencil centred around a lower X value
      fval.m = operator()(bx.jxm,bx.jy,bx.jz);
      fval.c = operator()(bx.jx,bx.jy,bx.jz);
      fval.p = operator()(bx.jxp,bx.jy,bx.jz);
      fval.p2 = operator()(bx.jx2p,bx.jy,bx.jz);
      fval.p3 = operator()(bx.jx+3,bx.jy,bx.jz);
      fval.p4 = operator()(bx.jx+4,bx.jy,bx.jz);
      
    }else if(location == CELL_XLOW) {
      // Stencil centred around a cell centre
      fval.m = operator()(bx.jx,bx.jy,bx.jz);
      fval.c = operator()(bx.jxp,bx.jy,bx.jz);
      fval.p = operator()(bx.jx2p,bx.jy,bx.jz);
      fval.p2 = operator()(bx.jx+3,bx.jy,bx.jz);
      fval.p3 = operator()(bx.jx+4,bx.jy,bx.jz);
      fval.p4 = operator()(bx.jx+5,bx.jy,bx.jz);
    }
    // Shifted in one direction -> shift in another
    // Could produce warning
  }
  else {
    // No shift in the z direction
    fval.m = operator()(bx.jxm,bx.jy,bx.jz);
    fval.c = operator()(bx.jx,bx.jy,bx.jz);
    fval.p = operator()(bx.jxp,bx.jy,bx.jz);
    fval.p2 = operator()(bx.jx2p,bx.jy,bx.jz);
    fval.p3 = operator()(bx.jx+3,bx.jy,bx.jz);
    fval.p4 = operator()(bx.jx+4,bx.jy,bx.jz);
  }
}

void Field3D::setXStencil(backward_stencil &fval, const bindex &bx, CELL_LOC loc) const
{
  fval.jx = bx.jx;
  fval.jy = bx.jy;
  fval.jz = bx.jz;

  ASSERT1(isAllocated());

  if(mesh->StaggerGrids && (loc != CELL_DEFAULT) && (loc != location)) {
    // Non-centred stencil

    if((location == CELL_CENTRE) && (loc == CELL_XLOW)) {
      // Producing a stencil centred around a lower X value
      fval.p = operator()(bx.jx,bx.jy,bx.jz);
      fval.c = operator()(bx.jxm,bx.jy,bx.jz);
      fval.m = operator()(bx.jx2m,bx.jy,bx.jz);
      fval.m2 = operator()(bx.jx-3,bx.jy,bx.jz);
      fval.m3 = operator()(bx.jx-4,bx.jy,bx.jz);
      fval.m4 = operator()(bx.jx-5,bx.jy,bx.jz);
      
    }else if(location == CELL_XLOW) {
      // Stencil centred around a cell centre
      fval.p = operator()(bx.jxp,bx.jy,bx.jz);
      fval.c = operator()(bx.jx,bx.jy,bx.jz);
      fval.m = operator()(bx.jxm,bx.jy,bx.jz);
      fval.m2 = operator()(bx.jx2m,bx.jy,bx.jz);
      fval.m3 = operator()(bx.jx-3,bx.jy,bx.jz);
      fval.m4 = operator()(bx.jx-4,bx.jy,bx.jz);
    }
    // Shifted in one direction -> shift in another
    // Could produce warning
  }
  else {
    // No shift in the z direction
    fval.p = operator()(bx.jxp,bx.jy,bx.jz);
    fval.c = operator()(bx.jx,bx.jy,bx.jz);
    fval.m = operator()(bx.jxm,bx.jy,bx.jz);
    fval.m2 = operator()(bx.jx2m,bx.jy,bx.jz);
    fval.m3 = operator()(bx.jx-3,bx.jy,bx.jz);
    fval.m4 = operator()(bx.jx-4,bx.jy,bx.jz);
  }
}

void Field3D::setYStencil(stencil &fval, const bindex &bx, CELL_LOC loc) const
{
  fval.jx = bx.jx;
  fval.jy = bx.jy;
  fval.jz = bx.jz;

  ASSERT0(isAllocated());
  
  fval.c = (*this)(bx.jx,bx.jy,bx.jz);
  fval.p = yup()(bx.jx,bx.jyp,bx.jz);
  fval.m = ydown()(bx.jx,bx.jym,bx.jz);
  fval.pp = 0.0;
  fval.mm = 0.0;

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

void Field3D::setYStencil(forward_stencil &fval, const bindex &bx, CELL_LOC loc) const
{
  fval.jx = bx.jx;
  fval.jy = bx.jy;
  fval.jz = bx.jz;

  ASSERT0(isAllocated());
  
  if(mesh->StaggerGrids && (loc != CELL_DEFAULT) && (loc != location)) {
    // Non-centred stencil

    if((location == CELL_CENTRE) && (loc == CELL_YLOW)) {
      // Producing a stencil centred around a lower Y value
      fval.m = operator()(bx.jx,bx.jym,bx.jz);
      fval.c = operator()(bx.jx,bx.jy,bx.jz);
      fval.p = operator()(bx.jx,bx.jyp,bx.jz);
      fval.p2 = operator()(bx.jx,bx.jy2p,bx.jz);
      fval.p3 = operator()(bx.jx,bx.jy+3,bx.jz);
      fval.p4 = operator()(bx.jx,bx.jy+4,bx.jz);
    }else if(location == CELL_YLOW) {
      // Stencil centred around a cell centre
      fval.m = operator()(bx.jx,bx.jy,bx.jz);
      fval.c = operator()(bx.jx,bx.jyp,bx.jz);
      fval.p = operator()(bx.jx,bx.jy2p,bx.jz);
      fval.p2 = operator()(bx.jx,bx.jy+3,bx.jz);
      fval.p3 = operator()(bx.jx,bx.jy+4,bx.jz);
      fval.p4 = operator()(bx.jx,bx.jy+5,bx.jz);
    }
    // Shifted in one direction -> shift in another
    // Could produce warning
  }
  else {
    fval.m = operator()(bx.jx,bx.jym,bx.jz);
    fval.c = operator()(bx.jx,bx.jy,bx.jz);
    fval.p = operator()(bx.jx,bx.jyp,bx.jz);
    fval.p2 = operator()(bx.jx,bx.jy2p,bx.jz);
    fval.p3 = operator()(bx.jx,bx.jy+3,bx.jz);
    fval.p4 = operator()(bx.jx,bx.jy+4,bx.jz);
  }
}

void Field3D::setYStencil(backward_stencil &fval, const bindex &bx, CELL_LOC loc) const {
  fval.jx = bx.jx;
  fval.jy = bx.jy;
  fval.jz = bx.jz;

  ASSERT0(isAllocated());

  if(mesh->StaggerGrids && (loc != CELL_DEFAULT) && (loc != location)) {
    // Non-centred stencil

    if((location == CELL_CENTRE) && (loc == CELL_YLOW)) {
      // Producing a stencil centred around a lower Y value
      fval.p = operator()(bx.jx,bx.jy,bx.jz);
      fval.c = operator()(bx.jx,bx.jym,bx.jz);
      fval.m = operator()(bx.jx,bx.jy2m,bx.jz);
      fval.m2 = operator()(bx.jx,bx.jy+3,bx.jz);
      fval.m3 = operator()(bx.jx,bx.jy+4,bx.jz);
      fval.m4 = operator()(bx.jx,bx.jy+5,bx.jz);
    }else if(location == CELL_YLOW) {
      // Stencil centred around a cell centre
      fval.p = operator()(bx.jx,bx.jyp,bx.jz);
      fval.c = operator()(bx.jx,bx.jy,bx.jz);
      fval.m = operator()(bx.jx,bx.jym,bx.jz);
      fval.m2 = operator()(bx.jx,bx.jy2m,bx.jz);
      fval.m3 = operator()(bx.jx,bx.jy+3,bx.jz);
      fval.m4 = operator()(bx.jx,bx.jy+4,bx.jz);
    }
    // Shifted in one direction -> shift in another
    // Could produce warning
  }
  else {
    fval.p = operator()(bx.jx,bx.jyp,bx.jz);
    fval.c = operator()(bx.jx,bx.jy,bx.jz);
    fval.m = operator()(bx.jx,bx.jym,bx.jz);
    fval.m2 = operator()(bx.jx,bx.jy2m,bx.jz);
    fval.m3 = operator()(bx.jx,bx.jy+3,bx.jz);
    fval.m4 = operator()(bx.jx,bx.jy+4,bx.jz);
  }
}

void Field3D::setZStencil(stencil &fval, const bindex &bx, CELL_LOC loc) const {
  fval.jx = bx.jx;
  fval.jy = bx.jy;
  fval.jz = bx.jz;

  ASSERT0(isAllocated());

  fval.c = operator()(bx.jx,bx.jy,bx.jz);

  fval.p = operator()(bx.jx,bx.jy,bx.jzp);
  fval.m = operator()(bx.jx,bx.jy,bx.jzm);
  fval.pp = operator()(bx.jx,bx.jy,bx.jz2p);
  fval.mm = operator()(bx.jx,bx.jy,bx.jz2m);

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

///////////////////// FieldData VIRTUAL FUNCTIONS //////////

int Field3D::getData(int x, int y, int z, void *vptr) const {

  // Check data set
  ASSERT0(isAllocated());

#if CHECK > 2
  // check ranges
  if((x < 0) || (x >= nx) || (y < 0) || (y >= ny) || (z < 0) || (z >= nz))
    throw BoutException("Field3D: getData (%d,%d,%d) out of bounds\n", x, y, z);
#endif
  
  BoutReal *ptr = (BoutReal*) vptr;
  *ptr = operator()(x,y,z);
  
  return sizeof(BoutReal);
}

int Field3D::getData(int x, int y, int z, BoutReal *rptr) const {
  ASSERT0(isAllocated());
  
#if CHECK > 2
  // check ranges
  if((x < 0) || (x >= nx) || (y < 0) || (y >= ny) || (z < 0) || (z >= nz))
    throw BoutException("Field3D: getData (%d,%d,%d) out of bounds\n", x, y, z);
#endif

  *rptr = operator()(x,y,z);
  return 1;
}

int Field3D::setData(int x, int y, int z, void *vptr) {
  allocate();
  
#if CHECK > 2
  // check ranges
  if((x < 0) || (x >= nx) || (y < 0) || (y >= ny) || (z < 0) || (z >= nz))
    throw BoutException("Field3D: setData (%d,%d,%d) out of bounds\n", x, y, z);
#endif
  BoutReal *ptr = (BoutReal*) vptr;
  operator()(x,y,z) = *ptr;
  
  return sizeof(BoutReal);
}

int Field3D::setData(int x, int y, int z, BoutReal *rptr) {
  allocate();
  
#if CHECK > 2
  // check ranges
  if((x < 0) || (x >= nx) || (y < 0) || (y >= ny) || (z < 0) || (z >= nz))
    throw BoutException("Field3D: setData (%d,%d,%d) out of bounds\n", x, y, z);
#endif

  operator()(x,y,z) = *rptr;
  return 1;
}

///////////////////// BOUNDARY CONDITIONS //////////////////

void Field3D::setBackground(const Field2D &f2d) {
  background = &f2d;
}

void Field3D::applyBoundary(bool init) {
  MsgStackItem trace("Field3D::applyBoundary()");

  ASSERT1(isAllocated());
  
  if(background != NULL) {
    // Apply boundary to the total of this and background
    
    Field3D tot = *this + (*background);
    tot.applyBoundary(init);
    *this = tot - (*background);
  } else {
    // Apply boundary to this field
    for(const auto& bndry : bndry_op)
      if ( !bndry->apply_to_ddt || init) // Always apply to the values when initialising fields, otherwise apply only if wanted
        bndry->apply(*this);
  }
  
  if (init) {
    // Set the corners to zero. ddt vanishes for the corners, so only need to be set once
    for(int jx=0;jx<fieldmesh->xstart;jx++) {
      for(int jy=0;jy<fieldmesh->ystart;jy++) {
        for(int jz=0;jz<nz;jz++)
          operator()(jx,jy,jz) = 0.;
      }
      for(int jy=fieldmesh->yend+1;jy<fieldmesh->ngy;jy++) {
        for(int jz=0;jz<nz;jz++)
          operator()(jx,jy,jz) = 0.;
      }
    }
    for(int jx=fieldmesh->xend+1;jx<fieldmesh->ngx;jx++) {
      for(int jy=0;jy<fieldmesh->ystart;jy++) {
        for(int jz=0;jz<nz;jz++)
          operator()(jx,jy,jz) = 0.;
      }
      for(int jy=fieldmesh->yend+1;jy<fieldmesh->ngy;jy++) {
        for(int jz=0;jz<nz;jz++)
          operator()(jx,jy,jz) = 0.;
      }
    }
  }
}

//JMAD
void Field3D::applyBoundary(BoutReal t) {
  MsgStackItem trace("Field3D::applyBoundary()");
  
#ifdef CHECK
  if(!boundaryIsSet)
    output << "WARNING: Call to Field3D::applyBoundary(), but no boundary set." << endl;
#endif

  ASSERT1(isAllocated())

  if(background != NULL) {
    // Apply boundary to the total of this and background

    Field3D tot = *this + (*background);
    tot.applyBoundary();
    *this = tot - (*background);
  }else {
    // Apply boundary to this field
    for(const auto& bndry : bndry_op)
      bndry->apply(*this,t);
  }

  // Set the corners to zero
  for(int jx=0;jx<fieldmesh->xstart;jx++) {
    for(int jy=0;jy<fieldmesh->ystart;jy++) {
      for(int jz=0;jz<nz;jz++)
        operator()(jx,jy,jz) = 0.;
    }
    for(int jy=fieldmesh->yend+1;jy<fieldmesh->ngy;jy++) {
      for(int jz=0;jz<nz;jz++)
        operator()(jx,jy,jz) = 0.;
    }
  }
  for(int jx=fieldmesh->xend+1;jx<fieldmesh->ngx;jx++) {
    for(int jy=0;jy<fieldmesh->ystart;jy++) {
      for(int jz=0;jz<nz;jz++)
        operator()(jx,jy,jz) = 0.;
    }
    for(int jy=fieldmesh->yend+1;jy<fieldmesh->ngy;jy++) {
      for(int jz=0;jz<nz;jz++)
        operator()(jx,jy,jz) = 0.;
    }
  }
}

void Field3D::applyBoundary(const string &condition) {
  MsgStackItem trace("Field3D::applyBoundary(condition)");
  
  ASSERT1(isAllocated());
  
  if(background != NULL) {
    // Apply boundary to the total of this and background
    
    Field3D tot = *this + (*background);
    tot.applyBoundary(condition);
    *this = tot - (*background);
    return;
  }

  /// Get the boundary factory (singleton)
  BoundaryFactory *bfact = BoundaryFactory::getInstance();
  
  /// Loop over the mesh boundary regions
  for(const auto& reg : fieldmesh->getBoundaries()) {
    BoundaryOp* op = bfact->create(condition, reg);
    op->apply(*this);
    delete op;
  }
  
  // Set the corners to zero
  for(int jx=0;jx<fieldmesh->xstart;jx++) {
    for(int jy=0;jy<fieldmesh->ystart;jy++) {
      for(int jz=0;jz<nz;jz++)
        operator()(jx,jy,jz) = 0.;
    }
    for(int jy=fieldmesh->yend+1;jy<fieldmesh->ngy;jy++) {
      for(int jz=0;jz<nz;jz++)
        operator()(jx,jy,jz) = 0.;
    }
  }
  for(int jx=fieldmesh->xend+1;jx<fieldmesh->ngx;jx++) {
    for(int jy=0;jy<fieldmesh->ystart;jy++) {
      for(int jz=0;jz<nz;jz++)
        operator()(jx,jy,jz) = 0.;
    }
    for(int jy=fieldmesh->yend+1;jy<fieldmesh->ngy;jy++) {
      for(int jz=0;jz<nz;jz++)
        operator()(jx,jy,jz) = 0.;
    }
  }
}

void Field3D::applyBoundary(const string &region, const string &condition) {
  ASSERT1(isAllocated());

  /// Get the boundary factory (singleton)
  BoundaryFactory *bfact = BoundaryFactory::getInstance();
  
  /// Loop over the mesh boundary regions
  for(const auto& reg : fieldmesh->getBoundaries()) {
    if(reg->label.compare(region) == 0) {
      BoundaryOp* op = bfact->create(condition, reg);
      op->apply(*this);
      delete op;
      break;
    }
  }
  
  // Set the corners to zero
  for(int jx=0;jx<fieldmesh->xstart;jx++) {
    for(int jy=0;jy<fieldmesh->ystart;jy++) {
      for(int jz=0;jz<nz;jz++)
        operator()(jx,jy,jz) = 0.;
    }
    for(int jy=fieldmesh->yend+1;jy<fieldmesh->ngy;jy++) {
      for(int jz=0;jz<nz;jz++)
        operator()(jx,jy,jz) = 0.;
    }
  }
  for(int jx=fieldmesh->xend+1;jx<fieldmesh->ngx;jx++) {
    for(int jy=0;jy<fieldmesh->ystart;jy++) {
      for(int jz=0;jz<nz;jz++)
        operator()(jx,jy,jz) = 0.;
    }
    for(int jy=fieldmesh->yend+1;jy<fieldmesh->ngy;jy++) {
      for(int jz=0;jz<nz;jz++)
        operator()(jx,jy,jz) = 0.;
    }
  }
}

void Field3D::applyTDerivBoundary() {
  MsgStackItem trace("Field3D::applyTDerivBoundary()");
  
  ASSERT1(isAllocated());
  ASSERT1(deriv != NULL);
  ASSERT1(deriv->isAllocated());
  
  if(background != NULL)
    *this += *background;
    
  for(const auto& bndry : bndry_op)
    bndry->apply_ddt(*this);
  
  if(background != NULL)
    *this -= *background;

  // Set the corners to zero
  for(int jx=0;jx<fieldmesh->xstart;jx++) {
    for(int jy=0;jy<fieldmesh->ystart;jy++) {
      for(int jz=0;jz<nz;jz++)
        (*deriv)(jx,jy,jz) = 0.;
    }
    for(int jy=fieldmesh->yend+1;jy<fieldmesh->ngy;jy++) {
      for(int jz=0;jz<nz;jz++)
        (*deriv)(jx,jy,jz) = 0.;
    }
  }
  for(int jx=fieldmesh->xend+1;jx<fieldmesh->ngx;jx++) {
    for(int jy=0;jy<fieldmesh->ystart;jy++) {
      for(int jz=0;jz<nz;jz++)
        (*deriv)(jx,jy,jz) = 0.;
    }
    for(int jy=fieldmesh->yend+1;jy<fieldmesh->ngy;jy++) {
      for(int jz=0;jz<nz;jz++)
        (*deriv)(jx,jy,jz) = 0.;
    }
  }
}

void Field3D::setBoundaryTo(const Field3D &f3d) {
  MsgStackItem trace("Field3D::setBoundary(const Field3D&)");
  
  allocate(); // Make sure data allocated

  ASSERT1(f3d.isAllocated());

  /// Loop over boundary regions
  for(const auto& reg : fieldmesh->getBoundaries()) {
    /// Loop within each region
    for(reg->first(); !reg->isDone(); reg->next())
      for(int z=0;z<nz;z++)
        (*this)(reg->x,reg->y,z) = f3d(reg->x,reg->y,z);
  }
}

/***************************************************************
 *               NON-MEMBER OVERLOADED OPERATORS
 ***************************************************************/


const Field3D operator-(const Field3D &f) {
  return -1.0*f;
}

#define F3D_OP_FPERP(op)                     	                          \
  const FieldPerp operator op(const Field3D &lhs, const FieldPerp &rhs) { \
    FieldPerp result;                                                     \
    result.allocate();                                                    \
    result.setIndex(rhs.getIndex());                                      \
    for(auto i : rhs)                                                     \
      result[i] = lhs[i] op rhs[i];                                       \
    return result;                                                        \
  }

F3D_OP_FPERP(+);
F3D_OP_FPERP(-);
F3D_OP_FPERP(/);
F3D_OP_FPERP(*);

#define F3D_OP_FIELD(op, ftype)                                     \
  const Field3D operator op(const Field3D &lhs, const ftype &rhs) { \
    Field3D result;                                                 \
    result.allocate();                                              \
    for(auto i : lhs)                                               \
      result[i] = lhs[i] op rhs[i];                                 \
    result.setLocation( lhs.getLocation() );                        \
    return result;                                                  \
  }

F3D_OP_FIELD(+, Field3D);   // Field3D + Field3D
F3D_OP_FIELD(-, Field3D);   // Field3D - Field3D
F3D_OP_FIELD(*, Field3D);   // Field3D * Field3D
F3D_OP_FIELD(/, Field3D);   // Field3D / Field3D

F3D_OP_FIELD(+, Field2D);   // Field3D + Field2D
F3D_OP_FIELD(-, Field2D);   // Field3D - Field2D
F3D_OP_FIELD(*, Field2D);   // Field3D * Field2D
F3D_OP_FIELD(/, Field2D);   // Field3D / Field2D

#define F3D_OP_REAL(op)                                         \
  const Field3D operator op(const Field3D &lhs, BoutReal rhs) { \
    Field3D result;                                             \
    result.allocate();                                          \
    for(auto i : lhs)                                           \
      result[i] = lhs[i] op rhs;                                \
    result.setLocation( lhs.getLocation() );                    \
    return result;                                              \
  }

F3D_OP_REAL(+); // Field3D + BoutReal
F3D_OP_REAL(-); // Field3D - BoutReal
F3D_OP_REAL(*); // Field3D * BoutReal
F3D_OP_REAL(/); // Field3D / BoutReal

#define REAL_OP_F3D(op)                                         \
  const Field3D operator op(BoutReal lhs, const Field3D &rhs) { \
    Field3D result;                                             \
    result.allocate();                                          \
    for(auto i : rhs)                                           \
      result[i] = lhs op rhs[i];                                \
    result.setLocation( rhs.getLocation() );                    \
    return result;                                              \
  }

REAL_OP_F3D(+); // BoutReal + Field3D
REAL_OP_F3D(-); // BoutReal - Field3D
REAL_OP_F3D(*); // BoutReal * Field3D
REAL_OP_F3D(/); // BoutReal / Field3D

//////////////// NON-MEMBER FUNCTIONS //////////////////

Field3D pow(const Field3D &lhs, const Field3D &rhs) {
  MsgStackItem trace("pow(Field3D, Field3D)");

  if(mesh->StaggerGrids && (lhs.getLocation() != rhs.getLocation())) {
    // Interpolate and call again
    return pow(lhs, interp_to(rhs, lhs.getLocation()));
  }
  
  Field3D result;
  result.allocate();

  // Iterate over indices
  for(auto i : result) {
    result[i] = ::pow(lhs[i], rhs[i]);
    ASSERT2( ::finite( result[i] ) );
  }
  
  result.setLocation( lhs.getLocation() );
  
  return result;
}

Field3D pow(const Field3D &lhs, const Field2D &rhs) {
  MsgStackItem trace("pow(Field3D, Field2D)");
  
  Field3D result;
  result.allocate();

  // Iterate over indices
  for(auto i : result) {
    result[i] = ::pow(lhs[i], rhs[i]);
    ASSERT2( ::finite( result[i] ) );
  }

  result.setLocation( lhs.getLocation() );
  
  return result;
}

Field3D pow(const Field3D &lhs, const FieldPerp &rhs) {
  MsgStackItem trace("pow(Field3D, FieldPerp)");
  
  Field3D result;
  result.allocate();

  // Iterate over indices
  for(auto i : result) {
    result[i] = ::pow(lhs[i], rhs[i]);
    ASSERT2( ::finite( result[i] ) );
  }

  result.setLocation( lhs.getLocation() );
  return result;
}

Field3D pow(const Field3D &f, BoutReal rhs) {
  Field3D result;
  result.allocate();
  for(auto i : result)
    result[i] = ::pow(f[i], rhs);
  
  result.setLocation( f.getLocation() );
  return result;
}

Field3D pow(BoutReal lhs, const Field3D &rhs) {
  Field3D result;
  result.allocate();
  for(auto i : result)
    result[i] = ::pow(lhs, rhs[i]);
  
  result.setLocation( rhs.getLocation() );
  return result;
}

BoutReal min(const Field3D &f, bool allpe) {
#ifdef CHECK
  if(!f.isAllocated())
    throw BoutException("Field3D: min() method on empty data");

  if(allpe) {
    msg_stack.push("Field3D::Min() over all PEs");
  }else
    msg_stack.push("Field3D::Min()");
#endif

  BoutReal result = f(mesh->xstart,mesh->ystart,0);
  
  for(int i=mesh->xstart; i<=mesh->xend; i++)
    for(int j=mesh->ystart; j<=mesh->yend; j++)
      for(int k=0;k<mesh->ngz-1;k++)
        if(f(i,j,k) < result)
          result = f(i,j,k);
  
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

BoutReal max(const Field3D &f, bool allpe) {
#ifdef CHECK
  if(!f.isAllocated())
    throw BoutException("Field3D: max() method on empty data");
  if(allpe) {
    msg_stack.push("Field3D::Max() over all PEs");
  }else
    msg_stack.push("Field3D::Max()");
#endif
  
  BoutReal result = f(mesh->xstart,mesh->ystart,0);
  
  for(int i=mesh->xstart; i<=mesh->xend; i++)
    for(int j=mesh->ystart; j<=mesh->yend; j++)
      for(int k=0;k<mesh->ngz-1;k++)
        if(f(i,j,k) > result)
          result = f(i,j,k);
  
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

/////////////////////////////////////////////////////////////////////
// Friend functions

#define F3D_FUNC(name, func)                               \
  const Field3D name(const Field3D &f) {                   \
    msg_stack.push(#name "(Field3D)");                     \
    /* Check if the input is allocated */                  \
    ASSERT1(f.isAllocated());                              \
    /* Define and allocate the output result */            \
    Field3D result;                                        \
    result.allocate();                                     \
    /* Loop over domain */                                 \
    for(auto d : result) {                                 \
      result[d] = func(f[d]);                              \
      /* If checking is set to 3 or higher, test result */ \
      ASSERT3(finite(result[d]));                          \
    }                                                      \
    result.setLocation(f.getLocation());                   \
    msg_stack.pop();                                       \
    return result;                                         \
  }

const Field3D SQ(const Field3D &f) {
  return f * f;
}

F3D_FUNC(sqrt, ::sqrt);
F3D_FUNC(abs, ::fabs);

F3D_FUNC(exp, ::exp);
F3D_FUNC(log, ::log);

F3D_FUNC(sin, ::sin);
F3D_FUNC(cos, ::cos);
F3D_FUNC(tan, ::tan);

F3D_FUNC(sinh, ::sinh);
F3D_FUNC(cosh, ::cosh);
F3D_FUNC(tanh, ::tanh);

const Field3D filter(const Field3D &var, int N0) {
  MsgStackItem trace("filter(Field3D, int)");
  
  ASSERT1(var.isAllocated());
  
  int ncz = mesh->ngz-1;
  Array<dcomplex> f(ncz/2 + 1);
  
  Field3D result;
  result.allocate();
  
  for(int jx=0;jx<mesh->ngx;jx++) {
    for(int jy=0;jy<mesh->ngy;jy++) {

      rfft(&(var(jx, jy, 0)), ncz, f.begin()); // Forward FFT

      for(int jz=0;jz<=ncz/2;jz++) {
	
	if(jz != N0) {
	  // Zero this component
	  f[jz] = 0.0;
	}
      }

      irfft(f.begin(), ncz, &(result(jx, jy, 0))); // Reverse FFT

      result(jx,jy,ncz) = result(jx,jy,0);
    }
  }
  
#ifdef TRACK
  result.name = "filter("+var.name+")";
#endif
  
  result.setLocation(var.getLocation());

  return result;
}

// Fourier filter in z
const Field3D lowPass(const Field3D &var, int zmax) {
  
  msg_stack.push("lowPass(Field3D, %d)", zmax);

  ASSERT1(var.isAllocated());
  
  int ncz = mesh->ngz-1;
  
  // Create an array 
  Array<dcomplex> f(ncz/2 + 1);
  
  if((zmax >= ncz/2) || (zmax < 0)) {
    // Removing nothing
    return var;
  }

  Field3D result;
  result.allocate();
  
  for(int jx=0;jx<mesh->ngx;jx++) {
    for(int jy=0;jy<mesh->ngy;jy++) {
      // Take FFT in the Z direction
      rfft(&(var(jx,jy,0)), ncz, f.begin());
      
      // Filter in z
      for(int jz=zmax+1;jz<=ncz/2;jz++)
	f[jz] = 0.0;

      irfft(f.begin(), ncz, &(result(jx,jy,0))); // Reverse FFT
      result(jx,jy,ncz) = result(jx,jy,0);
    }
  }
  
  result.setLocation(var.getLocation());

  msg_stack.pop();
  
  return result;
}

// Fourier filter in z with zmin
const Field3D lowPass(const Field3D &var, int zmax, int zmin) {

#ifdef CHECK
  msg_stack.push("lowPass(Field3D, %d, %d)", zmax, zmin);
#endif

  ASSERT1(var.isAllocated());

  int ncz = mesh->ngz-1;
  Array<dcomplex> f(ncz/2 + 1);
 
  if(((zmax >= ncz/2) || (zmax < 0)) && (zmin < 0)) {
    // Removing nothing
    return var;
  }

  Field3D result;
  result.allocate();
  
  for(int jx=0;jx<mesh->ngx;jx++) {
    for(int jy=0;jy<mesh->ngy;jy++) {
      // Take FFT in the Z direction
      rfft(&(var(jx,jy,0)), ncz, f.begin());
      
      // Filter in z
      for(int jz=zmax+1;jz<=ncz/2;jz++)
	f[jz] = 0.0;

      // Filter zonal mode
      if(zmin==0) {
	f[0] = 0.0;
      }
      irfft(f.begin(), ncz, &(result(jx,jy,0))); // Reverse FFT
      result(jx,jy,ncz) = result(jx,jy,0);
    }
  }
  
  result.setLocation(var.getLocation());
  
#ifdef CHECK
  msg_stack.pop();
#endif
  
  return result;
}

bool finite(const Field3D &f) {
  MsgStackItem trace("finite( Field3D )");
  
  if(!f.isAllocated()) {
    return false;
  }
  
  for(auto d : f)
    if(!finite(f[d]))
      return false;
  
  return true;
}

#ifdef CHECK
/// Check if the data is valid
void checkData(const Field3D &f)  {
  if(!f.isAllocated())
    throw BoutException("Field3D: Operation on empty data\n");
  
  for(auto d : f) {
    if( (d.x < mesh->xstart) or (d.x > mesh->xend) or (d.y < mesh->ystart) or (d.y > mesh->yend) or (d.z >= mesh->ngz-1))
      continue; // Exclude boundary cells
    
    if(!finite(f[d]))
      throw BoutException("Field3D: Operation on non-finite data at [%d][%d][%d]\n", d.x, d.y, d.z);
  }
}
#endif

const Field3D copy(const Field3D &f) {
  Field3D result = f;
  result.allocate();
  return result;
}

const Field3D floor(const Field3D &var, BoutReal f) {
  Field3D result = copy(var);
  
  for(auto d : result)
    if(result[d] < f)
      result[d] = f;
  
  return result;
}

Field2D DC(const Field3D &f) {
  MsgStackItem trace("DC(Field3D)");
  
  Field2D result;
  result.allocate();

  for(int i=0;i<mesh->ngx;i++)
    for(int j=0;j<mesh->ngy;j++) {
      result(i,j) = 0.0;
      for(int k=0;k<mesh->ngz-1;k++)
	result(i,j) += f(i,j,k);
      result(i,j) /= (mesh->ngz-1);
    }
  
  return result;
}

