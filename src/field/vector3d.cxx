/**************************************************************************
 * Class for 3D vectors. Built on the Field3D class,
 * all operators relating to vectors are here (none in Field classes)
 *
 * NOTE: COMPONENTS STORED ARE COVARIANT
 *
 * B.Dudson, October 2007
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

#include <vector3d.hxx>
#include <boundary_op.hxx>
#include <boutexception.hxx>

Vector3D::Vector3D() : covariant(true), deriv(NULL) { }

Vector3D::Vector3D(const Vector3D &f) : covariant(f.covariant), deriv(NULL) {
  x = f.x;
  y = f.y;
  z = f.z;
}

Vector3D::~Vector3D() {
  if(deriv != NULL) {
    // The ddt of the components (x.ddt) point to the same place as ddt.x
    // only delete once
    x.deriv = NULL;
    y.deriv = NULL;
    z.deriv = NULL;
    
    // Now delete them as part of the deriv vector
    delete deriv;
  }
}

void Vector3D::toCovariant() {  
  if(!covariant) {
    Field3D gx, gy, gz;

    // multiply by g_{ij}
    gx = x*mesh->g_11 + mesh->g_12*y + mesh->g_13*z;
    gy = y*mesh->g_22 + mesh->g_12*x + mesh->g_23*z;
    gz = z*mesh->g_33 + mesh->g_13*x + mesh->g_23*y;

    x = gx;
    y = gy;
    z = gz;
    
    covariant = true;
  }
}
void Vector3D::toContravariant()
{  
  if(covariant) {
    // multiply by g^{ij}
    
    Field3D gx, gy, gz;

    gx = x*mesh->g11 + mesh->g12*y + mesh->g13*z;
    gy = y*mesh->g22 + mesh->g12*x + mesh->g23*z;
    gz = z*mesh->g33 + mesh->g13*x + mesh->g23*y;

    x = gx;
    y = gy;
    z = gz;
    
    covariant = false;
  }
}

Vector3D* Vector3D::timeDeriv()
{
  if(deriv == NULL) {
    deriv = new Vector3D();
    
    // Check if the components have a time-derivative
    // Need to make sure that ddt(v.x) = ddt(v).x
    
    if(x.deriv != NULL) {
      // already set. Copy across then delete
      deriv->x = *(x.deriv);
      delete x.deriv;
    }
    if(y.deriv != NULL) {
      deriv->y = *(y.deriv);
      delete y.deriv;
    }
    if(z.deriv != NULL) {
      deriv->z = *(z.deriv);
      delete z.deriv;
    }
    // Set the component time-derivatives
    x.deriv = &(deriv->x);
    y.deriv = &(deriv->y);
    z.deriv = &(deriv->z);
  }
  return deriv;
}

/***************************************************************
 *                         OPERATORS 
 ***************************************************************/

/////////////////// ASSIGNMENT ////////////////////

Vector3D & Vector3D::operator=(const Vector3D &rhs)
{
  x = rhs.x;
  y = rhs.y;
  z = rhs.z;

  covariant = rhs.covariant;

  return *this;
}

Vector3D & Vector3D::operator=(const Vector2D &rhs)
{
  x = rhs.x;
  y = rhs.y;
  z = rhs.z;
  
  covariant = rhs.covariant;

  return *this;
}

BoutReal Vector3D::operator=(const BoutReal val)
{
  x = val;
  y = val;
  z = val;
  
  return val;
}

////////////////// ADDITION //////////////////////

Vector3D & Vector3D::operator+=(const Vector3D &rhs)
{
  // Make sure they're of the same type (co/contra-variant)
  if(rhs.covariant) {
    toCovariant();
  }else {
    toContravariant();
  }

  x += rhs.x;
  y += rhs.y;
  z += rhs.z;

  return *this;
}

Vector3D & Vector3D::operator+=(const Vector2D &rhs)
{
  if(rhs.covariant) {
    toCovariant();
  }else {
    toContravariant();
  }

  x += rhs.x;
  y += rhs.y;
  z += rhs.z;

  return *this;
}

///////////////// SUBTRACTION ////////////////////

const Vector3D Vector3D::operator-() const
{
  Vector3D result = *this;

  result.x *= -1.0;
  result.y *= -1.0;
  result.z *= -1.0;

  return result;
}

Vector3D & Vector3D::operator-=(const Vector3D &rhs)
{
  if(rhs.covariant) {
    toCovariant();
  }else {
    toContravariant();
  }
  
  x -= rhs.x;
  y -= rhs.y;
  z -= rhs.z;

  return *this;
}

Vector3D & Vector3D::operator-=(const Vector2D &rhs)
{
  if(rhs.covariant) {
    toCovariant();
  }else {
    toContravariant();
  }
  
  x -= rhs.x;
  y -= rhs.y;
  z -= rhs.z;

  return *this;
}

//////////////// MULTIPLICATION //////////////////

Vector3D & Vector3D::operator*=(const BoutReal rhs)
{
  x *= rhs;
  y *= rhs;
  z *= rhs;
  
  return *this;
}

Vector3D & Vector3D::operator*=(const Field2D &rhs)
{
  x *= rhs;
  y *= rhs;
  z *= rhs;
  
  return *this;
}

Vector3D & Vector3D::operator*=(const Field3D &rhs)
{
  x *= rhs;
  y *= rhs;
  z *= rhs;

  return *this;
}

/////////////////// DIVISION /////////////////////

Vector3D & Vector3D::operator/=(const BoutReal rhs)
{
  x /= rhs;
  y /= rhs;
  z /= rhs;
  
  return *this;
}

Vector3D & Vector3D::operator/=(const Field2D &rhs)
{
  x /= rhs;
  y /= rhs;
  z /= rhs;

  return *this;
}

Vector3D & Vector3D::operator/=(const Field3D &rhs)
{
  x /= rhs;
  y /= rhs;
  z /= rhs;

  return *this;
}

///////////////// CROSS PRODUCT //////////////////

Vector3D & Vector3D::operator^=(const Vector3D &rhs)
{
  Vector3D result;

  // Make sure both vector components are covariant
  Vector3D rco = rhs;
  rco.toCovariant();
  toCovariant();

  // calculate contravariant components of cross-product
  result.x = (y*rco.z - z*rco.y)/mesh->J;
  result.y = (z*rco.x - x*rco.z)/mesh->J;
  result.z = (x*rco.y - y*rco.x)/mesh->J;
  result.covariant = false;

  *this = result;

  return *this;
}

Vector3D & Vector3D::operator^=(const Vector2D &rhs)
{
  Vector3D result;
  
  // Make sure both vector components are covariant
  Vector2D rco = rhs;
  rco.toCovariant();
  toCovariant();
  
  // calculate contravariant components of cross-product
  result.x = (y*rco.z - z*rco.y)/mesh->J;
  result.y = (z*rco.x - x*rco.z)/mesh->J;
  result.z = (x*rco.y - y*rco.x)/mesh->J;
  result.covariant = false;

  *this = result;

  return *this;
}

/***************************************************************
 *                      BINARY OPERATORS 
 ***************************************************************/

////////////////// ADDITION //////////////////////

const Vector3D Vector3D::operator+(const Vector3D &rhs) const
{
  Vector3D result = *this;
  result += rhs;
  return result;
}

const Vector3D Vector3D::operator+(const Vector2D &rhs) const
{
  Vector3D result = *this;
  result += rhs;
  return result;
}

///////////////// SUBTRACTION ////////////////////

const Vector3D Vector3D::operator-(const Vector3D &rhs) const
{
  Vector3D result = *this;
  result -= rhs;
  return result;
}

const Vector3D Vector3D::operator-(const Vector2D &rhs) const
{
  Vector3D result = *this;
  result -= rhs;
  return result;
}

//////////////// MULTIPLICATION //////////////////

const Vector3D Vector3D::operator*(const BoutReal rhs) const
{
  Vector3D result = *this;
  result *= rhs;
  return result;
}

const Vector3D Vector3D::operator*(const Field2D &rhs) const
{
  Vector3D result = *this;
  result *= rhs;
  return result;
}

const Vector3D Vector3D::operator*(const Field3D &rhs) const
{
  Vector3D result = *this;
  result *= rhs;
  return result;
}

/////////////////// DIVISION /////////////////////

const Vector3D Vector3D::operator/(const BoutReal rhs) const
{
  Vector3D result = *this;
  result /= rhs;
  return result;
}

const Vector3D Vector3D::operator/(const Field2D &rhs) const
{
  Vector3D result = *this;
  result /= rhs;
  return result;
}

const Vector3D Vector3D::operator/(const Field3D &rhs) const
{
  Vector3D result = *this;
  result /= rhs;
  return result;
}

////////////////// DOT PRODUCT ///////////////////

const Field3D Vector3D::operator*(const Vector3D &rhs) const
{
  Field3D result;

  if(rhs.covariant ^ covariant) {
    // Both different - just multiply components
    result = x*rhs.x + y*rhs.y + z*rhs.z;
  }else {
    // Both are covariant or contravariant
    if(covariant) {
      // Both covariant
      result = x*rhs.x*mesh->g11 + y*rhs.y*mesh->g22 + z*rhs.z*mesh->g33;
      result += (x*rhs.y + y*rhs.x)*mesh->g12
	+ (x*rhs.z + z*rhs.x)*mesh->g13
	+ (y*rhs.z + z*rhs.y)*mesh->g23;
    }else {
      // Both contravariant
      result = x*rhs.x*mesh->g_11 + y*rhs.y*mesh->g_22 + z*rhs.z*mesh->g_33;
      result += (x*rhs.y + y*rhs.x)*mesh->g_12
	+ (x*rhs.z + z*rhs.x)*mesh->g_13
	+ (y*rhs.z + z*rhs.y)*mesh->g_23;
    }
  }
  
  return result;
}

const Field3D Vector3D::operator*(const Vector2D &rhs) const
{
  Field3D result;

  if(rhs.covariant ^ covariant) {
    // Both different - just multiply components
    result = x*rhs.x + y*rhs.y + z*rhs.z;
  }else {
    // Both are covariant or contravariant
    if(covariant) {
      // Both covariant
      result = x*rhs.x*mesh->g11 + y*rhs.y*mesh->g22 + z*rhs.z*mesh->g33;
      result += (x*rhs.y + y*rhs.x)*mesh->g12
	+ (x*rhs.z + z*rhs.x)*mesh->g13
	+ (y*rhs.z + z*rhs.y)*mesh->g23;
    }else {
      // Both contravariant
      result = x*rhs.x*mesh->g_11 + y*rhs.y*mesh->g_22 + z*rhs.z*mesh->g_33;
      result += (x*rhs.y + y*rhs.x)*mesh->g_12
	+ (x*rhs.z + z*rhs.x)*mesh->g_13
	+ (y*rhs.z + z*rhs.y)*mesh->g_23;
    }
  }

  return result;
}
 
///////////////// CROSS PRODUCT //////////////////

const Vector3D Vector3D::operator^(const Vector3D &rhs) const
{
  Vector3D result = *this;
  
  result ^= rhs;

  return result;
}

const Vector3D Vector3D::operator^(const Vector2D &rhs) const
{
  Vector3D result = *this;
  
  result ^= rhs;

  return result;
}

/***************************************************************
 *       Set variable location for staggered meshes
 ***************************************************************/

void Vector3D::setLocation(CELL_LOC loc) {
  if(loc == CELL_VSHIFT) {
    x.setLocation(CELL_XLOW);
    y.setLocation(CELL_YLOW);
    z.setLocation(CELL_ZLOW);
  }else {
    x.setLocation(loc);
    y.setLocation(loc);
    z.setLocation(loc);
  }
}

/***************************************************************
 *               Z SHIFTING
 ***************************************************************/

const Vector3D Vector3D::shiftZ(const BoutReal zangle) const {
  Vector3D result;

  result.covariant = covariant;

  result.x = x.shiftZ(zangle);
  result.y = y.shiftZ(zangle);
  result.z = z.shiftZ(zangle);

  return result;
}

/***************************************************************
 *               NON-MEMBER OVERLOADED OPERATORS
 ***************************************************************/

const Vector3D operator*(const BoutReal lhs, const Vector3D &rhs) {
  return(rhs * lhs);
}

const Vector3D operator*(const Field2D &lhs, const Vector3D &rhs) {
  return(rhs * lhs);
}

const Vector3D operator*(const Field3D &lhs, const Vector3D &rhs)
{
  return(rhs * lhs);
}

/***************************************************************
 *               NON-MEMBER FUNCTIONS
 ***************************************************************/

// Return the magnitude of a vector
const Field3D abs(const Vector3D &v)
{
  return sqrt(v*v);
}

/***************************************************************
 *               FieldData VIRTUAL FUNCTIONS
 ***************************************************************/

int Vector3D::getData(int jx, int jy, int jz, void *vptr) const
{
#ifdef CHECK
  // check ranges
  if((jx < 0) || (jx >= mesh->ngx) || (jy < 0) || (jy >= mesh->ngy) || (jz < 0) || (jz >= mesh->ngz)) {
    throw BoutException("Vector3D: getData (%d,%d,%d) out of bounds\n", jx, jy, jz);
  }
#endif
  BoutReal *ptr = (BoutReal*) vptr;
  *ptr = x[jx][jy][jz]; ptr++;
  *ptr = y[jx][jy][jz]; ptr++;
  *ptr = z[jx][jy][jz];
  
  return 3*sizeof(BoutReal);
}

int Vector3D::getData(int jx, int jy, int jz, BoutReal *rptr) const
{
#ifdef CHECK
  // check ranges
  if((jx < 0) || (jx >= mesh->ngx) || (jy < 0) || (jy > mesh->ngy) || (jz < 0) || (jz >= mesh->ngz)) {
    throw BoutException("Vector3D: getData (%d,%d,%d) out of bounds\n", jx, jy, jz);
  }
#endif

  *rptr = x[jx][jy][jz]; rptr++;
  *rptr = y[jx][jy][jz]; rptr++;
  *rptr = z[jx][jy][jz];
  
  return 3;
}

int Vector3D::setData(int jx, int jy, int jz, void *vptr)
{
#ifdef CHECK
  // check ranges
  if((jx < 0) || (jx >= mesh->ngx) || (jy < 0) || (jy >= mesh->ngy) || (jz < 0) || (jz >= mesh->ngz)) {
    throw BoutException("Vector3D: setData (%d,%d,%d) out of bounds\n", jx, jy, jz);
  }
#endif
  BoutReal *rptr = (BoutReal*) vptr;
  x[jx][jy][jz] = *rptr; rptr++;
  y[jx][jy][jz] = *rptr; rptr++;
  z[jx][jy][jz] = *rptr;

  return 3*sizeof(BoutReal);
}

int Vector3D::setData(int jx, int jy, int jz, BoutReal *rptr)
{
#ifdef CHECK
  // check ranges
  if((jx < 0) || (jx >= mesh->ngx) || (jy < 0) || (jy >= mesh->ngy) || (jz < 0) || (jz >= mesh->ngz)) {
    throw BoutException("Vector3D: setData (%d,%d,%d) out of bounds\n", jx, jy, jz);
  }
#endif

  x[jx][jy][jz] = *rptr; rptr++;
  y[jx][jy][jz] = *rptr; rptr++;
  z[jx][jy][jz] = *rptr;
  
  return 3;
}

///////////////////// BOUNDARY CONDITIONS //////////////////

void Vector3D::applyBoundary()
{
  for(vector<BoundaryOp*>::iterator it = bndry_op.begin(); it != bndry_op.end(); it++)
    (*it)->apply(*this);
}

void Vector3D::applyTDerivBoundary()
{
  for(vector<BoundaryOp*>::iterator it = bndry_op.begin(); it != bndry_op.end(); it++)
    (*it)->apply_ddt(*this);
}
