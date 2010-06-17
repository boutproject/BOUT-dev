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

#include "globals.h"

#include "vector3d.h"

Vector3D::Vector3D()
{
  covariant = true;
}

Vector3D::Vector3D(const Vector3D &f)
{
  *this = f;
}

Vector3D::~Vector3D()
{

}

void Vector3D::to_covariant()
{  
  if(!covariant) {
    Field3D gx, gy, gz;

    // multiply by g_{ij}
    gx = g_11*x + g_12*y + g_13*z;
    gy = g_12*x + g_22*y + g_23*z;
    gz = g_13*x + g_23*y + g_33*z;

    x = gx;
    y = gy;
    z = gz;
    
    covariant = true;
  }
}
void Vector3D::to_contravariant()
{  
  if(covariant) {
    // multiply by g^{ij}
    
    Field3D gx, gy, gz;

    gx = g11*x + g12*y + g13*z;
    gy = g12*x + g22*y + g23*z;
    gz = g13*x + g23*y + g33*z;

    x = gx;
    y = gy;
    z = gz;
    
    covariant = false;
  }
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

real Vector3D::operator=(const real val)
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
    to_covariant();
  }else {
    to_contravariant();
  }

  x += rhs.x;
  y += rhs.y;
  z += rhs.z;

  return *this;
}

Vector3D & Vector3D::operator+=(const Vector2D &rhs)
{
  if(rhs.covariant) {
    to_covariant();
  }else {
    to_contravariant();
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
    to_covariant();
  }else {
    to_contravariant();
  }
  
  x -= rhs.x;
  y -= rhs.y;
  z -= rhs.z;

  return *this;
}

Vector3D & Vector3D::operator-=(const Vector2D &rhs)
{
  if(rhs.covariant) {
    to_covariant();
  }else {
    to_contravariant();
  }
  
  x -= rhs.x;
  y -= rhs.y;
  z -= rhs.z;

  return *this;
}

//////////////// MULTIPLICATION //////////////////

Vector3D & Vector3D::operator*=(const real rhs)
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

Vector3D & Vector3D::operator/=(const real rhs)
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
  rco.to_covariant();
  to_covariant();

  // calculate contravariant components of cross-product
  result.x = (y*rco.z - z*rco.y)/J;
  result.y = (z*rco.x - x*rco.z)/J;
  result.z = (x*rco.y - y*rco.x)/J;
  result.covariant = false;

  *this = result;

  return *this;
}

Vector3D & Vector3D::operator^=(const Vector2D &rhs)
{
  Vector3D result;
  
  // Make sure both vector components are covariant
  Vector2D rco = rhs;
  rco.to_covariant();
  to_covariant();
  
  // calculate contravariant components of cross-product
  result.x = (y*rco.z - z*rco.y)/J;
  result.y = (z*rco.x - x*rco.z)/J;
  result.z = (x*rco.y - y*rco.x)/J;
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

const Vector3D Vector3D::operator*(const real rhs) const
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

const Vector3D Vector3D::operator/(const real rhs) const
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
      result = x*rhs.x*g11 + y*rhs.y*g22 + z*rhs.z*g33;
      result += (x*rhs.y + y*rhs.x)*g12
	+ (x*rhs.z + z*rhs.x)*g13
	+ (y*rhs.z + z*rhs.y)*g23;
    }else {
      // Both contravariant
      result = x*rhs.x*g_11 + y*rhs.y*g_22 + z*rhs.z*g_33;
      result += (x*rhs.y + y*rhs.x)*g_12
	+ (x*rhs.z + z*rhs.x)*g_13
	+ (y*rhs.z + z*rhs.y)*g_23;
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
      result = x*rhs.x*g11 + y*rhs.y*g22 + z*rhs.z*g33;
      result += (x*rhs.y + y*rhs.x)*g12
	+ (x*rhs.z + z*rhs.x)*g13
	+ (y*rhs.z + z*rhs.y)*g23;
    }else {
      // Both contravariant
      result = x*rhs.x*g_11 + y*rhs.y*g_22 + z*rhs.z*g_33;
      result += (x*rhs.y + y*rhs.x)*g_12
	+ (x*rhs.z + z*rhs.x)*g_13
	+ (y*rhs.z + z*rhs.y)*g_23;
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
 *               Z SHIFTING
 ***************************************************************/

const Vector3D Vector3D::ShiftZ(const real zangle) const
{
  Vector3D result;

  result.covariant = covariant;

  result.x = x.ShiftZ(zangle);
  result.y = y.ShiftZ(zangle);
  result.z = z.ShiftZ(zangle);

  return result;
}

/***************************************************************
 *               NON-MEMBER OVERLOADED OPERATORS
 ***************************************************************/

const Vector3D operator*(const real lhs, const Vector3D &rhs)
{
  return(rhs * lhs);
}

const Vector3D operator*(const Field2D &lhs, const Vector3D &rhs)
{
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
  if((jx < 0) || (jx > ncx) || (jy < 0) || (jy > ncy) || (jz < 0) || (jz >= ncz)) {
    output.write("Vector3D: getData (%d,%d,%d) out of bounds\n", jx, jy, jz);
    exit(1);
  }
#endif
  real *ptr = (real*) vptr;
  *ptr = x[jx][jy][jz]; ptr++;
  *ptr = y[jx][jy][jz]; ptr++;
  *ptr = z[jx][jy][jz];
  
  return 3*sizeof(real);
}

int Vector3D::getData(int jx, int jy, int jz, real *rptr) const
{
#ifdef CHECK
  // check ranges
  if((jx < 0) || (jx > ncx) || (jy < 0) || (jy > ncy) || (jz < 0) || (jz >= ncz)) {
    output.write("Vector3D: getData (%d,%d,%d) out of bounds\n", jx, jy, jz);
    exit(1);
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
  if((jx < 0) || (jx > ncx) || (jy < 0) || (jy > ncy) || (jz < 0) || (jz >= ncz)) {
    output.write("Vector3D: setData (%d,%d,%d) out of bounds\n", jx, jy, jz);
    exit(1);
  }
#endif
  real *rptr = (real*) vptr;
  x[jx][jy][jz] = *rptr; rptr++;
  y[jx][jy][jz] = *rptr; rptr++;
  z[jx][jy][jz] = *rptr;

  return 3*sizeof(real);
}

int Vector3D::setData(int jx, int jy, int jz, real *rptr)
{
#ifdef CHECK
  // check ranges
  if((jx < 0) || (jx > ncx) || (jy < 0) || (jy > ncy) || (jz < 0) || (jz >= ncz)) {
    output.write("Vector3D: setData (%d,%d,%d) out of bounds\n", jx, jy, jz);
    exit(1);
  }
#endif

  x[jx][jy][jz] = *rptr; rptr++;
  y[jx][jy][jz] = *rptr; rptr++;
  z[jx][jy][jz] = *rptr;
  
  return 3;
}

