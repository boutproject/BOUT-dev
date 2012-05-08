/**************************************************************************
 * Class for 2D vectors. Built on the Field2D class,
 * all operators relating to vectors are here (none in Field classes)
 *
 * As with Field2D, Vector2D are constant in z (toroidal angle) 
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

#include <vector2d.hxx>
#include <boundary_op.hxx>
#include <output.hxx>

Vector2D::Vector2D() : covariant(true), deriv(NULL) { }

Vector2D::Vector2D(const Vector2D &f) : x(f.x), y(f.y), z(f.z), covariant(f.covariant), deriv(NULL) { }

Vector2D::~Vector2D() {
  if(deriv != NULL) {
    // The ddt of the components (x.ddt) point to the same place as ddt.x
    // only delete once
    x.deriv = NULL;
    y.deriv = NULL;
    z.deriv = NULL;
    
    // Now delete them as part of the ddt vector
    delete deriv;
  }
}

void Vector2D::toCovariant() {  
  if(!covariant) {
    Field2D gx, gy, gz;

    // multiply by g_{ij}
    gx = mesh->g_11*x + mesh->g_12*y + mesh->g_13*z;
    gy = mesh->g_12*x + mesh->g_22*y + mesh->g_23*z;
    gz = mesh->g_13*x + mesh->g_23*y + mesh->g_33*z;

    x = gx;
    y = gy;
    z = gz;
    
    covariant = true;
  }
}

void Vector2D::toContravariant() {  
  if(covariant) {
    // multiply by g^{ij}
    
    Field2D gx, gy, gz;

    gx = mesh->g11*x + mesh->g12*y + mesh->g13*z;
    gy = mesh->g12*x + mesh->g22*y + mesh->g23*z;
    gz = mesh->g13*x + mesh->g23*y + mesh->g33*z;

    x = gx;
    y = gy;
    z = gz;
    
    covariant = false;
  }
}

Vector2D* Vector2D::timeDeriv() {
  if(deriv == NULL) {
    deriv = new Vector2D();
    
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

Vector2D & Vector2D::operator=(const Vector2D &rhs) {
  x = rhs.x;
  y = rhs.y;
  z = rhs.z;

  covariant = rhs.covariant;

  return *this;
}

BoutReal Vector2D::operator=(const BoutReal val) {
  x = val;
  y = val;
  z = val;

  return val;
}

////////////////// ADDITION //////////////////////

Vector2D & Vector2D::operator+=(const Vector2D &rhs) {
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

const Vector2D Vector2D::operator-() const {
  Vector2D result = *this;

  result.x *= -1.0;
  result.y *= -1.0;
  result.z *= -1.0;

  return result;
}

Vector2D & Vector2D::operator-=(const Vector2D &rhs) {
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

Vector2D & Vector2D::operator*=(const BoutReal rhs) {
  x *= rhs;
  y *= rhs;
  z *= rhs;
  
  return *this;
}

Vector2D & Vector2D::operator*=(const Field2D &rhs) {
  x *= rhs;
  y *= rhs;
  z *= rhs;
  
  return *this;
}

/////////////////// DIVISION /////////////////////

Vector2D & Vector2D::operator/=(const BoutReal rhs) {
  x /= rhs;
  y /= rhs;
  z /= rhs;
  
  return *this;
}

Vector2D & Vector2D::operator/=(const Field2D &rhs) {
  x /= rhs;
  y /= rhs;
  z /= rhs;

  return *this;
}

///////////////// CROSS PRODUCT //////////////////

Vector2D & Vector2D::operator^=(const Vector2D &rhs) {
  Vector2D result;

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

const Vector2D Vector2D::operator+(const Vector2D &rhs) const {
  Vector2D result = *this;
  result += rhs;
  return result;
}

const Vector3D Vector2D::operator+(const Vector3D &rhs) const {
  return rhs+(*this);
}

///////////////// SUBTRACTION ////////////////////

const Vector2D Vector2D::operator-(const Vector2D &rhs) const {
  Vector2D result = *this;
  result -= rhs;
  return result;
}

const Vector3D Vector2D::operator-(const Vector3D &rhs) const {
  Vector3D result;
  result = *this;
  result -= rhs;
  return result;
}

/////////////// MULTIPLICATION //////////////////

const Vector2D Vector2D::operator*(const BoutReal rhs) const {
  Vector2D result = *this;
  result *= rhs;
  return result;
}

const Vector2D Vector2D::operator*(const Field2D &rhs) const {
  Vector2D result = *this;
  result *= rhs;
  return result;
}

const Vector3D Vector2D::operator*(const Field3D &rhs) const {
  Vector3D result;
  result = *this;
  result *= rhs;
  return result;
}

/////////////////// DIVISION /////////////////////

const Vector2D Vector2D::operator/(const BoutReal rhs) const {
  Vector2D result = *this;
  result /= rhs;
  return result;
}

const Vector2D Vector2D::operator/(const Field2D &rhs) const {
  Vector2D result = *this;
  result /= rhs;
  return result;
}

const Vector3D Vector2D::operator/(const Field3D &rhs) const {
  Vector3D result;
  result = *this;
  result /= rhs;
  return result;
}

////////////////// DOT PRODUCT ///////////////////

const Field2D Vector2D::operator*(const Vector2D &rhs) const {
  Field2D result;

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

const Field3D Vector2D::operator*(const Vector3D &rhs) const {
  return rhs*(*this);
}

///////////////// CROSS PRODUCT //////////////////

const Vector2D Vector2D::operator^(const Vector2D &rhs) const {
  Vector2D result = *this;
  result ^= rhs;
  return result;
}

const Vector3D Vector2D::operator^(const Vector3D &rhs) const {
  return -1.0*rhs^(*this);
}

/***************************************************************
 *               NON-MEMBER OVERLOADED OPERATORS
 ***************************************************************/

const Vector2D operator*(const BoutReal lhs, const Vector2D &rhs) {
  return rhs*lhs;
}

const Vector2D operator*(const Field2D &lhs, const Vector2D &rhs) {
  return rhs*lhs;
}

const Vector3D operator*(const Field3D &lhs, const Vector2D &rhs) {
  return rhs*lhs;
}

/***************************************************************
 *               NON-MEMBER FUNCTIONS
 ***************************************************************/

// Return the magnitude of a vector
const Field2D abs(const Vector2D &v) {
  return sqrt(v*v);
}

/***************************************************************
 *               FieldData VIRTUAL FUNCTIONS
 ***************************************************************/

int Vector2D::getData(int jx, int jy, int jz, void *vptr) const {
#ifdef CHECK
  // check ranges
  if((jx < 0) || (jx >= mesh->ngx) || (jy < 0) || (jy >= mesh->ngy) || (jz < 0) || (jz >= mesh->ngz)) {
    output.write("Vector2D: getData (%d,%d,%d) out of bounds\n", jx, jy, jz);
    exit(1);
  }
#endif
  BoutReal *ptr = (BoutReal*) vptr;
  *ptr = x[jx][jy]; ptr++;
  *ptr = y[jx][jy]; ptr++;
  *ptr = z[jx][jy];
  
  return 3*sizeof(BoutReal);
}

int Vector2D::getData(int jx, int jy, int jz, BoutReal *rptr) const
{
#ifdef CHECK
  // check ranges
  if((jx < 0) || (jx >= mesh->ngx) || (jy < 0) || (jy >= mesh->ngy) || (jz < 0) || (jz >= mesh->ngz)) {
    output.write("Vector2D: getData (%d,%d,%d) out of bounds\n", jx, jy, jz);
    exit(1);
  }
#endif

  *rptr = x[jx][jy]; rptr++;
  *rptr = y[jx][jy]; rptr++;
  *rptr = z[jx][jy];
  
  return 3;
}

int Vector2D::setData(int jx, int jy, int jz, void *vptr)
{
#ifdef CHECK
  // check ranges
  if((jx < 0) || (jx >= mesh->ngx) || (jy < 0) || (jy >= mesh->ngy) || (jz < 0) || (jz >= mesh->ngz)) {
    output.write("Vector2D: setData (%d,%d,%d) out of bounds\n", jx, jy, jz);
    exit(1);
  }
#endif
  BoutReal *rptr = (BoutReal*) vptr;
  x[jx][jy] = *rptr; rptr++;
  y[jx][jy] = *rptr; rptr++;
  z[jx][jy] = *rptr;

  return 3*sizeof(BoutReal);
}

int Vector2D::setData(int jx, int jy, int jz, BoutReal *rptr)
{
#ifdef CHECK
  // check ranges
  if((jx < 0) || (jx >= mesh->ngx) || (jy < 0) || (jy >= mesh->ngy) || (jz < 0) || (jz >= mesh->ngz)) {
    output.write("Vector2D: setData (%d,%d,%d) out of bounds\n", jx, jy, jz);
    exit(1);
  }
#endif

  x[jx][jy] = *rptr; rptr++;
  y[jx][jy] = *rptr; rptr++;
  z[jx][jy] = *rptr;
  
  return 3;
}

///////////////////// BOUNDARY CONDITIONS //////////////////

void Vector2D::applyBoundary()
{
  for(vector<BoundaryOp*>::iterator it = bndry_op.begin(); it != bndry_op.end(); it++)
    (*it)->apply(*this);
}

void Vector2D::applyTDerivBoundary()
{
  for(vector<BoundaryOp*>::iterator it = bndry_op.begin(); it != bndry_op.end(); it++)
    (*it)->apply_ddt(*this);
}

