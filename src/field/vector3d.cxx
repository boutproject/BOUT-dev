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
#include <bout/assert.hxx>

Vector3D::Vector3D(Mesh *localmesh)
    : x(localmesh), y(localmesh), z(localmesh), covariant(true), deriv(nullptr), location(CELL_CENTRE) {}

Vector3D::Vector3D(const Vector3D &f)
    : x(f.x), y(f.y), z(f.y), covariant(f.covariant), deriv(nullptr), location(CELL_CENTRE) {}

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
    Mesh *localmesh = x.getMesh();
    Field3D gx(localmesh), gy(localmesh), gz(localmesh);

    Coordinates *metric = localmesh->coordinates();

    // multiply by g_{ij}
    gx = x*metric->g_11 + metric->g_12*y + metric->g_13*z;
    gy = y*metric->g_22 + metric->g_12*x + metric->g_23*z;
    gz = z*metric->g_33 + metric->g_13*x + metric->g_23*y;

    x = gx;
    y = gy;
    z = gz;
    
    covariant = true;
  }
}
void Vector3D::toContravariant() {  
  if(covariant) {
    // multiply by g^{ij}
    Mesh *localmesh = x.getMesh();
    Field3D gx(localmesh), gy(localmesh), gz(localmesh);

    Coordinates *metric = localmesh->coordinates();

    gx = x*metric->g11 + metric->g12*y + metric->g13*z;
    gy = y*metric->g22 + metric->g12*x + metric->g23*z;
    gz = z*metric->g33 + metric->g13*x + metric->g23*y;

    x = gx;
    y = gy;
    z = gz;
    
    covariant = false;
  }
}

Vector3D* Vector3D::timeDeriv() {
  if(deriv == NULL) {
    deriv = new Vector3D(x.getMesh());

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

Vector3D & Vector3D::operator=(const Vector3D &rhs) {
  x = rhs.x;
  y = rhs.y;
  z = rhs.z;

  covariant = rhs.covariant;

  return *this;
}

Vector3D & Vector3D::operator=(const Vector2D &rhs) {
  x = rhs.x;
  y = rhs.y;
  z = rhs.z;
  
  covariant = rhs.covariant;

  return *this;
}

Vector3D & Vector3D::operator=(const BoutReal val)
{
  x = val;
  y = val;
  z = val;
  
  return *this;
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

Vector3D & Vector3D::operator^=(const Vector3D &rhs) {
  Mesh *localmesh = x.getMesh();
  Vector3D result(localmesh);

  // Make sure both vector components are covariant
  Vector3D rco = rhs;
  rco.toCovariant();
  toCovariant();

  Coordinates *metric = localmesh->coordinates();

  // calculate contravariant components of cross-product
  result.x = (y*rco.z - z*rco.y)/metric->J;
  result.y = (z*rco.x - x*rco.z)/metric->J;
  result.z = (x*rco.y - y*rco.x)/metric->J;
  result.covariant = false;

  *this = result;

  return *this;
}

Vector3D & Vector3D::operator^=(const Vector2D &rhs) {
  Mesh *localmesh = x.getMesh();
  Vector3D result(localmesh);

  // Make sure both vector components are covariant
  Vector2D rco = rhs;
  rco.toCovariant();
  toCovariant();

  Coordinates *metric = localmesh->coordinates();

  // calculate contravariant components of cross-product
  result.x = (y*rco.z - z*rco.y)/metric->J;
  result.y = (z*rco.x - x*rco.z)/metric->J;
  result.z = (x*rco.y - y*rco.x)/metric->J;
  result.covariant = false;

  *this = result;

  return *this;
}

/***************************************************************
 *                      BINARY OPERATORS 
 ***************************************************************/

////////////////// ADDITION //////////////////////

const Vector3D Vector3D::operator+(const Vector3D &rhs) const {
  Vector3D result = *this;
  result += rhs;
  return result;
}

const Vector3D Vector3D::operator+(const Vector2D &rhs) const {
  Vector3D result = *this;
  result += rhs;
  return result;
}

///////////////// SUBTRACTION ////////////////////

const Vector3D Vector3D::operator-(const Vector3D &rhs) const {
  Vector3D result = *this;
  result -= rhs;
  return result;
}

const Vector3D Vector3D::operator-(const Vector2D &rhs) const {
  Vector3D result = *this;
  result -= rhs;
  return result;
}

//////////////// MULTIPLICATION //////////////////

const Vector3D Vector3D::operator*(const BoutReal rhs) const {
  Vector3D result = *this;
  result *= rhs;
  return result;
}

const Vector3D Vector3D::operator*(const Field2D &rhs) const {
  Vector3D result = *this;
  result *= rhs;
  return result;
}

const Vector3D Vector3D::operator*(const Field3D &rhs) const {
  Vector3D result = *this;
  result *= rhs;
  return result;
}

/////////////////// DIVISION /////////////////////

const Vector3D Vector3D::operator/(const BoutReal rhs) const {
  Vector3D result = *this;
  result /= rhs;
  return result;
}

const Vector3D Vector3D::operator/(const Field2D &rhs) const {
  Vector3D result = *this;
  result /= rhs;
  return result;
}

const Vector3D Vector3D::operator/(const Field3D &rhs) const {
  Vector3D result = *this;
  result /= rhs;
  return result;
}

////////////////// DOT PRODUCT ///////////////////

const Field3D Vector3D::operator*(const Vector3D &rhs) const {
  Field3D result(x.getMesh());

  if(rhs.covariant ^ covariant) {
    // Both different - just multiply components
    result = x*rhs.x + y*rhs.y + z*rhs.z;
  }else {
    // Both are covariant or contravariant

    Coordinates *metric = mesh->coordinates();
    
    if(covariant) {
      // Both covariant
      result = x*rhs.x*metric->g11 + y*rhs.y*metric->g22 + z*rhs.z*metric->g33;
      result += (x*rhs.y + y*rhs.x)*metric->g12
	+ (x*rhs.z + z*rhs.x)*metric->g13
	+ (y*rhs.z + z*rhs.y)*metric->g23;
    }else {
      // Both contravariant
      result = x*rhs.x*metric->g_11 + y*rhs.y*metric->g_22 + z*rhs.z*metric->g_33;
      result += (x*rhs.y + y*rhs.x)*metric->g_12
	+ (x*rhs.z + z*rhs.x)*metric->g_13
	+ (y*rhs.z + z*rhs.y)*metric->g_23;
    }
  }
  
  return result;
}

const Field3D Vector3D::operator*(const Vector2D &rhs) const
{
  Field3D result(x.getMesh());

  if(rhs.covariant ^ covariant) {
    // Both different - just multiply components
    result = x*rhs.x + y*rhs.y + z*rhs.z;
  }else {
    // Both are covariant or contravariant

    Coordinates *metric = x.getMesh()->coordinates();
    if(covariant) {
      // Both covariant
      result = x*rhs.x*metric->g11 + y*rhs.y*metric->g22 + z*rhs.z*metric->g33;
      result += (x*rhs.y + y*rhs.x)*metric->g12
	+ (x*rhs.z + z*rhs.x)*metric->g13
	+ (y*rhs.z + z*rhs.y)*metric->g23;
    }else {
      // Both contravariant
      result = x*rhs.x*metric->g_11 + y*rhs.y*metric->g_22 + z*rhs.z*metric->g_33;
      result += (x*rhs.y + y*rhs.x)*metric->g_12
	+ (x*rhs.z + z*rhs.x)*metric->g_13
	+ (y*rhs.z + z*rhs.y)*metric->g_23;
    }
  }

  return result;
}
 
///////////////// CROSS PRODUCT //////////////////

const Vector3D Vector3D::operator^(const Vector3D &rhs) const {
  Vector3D result = *this;
  
  result ^= rhs;

  return result;
}

const Vector3D Vector3D::operator^(const Vector2D &rhs) const {
  Vector3D result = *this;
  
  result ^= rhs;

  return result;
}

/***************************************************************
 *       Get/set variable location for staggered meshes
 ***************************************************************/

CELL_LOC Vector3D::getLocation() const {

  ASSERT1(((location == CELL_VSHIFT) && (x.getLocation() == CELL_XLOW) &&
           (y.getLocation() == CELL_YLOW) && (z.getLocation() == CELL_ZLOW)) ||
          ((location == x.getLocation()) && (location == y.getLocation()) &&
           (location == z.getLocation())));

  return location;
}

void Vector3D::setLocation(CELL_LOC loc) {
  location = loc;
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
const Field3D abs(const Vector3D &v) {
  return sqrt(v*v);
}

/***************************************************************
 *               FieldData VIRTUAL FUNCTIONS
 ***************************************************************/

////////////////////////////////////////////////////////////
// Visitor pattern support
void Vector3D::accept(FieldVisitor &v) {
  v.accept(*this);
}

///////////////////// BOUNDARY CONDITIONS //////////////////

void Vector3D::applyBoundary(bool init)
{
  for(const auto& bndry : bndry_op)
    if ( !bndry->apply_to_ddt || init) // Always apply to the values when initialising fields, otherwise apply only if wanted
      bndry->apply(*this);
}

void Vector3D::applyTDerivBoundary()
{
  for(const auto& bndry : bndry_op)
    bndry->apply_ddt(*this);
}
