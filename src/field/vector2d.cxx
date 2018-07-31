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
#include <boutexception.hxx>
#include <interpolation.hxx>

Vector2D::Vector2D(Mesh *localmesh)
    : x(localmesh), y(localmesh), z(localmesh), covariant(true), deriv(nullptr), location(CELL_CENTRE) {}

Vector2D::Vector2D(const Vector2D &f)
    : x(f.x), y(f.y), z(f.z), covariant(f.covariant), deriv(nullptr), location(CELL_CENTRE) {}

Vector2D::~Vector2D() {
  if (deriv != nullptr) {
    // The ddt of the components (x.ddt) point to the same place as ddt.x
    // only delete once
    x.deriv = nullptr;
    y.deriv = nullptr;
    z.deriv = nullptr;

    // Now delete them as part of the ddt vector
    delete deriv;
  }
}

void Vector2D::toCovariant() {  
  if(!covariant) {
    Mesh *localmesh = x.getMesh();
    Field2D gx(localmesh), gy(localmesh), gz(localmesh);

    Coordinates *metric_x, *metric_y, *metric_z;
    if (location == CELL_VSHIFT) {
      metric_x = localmesh->coordinates(CELL_XLOW);
      metric_y = localmesh->coordinates(CELL_YLOW);
      metric_z = localmesh->coordinates(CELL_ZLOW);
    } else {
      metric_x = localmesh->coordinates(location);
      metric_y = localmesh->coordinates(location);
      metric_z = localmesh->coordinates(location);
    }

    // multiply by g_{ij}
    gx = x*metric_x->g_11 + metric_x->g_12*interp_to(y, x.getLocation()) + metric_x->g_13*interp_to(z, x.getLocation());
    gy = y*metric_y->g_22 + metric_y->g_12*interp_to(x, y.getLocation()) + metric_y->g_23*interp_to(z, y.getLocation());
    gz = z*metric_z->g_33 + metric_z->g_13*interp_to(x, z.getLocation()) + metric_z->g_23*interp_to(y, z.getLocation());

    x = gx;
    y = gy;
    z = gz;
    
    covariant = true;
  }
}

void Vector2D::toContravariant() {  
  if(covariant) {
    // multiply by g^{ij}
    Mesh *localmesh = x.getMesh();
    Field2D gx(localmesh), gy(localmesh), gz(localmesh);

    Coordinates *metric_x, *metric_y, *metric_z;
    if (location == CELL_VSHIFT) {
      metric_x = localmesh->coordinates(CELL_XLOW);
      metric_y = localmesh->coordinates(CELL_YLOW);
      metric_z = localmesh->coordinates(CELL_ZLOW);
    } else {
      metric_x = localmesh->coordinates(location);
      metric_y = localmesh->coordinates(location);
      metric_z = localmesh->coordinates(location);
    }

    // multiply by g_{ij}
    gx = x*metric_x->g11 + metric_x->g12*interp_to(y, x.getLocation()) + metric_x->g13*interp_to(z, x.getLocation());
    gy = y*metric_y->g22 + metric_y->g12*interp_to(x, y.getLocation()) + metric_y->g23*interp_to(z, y.getLocation());
    gz = z*metric_z->g33 + metric_z->g13*interp_to(x, z.getLocation()) + metric_z->g23*interp_to(y, z.getLocation());

    x = gx;
    y = gy;
    z = gz;
    
    covariant = false;
  }
}

Vector2D* Vector2D::timeDeriv() {
  if (deriv == nullptr) {
    deriv = new Vector2D(x.getMesh());

    // Check if the components have a time-derivative
    // Need to make sure that ddt(v.x) = ddt(v).x

    if (x.deriv != nullptr) {
      // already set. Copy across then delete
      deriv->x = *(x.deriv);
      delete x.deriv;
    }
    if (y.deriv != nullptr) {
      deriv->y = *(y.deriv);
      delete y.deriv;
    }
    if (z.deriv != nullptr) {
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

Vector2D & Vector2D::operator=(const BoutReal val) {
  x = val;
  y = val;
  z = val;

  return *this;
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

// cross product implementation in vector3d.cxx
Vector2D & Vector2D::operator^=(const Vector2D &rhs) {
  *this = cross(*this, rhs);
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
  Vector3D result(x.getMesh());
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
  Vector3D result(x.getMesh());
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
  Vector3D result(x.getMesh());
  result = *this;
  result /= rhs;
  return result;
}

////////////////// DOT PRODUCT ///////////////////

const Field2D Vector2D::operator*(const Vector2D &rhs) const {
  ASSERT2(location == rhs.getLocation());

  Mesh *localmesh = x.getMesh();
  Field2D result(localmesh);

  if(rhs.covariant ^ covariant) {
    // Both different - just multiply components
    result = x*rhs.x + y*rhs.y + z*rhs.z;
  }else {
    // Both are covariant or contravariant
    Coordinates *metric = localmesh->coordinates(location);

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

const Field3D Vector2D::operator*(const Vector3D &rhs) const {
  return rhs*(*this);
}

///////////////// CROSS PRODUCT //////////////////

const Vector2D Vector2D::operator^(const Vector2D &rhs) const {
  return cross(*this,rhs);
}

const Vector3D Vector2D::operator^(const Vector3D &rhs) const {
  return cross(*this,rhs);
}

/***************************************************************
 *       Get/set variable location for staggered meshes
 ***************************************************************/

CELL_LOC Vector2D::getLocation() const {

  if (location == CELL_VSHIFT) {
    ASSERT1((x.getLocation() == CELL_XLOW) && (y.getLocation() == CELL_YLOW) &&
            (z.getLocation() == CELL_ZLOW));
  } else {
    ASSERT1((location == x.getLocation()) && (location == y.getLocation()) &&
            (location == z.getLocation()));
  }

  return location;
}

void Vector2D::setLocation(CELL_LOC loc) {
  location = loc;
  if(loc == CELL_VSHIFT) {
    x.setLocation(CELL_XLOW);
    y.setLocation(CELL_YLOW);
    z.setLocation(CELL_ZLOW);
  } else {
    x.setLocation(loc);
    y.setLocation(loc);
    z.setLocation(loc);
  }
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
const Field2D abs(const Vector2D &v, REGION region) {
  return sqrt(v*v, region);
}

/***************************************************************
 *               FieldData VIRTUAL FUNCTIONS
 ***************************************************************/

////////////////////////////////////////////////////////////
// Visitor pattern support

void Vector2D::accept(FieldVisitor &v) {
  v.accept(*this);
}

///////////////////// BOUNDARY CONDITIONS //////////////////

void Vector2D::applyBoundary(bool init)
{
  for(const auto& bndry : bndry_op)
    if ( !bndry->apply_to_ddt || init) // Always apply to the values when initialising fields, otherwise apply only if wanted
      bndry->apply(*this);
}

void Vector2D::applyTDerivBoundary()
{
  for(const auto& bndry : bndry_op)
    bndry->apply_ddt(*this);
}

