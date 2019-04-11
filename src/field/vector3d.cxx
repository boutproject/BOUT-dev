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
#include <bout/scorepwrapper.hxx>
#include <interpolation.hxx>

Vector3D::Vector3D(Mesh *localmesh)
    : x(localmesh), y(localmesh), z(localmesh), covariant(true), deriv(nullptr), location(CELL_CENTRE) {}

Vector3D::Vector3D(const Vector3D &f)
    : x(f.x), y(f.y), z(f.z), covariant(f.covariant), deriv(nullptr),
      location(f.getLocation()) {}

Vector3D::~Vector3D() {
  if (deriv != nullptr) {
    // The ddt of the components (x.ddt) point to the same place as ddt.x
    // only delete once
    x.deriv = nullptr;
    y.deriv = nullptr;
    z.deriv = nullptr;

    // Now delete them as part of the deriv vector
    delete deriv;
  }
}

void Vector3D::toCovariant() {
  SCOREP0();  
  if(!covariant) {
    Mesh *localmesh = x.getMesh();

    if (location == CELL_VSHIFT) {
      Coordinates *metric_x, *metric_y, *metric_z;
      metric_x = localmesh->getCoordinates(CELL_XLOW);
      metric_y = localmesh->getCoordinates(CELL_YLOW);
      metric_z = localmesh->getCoordinates(CELL_ZLOW);

      // Fields at different locations so we need to interpolate
      // Note : Could reduce peak memory requirement here by just
      // dealing with the three components seperately. This would
      // require the use of temporary fields to hold the intermediate
      // result so would likely only reduce memory usage by one field
      const auto y_at_x = interp_to(y, x.getLocation());
      const auto z_at_x = interp_to(z, x.getLocation());
      const auto x_at_y = interp_to(x, y.getLocation());
      const auto z_at_y = interp_to(z, y.getLocation());
      const auto x_at_z = interp_to(x, z.getLocation());
      const auto y_at_z = interp_to(y, z.getLocation());

      // multiply by g_{ij}
      BOUT_FOR(i, localmesh->getRegion3D("RGN_ALL")){
        x[i] = metric_x->g_11[i]*x[i] + metric_x->g_12[i]*y_at_x[i] + metric_x->g_13[i]*z_at_x[i];
        y[i] = metric_y->g_22[i]*y[i] + metric_y->g_12[i]*x_at_y[i] + metric_y->g_23[i]*z_at_y[i];
        z[i] = metric_z->g_33[i]*z[i] + metric_z->g_13[i]*x_at_z[i] + metric_z->g_23[i]*y_at_z[i];
      };
    } else {
      const auto metric = localmesh->getCoordinates(location);

      // Need to use temporary arrays to store result
      Field3D gx{emptyFrom(x)}, gy{emptyFrom(y)}, gz{emptyFrom(z)};

      BOUT_FOR(i, localmesh->getRegion3D("RGN_ALL")){
        gx[i] = metric->g_11[i]*x[i] + metric->g_12[i]*y[i] + metric->g_13[i]*z[i];
        gy[i] = metric->g_22[i]*y[i] + metric->g_12[i]*x[i] + metric->g_23[i]*z[i];
        gz[i] = metric->g_33[i]*z[i] + metric->g_13[i]*x[i] + metric->g_23[i]*y[i];
      };

      x = gx;
      y = gy;
      z = gz;
    }

    covariant = true;
  }
}
void Vector3D::toContravariant() {  
  SCOREP0();
  if(covariant) {
    // multiply by g^{ij}
    Mesh *localmesh = x.getMesh();

    if (location == CELL_VSHIFT) {
      Coordinates *metric_x, *metric_y, *metric_z;
    
      metric_x = localmesh->getCoordinates(CELL_XLOW);
      metric_y = localmesh->getCoordinates(CELL_YLOW);
      metric_z = localmesh->getCoordinates(CELL_ZLOW);

      // Fields at different locations so we need to interpolate
      // Note : Could reduce peak memory requirement here by just
      // dealing with the three components seperately. This would
      // require the use of temporary fields to hold the intermediate
      // result so would likely only reduce memory usage by one field
      const auto y_at_x = interp_to(y, x.getLocation());
      const auto z_at_x = interp_to(z, x.getLocation());
      const auto x_at_y = interp_to(x, y.getLocation());
      const auto z_at_y = interp_to(z, y.getLocation());
      const auto x_at_z = interp_to(x, z.getLocation());
      const auto y_at_z = interp_to(y, z.getLocation());

      // multiply by g_{ij}
      BOUT_FOR(i, localmesh->getRegion3D("RGN_ALL")){
        x[i] = metric_x->g11[i]*x[i] + metric_x->g12[i]*y_at_x[i] + metric_x->g13[i]*z_at_x[i];
        y[i] = metric_y->g22[i]*y[i] + metric_y->g12[i]*x_at_y[i] + metric_y->g23[i]*z_at_y[i];
        z[i] = metric_z->g33[i]*z[i] + metric_z->g13[i]*x_at_z[i] + metric_z->g23[i]*y_at_z[i];
      };

    } else {
      const auto metric = localmesh->getCoordinates(location);

      // Need to use temporary arrays to store result
      Field3D gx{emptyFrom(x)}, gy{emptyFrom(y)}, gz{emptyFrom(z)};

      BOUT_FOR(i, localmesh->getRegion3D("RGN_ALL")){
        gx[i] = metric->g11[i]*x[i] + metric->g12[i]*y[i] + metric->g13[i]*z[i];
        gy[i] = metric->g22[i]*y[i] + metric->g12[i]*x[i] + metric->g23[i]*z[i];
        gz[i] = metric->g33[i]*z[i] + metric->g13[i]*x[i] + metric->g23[i]*y[i];
      };

      x = gx;
      y = gy;
      z = gz;
    }
    
    covariant = false;
  }
}

Vector3D* Vector3D::timeDeriv() {
  if (deriv == nullptr) {
    deriv = new Vector3D(x.getMesh());

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
  SCOREP0();
  x = rhs.x;
  y = rhs.y;
  z = rhs.z;

  covariant = rhs.covariant;

  setLocation(rhs.getLocation());
  return *this;
}

Vector3D & Vector3D::operator=(const Vector2D &rhs) {
  SCOREP0();  
  x = rhs.x;
  y = rhs.y;
  z = rhs.z;
  
  covariant = rhs.covariant;

  setLocation(rhs.getLocation());

  return *this;
}

Vector3D & Vector3D::operator=(const BoutReal val)
{
  SCOREP0();
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

#define CROSS(v0, v1, v2)                                               \
                                                                        \
  const v0 cross(const v1 &lhs, const v2 &rhs) {                        \
    ASSERT1(lhs.getLocation() == rhs.getLocation());                    \
    Mesh *localmesh = lhs.x.getMesh();                                  \
    v0 result(localmesh);                                               \
                                                                        \
    /* Make sure both vector components are covariant */                \
    v2 rco = rhs;                                                       \
    rco.toCovariant();                                                  \
    v1 lco = lhs;                                                       \
    lco.toCovariant();                                                  \
                                                                        \
    Coordinates *metric = localmesh->getCoordinates(lhs.getLocation());    \
                                                                        \
    /* calculate contravariant components of cross-product */           \
    result.x = (lco.y * rco.z - lco.z * rco.y) / metric->J;             \
    result.y = (lco.z * rco.x - lco.x * rco.z) / metric->J;             \
    result.z = (lco.x * rco.y - lco.y * rco.x) / metric->J;             \
    result.covariant = false;                                           \
                                                                        \
    return result;                                                      \
  };                                                                    \


CROSS(Vector3D, Vector3D, Vector3D);
CROSS(Vector3D, Vector3D, Vector2D);
CROSS(Vector3D, Vector2D, Vector3D);
CROSS(Vector2D, Vector2D, Vector2D);

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
  Mesh* mesh = x.getMesh();

  Field3D result{emptyFrom(x)};
  ASSERT2(location == rhs.getLocation())

  if(rhs.covariant ^ covariant) {
    // Both different - just multiply components
    result = x*rhs.x + y*rhs.y + z*rhs.z;
  }else {
    // Both are covariant or contravariant

    Coordinates *metric = mesh->getCoordinates(location);
    
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
  ASSERT2(location == rhs.getLocation());

  Field3D result{emptyFrom(x)};

  if(rhs.covariant ^ covariant) {
    // Both different - just multiply components
    result = x*rhs.x + y*rhs.y + z*rhs.z;
  }else {
    // Both are covariant or contravariant

    Coordinates *metric = x.getCoordinates(location);
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
 
/***************************************************************
 *       Get/set variable location for staggered meshes
 ***************************************************************/

CELL_LOC Vector3D::getLocation() const {

  if (location == CELL_VSHIFT) {
    ASSERT1((x.getLocation() == CELL_XLOW) && (y.getLocation() == CELL_YLOW) &&
            (z.getLocation() == CELL_ZLOW));
  } else {
    ASSERT1((location == x.getLocation()) && (location == y.getLocation()) &&
            (location == z.getLocation()));
  }

  return location;
}

void Vector3D::setLocation(CELL_LOC loc) {
  SCOREP0();  
  TRACE("Vector3D::setLocation");
  if (loc == CELL_DEFAULT) {
    loc = CELL_CENTRE;
  }

  if (x.getMesh()->StaggerGrids) {
    if (loc == CELL_VSHIFT) {
      x.setLocation(CELL_XLOW);
      y.setLocation(CELL_YLOW);
      z.setLocation(CELL_ZLOW);
    } else {
      x.setLocation(loc);
      y.setLocation(loc);
      z.setLocation(loc);
    }
  } else {
#if CHECK > 0
    if (loc != CELL_CENTRE) {
      throw BoutException("Vector3D: Trying to set off-centre location on "
                          "non-staggered grid\n"
                          "         Did you mean to enable staggered grids?");
    }
#endif
  }

  location = loc;
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
const Field3D abs(const Vector3D &v, REGION region) {
  return sqrt(v*v, region);
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
