/**************************************************************************
 * Operators on vector objects
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
#include <vecops.hxx>
#include <derivs.hxx>
#include <msg_stack.hxx>
#include <unused.hxx>

/**************************************************************************
 * Gradient operators
 **************************************************************************/

const Vector2D Grad(const Field2D &f, CELL_LOC outloc) {
  TRACE("Grad( Field2D )");

  CELL_LOC outloc_x, outloc_y, outloc_z;
  if (outloc == CELL_VSHIFT) {
    outloc_x = CELL_XLOW;
    outloc_y = CELL_YLOW;
    outloc_z = CELL_ZLOW;
  } else {
    outloc_x = outloc_y = outloc_z = outloc;
  }

  Vector2D result(f.getMesh());

  result.x = DDX(f, outloc_x);
  result.y = DDY(f, outloc_y);
  result.z = DDZ(f, outloc_z);

  if (outloc == CELL_DEFAULT) {
    result.setLocation(result.x.getLocation());
  } else {
    result.setLocation(outloc);
  }

  result.covariant = true;

  return result;
}

const Vector3D Grad(const Field3D &f, CELL_LOC outloc_x, CELL_LOC outloc_y,
                    CELL_LOC outloc_z) {
  // Note no Vector2D equivalent to this three location overload
  TRACE("Grad( Field3D )");

  ASSERT1((outloc_x == outloc_y && outloc_x == outloc_z) ||
          (outloc_x == CELL_XLOW && outloc_y == CELL_YLOW &&
           outloc_z == CELL_ZLOW)); // CELL_VSHIFT

  CELL_LOC outloc =
      (outloc_x == outloc_y && outloc_x == outloc_z) ? outloc_x : CELL_VSHIFT;
  return Grad(f, outloc);
}

const Vector3D Grad(const Field3D &f, CELL_LOC outloc) {
  TRACE("Grad( Field3D )");

  CELL_LOC outloc_x, outloc_y, outloc_z;
  if (outloc == CELL_VSHIFT) {
    outloc_x = CELL_XLOW;
    outloc_y = CELL_YLOW;
    outloc_z = CELL_ZLOW;
  } else {
    outloc_x = outloc_y = outloc_z = outloc;
  }

  Vector3D result(f.getMesh());

  result.x = DDX(f, outloc_x);
  result.y = DDY(f, outloc_y);
  result.z = DDZ(f, outloc_z);

  if (outloc == CELL_DEFAULT) {
    result.setLocation(result.x.getLocation());
  } else {
    result.setLocation(outloc);
  }

  result.covariant = true;

  return result;
}

const Vector3D Grad_perp(const Field3D &f, CELL_LOC outloc_x,
                         MAYBE_UNUSED(CELL_LOC outloc_y),
                         MAYBE_UNUSED(CELL_LOC outloc_z)) {
  TRACE("Grad_perp( Field3D )");
  ASSERT1(outloc_x == outloc_y && outloc_x == outloc_z);
  ASSERT1(outloc_x == CELL_DEFAULT || outloc_x == f.getLocation());
  return Grad_perp(f, outloc_x);
}

const Vector3D Grad_perp(const Field3D &f, CELL_LOC outloc) {
  TRACE("Grad_perp( Field3D )");

  ASSERT1(outloc == CELL_DEFAULT || outloc == f.getLocation());

  Coordinates *metric = f.getCoordinates(outloc);

  Vector3D result(f.getMesh());

  result.x = DDX(f, outloc) - metric->g_12 * DDY(f, outloc) / SQ(metric->J * metric->Bxy);
  result.y = 0.0;
  result.z = DDZ(f, outloc) - metric->g_23 * DDY(f, outloc) / SQ(metric->J * metric->Bxy);

  result.setLocation(result.x.getLocation());

  result.covariant = true;

  return result;
}

/**************************************************************************
 * Divergence operators
 **************************************************************************/

const Field2D Div(const Vector2D &v, CELL_LOC outloc) {
  TRACE("Div( Vector2D )");

  if (outloc == CELL_DEFAULT) {
    outloc = v.getLocation();
  }

  ASSERT1(outloc != CELL_VSHIFT);

  Mesh *localmesh = v.x.getMesh();
  Field2D result(localmesh);

  Coordinates *metric = localmesh->getCoordinates(outloc);

  // get contravariant components of v
  Vector2D vcn = v;
  vcn.toContravariant();
  
  result = DDX(metric->J*vcn.x, outloc);
  result += DDY(metric->J*vcn.y, outloc);
  result += DDZ(metric->J*vcn.z, outloc);
  result /= metric->J;

  return result;
}

const Field3D Div(const Vector3D &v, CELL_LOC outloc) {
  TRACE("Div( Vector3D )");

  if (outloc == CELL_DEFAULT) {
    outloc = v.getLocation();
  }

  // This also catches the combination of v at VSHIFT and outloc at DEFAULT
  ASSERT1(outloc != CELL_VSHIFT);

  Mesh *localmesh = v.x.getMesh();
  Field3D result(localmesh);

  Coordinates *metric = localmesh->getCoordinates(outloc);

  // get contravariant components of v
  Vector3D vcn = v;
  vcn.toContravariant();

  result = DDX(vcn.x.getCoordinates()->J * vcn.x, outloc);
  result += DDY(vcn.y.getCoordinates()->J * vcn.y, outloc);
  result += DDZ(vcn.z.getCoordinates()->J * vcn.z, outloc);
  result /= metric->J;

  return result;
}

/**************************************************************************
 * Divergence operators for flux methods
 **************************************************************************/

const Field2D Div(const Vector2D &v, const Field2D &f, CELL_LOC outloc) {
  TRACE("Div( Vector2D, Field2D )");

  if (outloc == CELL_DEFAULT) {
    outloc = v.getLocation();
  }

  ASSERT1(outloc != CELL_VSHIFT);

  Mesh *localmesh = f.getMesh();

  Coordinates *metric = localmesh->getCoordinates(outloc);

  // get contravariant components of v
  Vector2D vcn = v;
  vcn.toContravariant();

  Field2D result(localmesh);
  result = FDDX(vcn.x.getCoordinates()->J * vcn.x, f, outloc);
  result += FDDY(vcn.y.getCoordinates()->J * vcn.y, f, outloc);
  result += FDDZ(vcn.z.getCoordinates()->J * vcn.z, f, outloc);
  result /= metric->J;

  return result;
}

const Field3D Div(const Vector3D &v, const Field3D &f, DIFF_METHOD method,
                  CELL_LOC outloc) {
  TRACE("Div( Vector3D, Field3D )");

  if (outloc == CELL_DEFAULT) {
    outloc = v.getLocation();
  }
  ASSERT1(outloc != CELL_VSHIFT);

  Mesh *localmesh = f.getMesh();

  Coordinates *metric = localmesh->getCoordinates(outloc);

  // get contravariant components of v
  Vector3D vcn = v;
  vcn.toContravariant();

  Field3D result(localmesh);
  result = FDDX(vcn.x.getCoordinates()->J * vcn.x, f, outloc, method);
  result += FDDY(vcn.y.getCoordinates()->J * vcn.y, f, outloc, method);
  result += FDDZ(vcn.z.getCoordinates()->J * vcn.z, f, outloc, method);
  result /= metric->J;

  return result;
}

const Field3D Div(const Vector3D &v, const Field3D &f, CELL_LOC outloc, DIFF_METHOD method) {
  TRACE("Div( Vector3D, Field3D)");
  return Div(v, f, method, outloc);
}

const Field3D Div(const Vector3D &v, const Field3D &f) {
  TRACE("Div( Vector3D, Field3D)");
  return Div(v, f, DIFF_DEFAULT, CELL_DEFAULT);
}

/**************************************************************************
 * Curl operators
 **************************************************************************/

const Vector2D Curl(const Vector2D &v) {

  TRACE("Curl( Vector2D )");

  ASSERT1(v.getLocation() != CELL_VSHIFT);
  Mesh *localmesh = v.x.getMesh();
  auto metric = v.x.getCoordinates();

  // Get covariant components of v
  Vector2D vco = v;
  vco.toCovariant();

  // get components (curl(v))^j
  Vector2D result(localmesh);
  result.x = (DDY(vco.z) - DDZ(vco.y)) / metric->J;
  result.y = (DDZ(vco.x) - DDX(vco.z)) / metric->J;
  result.z = (DDX(vco.y) - DDY(vco.x)) / metric->J;

  /// Coordinate torsion
  result.z -= metric->ShiftTorsion * vco.z / metric->J;

  result.setLocation(v.getLocation());

  result.covariant = false; // result is contravariant

  return result;
}

const Vector3D Curl(const Vector3D &v) {
  TRACE("Curl( Vector3D )");

  ASSERT1(v.getLocation() != CELL_VSHIFT);

  Mesh *localmesh = v.x.getMesh();
  auto metric = v.x.getCoordinates();

  // Get covariant components of v
  Vector3D vco = v;
  vco.toCovariant();

  // get components (curl(v))^j
  Vector3D result(localmesh);
  result.x = (DDY(vco.z) - DDZ(vco.y)) / metric->J;
  result.y = (DDZ(vco.x) - DDX(vco.z)) / metric->J;
  result.z = (DDX(vco.y) - DDY(vco.x)) / metric->J;

  // Coordinate torsion
  result.z -= metric->ShiftTorsion * vco.z / metric->J;

  result.setLocation(v.getLocation());

  result.covariant = false; // result is contravariant

  return result;
}

/**************************************************************************
 * Upwinding operators
 **************************************************************************/

const Field2D V_dot_Grad(const Vector2D &v, const Field2D &f) {
  TRACE("V_dot_Grad( Vector2D , Field2D )");

  Field2D result(f.getMesh());

  // Get contravariant components of v
  Vector2D vcn = v;
  vcn.toContravariant();

  result = VDDX(vcn.x, f) + VDDY(vcn.y, f) + VDDZ(vcn.z, f);

  return result;
}

const Field3D V_dot_Grad(const Vector2D &v, const Field3D &f) {
  TRACE("V_dot_Grad( Vector2D , Field3D )");

  Field3D result(f.getMesh());

  // Get contravariant components of v
  Vector2D vcn = v;
  vcn.toContravariant();

  result = VDDX(vcn.x, f) + VDDY(vcn.y, f) + VDDZ(vcn.z, f);

  return result;
}

const Field3D V_dot_Grad(const Vector3D &v, const Field2D &f) {
  TRACE("V_dot_Grad( Vector3D , Field2D )");

  Field3D result(f.getMesh());

  // Get contravariant components of v
  Vector3D vcn = v;
  vcn.toContravariant();

  result = VDDX(vcn.x, f) + VDDY(vcn.y, f) + VDDZ(vcn.z, f);

  return result;
}

const Field3D V_dot_Grad(const Vector3D &v, const Field3D &f) {
  TRACE("V_dot_Grad( Vector3D , Field3D )");

  Field3D result(f.getMesh());

  // Get contravariant components of v
  Vector3D vcn = v;
  vcn.toContravariant();
  
  result = VDDX(vcn.x, f) + VDDY(vcn.y, f) + VDDZ(vcn.z, f);

  return result;
}

const Vector2D V_dot_Grad(const Vector2D &v, const Vector2D &a) {
  TRACE("V_dot_Grad( Vector2D , Vector2D )");

  ASSERT1(v.getLocation() == a.getLocation());
  ASSERT1(v.getLocation() != CELL_VSHIFT);

  Vector2D result{v.x.getMesh()};

  auto metric = v.x.getCoordinates();

  Vector2D vcn = v;
  vcn.toContravariant();

  if (a.covariant) {

    result.x = VDDX(vcn.x, a.x) + VDDY(vcn.y, a.x) + VDDZ(vcn.z, a.x);
    result.x -= vcn.x * (metric->G1_11 * a.x + metric->G2_11 * a.y + metric->G3_11 * a.z);
    result.x -= vcn.y * (metric->G1_12 * a.x + metric->G2_12 * a.y + metric->G3_12 * a.z);
    result.x -= vcn.z * (metric->G1_13 * a.x + metric->G2_13 * a.y + metric->G3_13 * a.z);

    result.y = VDDX(vcn.x, a.y) + VDDY(vcn.y, a.y) + VDDZ(vcn.z, a.y);
    result.y -= vcn.x * (metric->G1_12 * a.x + metric->G2_12 * a.y + metric->G3_12 * a.z);
    result.y -= vcn.y * (metric->G1_22 * a.x + metric->G2_22 * a.y + metric->G3_22 * a.z);
    result.y -= vcn.z * (metric->G1_23 * a.x + metric->G2_23 * a.y + metric->G3_23 * a.z);

    result.z = VDDX(vcn.x, a.z) + VDDY(vcn.y, a.z) + VDDZ(vcn.z, a.z);
    result.z -= vcn.x * (metric->G1_13 * a.x + metric->G2_13 * a.y + metric->G3_13 * a.z);
    result.z -= vcn.y * (metric->G1_23 * a.x + metric->G2_23 * a.y + metric->G3_23 * a.z);
    result.z -= vcn.z * (metric->G1_33 * a.x + metric->G2_33 * a.y + metric->G3_33 * a.z);

    result.covariant = true;
  } else {

    result.x = VDDX(vcn.x, a.x) + VDDY(vcn.y, a.x) + VDDZ(vcn.z, a.x);
    result.x += vcn.x * (metric->G1_11 * a.x + metric->G1_12 * a.y + metric->G1_13 * a.z);
    result.x += vcn.y * (metric->G1_12 * a.x + metric->G1_22 * a.y + metric->G1_23 * a.z);
    result.x += vcn.z * (metric->G1_13 * a.x + metric->G1_23 * a.y + metric->G1_33 * a.z);

    result.y = VDDX(vcn.x, a.y) + VDDY(vcn.y, a.y) + VDDZ(vcn.z, a.y);
    result.y += vcn.x * (metric->G2_11 * a.x + metric->G2_12 * a.y + metric->G2_13 * a.z);
    result.y += vcn.y * (metric->G2_12 * a.x + metric->G2_22 * a.y + metric->G2_23 * a.z);
    result.y += vcn.z * (metric->G2_13 * a.x + metric->G2_23 * a.y + metric->G2_33 * a.z);

    result.z = VDDX(vcn.x, a.z) + VDDY(vcn.y, a.z) + VDDZ(vcn.z, a.z);
    result.z += vcn.x * (metric->G3_11 * a.x + metric->G3_12 * a.y + metric->G3_13 * a.z);
    result.z += vcn.y * (metric->G3_12 * a.x + metric->G3_22 * a.y + metric->G3_23 * a.z);
    result.z += vcn.z * (metric->G3_13 * a.x + metric->G3_23 * a.y + metric->G3_33 * a.z);

    result.covariant = false;
  }

  result.setLocation(v.getLocation());

  return result;
}

const Vector3D V_dot_Grad(const Vector2D &v, const Vector3D &a) {
  TRACE("V_dot_Grad( Vector2D , Vector3D )");

  ASSERT1(v.getLocation() == a.getLocation());
  ASSERT1(v.getLocation() != CELL_VSHIFT);

  Vector3D result{v.x.getMesh()};

  auto metric = v.x.getCoordinates();

  Vector2D vcn = v;
  vcn.toContravariant();

  if (a.covariant) {
    result.x = VDDX(vcn.x, a.x) + VDDY(vcn.y, a.x) + VDDZ(vcn.z, a.x);
    result.x -= vcn.x * (metric->G1_11 * a.x + metric->G2_11 * a.y + metric->G3_11 * a.z);
    result.x -= vcn.y * (metric->G1_12 * a.x + metric->G2_12 * a.y + metric->G3_12 * a.z);
    result.x -= vcn.z * (metric->G1_13 * a.x + metric->G2_13 * a.y + metric->G3_13 * a.z);

    result.y = VDDX(vcn.x, a.y) + VDDY(vcn.y, a.y) + VDDZ(vcn.z, a.y);
    result.y -= vcn.x * (metric->G1_12 * a.x + metric->G2_12 * a.y + metric->G3_12 * a.z);
    result.y -= vcn.y * (metric->G1_22 * a.x + metric->G2_22 * a.y + metric->G3_22 * a.z);
    result.y -= vcn.z * (metric->G1_23 * a.x + metric->G2_23 * a.y + metric->G3_23 * a.z);

    result.z = VDDX(vcn.x, a.z) + VDDY(vcn.y, a.z) + VDDZ(vcn.z, a.z);
    result.z -= vcn.x * (metric->G1_13 * a.x + metric->G2_13 * a.y + metric->G3_13 * a.z);
    result.z -= vcn.y * (metric->G1_23 * a.x + metric->G2_23 * a.y + metric->G3_23 * a.z);
    result.z -= vcn.z * (metric->G1_33 * a.x + metric->G2_33 * a.y + metric->G3_33 * a.z);

    result.covariant = true;
  } else {
    result.x = VDDX(vcn.x, a.x) + VDDY(vcn.y, a.x) + VDDZ(vcn.z, a.x);
    result.x += vcn.x * (metric->G1_11 * a.x + metric->G1_12 * a.y + metric->G1_13 * a.z);
    result.x += vcn.y * (metric->G1_12 * a.x + metric->G1_22 * a.y + metric->G1_23 * a.z);
    result.x += vcn.z * (metric->G1_13 * a.x + metric->G1_23 * a.y + metric->G1_33 * a.z);

    result.y = VDDX(vcn.x, a.y) + VDDY(vcn.y, a.y) + VDDZ(vcn.z, a.y);
    result.y += vcn.x * (metric->G2_11 * a.x + metric->G2_12 * a.y + metric->G2_13 * a.z);
    result.y += vcn.y * (metric->G2_12 * a.x + metric->G2_22 * a.y + metric->G2_23 * a.z);
    result.y += vcn.z * (metric->G2_13 * a.x + metric->G2_23 * a.y + metric->G2_33 * a.z);

    result.z = VDDX(vcn.x, a.z) + VDDY(vcn.y, a.z) + VDDZ(vcn.z, a.z);
    result.z += vcn.x * (metric->G3_11 * a.x + metric->G3_12 * a.y + metric->G3_13 * a.z);
    result.z += vcn.y * (metric->G3_12 * a.x + metric->G3_22 * a.y + metric->G3_23 * a.z);
    result.z += vcn.z * (metric->G3_13 * a.x + metric->G3_23 * a.y + metric->G3_33 * a.z);

    result.covariant = false;
  }

  result.setLocation(v.getLocation());

  return result;
}

const Vector3D V_dot_Grad(const Vector3D &v, const Vector2D &a) {
  TRACE("V_dot_Grad( Vector3D , Vector2D )");

  ASSERT1(v.getLocation() == a.getLocation());
  ASSERT1(v.getLocation() != CELL_VSHIFT);

  Vector3D result{v.x.getMesh()};

  auto metric = v.x.getCoordinates();

  Vector3D vcn = v;
  vcn.toContravariant();

  if (a.covariant) {
    result.x = VDDX(vcn.x, a.x) + VDDY(vcn.y, a.x) + VDDZ(vcn.z, a.x);
    result.x -= vcn.x * (metric->G1_11 * a.x + metric->G2_11 * a.y + metric->G3_11 * a.z);
    result.x -= vcn.y * (metric->G1_12 * a.x + metric->G2_12 * a.y + metric->G3_12 * a.z);
    result.x -= vcn.z * (metric->G1_13 * a.x + metric->G2_13 * a.y + metric->G3_13 * a.z);

    result.y = VDDX(vcn.x, a.y) + VDDY(vcn.y, a.y) + VDDZ(vcn.z, a.y);
    result.y -= vcn.x * (metric->G1_12 * a.x + metric->G2_12 * a.y + metric->G3_12 * a.z);
    result.y -= vcn.y * (metric->G1_22 * a.x + metric->G2_22 * a.y + metric->G3_22 * a.z);
    result.y -= vcn.z * (metric->G1_23 * a.x + metric->G2_23 * a.y + metric->G3_23 * a.z);

    result.z = VDDX(vcn.x, a.z) + VDDY(vcn.y, a.z) + VDDZ(vcn.z, a.z);
    result.z -= vcn.x * (metric->G1_13 * a.x + metric->G2_13 * a.y + metric->G3_13 * a.z);
    result.z -= vcn.y * (metric->G1_23 * a.x + metric->G2_23 * a.y + metric->G3_23 * a.z);
    result.z -= vcn.z * (metric->G1_33 * a.x + metric->G2_33 * a.y + metric->G3_33 * a.z);

    result.covariant = true;
  } else {
    result.x = VDDX(vcn.x, a.x) + VDDY(vcn.y, a.x) + VDDZ(vcn.z, a.x);
    result.x += vcn.x * (metric->G1_11 * a.x + metric->G1_12 * a.y + metric->G1_13 * a.z);
    result.x += vcn.y * (metric->G1_12 * a.x + metric->G1_22 * a.y + metric->G1_23 * a.z);
    result.x += vcn.z * (metric->G1_13 * a.x + metric->G1_23 * a.y + metric->G1_33 * a.z);

    result.y = VDDX(vcn.x, a.y) + VDDY(vcn.y, a.y) + VDDZ(vcn.z, a.y);
    result.y += vcn.x * (metric->G2_11 * a.x + metric->G2_12 * a.y + metric->G2_13 * a.z);
    result.y += vcn.y * (metric->G2_12 * a.x + metric->G2_22 * a.y + metric->G2_23 * a.z);
    result.y += vcn.z * (metric->G2_13 * a.x + metric->G2_23 * a.y + metric->G2_33 * a.z);

    result.z = VDDX(vcn.x, a.z) + VDDY(vcn.y, a.z) + VDDZ(vcn.z, a.z);
    result.z += vcn.x * (metric->G3_11 * a.x + metric->G3_12 * a.y + metric->G3_13 * a.z);
    result.z += vcn.y * (metric->G3_12 * a.x + metric->G3_22 * a.y + metric->G3_23 * a.z);
    result.z += vcn.z * (metric->G3_13 * a.x + metric->G3_23 * a.y + metric->G3_33 * a.z);

    result.covariant = false;
  }

  result.setLocation(v.getLocation());

  return result;
}

const Vector3D V_dot_Grad(const Vector3D &v, const Vector3D &a) {
  TRACE("V_dot_Grad( Vector3D , Vector3D )");

  ASSERT1(v.getLocation() == a.getLocation());
  ASSERT1(v.getLocation() != CELL_VSHIFT);

  Vector3D result{v.x.getMesh()};

  auto metric = v.x.getCoordinates();

  Vector3D vcn = v;
  vcn.toContravariant();

  if (a.covariant) {
    result.x = VDDX(vcn.x, a.x) + VDDY(vcn.y, a.x) + VDDZ(vcn.z, a.x);
    result.x -= vcn.x * (metric->G1_11 * a.x + metric->G2_11 * a.y + metric->G3_11 * a.z);
    result.x -= vcn.y * (metric->G1_12 * a.x + metric->G2_12 * a.y + metric->G3_12 * a.z);
    result.x -= vcn.z * (metric->G1_13 * a.x + metric->G2_13 * a.y + metric->G3_13 * a.z);

    result.y = VDDX(vcn.x, a.y) + VDDY(vcn.y, a.y) + VDDZ(vcn.z, a.y);
    result.y -= vcn.x * (metric->G1_12 * a.x + metric->G2_12 * a.y + metric->G3_12 * a.z);
    result.y -= vcn.y * (metric->G1_22 * a.x + metric->G2_22 * a.y + metric->G3_22 * a.z);
    result.y -= vcn.z * (metric->G1_23 * a.x + metric->G2_23 * a.y + metric->G3_23 * a.z);

    result.z = VDDX(vcn.x, a.z) + VDDY(vcn.y, a.z) + VDDZ(vcn.z, a.z);
    result.z -= vcn.x * (metric->G1_13 * a.x + metric->G2_13 * a.y + metric->G3_13 * a.z);
    result.z -= vcn.y * (metric->G1_23 * a.x + metric->G2_23 * a.y + metric->G3_23 * a.z);
    result.z -= vcn.z * (metric->G1_33 * a.x + metric->G2_33 * a.y + metric->G3_33 * a.z);

    result.covariant = true;
  } else {
    result.x = VDDX(vcn.x, a.x) + VDDY(vcn.y, a.x) + VDDZ(vcn.z, a.x);
    result.x += vcn.x * (metric->G1_11 * a.x + metric->G1_12 * a.y + metric->G1_13 * a.z);
    result.x += vcn.y * (metric->G1_12 * a.x + metric->G1_22 * a.y + metric->G1_23 * a.z);
    result.x += vcn.z * (metric->G1_13 * a.x + metric->G1_23 * a.y + metric->G1_33 * a.z);

    result.y = VDDX(vcn.x, a.y) + VDDY(vcn.y, a.y) + VDDZ(vcn.z, a.y);
    result.y += vcn.x * (metric->G2_11 * a.x + metric->G2_12 * a.y + metric->G2_13 * a.z);
    result.y += vcn.y * (metric->G2_12 * a.x + metric->G2_22 * a.y + metric->G2_23 * a.z);
    result.y += vcn.z * (metric->G2_13 * a.x + metric->G2_23 * a.y + metric->G2_33 * a.z);

    result.z = VDDX(vcn.x, a.z) + VDDY(vcn.y, a.z) + VDDZ(vcn.z, a.z);
    result.z += vcn.x * (metric->G3_11 * a.x + metric->G3_12 * a.y + metric->G3_13 * a.z);
    result.z += vcn.y * (metric->G3_12 * a.x + metric->G3_22 * a.y + metric->G3_23 * a.z);
    result.z += vcn.z * (metric->G3_13 * a.x + metric->G3_23 * a.y + metric->G3_33 * a.z);

    result.covariant = false;
  }

  result.setLocation(v.getLocation());

  return result;
}
