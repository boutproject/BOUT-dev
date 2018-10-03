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
  Vector2D result(f.getMesh());

  TRACE("Grad( Field2D )");
  
  if (outloc == CELL_DEFAULT) {
    outloc = f.getLocation();
  }

  result.x = DDX(f, outloc);
  result.y = DDY(f, outloc);
  result.z = DDZ(f, outloc);

  result.covariant = true;

  return result;
}

const Vector3D Grad(const Field3D &f, CELL_LOC outloc_x, CELL_LOC outloc_y,
                    CELL_LOC outloc_z) {
  Vector3D result(f.getMesh());

  TRACE("Grad( Field3D )");

  if (outloc_x == CELL_DEFAULT) {
    outloc_x = f.getLocation();
  }
  if (outloc_y == CELL_DEFAULT) {
    outloc_y = f.getLocation();
  }
  if (outloc_z == CELL_DEFAULT) {
    outloc_z = f.getLocation();
  }

  result.x = DDX(f, outloc_x);
  result.y = DDY(f, outloc_y);
  result.z = DDZ(f, outloc_z);

  result.covariant = true;

  return result;
}

const Vector3D Grad(const Field3D &f, CELL_LOC outloc) {
  if (outloc == CELL_VSHIFT) {
    return Grad(f, CELL_XLOW, CELL_YLOW, CELL_ZLOW);
  }

  if (outloc == CELL_DEFAULT) {
    outloc = f.getLocation();
  }

  return Grad(f, outloc, outloc, outloc);
}

const Vector3D Grad_perp(const Field3D &f, CELL_LOC outloc_x, CELL_LOC UNUSED(outloc_y),
                         CELL_LOC outloc_z) {
  Vector3D result(f.getMesh());

  TRACE("Grad_perp( Field3D )");

  if (outloc_x == CELL_DEFAULT) {
    outloc_x = f.getLocation();
  }
  if (outloc_z == CELL_DEFAULT) {
    outloc_z = f.getLocation();
  }

  Coordinates *metric_x = mesh->coordinates(outloc_x);
  Coordinates *metric_z = mesh->coordinates(outloc_z);

  result.x = DDX(f, outloc_x) -
             metric_x->g_12 * DDY(f, outloc_x) / SQ(metric_x->J * metric_x->Bxy);
  result.y = 0.0;
  result.z = DDZ(f, outloc_z) -
             metric_z->g_23 * DDY(f, outloc_z) / SQ(metric_z->J * metric_z->Bxy);

  result.covariant = true;

  return result;
}

/**************************************************************************
 * Divergence operators
 **************************************************************************/

const Field2D Div(const Vector2D &v, CELL_LOC outloc) {
  TRACE("Div( Vector2D )");

  Mesh *localmesh = v.x.getMesh();
  Field2D result(localmesh);

  if (outloc == CELL_DEFAULT) {
    outloc = v.getLocation();
  }

  Coordinates *metric = localmesh->coordinates(outloc);

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

  Mesh *localmesh = v.x.getMesh();
  Field3D result(localmesh);

  if (outloc == CELL_DEFAULT) {
    outloc = v.getLocation();
  }

  Coordinates *metric = localmesh->coordinates(outloc);

  // get contravariant components of v
  Vector3D vcn = v;
  vcn.toContravariant();

  result = DDX(metric->J * vcn.x, outloc);
  result += DDY(metric->J * vcn.y, outloc);
  result += DDZ(metric->J * vcn.z, outloc);
  result /= metric->J;

  return result;
}

/**************************************************************************
 * Divergence operators for flux methods
 **************************************************************************/

const Field2D Div(const Vector2D &v, const Field2D &f, CELL_LOC outloc) {
  TRACE("Div( Vector2D, Field2D )");

  ASSERT1(v.getLocation() == f.getLocation());

  Mesh *localmesh = f.getMesh();

  if (outloc == CELL_DEFAULT) {
    outloc = v.getLocation();
  }

  Coordinates *metric = localmesh->coordinates(outloc);

  // get contravariant components of v
  Vector2D vcn = v;
  vcn.toContravariant();

  Field2D result(localmesh);
  result = FDDX(metric->J * vcn.x, f, outloc);
  result += FDDY(metric->J * vcn.y, f, outloc);
  result += FDDZ(metric->J * vcn.z, f, outloc);
  result /= metric->J;

  return result;
}

const Field3D Div(const Vector3D &v, const Field3D &f, DIFF_METHOD method,
                  CELL_LOC outloc) {
  TRACE("Div( Vector3D, Field3D )");

  Mesh *localmesh = f.getMesh();
  Field3D result(localmesh);

  if (outloc == CELL_DEFAULT) {
    outloc = v.getLocation();
  }

  Coordinates *metric = localmesh->coordinates(outloc);

  // get contravariant components of v
  Vector3D vcn = v;
  vcn.toContravariant();

  result = FDDX(metric->J * vcn.x, f, outloc, method);
  result += FDDY(metric->J * vcn.y, f, outloc, method);
  result += FDDZ(metric->J * vcn.z, f, outloc, method);
  result /= metric->J;

  return result;
}

const Field3D Div(const Vector3D &v, const Field3D &f, CELL_LOC outloc, DIFF_METHOD method) {
  return Div(v, f, method, outloc);
}

const Field3D Div(const Vector3D &v, const Field3D &f) {
  return Div(v, f, DIFF_DEFAULT, CELL_DEFAULT);
}

/**************************************************************************
 * Curl operators
 **************************************************************************/

const Vector2D Curl(const Vector2D &v, CELL_LOC outloc) {

  TRACE("Curl( Vector2D )");

  Mesh *localmesh = v.x.getMesh();

  if (outloc == CELL_DEFAULT) {
    outloc = v.getLocation();
  }

  Coordinates *metric = localmesh->coordinates(outloc);

  // Get covariant components of v
  Vector2D vco = v;
  vco.toCovariant();

  // get components (curl(v))^j
  Vector2D result(localmesh);
  result.x = (DDY(vco.z, outloc) - DDZ(vco.y, outloc)) / metric->J;
  result.y = (DDZ(vco.x, outloc) - DDX(vco.z, outloc)) / metric->J;
  result.z = (DDX(vco.y, outloc) - DDY(vco.x, outloc)) / metric->J;

  /// Coordinate torsion
  result.z -= metric->ShiftTorsion * vco.z / metric->J;

  result.covariant = false; // result is contravariant

  return result;
}

const Vector3D Curl(const Vector3D &v, CELL_LOC outloc_x, CELL_LOC outloc_y,
                    CELL_LOC outloc_z) {

  TRACE("Curl( Vector3D )");

  Mesh *localmesh = v.x.getMesh();

  if (outloc_x == CELL_DEFAULT) {
    outloc_x = v.getLocation();
  }
  if (outloc_y == CELL_DEFAULT) {
    outloc_y = v.getLocation();
  }
  if (outloc_z == CELL_DEFAULT) {
    outloc_z = v.getLocation();
  }

  Coordinates* metric_z = localmesh->coordinates(outloc_z);

  // Get covariant components of v
  Vector3D vco = v;
  vco.toCovariant();

  // get components (curl(v))^j
  Vector3D result(localmesh);
  result.x = (DDY(vco.z, outloc_x) - DDZ(vco.y, outloc_x))/localmesh->coordinates(outloc_x)->J;
  result.y = (DDZ(vco.x, outloc_y) - DDX(vco.z, outloc_y))/localmesh->coordinates(outloc_y)->J;
  result.z = (DDX(vco.y, outloc_z) - DDY(vco.x, outloc_z))/metric_z->J;

  // Coordinate torsion
  result.z -= metric_z->ShiftTorsion*vco.z / metric_z->J;

  result.covariant = false; // result is contravariant

  return result;
}

const Vector3D Curl(const Vector3D &v, CELL_LOC outloc) {
  if (outloc == CELL_VSHIFT) {
    return Curl(v, CELL_XLOW, CELL_YLOW, CELL_ZLOW);
  }

  if (outloc == CELL_DEFAULT) {
    outloc = v.getLocation();
  }

  return Curl(v, outloc, outloc, outloc);
}

/**************************************************************************
 * Upwinding operators
 **************************************************************************/

const Field2D V_dot_Grad(const Vector2D &v, const Field2D &f) {
  Field2D result(f.getMesh());

  TRACE("V_dot_Grad( Vector2D , Field2D )");

  // Get contravariant components of v
  Vector2D vcn = v;
  vcn.toContravariant();

  result = VDDX(vcn.x, f) + VDDY(vcn.y, f) + VDDZ(vcn.z, f);

  return result;
}

const Field3D V_dot_Grad(const Vector2D &v, const Field3D &f) {
  Field3D result(f.getMesh());

  TRACE("V_dot_Grad( Vector2D , Field3D )");

  // Get contravariant components of v
  Vector2D vcn = v;
  vcn.toContravariant();

  result = VDDX(vcn.x, f) + VDDY(vcn.y, f) + VDDZ(vcn.z, f);

  return result;
}

const Field3D V_dot_Grad(const Vector3D &v, const Field2D &f) {
  Field3D result(f.getMesh());

  TRACE("V_dot_Grad( Vector3D , Field2D )");

  // Get contravariant components of v
  Vector3D vcn = v;
  vcn.toContravariant();

  result = VDDX(vcn.x, f) + VDDY(vcn.y, f) + VDDZ(vcn.z, f);

  return result;
}

const Field3D V_dot_Grad(const Vector3D &v, const Field3D &f) {
  Field3D result(f.getMesh());

  TRACE("V_dot_Grad( Vector3D , Field3D )");

  // Get contravariant components of v
  Vector3D vcn = v;
  vcn.toContravariant();
  
  result = VDDX(vcn.x, f) + VDDY(vcn.y, f) + VDDZ(vcn.z, f);

  return result;
}

const Vector2D V_dot_Grad(const Vector2D &v, const Vector2D &a, const CELL_LOC outloc) {
  TRACE("V_dot_Grad( Vector2D , Vector2D )");

  Mesh *localmesh = v.x.getMesh();
  Vector2D result(localmesh);

  Coordinates *metric = localmesh->coordinates(outloc);

  Vector2D vcn = v;
  vcn.toContravariant();

  if(a.covariant) {
    
    result.x = VDDX(vcn.x, a.x, outloc) + VDDY(vcn.y, a.x, outloc) + VDDZ(vcn.z, a.x, outloc);
    result.x -= vcn.x*(metric->G1_11*a.x + metric->G2_11*a.y + metric->G3_11*a.z);
    result.x -= vcn.y*(metric->G1_12*a.x + metric->G2_12*a.y + metric->G3_12*a.z);
    result.x -= vcn.z*(metric->G1_13*a.x + metric->G2_13*a.y + metric->G3_13*a.z);

    result.y = VDDX(vcn.x, a.y, outloc) + VDDY(vcn.y, a.y, outloc) + VDDZ(vcn.z, a.y, outloc);
    result.y -= vcn.x*(metric->G1_12*a.x + metric->G2_12*a.y + metric->G3_12*a.z);
    result.y -= vcn.y*(metric->G1_22*a.x + metric->G2_22*a.y + metric->G3_22*a.z);
    result.y -= vcn.z*(metric->G1_23*a.x + metric->G2_23*a.y + metric->G3_23*a.z);

    result.z = VDDX(vcn.x, a.z, outloc) + VDDY(vcn.y, a.z, outloc) + VDDZ(vcn.z, a.z, outloc);
    result.z -= vcn.x*(metric->G1_13*a.x + metric->G2_13*a.y + metric->G3_13*a.z);
    result.z -= vcn.y*(metric->G1_23*a.x + metric->G2_23*a.y + metric->G3_23*a.z);
    result.z -= vcn.z*(metric->G1_33*a.x + metric->G2_33*a.y + metric->G3_33*a.z);

    result.covariant = true;
  }else {
    
    result.x = VDDX(vcn.x, a.x, outloc) + VDDY(vcn.y, a.x, outloc) + VDDZ(vcn.z, a.x, outloc);
    result.x += vcn.x*(metric->G1_11*a.x + metric->G1_12*a.y + metric->G1_13*a.z);
    result.x += vcn.y*(metric->G1_12*a.x + metric->G1_22*a.y + metric->G1_23*a.z);
    result.x += vcn.z*(metric->G1_13*a.x + metric->G1_23*a.y + metric->G1_33*a.z);

    result.y = VDDX(vcn.x, a.y, outloc) + VDDY(vcn.y, a.y, outloc) + VDDZ(vcn.z, a.y, outloc);
    result.y += vcn.x*(metric->G2_11*a.x + metric->G2_12*a.y + metric->G2_13*a.z);
    result.y += vcn.y*(metric->G2_12*a.x + metric->G2_22*a.y + metric->G2_23*a.z);
    result.y += vcn.z*(metric->G2_13*a.x + metric->G2_23*a.y + metric->G2_33*a.z);

    result.z = VDDX(vcn.x, a.z, outloc) + VDDY(vcn.y, a.z, outloc) + VDDZ(vcn.z, a.z, outloc);
    result.z += vcn.x*(metric->G3_11*a.x + metric->G3_12*a.y + metric->G3_13*a.z);
    result.z += vcn.y*(metric->G3_12*a.x + metric->G3_22*a.y + metric->G3_23*a.z);
    result.z += vcn.z*(metric->G3_13*a.x + metric->G3_23*a.y + metric->G3_33*a.z);

    result.covariant = false;
  }

  return result;
}

const Vector3D V_dot_Grad(const Vector2D &v, const Vector3D &a, const CELL_LOC outloc) {
  Mesh *localmesh = v.x.getMesh();
  Vector3D result(localmesh);

  TRACE("V_dot_Grad( Vector2D , Vector3D )");

  Coordinates *metric = localmesh->coordinates(outloc);

  Vector2D vcn = v;
  vcn.toContravariant();

  if(a.covariant) {
    result.x = VDDX(vcn.x, a.x, outloc) + VDDY(vcn.y, a.x, outloc) + VDDZ(vcn.z, a.x, outloc);
    result.x -= vcn.x*(metric->G1_11*a.x + metric->G2_11*a.y + metric->G3_11*a.z);
    result.x -= vcn.y*(metric->G1_12*a.x + metric->G2_12*a.y + metric->G3_12*a.z);
    result.x -= vcn.z*(metric->G1_13*a.x + metric->G2_13*a.y + metric->G3_13*a.z);

    result.y = VDDX(vcn.x, a.y, outloc) + VDDY(vcn.y, a.y, outloc) + VDDZ(vcn.z, a.y, outloc);
    result.y -= vcn.x*(metric->G1_12*a.x + metric->G2_12*a.y + metric->G3_12*a.z);
    result.y -= vcn.y*(metric->G1_22*a.x + metric->G2_22*a.y + metric->G3_22*a.z);
    result.y -= vcn.z*(metric->G1_23*a.x + metric->G2_23*a.y + metric->G3_23*a.z);

    result.z = VDDX(vcn.x, a.z, outloc) + VDDY(vcn.y, a.z, outloc) + VDDZ(vcn.z, a.z, outloc);
    result.z -= vcn.x*(metric->G1_13*a.x + metric->G2_13*a.y + metric->G3_13*a.z);
    result.z -= vcn.y*(metric->G1_23*a.x + metric->G2_23*a.y + metric->G3_23*a.z);
    result.z -= vcn.z*(metric->G1_33*a.x + metric->G2_33*a.y + metric->G3_33*a.z);

    result.covariant = true;
  }else {
    result.x = VDDX(vcn.x, a.x, outloc) + VDDY(vcn.y, a.x, outloc) + VDDZ(vcn.z, a.x, outloc);
    result.x += vcn.x*(metric->G1_11*a.x + metric->G1_12*a.y + metric->G1_13*a.z);
    result.x += vcn.y*(metric->G1_12*a.x + metric->G1_22*a.y + metric->G1_23*a.z);
    result.x += vcn.z*(metric->G1_13*a.x + metric->G1_23*a.y + metric->G1_33*a.z);

    result.y = VDDX(vcn.x, a.y, outloc) + VDDY(vcn.y, a.y, outloc) + VDDZ(vcn.z, a.y, outloc);
    result.y += vcn.x*(metric->G2_11*a.x + metric->G2_12*a.y + metric->G2_13*a.z);
    result.y += vcn.y*(metric->G2_12*a.x + metric->G2_22*a.y + metric->G2_23*a.z);
    result.y += vcn.z*(metric->G2_13*a.x + metric->G2_23*a.y + metric->G2_33*a.z);

    result.z = VDDX(vcn.x, a.z, outloc) + VDDY(vcn.y, a.z, outloc) + VDDZ(vcn.z, a.z, outloc);
    result.z += vcn.x*(metric->G3_11*a.x + metric->G3_12*a.y + metric->G3_13*a.z);
    result.z += vcn.y*(metric->G3_12*a.x + metric->G3_22*a.y + metric->G3_23*a.z);
    result.z += vcn.z*(metric->G3_13*a.x + metric->G3_23*a.y + metric->G3_33*a.z);

    result.covariant = false;
  }

  return result;
}

const Vector3D V_dot_Grad(const Vector3D &v, const Vector2D &a, const CELL_LOC outloc) {
  Mesh *localmesh = v.x.getMesh();
  Vector3D result(localmesh);

  TRACE("V_dot_Grad( Vector3D , Vector2D )");

  Coordinates *metric = localmesh->coordinates(outloc);

  Vector3D vcn = v;
  vcn.toContravariant();

  if(a.covariant) {
    result.x = VDDX(vcn.x, a.x, outloc) + VDDY(vcn.y, a.x, outloc) + VDDZ(vcn.z, a.x, outloc);
    result.x -= vcn.x*(metric->G1_11*a.x + metric->G2_11*a.y + metric->G3_11*a.z);
    result.x -= vcn.y*(metric->G1_12*a.x + metric->G2_12*a.y + metric->G3_12*a.z);
    result.x -= vcn.z*(metric->G1_13*a.x + metric->G2_13*a.y + metric->G3_13*a.z);

    result.y = VDDX(vcn.x, a.y, outloc) + VDDY(vcn.y, a.y, outloc) + VDDZ(vcn.z, a.y, outloc);
    result.y -= vcn.x*(metric->G1_12*a.x + metric->G2_12*a.y + metric->G3_12*a.z);
    result.y -= vcn.y*(metric->G1_22*a.x + metric->G2_22*a.y + metric->G3_22*a.z);
    result.y -= vcn.z*(metric->G1_23*a.x + metric->G2_23*a.y + metric->G3_23*a.z);

    result.z = VDDX(vcn.x, a.z, outloc) + VDDY(vcn.y, a.z, outloc) + VDDZ(vcn.z, a.z, outloc);
    result.z -= vcn.x*(metric->G1_13*a.x + metric->G2_13*a.y + metric->G3_13*a.z);
    result.z -= vcn.y*(metric->G1_23*a.x + metric->G2_23*a.y + metric->G3_23*a.z);
    result.z -= vcn.z*(metric->G1_33*a.x + metric->G2_33*a.y + metric->G3_33*a.z);

    result.covariant = true;
  }else {
    result.x = VDDX(vcn.x, a.x, outloc) + VDDY(vcn.y, a.x, outloc) + VDDZ(vcn.z, a.x, outloc);
    result.x += vcn.x*(metric->G1_11*a.x + metric->G1_12*a.y + metric->G1_13*a.z);
    result.x += vcn.y*(metric->G1_12*a.x + metric->G1_22*a.y + metric->G1_23*a.z);
    result.x += vcn.z*(metric->G1_13*a.x + metric->G1_23*a.y + metric->G1_33*a.z);

    result.y = VDDX(vcn.x, a.y, outloc) + VDDY(vcn.y, a.y, outloc) + VDDZ(vcn.z, a.y, outloc);
    result.y += vcn.x*(metric->G2_11*a.x + metric->G2_12*a.y + metric->G2_13*a.z);
    result.y += vcn.y*(metric->G2_12*a.x + metric->G2_22*a.y + metric->G2_23*a.z);
    result.y += vcn.z*(metric->G2_13*a.x + metric->G2_23*a.y + metric->G2_33*a.z);

    result.z = VDDX(vcn.x, a.z, outloc) + VDDY(vcn.y, a.z, outloc) + VDDZ(vcn.z, a.z, outloc);
    result.z += vcn.x*(metric->G3_11*a.x + metric->G3_12*a.y + metric->G3_13*a.z);
    result.z += vcn.y*(metric->G3_12*a.x + metric->G3_22*a.y + metric->G3_23*a.z);
    result.z += vcn.z*(metric->G3_13*a.x + metric->G3_23*a.y + metric->G3_33*a.z);

    result.covariant = false;
  }

  return result;
}

const Vector3D V_dot_Grad(const Vector3D &v, const Vector3D &a, const CELL_LOC outloc) {
  Mesh *localmesh = v.x.getMesh();
  Vector3D result(localmesh);

  TRACE("V_dot_Grad( Vector3D , Vector3D )");

  Coordinates *metric = localmesh->coordinates(outloc);

  Vector3D vcn = v;
  vcn.toContravariant();

  if(a.covariant) {
    result.x = VDDX(vcn.x, a.x, outloc) + VDDY(vcn.y, a.x, outloc) + VDDZ(vcn.z, a.x, outloc);
    result.x -= vcn.x*(metric->G1_11*a.x + metric->G2_11*a.y + metric->G3_11*a.z);
    result.x -= vcn.y*(metric->G1_12*a.x + metric->G2_12*a.y + metric->G3_12*a.z);
    result.x -= vcn.z*(metric->G1_13*a.x + metric->G2_13*a.y + metric->G3_13*a.z);

    result.y = VDDX(vcn.x, a.y, outloc) + VDDY(vcn.y, a.y, outloc) + VDDZ(vcn.z, a.y, outloc);
    result.y -= vcn.x*(metric->G1_12*a.x + metric->G2_12*a.y + metric->G3_12*a.z);
    result.y -= vcn.y*(metric->G1_22*a.x + metric->G2_22*a.y + metric->G3_22*a.z);
    result.y -= vcn.z*(metric->G1_23*a.x + metric->G2_23*a.y + metric->G3_23*a.z);

    result.z = VDDX(vcn.x, a.z, outloc) + VDDY(vcn.y, a.z, outloc) + VDDZ(vcn.z, a.z, outloc);
    result.z -= vcn.x*(metric->G1_13*a.x + metric->G2_13*a.y + metric->G3_13*a.z);
    result.z -= vcn.y*(metric->G1_23*a.x + metric->G2_23*a.y + metric->G3_23*a.z);
    result.z -= vcn.z*(metric->G1_33*a.x + metric->G2_33*a.y + metric->G3_33*a.z);

    result.covariant = true;
  }else {
    result.x = VDDX(vcn.x, a.x, outloc) + VDDY(vcn.y, a.x, outloc) + VDDZ(vcn.z, a.x, outloc);
    result.x += vcn.x*(metric->G1_11*a.x + metric->G1_12*a.y + metric->G1_13*a.z);
    result.x += vcn.y*(metric->G1_12*a.x + metric->G1_22*a.y + metric->G1_23*a.z);
    result.x += vcn.z*(metric->G1_13*a.x + metric->G1_23*a.y + metric->G1_33*a.z);

    result.y = VDDX(vcn.x, a.y, outloc) + VDDY(vcn.y, a.y, outloc) + VDDZ(vcn.z, a.y, outloc);
    result.y += vcn.x*(metric->G2_11*a.x + metric->G2_12*a.y + metric->G2_13*a.z);
    result.y += vcn.y*(metric->G2_12*a.x + metric->G2_22*a.y + metric->G2_23*a.z);
    result.y += vcn.z*(metric->G2_13*a.x + metric->G2_23*a.y + metric->G2_33*a.z);

    result.z = VDDX(vcn.x, a.z, outloc) + VDDY(vcn.y, a.z, outloc) + VDDZ(vcn.z, a.z, outloc);
    result.z += vcn.x*(metric->G3_11*a.x + metric->G3_12*a.y + metric->G3_13*a.z);
    result.z += vcn.y*(metric->G3_12*a.x + metric->G3_22*a.y + metric->G3_23*a.z);
    result.z += vcn.z*(metric->G3_13*a.x + metric->G3_23*a.y + metric->G3_33*a.z);

    result.covariant = false;
  }

  return result;
}


