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

#include <bout/mesh.hxx>
#include <bout/scorepwrapper.hxx>

#include <globals.hxx>
#include <vecops.hxx>
#include <derivs.hxx>
#include <msg_stack.hxx>
#include <unused.hxx>
#include <utils.hxx>

/**************************************************************************
 * Gradient operators
 **************************************************************************/

const Vector2D Grad(const Field2D& f, CELL_LOC outloc, const std::string& method) {
  TRACE("Grad( Field2D )");
  SCOREP0();
  CELL_LOC outloc_x, outloc_y, outloc_z;
  if (outloc == CELL_VSHIFT) {
    outloc_x = CELL_XLOW;
    outloc_y = CELL_YLOW;
    outloc_z = CELL_ZLOW;
  } else {
    outloc_x = outloc_y = outloc_z = outloc;
  }

  Vector2D result(f.getMesh());

  result.x = DDX(f, outloc_x, method);
  result.y = DDY(f, outloc_y, method);
  result.z = DDZ(f, outloc_z, method);

  if (outloc == CELL_DEFAULT) {
    result.setLocation(result.x.getLocation());
  } else {
    result.setLocation(outloc);
  }

  result.covariant = true;

  return result;
}

const Vector3D Grad(const Field3D &f, CELL_LOC outloc, const std::string& method) {
  TRACE("Grad( Field3D )");
  SCOREP0();
  CELL_LOC outloc_x, outloc_y, outloc_z;
  if (outloc == CELL_VSHIFT) {
    outloc_x = CELL_XLOW;
    outloc_y = CELL_YLOW;
    outloc_z = CELL_ZLOW;
  } else {
    outloc_x = outloc_y = outloc_z = outloc;
  }

  Vector3D result(f.getMesh());

  result.x = DDX(f, outloc_x, method);
  result.y = DDY(f, outloc_y, method);
  result.z = DDZ(f, outloc_z, method);

  if (outloc == CELL_DEFAULT) {
    result.setLocation(result.x.getLocation());
  } else {
    result.setLocation(outloc);
  }

  result.covariant = true;

  return result;
}

const Vector3D Grad_perp(const Field3D &f, CELL_LOC outloc, const std::string& method) {
  TRACE("Grad_perp( Field3D )");
  SCOREP0();
  ASSERT1(outloc == CELL_DEFAULT || outloc == f.getLocation());

  Coordinates *metric = f.getCoordinates(outloc);

  Vector3D result(f.getMesh());

  result.x = DDX(f, outloc, method)
             - metric->g_12 * DDY(f, outloc, method) / SQ(metric->J * metric->Bxy);
  result.y = 0.0;
  result.z = DDZ(f, outloc, method)
             - metric->g_23 * DDY(f, outloc, method) / SQ(metric->J * metric->Bxy);

  result.setLocation(result.x.getLocation());

  result.covariant = true;

  return result;
}

/**************************************************************************
 * Divergence operators
 **************************************************************************/

const Coordinates::metric_field_type Div(const Vector2D &v, CELL_LOC outloc, const std::string& method) {
  TRACE("Div( Vector2D )");
  SCOREP0();
  if (outloc == CELL_DEFAULT) {
    outloc = v.getLocation();
  }

  ASSERT1(outloc != CELL_VSHIFT);

  Mesh *localmesh = v.x.getMesh();

  Coordinates *metric = localmesh->getCoordinates(outloc);

  // get contravariant components of v
  Vector2D vcn = v;
  vcn.toContravariant();
  
  Coordinates::metric_field_type result = DDX(metric->J*vcn.x, outloc, method);
  result += DDY(metric->J*vcn.y, outloc, method);
  result += DDZ(metric->J*vcn.z, outloc, method);
  result /= metric->J;

  return result;
}

const Field3D Div(const Vector3D& v, CELL_LOC outloc, const std::string& method) {
  TRACE("Div( Vector3D )");
  SCOREP0();
  if (outloc == CELL_DEFAULT) {
    outloc = v.getLocation();
  }

  // This also catches the combination of v at VSHIFT and outloc at DEFAULT
  ASSERT1(outloc != CELL_VSHIFT);

  Mesh *localmesh = v.x.getMesh();

  Coordinates *metric = localmesh->getCoordinates(outloc);

  // get contravariant components of v
  Vector3D vcn = v;
  vcn.toContravariant();

  Field3D result = 0.0;
  
  if(!vcn.y.hasParallelSlices()){
    Field3D vcnJy = vcn.y.getCoordinates()->J * vcn.y;
    localmesh->communicate(vcnJy);
    result = DDY(vcnJy, outloc, method);
  }else{
    result = DDY(vcn.y.getCoordinates()->J * vcn.y, outloc, method);
  }

  result += DDX(vcn.x.getCoordinates()->J * vcn.x, outloc, method);
  result += DDZ(vcn.z.getCoordinates()->J * vcn.z, outloc, method);
  result /= metric->J;

  return result;
}

/**************************************************************************
 * Divergence operators for flux methods
 **************************************************************************/

const Coordinates::metric_field_type Div(const Vector2D &v, const Field2D &f, CELL_LOC outloc,
                                         const std::string& method) {
  TRACE("Div( Vector2D, Field2D )");
  SCOREP0();
  if (outloc == CELL_DEFAULT) {
    outloc = v.getLocation();
  }

  ASSERT1(outloc != CELL_VSHIFT);

  Mesh *localmesh = f.getMesh();

  Coordinates *metric = localmesh->getCoordinates(outloc);

  // get contravariant components of v
  Vector2D vcn = v;
  vcn.toContravariant();

  Coordinates::metric_field_type result = FDDX(vcn.x.getCoordinates()->J * vcn.x, f, outloc, method);
  result += FDDY(vcn.y.getCoordinates()->J * vcn.y, f, outloc, method);
  result += FDDZ(vcn.z.getCoordinates()->J * vcn.z, f, outloc, method);
  result /= metric->J;

  return result;
}

const Field3D Div(const Vector3D& v, const Field3D& f, CELL_LOC outloc,
                  const std::string& method) {
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

  Field3D result = FDDX(vcn.x.getCoordinates()->J * vcn.x, f, outloc, method);
  result += FDDY(vcn.y.getCoordinates()->J * vcn.y, f, outloc, method);
  result += FDDZ(vcn.z.getCoordinates()->J * vcn.z, f, outloc, method);
  result /= metric->J;

  return result;
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
  SCOREP0();
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
const Coordinates::metric_field_type V_dot_Grad(const Vector2D &v, const Field2D &f) {
  TRACE("V_dot_Grad( Vector2D , Field2D )");
  SCOREP0();

  // Get contravariant components of v
  auto vcn = v;
  vcn.toContravariant();
  
  return VDDX(vcn.x, f) + VDDY(vcn.y, f) + VDDZ(vcn.z, f);
}

const Field3D V_dot_Grad(const Vector2D &v, const Field3D &f) {
  TRACE("V_dot_Grad( Vector2D , Field3D )");
  SCOREP0();

  // Get contravariant components of v
  auto vcn = v;
  vcn.toContravariant();
  
  return VDDX(vcn.x, f) + VDDY(vcn.y, f) + VDDZ(vcn.z, f);
}

const Field3D V_dot_Grad(const Vector3D &v, const Field2D &f) {
  TRACE("V_dot_Grad( Vector3D , Field2D )");
  SCOREP0();

  // Get contravariant components of v
  auto vcn = v;
  vcn.toContravariant();
  
  return VDDX(vcn.x, f) + VDDY(vcn.y, f) + VDDZ(vcn.z, f);
}

const Field3D V_dot_Grad(const Vector3D &v, const Field3D &f) {
  TRACE("V_dot_Grad( Vector3D , Field3D )");
  SCOREP0();

  // Get contravariant components of v
  auto vcn = v;
  vcn.toContravariant();
  
  return VDDX(vcn.x, f) + VDDY(vcn.y, f) + VDDZ(vcn.z, f);
}

// Here R is the deduced return type based on a promoting
// operation (addition) between the two input types.
template<typename T, typename F, typename R = decltype(T{}+F{})>
R V_dot_Grad(const T &v, const F &a) {
  AUTO_TRACE();
  SCOREP0();
  ASSERT1(v.getLocation() == a.getLocation());
  ASSERT1(v.getLocation() != CELL_VSHIFT);

  // Note by default R will describe a const vector type. By using
  // the following form of declaring result we ignore the const
  // qualifier here but keep it on the return type in the function
  // signature.
  auto result = R{v.x.getMesh()};

  auto metric = v.x.getCoordinates();

  auto vcn = v;
  vcn.toContravariant();

   if (a.covariant) {
    result.x = VDDX(vcn.x, a.x) + VDDY(vcn.y, a.x) + VDDZ(vcn.z, a.x);
    BOUT_FOR(i, result.x.getRegion("RGN_ALL")) {
      result.x[i] -= vcn.x[i] * (metric->G1_11[i] * a.x[i] + metric->G2_11[i] * a.y[i] + metric->G3_11[i] * a.z[i]);
      result.x[i] -= vcn.y[i] * (metric->G1_12[i] * a.x[i] + metric->G2_12[i] * a.y[i] + metric->G3_12[i] * a.z[i]);
      result.x[i] -= vcn.z[i] * (metric->G1_13[i] * a.x[i] + metric->G2_13[i] * a.y[i] + metric->G3_13[i] * a.z[i]);
    }
    
    result.y = VDDX(vcn.x, a.y) + VDDY(vcn.y, a.y) + VDDZ(vcn.z, a.y);
    BOUT_FOR(i, result.y.getRegion("RGN_ALL")) {
      result.y[i] -= vcn.x[i] * (metric->G1_12[i] * a.x[i] + metric->G2_12[i] * a.y[i] + metric->G3_12[i] * a.z[i]);
      result.y[i] -= vcn.y[i] * (metric->G1_22[i] * a.x[i] + metric->G2_22[i] * a.y[i] + metric->G3_22[i] * a.z[i]);
      result.y[i] -= vcn.z[i] * (metric->G1_23[i] * a.x[i] + metric->G2_23[i] * a.y[i] + metric->G3_23[i] * a.z[i]);
    }
    
    result.z = VDDX(vcn.x, a.z) + VDDY(vcn.y, a.z) + VDDZ(vcn.z, a.z);
    BOUT_FOR(i, result.z.getRegion("RGN_ALL")) {
      result.z[i] -= vcn.x[i] * (metric->G1_13[i] * a.x[i] + metric->G2_13[i] * a.y[i] + metric->G3_13[i] * a.z[i]);
      result.z[i] -= vcn.y[i] * (metric->G1_23[i] * a.x[i] + metric->G2_23[i] * a.y[i] + metric->G3_23[i] * a.z[i]);
      result.z[i] -= vcn.z[i] * (metric->G1_33[i] * a.x[i] + metric->G2_33[i] * a.y[i] + metric->G3_33[i] * a.z[i]);
    }
    result.covariant = true;
  } else {
    result.x = VDDX(vcn.x, a.x) + VDDY(vcn.y, a.x) + VDDZ(vcn.z, a.x);
    BOUT_FOR(i, result.x.getRegion("RGN_ALL")) {    
      result.x[i] += vcn.x[i] * (metric->G1_11[i] * a.x[i] + metric->G1_12[i] * a.y[i] + metric->G1_13[i] * a.z[i]);
      result.x[i] += vcn.y[i] * (metric->G1_12[i] * a.x[i] + metric->G1_22[i] * a.y[i] + metric->G1_23[i] * a.z[i]);
      result.x[i] += vcn.z[i] * (metric->G1_13[i] * a.x[i] + metric->G1_23[i] * a.y[i] + metric->G1_33[i] * a.z[i]);
    }
    
    result.y = VDDX(vcn.x, a.y) + VDDY(vcn.y, a.y) + VDDZ(vcn.z, a.y);
    BOUT_FOR(i, result.y.getRegion("RGN_ALL")) {    
      result.y[i] += vcn.x[i] * (metric->G2_11[i] * a.x[i] + metric->G2_12[i] * a.y[i] + metric->G2_13[i] * a.z[i]);
      result.y[i] += vcn.y[i] * (metric->G2_12[i] * a.x[i] + metric->G2_22[i] * a.y[i] + metric->G2_23[i] * a.z[i]);
      result.y[i] += vcn.z[i] * (metric->G2_13[i] * a.x[i] + metric->G2_23[i] * a.y[i] + metric->G2_33[i] * a.z[i]);
    }
    
    result.z = VDDX(vcn.x, a.z) + VDDY(vcn.y, a.z) + VDDZ(vcn.z, a.z);
    BOUT_FOR(i, result.z.getRegion("RGN_ALL")) {
      result.z[i] += vcn.x[i] * (metric->G3_11[i] * a.x[i] + metric->G3_12[i] * a.y[i] + metric->G3_13[i] * a.z[i]);
      result.z[i] += vcn.y[i] * (metric->G3_12[i] * a.x[i] + metric->G3_22[i] * a.y[i] + metric->G3_23[i] * a.z[i]);
      result.z[i] += vcn.z[i] * (metric->G3_13[i] * a.x[i] + metric->G3_23[i] * a.y[i] + metric->G3_33[i] * a.z[i]);
    }
    
    result.covariant = false;
  }

  result.setLocation(v.getLocation());

  return result;
  
};

// Implement vector-vector operation in terms of templated routine above
const Vector2D V_dot_Grad(const Vector2D &v, const Vector2D &a) {
  return V_dot_Grad<Vector2D, Vector2D>(v, a);
}
const Vector3D V_dot_Grad(const Vector2D &v, const Vector3D &a) {
  return V_dot_Grad<Vector2D, Vector3D>(v, a);
}
const Vector3D V_dot_Grad(const Vector3D &v, const Vector2D &a) {
  return V_dot_Grad<Vector3D, Vector2D>(v, a);
}
const Vector3D V_dot_Grad(const Vector3D &v, const Vector3D &a) {
  return V_dot_Grad<Vector3D, Vector3D>(v, a);
}
