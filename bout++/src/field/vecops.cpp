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

#include "globals.h"
#include "vecops.h"
#include "derivs.h"

/**************************************************************************
 * Gradient operators
 **************************************************************************/

const Vector2D Grad(const Field2D &f, CELL_LOC outloc)
{
  Vector2D result;

#ifdef CHECK
  int msg_pos = msg_stack.push("Grad( Field2D )");
#endif

  result.x = DDX(f);
  result.y = DDY(f);
  result.z = DDZ(f);

  result.covariant = true;

#ifdef CHECK
  msg_stack.pop(msg_pos);
#endif

  return result;
}

const Vector3D Grad(const Field3D &f, CELL_LOC outloc)
{
  Vector3D result;

#ifdef CHECK
  int msg_pos = msg_stack.push("Grad( Field3D )");
#endif

  if(outloc == CELL_DEFAULT)
    outloc = f.getLocation();

  if(outloc == CELL_VSHIFT) {
    // Each vector component is shifted
    result.x = DDX(f, CELL_XLOW);
    result.y = DDY(f, CELL_YLOW);
    result.z = DDZ(f, CELL_ZLOW);
  }else {
    result.x = DDX(f);
    result.y = DDY(f);
    result.z = DDZ(f);
  }

  result.covariant = true;
  
#ifdef CHECK
  msg_stack.pop(msg_pos);
#endif

  return result;
}

/**************************************************************************
 * Divergence operators
 **************************************************************************/

const Field2D Div(const Vector2D &v, CELL_LOC outloc)
{
  Field2D result;

#ifdef CHECK
  int msg_pos = msg_stack.push("Div( Vector2D )");
#endif
  
  // get contravariant components of v
  Vector2D vcn = v;
  vcn.to_contravariant();
  
  result = DDX(J*vcn.x);
  result += DDY(J*vcn.y);
  result += DDZ(J*vcn.z);
  result /= J;

#ifdef CHECK
  msg_stack.pop(msg_pos);
#endif

  return result;
}

const Field3D Div(const Vector3D &v, CELL_LOC outloc)
{
  Field3D result;

#ifdef CHECK
  int msg_pos = msg_stack.push("Div( Vector3D )");
#endif

  if(outloc == CELL_DEFAULT) 
    outloc = CELL_CENTRE;

  // get contravariant components of v
  Vector3D vcn = v;
  vcn.to_contravariant();
  
  result = DDX(J*vcn.x, outloc);
  result += DDY(J*vcn.y, outloc);
  result += DDZ(J*vcn.z, outloc);
  result /= J;

#ifdef CHECK
  msg_stack.pop(msg_pos);
#endif

  return result;
}

/**************************************************************************
 * Curl operators
 **************************************************************************/

const Vector2D Curl(const Vector2D &v)
{
  Vector2D result;

#ifdef CHECK
  int msg_pos = msg_stack.push("Curl( Vector2D )");
#endif

  // Get covariant components of v
  Vector2D vco = v;
  vco.to_covariant();

  // get components (curl(v))^j

  result.x = (DDY(vco.z) - DDZ(vco.y))/J;
  result.y = (DDZ(vco.x) - DDX(vco.z))/J;
  result.z = (DDX(vco.y) - DDY(vco.x))/J;
  
  if(ShiftXderivs) {
    result.z -= ShiftTorsion*vco.z / J;
  }

  result.covariant = false; // result is contravariant
#ifdef CHECK
  msg_stack.pop(msg_pos);
#endif
  return result;
}

const Vector3D Curl(const Vector3D &v)
{
  Vector3D result;

#ifdef CHECK
  int msg_pos = msg_stack.push("Curl( Vector3D )");
#endif

  // Get covariant components of v
  Vector3D vco = v;
  vco.to_covariant();

  // get components (curl(v))^j

  result.x = (DDY(vco.z) - DDZ(vco.y))/J;
  result.y = (DDZ(vco.x) - DDX(vco.z))/J;
  result.z = (DDX(vco.y) - DDY(vco.x))/J;

  if(ShiftXderivs) {
    result.z -= ShiftTorsion*vco.z / J;
  }

  result.covariant = false; // result is contravariant

#ifdef CHECK
  msg_stack.pop(msg_pos);
#endif

  return result;
}

/**************************************************************************
 * Upwinding operators
 **************************************************************************/

const Field2D V_dot_Grad(const Vector2D &v, const Field2D &f)
{
  Field2D result;
  
#ifdef CHECK
  int msg_pos = msg_stack.push("V_dot_Grad( Vector2D , Field2D )");
#endif

  // Get contravariant components of v
  Vector2D vcn = v;
  vcn.to_contravariant();

  result = VDDX(vcn.x, f) + VDDY(vcn.y, f) + VDDZ(vcn.z, f);

#ifdef CHECK
  msg_stack.pop(msg_pos);
#endif

  return result;
}

const Field3D V_dot_Grad(const Vector2D &v, const Field3D &f)
{
  Field3D result;
  
#ifdef CHECK
  int msg_pos = msg_stack.push("V_dot_Grad( Vector2D , Field3D )");
#endif

  // Get contravariant components of v
  Vector2D vcn = v;
  vcn.to_contravariant();

  result = VDDX(vcn.x, f) + VDDY(vcn.y, f) + VDDZ(vcn.z, f);

#ifdef CHECK
  msg_stack.pop(msg_pos);
#endif

  return result;
}

const Field3D V_dot_Grad(const Vector3D &v, const Field2D &f)
{
  Field3D result;
  
#ifdef CHECK
  int msg_pos = msg_stack.push("V_dot_Grad( Vector3D , Field2D )");
#endif

  // Get contravariant components of v
  Vector3D vcn = v;
  vcn.to_contravariant();

  result = VDDX(vcn.x, f) + VDDY(vcn.y, f) + VDDZ(vcn.z, f);

#ifdef CHECK
  msg_stack.pop(msg_pos);
#endif

  return result;
}

const Field3D V_dot_Grad(const Vector3D &v, const Field3D &f)
{
  Field3D result;
  
#ifdef CHECK
  int msg_pos = msg_stack.push("V_dot_Grad( Vector3D , Field3D )");
#endif

  // Get contravariant components of v
  Vector3D vcn = v;
  vcn.to_contravariant();
  
  result = VDDX(vcn.x, f) + VDDY(vcn.y, f) + VDDZ(vcn.z, f);

#ifdef CHECK
  msg_stack.pop(msg_pos);
#endif

  return result;
}

const Vector2D V_dot_Grad(const Vector2D &v, const Vector2D &a)
{
  Vector2D result;

#ifdef CHECK
  int msg_pos = msg_stack.push("V_dot_Grad( Vector2D , Vector2D )");
#endif

  Vector2D vcn = v;
  vcn.to_contravariant();

  if(a.covariant) {
    
    result.x = VDDX(vcn.x, a.x) + VDDY(vcn.y, a.x) + VDDZ(vcn.z, a.x);
    result.x -= vcn.x*(G1_11*a.x + G2_11*a.y + G3_11*a.z);
    result.x -= vcn.y*(G1_12*a.x + G2_12*a.y);
    result.x -= vcn.z*(G1_13*a.x + G3_13*a.z);

    result.y = VDDX(vcn.x, a.y) + VDDY(vcn.y, a.y) + VDDZ(vcn.z, a.y);
    result.y -= vcn.x*(G1_12*a.x + G2_12*a.y);
    result.y -= vcn.y*(G1_22*a.x + G2_22*a.y + G3_22*a.z);
    result.y -= vcn.z*(G2_23*a.y + G3_23*a.z);

    result.z = VDDX(vcn.x, a.z) + VDDY(vcn.y, a.z) + VDDZ(vcn.z, a.z);
    result.z -= vcn.x*(G1_13*a.x + G3_13*a.z);
    result.z -= vcn.y*(G2_23*a.y + G3_23*a.z);
    result.z -= vcn.z*(G1_33*a.x + G2_33*a.y + G3_33*a.z);

    result.covariant = true;
  }else {
    
    result.x = VDDX(vcn.x, a.x) + VDDY(vcn.y, a.x) + VDDZ(vcn.z, a.x);
    result.x += vcn.x*(G1_11*a.x + G1_12*a.y + G1_13*a.z);
    result.x += vcn.y*(G1_12*a.x + G1_22*a.y);
    result.x += vcn.z*(G1_13*a.x + G1_33*a.z);
    
    result.y = VDDX(vcn.x, a.y) + VDDY(vcn.y, a.y) + VDDZ(vcn.z, a.y);
    result.y += vcn.x*(G2_11*a.x + G2_12*a.y);
    result.y += vcn.y*(G2_12*a.x + G2_22*a.y + G2_23*a.z);
    result.y += vcn.z*(G2_23*a.y + G2_33*a.z);
    
    result.z = VDDX(vcn.x, a.z) + VDDY(vcn.y, a.z) + VDDZ(vcn.z, a.z);
    result.z += vcn.x*(G3_11*a.x + G3_13*a.z);
    result.z += vcn.y*(G3_22*a.y + G3_23*a.z);
    result.z += vcn.z*(G3_13*a.x + G3_23*a.y + G3_33*a.z);

    result.covariant = false;
  }
  
#ifdef CHECK
  msg_stack.pop(msg_pos);
#endif

  return result;
}

const Vector3D V_dot_Grad(const Vector2D &v, const Vector3D &a)
{
  Vector3D result;

#ifdef CHECK
  int msg_pos = msg_stack.push("V_dot_Grad( Vector2D , Vector3D )");
#endif

  Vector2D vcn = v;
  vcn.to_contravariant();

  if(a.covariant) {
    result.x = VDDX(vcn.x, a.x) + VDDY(vcn.y, a.x) + VDDZ(vcn.z, a.x);
    result.x -= vcn.x*(G1_11*a.x + G2_11*a.y + G3_11*a.z);
    result.x -= vcn.y*(G1_12*a.x + G2_12*a.y);
    result.x -= vcn.z*(G1_13*a.x + G3_13*a.z);

    result.y = VDDX(vcn.x, a.y) + VDDY(vcn.y, a.y) + VDDZ(vcn.z, a.y);
    result.y -= vcn.x*(G1_12*a.x + G2_12*a.y);
    result.y -= vcn.y*(G1_22*a.x + G2_22*a.y + G3_22*a.z);
    result.y -= vcn.z*(G2_23*a.y + G3_23*a.z);

    result.z = VDDX(vcn.x, a.z) + VDDY(vcn.y, a.z) + VDDZ(vcn.z, a.z);
    result.z -= vcn.x*(G1_13*a.x + G3_13*a.z);
    result.z -= vcn.y*(G2_23*a.y + G3_23*a.z);
    result.z -= vcn.z*(G1_33*a.x + G2_33*a.y + G3_33*a.z);

    result.covariant = true;
  }else {
    result.x = VDDX(vcn.x, a.x) + VDDY(vcn.y, a.x) + VDDZ(vcn.z, a.x);
    result.x += vcn.x*(G1_11*a.x + G1_12*a.y + G1_13*a.z);
    result.x += vcn.y*(G1_12*a.x + G1_22*a.y);
    result.x += vcn.z*(G1_13*a.x + G1_33*a.z);

    result.y = VDDX(vcn.x, a.y) + VDDY(vcn.y, a.y) + VDDZ(vcn.z, a.y);
    result.y += vcn.x*(G2_11*a.x + G2_12*a.y);
    result.y += vcn.y*(G2_12*a.x + G2_22*a.y + G2_23*a.z);
    result.y += vcn.z*(G2_23*a.y + G2_33*a.z);

    result.z = VDDX(vcn.x, a.z) + VDDY(vcn.y, a.z) + VDDZ(vcn.z, a.z);
    result.z += vcn.x*(G3_11*a.x + G3_13*a.z);
    result.z += vcn.y*(G3_22*a.y + G3_23*a.z);
    result.z += vcn.z*(G3_13*a.x + G3_23*a.y + G3_33*a.z);

    result.covariant = false;
  }
  
#ifdef CHECK
  msg_stack.pop(msg_pos);
#endif

  return result;
}

const Vector3D V_dot_Grad(const Vector3D &v, const Vector2D &a)
{
  Vector3D result;
  
#ifdef CHECK
  int msg_pos = msg_stack.push("V_dot_Grad( Vector3D , Vector2D )");
#endif

  Vector3D vcn = v;
  vcn.to_contravariant();

  if(a.covariant) {
    result.x = VDDX(vcn.x, a.x) + VDDY(vcn.y, a.x) + VDDZ(vcn.z, a.x);
    result.x -= vcn.x*(G1_11*a.x + G2_11*a.y + G3_11*a.z);
    result.x -= vcn.y*(G1_12*a.x + G2_12*a.y);
    result.x -= vcn.z*(G1_13*a.x + G3_13*a.z);

    result.y = VDDX(vcn.x, a.y) + VDDY(vcn.y, a.y) + VDDZ(vcn.z, a.y);
    result.y -= vcn.x*(G1_12*a.x + G2_12*a.y);
    result.y -= vcn.y*(G1_22*a.x + G2_22*a.y + G3_22*a.z);
    result.y -= vcn.z*(G2_23*a.y + G3_23*a.z);

    result.z = VDDX(vcn.x, a.z) + VDDY(vcn.y, a.z) + VDDZ(vcn.z, a.z);
    result.z -= vcn.x*(G1_13*a.x + G3_13*a.z);
    result.z -= vcn.y*(G2_23*a.y + G3_23*a.z);
    result.z -= vcn.z*(G1_33*a.x + G2_33*a.y + G3_33*a.z);

    result.covariant = true;
  }else {
    result.x = VDDX(vcn.x, a.x) + VDDY(vcn.y, a.x) + VDDZ(vcn.z, a.x);
    result.x += vcn.x*(G1_11*a.x + G1_12*a.y + G1_13*a.z);
    result.x += vcn.y*(G1_12*a.x + G1_22*a.y);
    result.x += vcn.z*(G1_13*a.x + G1_33*a.z);

    result.y = VDDX(vcn.x, a.y) + VDDY(vcn.y, a.y) + VDDZ(vcn.z, a.y);
    result.y += vcn.x*(G2_11*a.x + G2_12*a.y);
    result.y += vcn.y*(G2_12*a.x + G2_22*a.y + G2_23*a.z);
    result.y += vcn.z*(G2_23*a.y + G2_33*a.z);

    result.z = VDDX(vcn.x, a.z) + VDDY(vcn.y, a.z) + VDDZ(vcn.z, a.z);
    result.z += vcn.x*(G3_11*a.x + G3_13*a.z);
    result.z += vcn.y*(G3_22*a.y + G3_23*a.z);
    result.z += vcn.z*(G3_13*a.x + G3_23*a.y + G3_33*a.z);

    result.covariant = false;
  }
  
#ifdef CHECK
  msg_stack.pop(msg_pos);
#endif

  return result;
}

const Vector3D V_dot_Grad(const Vector3D &v, const Vector3D &a)
{
  Vector3D result;

#ifdef CHECK
  int msg_pos = msg_stack.push("V_dot_Grad( Vector3D , Vector3D )");
#endif

  Vector3D vcn = v;
  vcn.to_contravariant();

  if(a.covariant) {
    result.x = VDDX(vcn.x, a.x) + VDDY(vcn.y, a.x) + VDDZ(vcn.z, a.x);
    result.x -= vcn.x*(G1_11*a.x + G2_11*a.y + G3_11*a.z);
    result.x -= vcn.y*(G1_12*a.x + G2_12*a.y);
    result.x -= vcn.z*(G1_13*a.x + G3_13*a.z);
    
    result.y = VDDX(vcn.x, a.y) + VDDY(vcn.y, a.y) + VDDZ(vcn.z, a.y);
    result.y -= vcn.x*(G1_12*a.x + G2_12*a.y);
    result.y -= vcn.y*(G1_22*a.x + G2_22*a.y + G3_22*a.z);
    result.y -= vcn.z*(G2_23*a.y + G3_23*a.z);

    result.z = VDDX(vcn.x, a.z) + VDDY(vcn.y, a.z) + VDDZ(vcn.z, a.z);
    result.z -= vcn.x*(G1_13*a.x + G3_13*a.z);
    result.z -= vcn.y*(G2_23*a.y + G3_23*a.z);
    result.z -= vcn.z*(G1_33*a.x + G2_33*a.y + G3_33*a.z);

    result.covariant = true;
  }else {
    result.x = VDDX(vcn.x, a.x) + VDDY(vcn.y, a.x) + VDDZ(vcn.z, a.x);
    result.x += vcn.x*(G1_11*a.x + G1_12*a.y + G1_13*a.z);
    result.x += vcn.y*(G1_12*a.x + G1_22*a.y);
    result.x += vcn.z*(G1_13*a.x + G1_33*a.z);
    
    result.y = VDDX(vcn.x, a.y) + VDDY(vcn.y, a.y) + VDDZ(vcn.z, a.y);
    result.y += vcn.x*(G2_11*a.x + G2_12*a.y);
    result.y += vcn.y*(G2_12*a.x + G2_22*a.y + G2_23*a.z);
    result.y += vcn.z*(G2_23*a.y + G2_33*a.z);
    
    result.z = VDDX(vcn.x, a.z) + VDDY(vcn.y, a.z) + VDDZ(vcn.z, a.z);
    result.z += vcn.x*(G3_11*a.x + G3_13*a.z);
    result.z += vcn.y*(G3_22*a.y + G3_23*a.z);
    result.z += vcn.z*(G3_13*a.x + G3_23*a.y + G3_33*a.z);

    result.covariant = false;
  }
  
#ifdef CHECK
  msg_stack.pop(msg_pos);
#endif

  return result;
}


