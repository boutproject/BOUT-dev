/**************************************************************************
 * Describes coordinate systems
 *
 * ChangeLog
 * =========
 * 
 * 2014-11-10 Ben Dudson <bd512@york.ac.uk>
 *    * Created by separating metric from Mesh
 *
 * 
 **************************************************************************
 * Copyright 2014 B.D.Dudson
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

class Coordinates;

#ifndef __COORDINATES_H__
#define __COORDINATES_H__

#define METRIC_FIELD_TYPE Field2D
//#define METRIC_FIELD_TYPE Field3D

#include "mesh.hxx"
#include "datafile.hxx"
#include "utils.hxx"
#include <bout_types.hxx>

/*!
 * Represents a coordinate system, and associated operators
 *
 * This is a container for a collection of metric tensor components
 */ 
class Coordinates {
public:
  using metric_type = METRIC_FIELD_TYPE;

  /// Standard constructor from input
  Coordinates(Mesh *mesh);

  /// Constructor interpolating from another Coordinates object
  Coordinates(Mesh *mesh, const CELL_LOC loc, const Coordinates* coords_in);
  
  ~Coordinates() {}
  
  /*!
   * Adds variables to the output file, for post-processing
   * 
   * Must be a better way so that Coordinates doesn't depend on Datafile
   */
  void outputVars(Datafile &file);

  metric_type dx, dy; ///< Mesh spacing in x and y
  BoutReal dz; ///< Mesh spacing in Z

  BoutReal zlength() const { return dz * nz; } ///< Length of the Z domain. Used for FFTs

  /// True if corrections for non-uniform mesh spacing should be included in operators
  bool non_uniform;
  metric_type d1_dx,
      d1_dy; ///< 2nd-order correction for non-uniform meshes d/di(1/dx) and d/di(1/dy)

  metric_type J; ///< Coordinate system Jacobian, so volume of cell is J*dx*dy*dz

  metric_type Bxy; ///< Magnitude of B = nabla z times nabla x

  /// Contravariant metric tensor (g^{ij})
  metric_type g11, g22, g33, g12, g13, g23;

  /// Covariant metric tensor
  metric_type g_11, g_22, g_33, g_12, g_13, g_23;

  /// Christoffel symbol of the second kind (connection coefficients)
  metric_type G1_11, G1_22, G1_33, G1_12, G1_13, G1_23;
  metric_type G2_11, G2_22, G2_33, G2_12, G2_13, G2_23;
  metric_type G3_11, G3_22, G3_33, G3_12, G3_13, G3_23;

  metric_type G1, G2, G3;

  metric_type
      ShiftTorsion; ///< d pitch angle / dx. Needed for vector differentials (Curl)

  metric_type IntShiftTorsion; ///< Integrated shear (I in BOUT notation)

  /// Calculate differential geometry quantities from the metric tensor
  int geometry();
  int calcCovariant(); ///< Inverts contravatiant metric to get covariant
  int calcContravariant(); ///< Invert covariant metric to get contravariant
  int jacobian(); ///< Calculate J and Bxy

  // Operators

  const metric_type DDX(const Field2D& f, CELL_LOC outloc = CELL_DEFAULT,
                        const std::string& method = "DEFAULT",
                        REGION region = RGN_NOBNDRY);
  const metric_type DDX(const Field2D& f, CELL_LOC outloc, DIFF_METHOD method,
                        REGION region = RGN_NOBNDRY) {
    return DDX(f, outloc, DIFF_METHOD_STRING(method), region);
  };

  const metric_type DDY(const Field2D& f, CELL_LOC outloc = CELL_DEFAULT,
                        const std::string& method = "DEFAULT",
                        REGION region = RGN_NOBNDRY);
  const metric_type DDY(const Field2D& f, CELL_LOC outloc, DIFF_METHOD method,
                        REGION region = RGN_NOBNDRY) {
    return DDY(f, outloc, DIFF_METHOD_STRING(method), region);
  };

  const metric_type DDZ(const Field2D& f, CELL_LOC outloc = CELL_DEFAULT,
                        const std::string& method = "DEFAULT",
                        REGION region = RGN_NOBNDRY);
  const metric_type DDZ(const Field2D& f, CELL_LOC outloc, DIFF_METHOD method,
                        REGION region = RGN_NOBNDRY) {
    return DDZ(f, outloc, DIFF_METHOD_STRING(method), region);
  };

  /// Gradient along magnetic field  b.Grad(f)
  const metric_type Grad_par(const Field2D& var, CELL_LOC outloc = CELL_DEFAULT,
                             const std::string& method = "DEFAULT");
  const metric_type Grad_par(const Field2D& var, CELL_LOC outloc, DIFF_METHOD method) {
    return Grad_par(var, outloc, DIFF_METHOD_STRING(method));
  };

  const Field3D Grad_par(const Field3D& var, CELL_LOC outloc = CELL_DEFAULT,
                         const std::string& method = "DEFAULT");
  const Field3D Grad_par(const Field3D& var, CELL_LOC outloc, DIFF_METHOD method) {
    return Grad_par(var, outloc, DIFF_METHOD_STRING(method));
  };

  /// Advection along magnetic field V*b.Grad(f)
  const metric_type Vpar_Grad_par(const Field2D& v, const Field2D& f,
                                  CELL_LOC outloc = CELL_DEFAULT,
                                  const std::string& method = "DEFAULT");
  const metric_type Vpar_Grad_par(const Field2D& v, const Field2D& f, CELL_LOC outloc,
                                  DIFF_METHOD method) {
    return Vpar_Grad_par(v, f, outloc, DIFF_METHOD_STRING(method));
  };

  const Field3D Vpar_Grad_par(const Field3D& v, const Field3D& f,
                              CELL_LOC outloc = CELL_DEFAULT,
                              const std::string& method = "DEFAULT");
  const Field3D Vpar_Grad_par(const Field3D& v, const Field3D& f, CELL_LOC outloc,
                              DIFF_METHOD method) {
    return Vpar_Grad_par(v, f, outloc, DIFF_METHOD_STRING(method));
  };

  /// Divergence along magnetic field  Div(b*f) = B.Grad(f/B)
  const metric_type Div_par(const Field2D& f, CELL_LOC outloc = CELL_DEFAULT,
                            const std::string& method = "DEFAULT");
  const metric_type Div_par(const Field2D& f, CELL_LOC outloc, DIFF_METHOD method) {
    return Div_par(f, outloc, DIFF_METHOD_STRING(method));
  };

  const Field3D Div_par(const Field3D& f, CELL_LOC outloc = CELL_DEFAULT,
                        const std::string& method = "DEFAULT");
  const Field3D Div_par(const Field3D& f, CELL_LOC outloc, DIFF_METHOD method) {
    return Div_par(f, outloc, DIFF_METHOD_STRING(method));
  };

  // Second derivative along magnetic field
  const metric_type Grad2_par2(const Field2D& f, CELL_LOC outloc = CELL_DEFAULT,
                               const std::string& method = "DEFAULT");
  const metric_type Grad2_par2(const Field2D& f, CELL_LOC outloc, DIFF_METHOD method) {
    return Grad2_par2(f, outloc, DIFF_METHOD_STRING(method));
  };

  const Field3D Grad2_par2(const Field3D& f, CELL_LOC outloc = CELL_DEFAULT,
                           const std::string& method = "DEFAULT");
  const Field3D Grad2_par2(const Field3D& f, CELL_LOC outloc, DIFF_METHOD method) {
    return Grad2_par2(f, outloc, DIFF_METHOD_STRING(method));
  };

  // Perpendicular Laplacian operator, using only X-Z derivatives
  // NOTE: This might be better bundled with the Laplacian inversion code
  // since it makes use of the same coefficients and FFT routines
  const metric_type Delp2(const Field2D& f, CELL_LOC outloc = CELL_DEFAULT);
  const Field3D Delp2(const Field3D &f, CELL_LOC outloc=CELL_DEFAULT);
  const FieldPerp Delp2(const FieldPerp &f, CELL_LOC outloc=CELL_DEFAULT);
  
  // Full parallel Laplacian operator on scalar field
  // Laplace_par(f) = Div( b (b dot Grad(f)) )
  const metric_type Laplace_par(const Field2D& f, CELL_LOC outloc = CELL_DEFAULT);
  const Field3D Laplace_par(const Field3D &f, CELL_LOC outloc=CELL_DEFAULT);
  
  // Full Laplacian operator on scalar field
  const metric_type Laplace(const Field2D& f, CELL_LOC outloc = CELL_DEFAULT);
  const Field3D Laplace(const Field3D &f, CELL_LOC outloc=CELL_DEFAULT);
  
private:
  int nz; // Size of mesh in Z. This is mesh->ngz-1
  Mesh * localmesh;
  CELL_LOC location;
};

/*
/// Standard coordinate system for tokamak simulations
class TokamakCoordinates : public Coordinates {
public:
  TokamakCoordinates(Mesh *mesh) : Coordinates(mesh) {
    
  }
private:
  
};
*/

#endif // __COORDINATES_H__
