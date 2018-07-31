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

#include "mesh.hxx"
#include "datafile.hxx"
#include "utils.hxx"
#include <bout_types.hxx>
#include <flexible.hxx>
/*!
 * Represents a coordinate system, and associated operators
 *
 * This is a container for a collection of metric tensor components
 */ 
class Coordinates {
public:
  /// Constructor
  Coordinates(Mesh *mesh);
  
  ~Coordinates() {}
  
  /*!
   * Adds variables to the output file, for post-processing
   * 
   * Must be a better way so that Coordinates doesn't depend on Datafile
   */
  void outputVars(Datafile &file);
  
  Flexible<Field2D> dx, dy; ///< Mesh spacing in x and y
  BoutReal dz; ///< Mesh spacing in Z

  BoutReal zlength() const { return dz * nz; } ///< Length of the Z domain. Used for FFTs

  /// True if corrections for non-uniform mesh spacing should be included in operators
  bool non_uniform;
  Flexible<Field2D> d1_dx, d1_dy;  ///< 2nd-order correction for non-uniform meshes d/di(1/dx) and d/di(1/dy)
  
  Flexible<Field2D> J; ///< Coordinate system Jacobian, so volume of cell is J*dx*dy*dz

  Flexible<Field2D> Bxy; ///< Magnitude of B = nabla z times nabla x
  
  /// Contravariant metric tensor (g^{ij})
  Flexible<Field2D> g11, g22, g33, g12, g13, g23;
  
  /// Covariant metric tensor
  Flexible<Field2D> g_11, g_22, g_33, g_12, g_13, g_23;
  
  /// Christoffel symbol of the second kind (connection coefficients)
  Flexible<Field2D> G1_11, G1_22, G1_33, G1_12, G1_13, G1_23;
  Flexible<Field2D> G2_11, G2_22, G2_33, G2_12, G2_13, G2_23;
  Flexible<Field2D> G3_11, G3_22, G3_33, G3_12, G3_13, G3_23;
  
  Flexible<Field2D> G1, G2, G3;
  
  Flexible<Field2D> ShiftTorsion; ///< d pitch angle / dx. Needed for vector differentials (Curl)

  Flexible<Field2D> IntShiftTorsion; ///< Integrated shear (I in BOUT notation)

  /// Calculate differential geometry quantities from the metric tensor
  int geometry();
  int calcCovariant(); ///< Inverts contravatiant metric to get covariant
  int calcContravariant(); ///< Invert covariant metric to get contravariant
  int jacobian(); ///< Calculate J and Bxy
  int setBxyBoundaries(); ///< If mesh->StaggerGrids==true, calculate CELL_XLOW and CELL_YLOW fields of Bxy and set their boundary cells

  // Operators

  const Field2D DDX(const Field2D &f, CELL_LOC outloc = CELL_DEFAULT,
                    DIFF_METHOD method = DIFF_DEFAULT,
                    REGION region = RGN_NOBNDRY);
  const Field2D DDY(const Field2D &f, CELL_LOC outloc = CELL_DEFAULT,
                    DIFF_METHOD method = DIFF_DEFAULT,
                    REGION region = RGN_NOBNDRY);
  const Field2D DDZ(const Field2D &f, CELL_LOC outloc = CELL_DEFAULT,
                    DIFF_METHOD method = DIFF_DEFAULT,
                    REGION region = RGN_NOBNDRY);
  
  /// Gradient along magnetic field  b.Grad(f)
  const Field2D Grad_par(const Field2D &var, CELL_LOC outloc=CELL_DEFAULT, DIFF_METHOD method=DIFF_DEFAULT);
  const Field3D Grad_par(const Field3D &var, CELL_LOC outloc=CELL_DEFAULT, DIFF_METHOD method=DIFF_DEFAULT);
  
  /// Advection along magnetic field V*b.Grad(f)
  const Field2D Vpar_Grad_par(const Field2D &v, const Field2D &f, CELL_LOC outloc=CELL_DEFAULT, DIFF_METHOD method=DIFF_DEFAULT);
  const Field3D Vpar_Grad_par(const Field3D &v, const Field3D &f, CELL_LOC outloc=CELL_DEFAULT, DIFF_METHOD method=DIFF_DEFAULT);
  
  /// Divergence along magnetic field  Div(b*f) = B.Grad(f/B)
  const Field2D Div_par(const Field2D &f, CELL_LOC outloc=CELL_DEFAULT, DIFF_METHOD method=DIFF_DEFAULT);
  const Field3D Div_par(const Field3D &f, CELL_LOC outloc=CELL_DEFAULT, DIFF_METHOD method=DIFF_DEFAULT);
  
  // Second derivative along magnetic field
  const Field2D Grad2_par2(const Field2D &f, CELL_LOC outloc=CELL_DEFAULT);
  const Field3D Grad2_par2(const Field3D &f, CELL_LOC outloc=CELL_DEFAULT);

  // Perpendicular Laplacian operator, using only X-Z derivatives
  // NOTE: This might be better bundled with the Laplacian inversion code
  // since it makes use of the same coefficients and FFT routines
  const Field2D Delp2(const Field2D &f);
  const Field3D Delp2(const Field3D &f);
  const FieldPerp Delp2(const FieldPerp &f);
  
  // Full parallel Laplacian operator on scalar field
  // Laplace_par(f) = Div( b (b dot Grad(f)) ) 
  const Field2D Laplace_par(const Field2D &f, CELL_LOC outloc = CELL_DEFAULT);
  const Field3D Laplace_par(const Field3D &f, CELL_LOC outloc = CELL_DEFAULT);
  
  // Full Laplacian operator on scalar field
  const Field2D Laplace(const Field2D &f, CELL_LOC outloc = CELL_DEFAULT);
  const Field3D Laplace(const Field3D &f, CELL_LOC outloc = CELL_DEFAULT);
  
private:
  int nz; // Size of mesh in Z. This is mesh->ngz-1
  Mesh * localmesh;
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
