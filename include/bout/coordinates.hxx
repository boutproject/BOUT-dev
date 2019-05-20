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

#ifndef __COORDINATES_H__
#define __COORDINATES_H__

#include "bout/paralleltransform.hxx"
#include "datafile.hxx"
#include "utils.hxx"
#include <bout_types.hxx>
#include "field2d.hxx"
#include "field3d.hxx"

class Mesh;

/*!
 * Represents a coordinate system, and associated operators
 *
 * This is a container for a collection of metric tensor components
 */ 
class Coordinates {
public:
  /// Standard constructor from input
  Coordinates(Mesh *mesh, Options* options = nullptr);

  /// Constructor interpolating from another Coordinates object
  /// By default attempts to read staggered Coordinates from grid data source,
  /// interpolating from CELL_CENTRE if not present. Set
  /// force_interpolate_from_centre argument to true to always interpolate
  /// (useful if CELL_CENTRE Coordinates have been changed, so reading from file
  /// would not be correct).
  Coordinates(Mesh *mesh, Options* options, const CELL_LOC loc, const Coordinates* coords_in,
      bool force_interpolate_from_centre=false);

  /// A constructor useful for testing purposes. To use it, inherit
  /// from Coordinates. If \p calculate_geometry is true (default),
  /// calculate the non-uniform variables, Christoffel symbols
  Coordinates(Mesh* mesh, Field2D dx, Field2D dy, BoutReal dz, Field2D J, Field2D Bxy,
              Field2D g11, Field2D g22, Field2D g33, Field2D g12, Field2D g13,
              Field2D g23, Field2D g_11, Field2D g_22, Field2D g_33, Field2D g_12,
              Field2D g_13, Field2D g_23, Field2D ShiftTorsion, Field2D IntShiftTorsion,
              bool calculate_geometry = true);

  Coordinates& operator=(Coordinates&&) = default;

  ~Coordinates() = default;

  /*!
   * Adds variables to the output file, for post-processing
   * 
   * Must be a better way so that Coordinates doesn't depend on Datafile
   */
  void outputVars(Datafile &file);
  
  Field2D dx, dy; ///< Mesh spacing in x and y
  BoutReal dz; ///< Mesh spacing in Z

  BoutReal zlength() const { return dz * nz; } ///< Length of the Z domain. Used for FFTs

  /// True if corrections for non-uniform mesh spacing should be included in operators
  bool non_uniform;
  Field2D d1_dx, d1_dy;  ///< 2nd-order correction for non-uniform meshes d/di(1/dx) and d/di(1/dy)
  
  Field2D J; ///< Coordinate system Jacobian, so volume of cell is J*dx*dy*dz

  Field2D Bxy; ///< Magnitude of B = nabla z times nabla x
  
  /// Contravariant metric tensor (g^{ij})
  Field2D g11, g22, g33, g12, g13, g23;
  
  /// Covariant metric tensor
  Field2D g_11, g_22, g_33, g_12, g_13, g_23;
  
  /// Christoffel symbol of the second kind (connection coefficients)
  Field2D G1_11, G1_22, G1_33, G1_12, G1_13, G1_23;
  Field2D G2_11, G2_22, G2_33, G2_12, G2_13, G2_23;
  Field2D G3_11, G3_22, G3_33, G3_12, G3_13, G3_23;
  
  Field2D G1, G2, G3;
  
  Field2D ShiftTorsion; ///< d pitch angle / dx. Needed for vector differentials (Curl)

  Field2D IntShiftTorsion; ///< Integrated shear (I in BOUT notation)

  /// Calculate differential geometry quantities from the metric tensor
  int geometry(bool recalculate_staggered = true,
      bool force_interpolate_from_centre = false);
  int calcCovariant(); ///< Inverts contravatiant metric to get covariant
  int calcContravariant(); ///< Invert covariant metric to get contravariant
  int jacobian(); ///< Calculate J and Bxy


  ///////////////////////////////////////////////////////////
  // Parallel transforms
  ///////////////////////////////////////////////////////////

  /// Set the parallel (y) transform for this mesh.
  /// Mostly useful for tests.
  void setParallelTransform(std::unique_ptr<ParallelTransform> pt) {
    transform = std::move(pt);
  }

  /// Return the parallel transform
  ParallelTransform& getParallelTransform() {
    ASSERT1(transform != nullptr);
    return *transform;
  }

  ///////////////////////////////////////////////////////////
  // Operators
  ///////////////////////////////////////////////////////////

  const Field2D DDX(const Field2D& f, CELL_LOC outloc = CELL_DEFAULT,
                    const std::string& method = "DEFAULT", REGION region = RGN_NOBNDRY);
  const Field2D DDX(const Field2D& f, CELL_LOC outloc, DIFF_METHOD method,
                    REGION region = RGN_NOBNDRY) {
    return DDX(f, outloc, toString(method), region);
  };

  const Field2D DDY(const Field2D& f, CELL_LOC outloc = CELL_DEFAULT,
                    const std::string& method = "DEFAULT", REGION region = RGN_NOBNDRY);
  const Field2D DDY(const Field2D& f, CELL_LOC outloc, DIFF_METHOD method,
                    REGION region = RGN_NOBNDRY) {
    return DDY(f, outloc, toString(method), region);
  };

  const Field2D DDZ(const Field2D& f, CELL_LOC outloc = CELL_DEFAULT,
                    const std::string& method = "DEFAULT", REGION region = RGN_NOBNDRY);
  const Field2D DDZ(const Field2D& f, CELL_LOC outloc, DIFF_METHOD method,
                    REGION region = RGN_NOBNDRY) {
    return DDZ(f, outloc, toString(method), region);
  };

  /// Gradient along magnetic field  b.Grad(f)
  const Field2D Grad_par(const Field2D& var, CELL_LOC outloc = CELL_DEFAULT,
                         const std::string& method = "DEFAULT");
  const Field2D Grad_par(const Field2D& var, CELL_LOC outloc, DIFF_METHOD method) {
    return Grad_par(var, outloc, toString(method));
  };

  const Field3D Grad_par(const Field3D& var, CELL_LOC outloc = CELL_DEFAULT,
                         const std::string& method = "DEFAULT");
  const Field3D Grad_par(const Field3D& var, CELL_LOC outloc, DIFF_METHOD method) {
    return Grad_par(var, outloc, toString(method));
  };

  /// Advection along magnetic field V*b.Grad(f)
  const Field2D Vpar_Grad_par(const Field2D& v, const Field2D& f,
                              CELL_LOC outloc = CELL_DEFAULT,
                              const std::string& method = "DEFAULT");
  const Field2D Vpar_Grad_par(const Field2D& v, const Field2D& f, CELL_LOC outloc,
                              DIFF_METHOD method) {
    return Vpar_Grad_par(v, f, outloc, toString(method));
  };

  const Field3D Vpar_Grad_par(const Field3D& v, const Field3D& f,
                              CELL_LOC outloc = CELL_DEFAULT,
                              const std::string& method = "DEFAULT");
  const Field3D Vpar_Grad_par(const Field3D& v, const Field3D& f, CELL_LOC outloc,
                              DIFF_METHOD method) {
    return Vpar_Grad_par(v, f, outloc, toString(method));
  };

  /// Divergence along magnetic field  Div(b*f) = B.Grad(f/B)
  const Field2D Div_par(const Field2D& f, CELL_LOC outloc = CELL_DEFAULT,
                        const std::string& method = "DEFAULT");
  const Field2D Div_par(const Field2D& f, CELL_LOC outloc, DIFF_METHOD method) {
    return Div_par(f, outloc, toString(method));
  };

  const Field3D Div_par(const Field3D& f, CELL_LOC outloc = CELL_DEFAULT,
                        const std::string& method = "DEFAULT");
  const Field3D Div_par(const Field3D& f, CELL_LOC outloc, DIFF_METHOD method) {
    return Div_par(f, outloc, toString(method));
  };

  // Second derivative along magnetic field
  const Field2D Grad2_par2(const Field2D& f, CELL_LOC outloc = CELL_DEFAULT,
                           const std::string& method = "DEFAULT");
  const Field2D Grad2_par2(const Field2D& f, CELL_LOC outloc, DIFF_METHOD method) {
    return Grad2_par2(f, outloc, toString(method));
  };

  const Field3D Grad2_par2(const Field3D& f, CELL_LOC outloc = CELL_DEFAULT,
                           const std::string& method = "DEFAULT");
  const Field3D Grad2_par2(const Field3D& f, CELL_LOC outloc, DIFF_METHOD method) {
    return Grad2_par2(f, outloc, toString(method));
  };

  // Perpendicular Laplacian operator, using only X-Z derivatives
  // NOTE: This might be better bundled with the Laplacian inversion code
  // since it makes use of the same coefficients and FFT routines
  const Field2D Delp2(const Field2D& f, CELL_LOC outloc = CELL_DEFAULT,
                      bool useFFT = true);
  const Field3D Delp2(const Field3D& f, CELL_LOC outloc = CELL_DEFAULT,
                      bool useFFT = true);
  const FieldPerp Delp2(const FieldPerp& f, CELL_LOC outloc = CELL_DEFAULT,
                        bool useFFT = true);

  // Full parallel Laplacian operator on scalar field
  // Laplace_par(f) = Div( b (b dot Grad(f)) ) 
  const Field2D Laplace_par(const Field2D &f, CELL_LOC outloc=CELL_DEFAULT);
  const Field3D Laplace_par(const Field3D &f, CELL_LOC outloc=CELL_DEFAULT);
  
  // Full Laplacian operator on scalar field
  const Field2D Laplace(const Field2D &f, CELL_LOC outloc=CELL_DEFAULT);
  const Field3D Laplace(const Field3D &f, CELL_LOC outloc=CELL_DEFAULT);

  // Full perpendicular Laplacian, in form of inverse of Laplacian operator in LaplaceXY solver
  const Field2D Laplace_perpXY(const Field2D &A, const Field2D &f);

private:
  int nz; // Size of mesh in Z. This is mesh->ngz-1
  Mesh * localmesh;
  CELL_LOC location;

  /// Handles calculation of yup and ydown
  std::unique_ptr<ParallelTransform> transform{nullptr};

  /// Set the parallel (y) transform from the options file.
  /// Used in the constructor to create the transform object.
  void setParallelTransform(Options* options);
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
