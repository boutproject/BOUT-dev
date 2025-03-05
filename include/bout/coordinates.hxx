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
 * Copyright 2014-2025 BOUT++ contributors
 *
 * Contact: Ben Dudson, dudson2@llnl.gov
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

#ifndef BOUT_COORDINATES_H
#define BOUT_COORDINATES_H

#include "bout/bout_types.hxx"
#include "bout/field2d.hxx"
#include "bout/field3d.hxx"
#include "bout/paralleltransform.hxx"

class Mesh;

/*!
 * Represents a coordinate system, and associated operators
 *
 * This is a container for a collection of metric tensor components
 */
class Coordinates {
public:
#if BOUT_USE_METRIC_3D
  using FieldMetric = Field3D;
#else
  using FieldMetric = Field2D;
#endif

  /// Standard constructor from input
  Coordinates(Mesh* mesh, Options* options = nullptr);

  /// Constructor interpolating from another Coordinates object
  /// By default attempts to read staggered Coordinates from grid data source,
  /// interpolating from CELL_CENTRE if not present. Set
  /// force_interpolate_from_centre argument to true to always interpolate
  /// (useful if CELL_CENTRE Coordinates have been changed, so reading from file
  /// would not be correct).
  Coordinates(Mesh* mesh, Options* options, const CELL_LOC loc,
              const Coordinates* coords_in, bool force_interpolate_from_centre = false);

  /// A constructor useful for testing purposes. To use it, inherit
  /// from Coordinates. If \p calculate_geometry is true (default),
  /// calculate the non-uniform variables, Christoffel symbols
  Coordinates(Mesh* mesh, FieldMetric dx, FieldMetric dy, FieldMetric dz, FieldMetric J,
              FieldMetric Bxy, FieldMetric g11, FieldMetric g22, FieldMetric g33,
              FieldMetric g12, FieldMetric g13, FieldMetric g23, FieldMetric g_11,
              FieldMetric g_22, FieldMetric g_33, FieldMetric g_12, FieldMetric g_13,
              FieldMetric g_23, FieldMetric ShiftTorsion, FieldMetric IntShiftTorsion);

  Coordinates& operator=(Coordinates&&) = default;

  ~Coordinates() = default;

  /// Add variables to \p output_options, for post-processing
  void outputVars(Options& output_options);

  FieldMetric dx, dy, dz; ///< Mesh spacing in x, y and z

  /// Length of the Z domain. Used for FFTs
  const Field2D& zlength() const;

  /// True if corrections for non-uniform mesh spacing should be included in operators
  bool non_uniform;
  /// 2nd-order correction for non-uniform meshes d/di(1/dx), d/di(1/dy) and d/di(1/dz)
  FieldMetric d1_dx, d1_dy, d1_dz;

  FieldMetric J; ///< Coordinate system Jacobian, so volume of cell is J*dx*dy*dz

  FieldMetric Bxy; ///< Magnitude of B = nabla z times nabla x

  /// Contravariant metric tensor (g^{ij})
  FieldMetric g11, g22, g33, g12, g13, g23;

  /// Covariant metric tensor
  FieldMetric g_11, g_22, g_33, g_12, g_13, g_23;

  /// Christoffel symbol of the second kind (connection coefficients)
  FieldMetric G1_11, G1_22, G1_33, G1_12, G1_13, G1_23;
  FieldMetric G2_11, G2_22, G2_33, G2_12, G2_13, G2_23;
  FieldMetric G3_11, G3_22, G3_33, G3_12, G3_13, G3_23;

  FieldMetric G1, G2, G3;

  /// d pitch angle / dx. Needed for vector differentials (Curl)
  FieldMetric ShiftTorsion;

  FieldMetric IntShiftTorsion; ///< Integrated shear (I in BOUT notation)

  /// Calculate differential geometry quantities from the metric tensor
  int geometry(bool recalculate_staggered = true,
               bool force_interpolate_from_centre = false);
  /// Invert contravatiant metric to get covariant components
  int calcCovariant(const std::string& region = "RGN_ALL");
  /// Invert covariant metric to get contravariant components
  int calcContravariant(const std::string& region = "RGN_ALL");
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

  FieldMetric DDX(const Field2D& f, CELL_LOC outloc = CELL_DEFAULT,
                  const std::string& method = "DEFAULT",
                  const std::string& region = "RGN_NOBNDRY");

  FieldMetric DDY(const Field2D& f, CELL_LOC outloc = CELL_DEFAULT,
                  const std::string& method = "DEFAULT",
                  const std::string& region = "RGN_NOBNDRY") const;

  FieldMetric DDZ(const Field2D& f, CELL_LOC outloc = CELL_DEFAULT,
                  const std::string& method = "DEFAULT",
                  const std::string& region = "RGN_NOBNDRY");

  Field3D DDX(const Field3D& f, CELL_LOC outloc = CELL_DEFAULT,
              const std::string& method = "DEFAULT",
              const std::string& region = "RGN_NOBNDRY");

  Field3D DDY(const Field3D& f, CELL_LOC outloc = CELL_DEFAULT,
              const std::string& method = "DEFAULT",
              const std::string& region = "RGN_NOBNDRY") const;

  Field3D DDZ(const Field3D& f, CELL_LOC outloc = CELL_DEFAULT,
              const std::string& method = "DEFAULT",
              const std::string& region = "RGN_NOBNDRY");

  /// Gradient along magnetic field  b.Grad(f)
  FieldMetric Grad_par(const Field2D& var, CELL_LOC outloc = CELL_DEFAULT,
                       const std::string& method = "DEFAULT");

  Field3D Grad_par(const Field3D& var, CELL_LOC outloc = CELL_DEFAULT,
                   const std::string& method = "DEFAULT");

  /// Advection along magnetic field V*b.Grad(f)
  FieldMetric Vpar_Grad_par(const Field2D& v, const Field2D& f,
                            CELL_LOC outloc = CELL_DEFAULT,
                            const std::string& method = "DEFAULT");

  Field3D Vpar_Grad_par(const Field3D& v, const Field3D& f,
                        CELL_LOC outloc = CELL_DEFAULT,
                        const std::string& method = "DEFAULT");

  /// Divergence along magnetic field  Div(b*f) = B.Grad(f/B)
  FieldMetric Div_par(const Field2D& f, CELL_LOC outloc = CELL_DEFAULT,
                      const std::string& method = "DEFAULT");

  Field3D Div_par(const Field3D& f, CELL_LOC outloc = CELL_DEFAULT,
                  const std::string& method = "DEFAULT");

  // Second derivative along magnetic field
  FieldMetric Grad2_par2(const Field2D& f, CELL_LOC outloc = CELL_DEFAULT,
                         const std::string& method = "DEFAULT");

  Field3D Grad2_par2(const Field3D& f, CELL_LOC outloc = CELL_DEFAULT,
                     const std::string& method = "DEFAULT");
  // Perpendicular Laplacian operator, using only X-Z derivatives
  // NOTE: This might be better bundled with the Laplacian inversion code
  // since it makes use of the same coefficients and FFT routines
  FieldMetric Delp2(const Field2D& f, CELL_LOC outloc = CELL_DEFAULT, bool useFFT = true);
  Field3D Delp2(const Field3D& f, CELL_LOC outloc = CELL_DEFAULT, bool useFFT = true);
  FieldPerp Delp2(const FieldPerp& f, CELL_LOC outloc = CELL_DEFAULT, bool useFFT = true);

  // Full parallel Laplacian operator on scalar field
  // Laplace_par(f) = Div( b (b dot Grad(f)) )
  FieldMetric Laplace_par(const Field2D& f, CELL_LOC outloc = CELL_DEFAULT);
  Field3D Laplace_par(const Field3D& f, CELL_LOC outloc = CELL_DEFAULT);

  // Full Laplacian operator on scalar field
  FieldMetric Laplace(const Field2D& f, CELL_LOC outloc = CELL_DEFAULT,
                      const std::string& dfdy_boundary_conditions = "free_o3",
                      const std::string& dfdy_dy_region = "");
  Field3D Laplace(const Field3D& f, CELL_LOC outloc = CELL_DEFAULT,
                  const std::string& dfdy_boundary_conditions = "free_o3",
                  const std::string& dfdy_dy_region = "");

  // Full perpendicular Laplacian, in form of inverse of Laplacian operator in LaplaceXY
  // solver
  Field2D Laplace_perpXY(const Field2D& A, const Field2D& f);

private:
  int nz; // Size of mesh in Z. This is mesh->ngz-1
  Mesh* localmesh;
  CELL_LOC location;

  /// Handles calculation of yup and ydown
  std::unique_ptr<ParallelTransform> transform{nullptr};

  /// Cache variable for `zlength`. Invalidated when
  /// `Coordinates::geometry` is called
  mutable std::unique_ptr<Field2D> zlength_cache{nullptr};

  /// Cache variable for Grad2_par2
  mutable std::map<std::string, std::unique_ptr<FieldMetric>> Grad2_par2_DDY_invSgCache;
  mutable std::unique_ptr<FieldMetric> invSgCache{nullptr};

  /// Set the parallel (y) transform from the options file.
  /// Used in the constructor to create the transform object.
  void setParallelTransform(Options* options);

  const FieldMetric& invSg() const;
  const FieldMetric& Grad2_par2_DDY_invSg(CELL_LOC outloc,
                                          const std::string& method) const;

  // check that covariant tensors are positive (if expected) and finite (always)
  void checkCovariant();
  // check that contravariant tensors are positive (if expected) and finite (always)
  void checkContravariant();
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

#endif // BOUT_COORDINATES_H
