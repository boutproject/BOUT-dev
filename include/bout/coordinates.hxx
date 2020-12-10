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

class Datafile;
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
  Coordinates(Mesh* mesh, FieldMetric dx, FieldMetric dy, FieldMetric dz, FieldMetric J,
              FieldMetric Bxy, FieldMetric g11, FieldMetric g22, FieldMetric g33,
              FieldMetric g12, FieldMetric g13, FieldMetric g23, FieldMetric g_11,
              FieldMetric g_22, FieldMetric g_33, FieldMetric g_12, FieldMetric g_13,
              FieldMetric g_23, FieldMetric ShiftTorsion, FieldMetric IntShiftTorsion,
              bool calculate_geometry = true);

  Coordinates& operator=(Coordinates&&) = default;

  ~Coordinates() = default;

  /*!
   * Adds variables to the output file, for post-processing
   * 
   * Must be a better way so that Coordinates doesn't depend on Datafile
   */
  void outputVars(Datafile &file);

  FieldMetric dx, dy, dz; ///< Mesh spacing in x and y
  // BoutReal dz; ///< Mesh spacing in Z

  Field2D zlength() const {
#if BOUT_USE_METRIC_3D
    Field2D result(0., localmesh);
    BOUT_FOR_SERIAL(i, dz.getRegion("RGN_ALL")) { result[i] += dz[i]; }
    return result;
#else
    return dz * nz;
#endif
  } ///< Length of the Z domain. Used for FFTs

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

  /// Return if the metrics are 3D
  // needs to be static for old gcc
  static constexpr bool is3D() {
    return bout::build::use_metric_3d;
  }

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

#ifdef DERIV_FUNC_REGION_ENUM_TO_STRING
#error This utility macro should not clash with another one
#else
#define DERIV_FUNC_REGION_ENUM_TO_STRING(func, ResultType, T)                          \
  [[deprecated(                                                                        \
      "Please use Coordinates::#func(const #T& f, "                                    \
      "CELL_LOC outloc = CELL_DEFAULT, const std::string& method = \"DEFAULT\", "      \
      "const std::string& region = \"RGN_ALL\") instead")]] inline ResultType          \
  func(const T& f, CELL_LOC outloc, const std::string& method, REGION region) {        \
    return func(f, outloc, method, toString(region));                                  \
  }                                                                                    \
  [[deprecated(                                                                        \
      "Please use Coordinates::#func(const #T& f, "                                    \
      "CELL_LOC outloc = CELL_DEFAULT, const std::string& method = \"DEFAULT\", "      \
      "const std::string& region = \"RGN_ALL\") instead")]] inline ResultType          \
  func(const T& f, CELL_LOC outloc, DIFF_METHOD method, REGION region = RGN_NOBNDRY) { \
    return func(f, outloc, toString(method), toString(region));                        \
  }
#endif

#ifdef GRAD_FUNC_REGION_ENUM_TO_STRING
#error This utility macro should not clash with another one
#else
#define GRAD_FUNC_REGION_ENUM_TO_STRING(func, ResultType, T)                      \
  [[deprecated(                                                                   \
      "Please use Coordinates::#func(const #T& f, "                               \
      "CELL_LOC outloc = CELL_DEFAULT, const std::string& method = \"DEFAULT\") " \
      "instead")]] inline ResultType                                              \
  func(const T& f, CELL_LOC outloc, DIFF_METHOD method) {                         \
    return func(f, outloc, toString(method));                                     \
  }
#endif

  FieldMetric DDX(const Field2D& f, CELL_LOC outloc = CELL_DEFAULT,
                        const std::string& method = "DEFAULT",
                        const std::string& region = "RGN_NOBNDRY");
  DERIV_FUNC_REGION_ENUM_TO_STRING(DDX, FieldMetric, Field2D);

  FieldMetric DDY(const Field2D& f, CELL_LOC outloc = CELL_DEFAULT,
                        const std::string& method = "DEFAULT",
                        const std::string& region = "RGN_NOBNDRY");
  DERIV_FUNC_REGION_ENUM_TO_STRING(DDY, FieldMetric, Field2D);

  FieldMetric DDZ(const Field2D& f, CELL_LOC outloc = CELL_DEFAULT,
                        const std::string& method = "DEFAULT",
                        const std::string& region = "RGN_NOBNDRY");
  DERIV_FUNC_REGION_ENUM_TO_STRING(DDZ, FieldMetric, Field2D);

  Field3D DDX(const Field3D& f, CELL_LOC outloc = CELL_DEFAULT,
              const std::string& method = "DEFAULT",
              const std::string& region = "RGN_NOBNDRY");

  Field3D DDY(const Field3D& f, CELL_LOC outloc = CELL_DEFAULT,
              const std::string& method = "DEFAULT",
              const std::string& region = "RGN_NOBNDRY");

  Field3D DDZ(const Field3D& f, CELL_LOC outloc = CELL_DEFAULT,
              const std::string& method = "DEFAULT",
              const std::string& region = "RGN_NOBNDRY");

  /// Gradient along magnetic field  b.Grad(f)
  FieldMetric Grad_par(const Field2D& var, CELL_LOC outloc = CELL_DEFAULT,
                             const std::string& method = "DEFAULT");
  GRAD_FUNC_REGION_ENUM_TO_STRING(Grad_par, FieldMetric, Field2D);

  Field3D Grad_par(const Field3D& var, CELL_LOC outloc = CELL_DEFAULT,
      const std::string& method = "DEFAULT");
  GRAD_FUNC_REGION_ENUM_TO_STRING(Grad_par, Field3D, Field3D);

  /// Advection along magnetic field V*b.Grad(f)
  FieldMetric Vpar_Grad_par(const Field2D& v, const Field2D& f,
                                  CELL_LOC outloc = CELL_DEFAULT,
                                  const std::string& method = "DEFAULT");
  [[deprecated(
      "Please use Coordinates::Vpar_Grad_par(const Field2D& v, "
      "const Field2D& f, CELL_LOC outloc = CELL_DEFAULT, "
      "const std::string& method = \"DEFAULT\") instead")]] inline FieldMetric
  Vpar_Grad_par(const Field2D& v, const Field2D& f, CELL_LOC outloc, DIFF_METHOD method) {
    return Vpar_Grad_par(v, f, outloc, toString(method));
  }

  Field3D Vpar_Grad_par(const Field3D& v, const Field3D& f,
      CELL_LOC outloc = CELL_DEFAULT, const std::string& method = "DEFAULT");
  [[deprecated("Please use Coordinates::Vpar_Grad_par(const Field3D& v, "
      "const Field3D& f, CELL_LOC outloc = CELL_DEFAULT, "
      "const std::string& method = \"DEFAULT\") instead")]]
  inline Field3D Vpar_Grad_par(const Field3D& v, const Field3D& f, CELL_LOC outloc,
      DIFF_METHOD method) {
    return Vpar_Grad_par(v, f, outloc, toString(method));
  }

  /// Divergence along magnetic field  Div(b*f) = B.Grad(f/B)
  FieldMetric Div_par(const Field2D& f, CELL_LOC outloc = CELL_DEFAULT,
                            const std::string& method = "DEFAULT");
  GRAD_FUNC_REGION_ENUM_TO_STRING(Div_par, FieldMetric, Field2D);

  Field3D Div_par(const Field3D& f, CELL_LOC outloc = CELL_DEFAULT,
      const std::string& method = "DEFAULT");
  GRAD_FUNC_REGION_ENUM_TO_STRING(Div_par, Field3D, Field3D);

  // Second derivative along magnetic field
  FieldMetric Grad2_par2(const Field2D& f, CELL_LOC outloc = CELL_DEFAULT,
                               const std::string& method = "DEFAULT");
  GRAD_FUNC_REGION_ENUM_TO_STRING(Grad2_par2, FieldMetric, Field2D);

  Field3D Grad2_par2(const Field3D& f, CELL_LOC outloc = CELL_DEFAULT,
      const std::string& method = "DEFAULT");
  GRAD_FUNC_REGION_ENUM_TO_STRING(Grad2_par2, Field3D, Field3D);

#undef DERIV_FUNC_REGION_ENUM_TO_STRING
#undef GRAD_FUNC_REGION_ENUM_TO_STRING

  // Perpendicular Laplacian operator, using only X-Z derivatives
  // NOTE: This might be better bundled with the Laplacian inversion code
  // since it makes use of the same coefficients and FFT routines
  FieldMetric Delp2(const Field2D& f, CELL_LOC outloc = CELL_DEFAULT,
                          bool useFFT = true);
  Field3D Delp2(const Field3D& f, CELL_LOC outloc = CELL_DEFAULT, bool useFFT = true);
  FieldPerp Delp2(const FieldPerp& f, CELL_LOC outloc = CELL_DEFAULT, bool useFFT = true);

  // Full parallel Laplacian operator on scalar field
  // Laplace_par(f) = Div( b (b dot Grad(f)) )
  FieldMetric Laplace_par(const Field2D& f, CELL_LOC outloc = CELL_DEFAULT);
  Field3D Laplace_par(const Field3D &f, CELL_LOC outloc=CELL_DEFAULT);
  
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
  Mesh * localmesh;
  CELL_LOC location;

  /// Handles calculation of yup and ydown
  std::unique_ptr<ParallelTransform> transform{nullptr};

  /// Set the parallel (y) transform from the options file.
  /// Used in the constructor to create the transform object.
  void setParallelTransform(Options* options);

  /// Alternative implementation for `fromFieldAligned` that can be
  /// called before the `Coordinates` constructor is finisehd
  Field3D maybeFromFieldAligned(const Field3D& f, const std::string& region = "RGN_ALL");

  Field2D maybeFromFieldAligned(const Field2D& f,
                                const std::string& UNUSED(region) = "RGN_ALL");

  /// A wrapper for index:DDY derivative that is able to tranform
  /// fields before the constructor is finished.
  FieldMetric indexDDY(const Field2D& f, CELL_LOC outloc = CELL_DEFAULT,
                       const std::string& method = "DEFAULT",
                       const std::string& region = "RGN_NOBNDRY");
  Field3D indexDDY(const Field3D& f, CELL_LOC outloc = CELL_DEFAULT,
                   const std::string& method = "DEFAULT",
                   const std::string& region = "RGN_NOBNDRY");
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
