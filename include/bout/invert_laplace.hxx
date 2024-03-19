/*!*************************************************************************
 * \file invert_laplace.hxx
 *
 * Perpendicular Laplacian inversion using FFT and Tridiagonal solver
 *
 * Equation solved is: \f$d*\nabla^2_\perp x + (1/c)\nabla_perp c\cdot\nabla_\perp x + a x = b\f$
 * 
 * Where a, c and d are functions of x and y only (not z)
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

class Laplacian;

#ifndef BOUT_LAPLACE_H
#define BOUT_LAPLACE_H

#include "bout/build_config.hxx"

#if BOUT_HAS_PETSC
#define PVEC_REAL_MPI_TYPE MPI_DOUBLE
#endif

#include "bout/field2d.hxx"
#include "bout/field3d.hxx"
#include "bout/fieldperp.hxx"
#include "bout/generic_factory.hxx"
#include "bout/monitor.hxx"
#include "bout/unused.hxx"
#include <bout/boutexception.hxx>

#include "bout/dcomplex.hxx"

class Solver;

constexpr auto LAPLACE_SPT = "spt";
constexpr auto LAPLACE_TRI = "tri";
constexpr auto LAPLACE_BAND = "band";
constexpr auto LAPLACE_PETSC = "petsc";
constexpr auto LAPLACE_PETSCAMG = "petscamg";
constexpr auto LAPLACE_PETSC3DAMG = "petsc3damg";
constexpr auto LAPLACE_HYPRE3D = "hypre3d";
constexpr auto LAPLACE_CYCLIC = "cyclic";
constexpr auto LAPLACE_MULTIGRID = "multigrid";
constexpr auto LAPLACE_NAULIN = "naulin";
constexpr auto LAPLACE_IPT = "ipt";
constexpr auto LAPLACE_PCR = "pcr";
constexpr auto LAPLACE_PCR_THOMAS = "pcr_thomas";

// Inversion flags for each boundary
/// Zero-gradient for DC (constant in Z) component. Default is zero value
constexpr int INVERT_DC_GRAD = 1;
/// Zero-gradient for AC (non-constant in Z) component. Default is zero value
constexpr int INVERT_AC_GRAD = 2;
/// Use zero-laplacian (decaying solution) to AC component
constexpr int INVERT_AC_LAP = 4;
/// Use symmetry to enforce either zero-value or zero-gradient
constexpr int INVERT_SYM = 8;
/// Set boundary to value
constexpr int INVERT_SET = 16;
/// Use input value in RHS boundary
constexpr int INVERT_RHS = 32;
/// Use zero-laplacian solution for DC component
constexpr int INVERT_DC_LAP = 64;
/// Only use one boundary point
constexpr int INVERT_BNDRY_ONE = 128;
constexpr int INVERT_DC_GRADPAR = 256;
constexpr int INVERT_DC_GRADPARINV = 512;
/// For use in cylindrical coordiate system.
constexpr int INVERT_IN_CYLINDER = 1024;

// Global flags
/// Zero the DC (constant in Z) component of the solution
constexpr int INVERT_ZERO_DC = 1;
/// Iterative method start from solution=0. Has no effect for direct solvers
constexpr int INVERT_START_NEW = 2;
/// Sets the width of the boundaries to 1
constexpr int INVERT_BOTH_BNDRY_ONE = 4;
/// Use band solver for 4th order in x
constexpr int INVERT_4TH_ORDER = 8;
/// Zero the kx=0, n = 0 component
constexpr int INVERT_KX_ZERO = 16;

/*
// Legacy flags, can be used in calls to setFlags()
// or in input option "flags"

  const int INVERT_DC_IN_GRAD  = 1;
  const int INVERT_AC_IN_GRAD  = 2;
  const int INVERT_DC_OUT_GRAD = 4;
  const int INVERT_AC_OUT_GRAD = 8;
  const int INVERT_ZERO_DC     = 16;
  const int INVERT_START_NEW   = 32;
  const int INVERT_BNDRY_ONE   = 64; // Sets the width of the boundary to 1
  const int INVERT_4TH_ORDER   = 128; // Use band solver for 4th order in x

  const int INVERT_AC_IN_LAP   = 256;
  const int INVERT_AC_OUT_LAP  = 512;

  const int INVERT_IN_SYM  =  1024; // Use symmetry to enforce either zero-value or zero-gradient
  const int INVERT_OUT_SYM =  2048; // Same for outer boundary
  const int INVERT_IN_SET  =  4096; // Set inner boundary
  const int INVERT_OUT_SET =  8192; // Set outer boundary
  const int INVERT_IN_RHS  = 16384; // Use input value in RHS at inner boundary
  const int INVERT_OUT_RHS = 32768; // Use input value in RHS at outer boundary
  const int INVERT_KX_ZERO = 65536; // Zero the kx=0, n = 0 component

  const int INVERT_DC_IN_LAP = 131072;

  const int INVERT_BNDRY_IN_ONE = 262144;
  const int INVERT_BNDRY_OUT_ONE = 524288;
  const int INVERT_DC_IN_GRADPAR = 1048576;
  const int INVERT_DC_IN_GRADPARINV = 2097152;
 */

class LaplaceFactory
    : public Factory<Laplacian, LaplaceFactory, Options*, CELL_LOC, Mesh*, Solver*> {
public:
  static constexpr auto type_name = "Laplacian";
  static constexpr auto section_name = "laplace";
  static constexpr auto option_name = "type";
  static constexpr auto default_type = LAPLACE_CYCLIC;

  ReturnType create(Options* options = nullptr, CELL_LOC loc = CELL_CENTRE,
                    Mesh* mesh = nullptr, Solver* solver = nullptr) {
    options = optionsOrDefaultSection(options);
    return Factory::create(getType(options), options, loc, mesh, solver);
  }
  ReturnType create(const std::string& type, Options* options) const {
    return Factory::create(type, options, CELL_CENTRE, nullptr, nullptr);
  }
};

/// Simpler name for Factory registration helper class
///
/// Usage:
///
///     #include <bout/laplacefactory.hxx>
///     namespace {
///     RegisterLaplace<MyLaplace> registerlaplacemine("mylaplace");
///     }
template <class DerivedType>
using RegisterLaplace = LaplaceFactory::RegisterInFactory<DerivedType>;

using RegisterUnavailableLaplace = LaplaceFactory::RegisterUnavailableInFactory;

class Options;
class Solver;

/// Base class for Laplacian inversion
class Laplacian {
public:
  Laplacian(Options* options = nullptr, const CELL_LOC loc = CELL_CENTRE,
            Mesh* mesh_in = nullptr, Solver* solver = nullptr);
  virtual ~Laplacian() = default;

  /// Set coefficients for inversion. Re-builds matrices if necessary
  virtual void setCoefA(const Field2D& val) = 0;
  virtual void setCoefA(const Field3D& val) { setCoefA(DC(val)); }
  virtual void setCoefA(BoutReal r) {
    Field2D f(r, localmesh);
    f.setLocation(location);
    setCoefA(f);
  }

  virtual void setCoefC(const Field2D& val) = 0;
  virtual void setCoefC(const Field3D& val) { setCoefC(DC(val)); }
  virtual void setCoefC(BoutReal r) {
    Field2D f(r, localmesh);
    f.setLocation(location);
    setCoefC(f);
  }

  virtual void setCoefC1(const Field2D& UNUSED(val)) {
    throw BoutException("setCoefC1 is not implemented for this Laplacian solver");
  }
  virtual void setCoefC1(const Field3D& val) { setCoefC1(DC(val)); }
  virtual void setCoefC1(BoutReal r) {
    Field2D f(r, localmesh);
    f.setLocation(location);
    setCoefC1(f);
  }

  virtual void setCoefC2(const Field2D& UNUSED(val)) {
    throw BoutException("setCoefC2 is not implemented for this Laplacian solver");
  }
  virtual void setCoefC2(const Field3D& val) { setCoefC2(DC(val)); }
  virtual void setCoefC2(BoutReal r) {
    Field2D f(r, localmesh);
    f.setLocation(location);
    setCoefC2(f);
  }

  virtual void setCoefD(const Field2D& val) = 0;
  virtual void setCoefD(const Field3D& val) { setCoefD(DC(val)); }
  virtual void setCoefD(BoutReal r) {
    Field2D f(r, localmesh);
    f.setLocation(location);
    setCoefD(f);
  }

  virtual void setCoefEx(const Field2D& val) = 0;
  virtual void setCoefEx(const Field3D& val) { setCoefEx(DC(val)); }
  virtual void setCoefEx(BoutReal r) {
    Field2D f(r, localmesh);
    f.setLocation(location);
    setCoefEx(f);
  }

  virtual void setCoefEz(const Field2D& val) = 0;
  virtual void setCoefEz(const Field3D& val) { setCoefEz(DC(val)); }
  virtual void setCoefEz(BoutReal r) {
    Field2D f(r, localmesh);
    f.setLocation(location);
    setCoefEz(f);
  }

  virtual void setGlobalFlags(int f) { global_flags = f; }
  virtual void setInnerBoundaryFlags(int f) { inner_boundary_flags = f; }
  virtual void setOuterBoundaryFlags(int f) { outer_boundary_flags = f; }

  /// Does this solver use Field3D coefficients (true) or only their DC component (false)
  virtual bool uses3DCoefs() const { return false; }

  virtual FieldPerp solve(const FieldPerp& b) = 0;
  virtual Field3D solve(const Field3D& b);
  virtual Field2D solve(const Field2D& b);

  virtual FieldPerp solve(const FieldPerp& b, const FieldPerp& UNUSED(x0)) {
    return solve(b);
  }
  virtual Field3D solve(const Field3D& b, const Field3D& x0);
  virtual Field2D solve(const Field2D& b, const Field2D& x0);

  /// Coefficients in tridiagonal inversion
  void tridagCoefs(int jx, int jy, int jz, dcomplex& a, dcomplex& b, dcomplex& c,
                   const Field2D* ccoef = nullptr, const Field2D* d = nullptr,
                   CELL_LOC loc = CELL_DEFAULT);

  /*!
   * Create a new Laplacian solver
   *
   * @param[in] opt  The options section to use. By default "laplace" will be used
   */
  static std::unique_ptr<Laplacian> create(Options* opts = nullptr,
                                           const CELL_LOC location = CELL_CENTRE,
                                           Mesh* mesh_in = nullptr,
                                           Solver* solver = nullptr) {
    return LaplaceFactory::getInstance().create(opts, location, mesh_in, solver);
  }
  static Laplacian* defaultInstance(); ///< Return pointer to global singleton

  static void cleanup(); ///< Frees all memory

  /// Add any output variables to \p output_options, for example,
  /// performance information, with optional name for the time
  /// dimension
  void outputVars(Options& output_options) const { outputVars(output_options, "t"); }
  virtual void outputVars([[maybe_unused]] Options& output_options,
                          [[maybe_unused]] const std::string& time_dimension) const {}

  /// Register performance monitor with \p solver, prefix output with
  /// `Options` section name
  void savePerformance(Solver& solver) { savePerformance(solver, getPerformanceName()); }
  /// Register performance monitor that is call every timestep with \p
  /// solver, prefix output with \p name. Call this function from your
  /// `PhysicsModel::init` to get time-dependent performance
  /// information, or call `outputVars` directly in non-model code to
  /// get the information at that point in time.
  ///
  /// To use this for a Laplacian implementation, override
  /// `outputVars(Options&, const std::string&)`, and add whatever
  /// information you would like to save to the `Options`
  /// argument. This will then be called automatically by
  /// `LaplacianMonitor`.
  ///
  /// See `LaplaceNaulin::outputVars` for an example.
  void savePerformance(Solver& solver, const std::string& name);

protected:
  bool async_send; ///< If true, use asyncronous send in parallel algorithms

  int maxmode; ///< The maximum Z mode to solve for

  bool low_mem;            ///< If true, reduce the amount of memory used
  bool all_terms;          ///< applies to Delp2 operator and laplacian inversion
  bool nonuniform;         ///< Non-uniform mesh correction
  bool include_yguards;    ///< solve in y-guard cells, default true.
  int extra_yguards_lower; ///< exclude some number of points at the lower boundary, useful for staggered grids or when boundary conditions make inversion redundant
  int extra_yguards_upper; ///< exclude some number of points at the upper boundary, useful for staggered grids or when boundary conditions make inversion redundant

  int global_flags;         ///< Default flags
  int inner_boundary_flags; ///< Flags to set inner boundary condition
  int outer_boundary_flags; ///< Flags to set outer boundary condition

  void tridagCoefs(int jx, int jy, BoutReal kwave, dcomplex& a, dcomplex& b, dcomplex& c,
                   const Field2D* ccoef = nullptr, const Field2D* d = nullptr,
                   CELL_LOC loc = CELL_DEFAULT) {
    tridagCoefs(jx, jy, kwave, a, b, c, ccoef, ccoef, d, loc);
  }
  void tridagCoefs(int jx, int jy, BoutReal kwave, dcomplex& a, dcomplex& b, dcomplex& c,
                   const Field2D* c1coef, const Field2D* c2coef, const Field2D* d,
                   CELL_LOC loc = CELL_DEFAULT);

  void tridagMatrix(dcomplex* avec, dcomplex* bvec, dcomplex* cvec, dcomplex* bk, int jy,
                    int kz, BoutReal kwave, int flags, int inner_boundary_flags,
                    int outer_boundary_flags, const Field2D* a, const Field2D* ccoef,
                    const Field2D* d, bool includeguards = true, bool zperiodic = true) {
    tridagMatrix(avec, bvec, cvec, bk, jy, kz, kwave, flags, inner_boundary_flags,
                 outer_boundary_flags, a, ccoef, ccoef, d, includeguards, zperiodic);
  }
  void tridagMatrix(dcomplex* avec, dcomplex* bvec, dcomplex* cvec, dcomplex* bk, int jy,
                    int kz, BoutReal kwave, int flags, int inner_boundary_flags,
                    int outer_boundary_flags, const Field2D* a, const Field2D* c1coef,
                    const Field2D* c2coef, const Field2D* d, bool includeguards = true,
                    bool zperiodic = true);
  CELL_LOC location;   ///< staggered grid location of this solver
  Mesh* localmesh;     ///< Mesh object for this solver
  Coordinates* coords; ///< Coordinates object, so we only have to call
                       ///  localmesh->getCoordinates(location) once

private:
  /// Singleton instance
  static std::unique_ptr<Laplacian> instance;
  /// Name for writing performance infomation; default taken from
  /// constructing `Options` section
  std::string performance_name;

  class LaplacianMonitor : public Monitor {
  public:
    LaplacianMonitor(Laplacian& owner) : laplacian(&owner) {}
    int call(Solver* solver, BoutReal time, int iter, int nout) override;
    void outputVars(Options& options, const std::string& time_dimension) override;

  private:
    Laplacian* laplacian{nullptr};
  };

  LaplacianMonitor monitor{*this};

public:
  /// Get name for writing performance information
  std::string getPerformanceName() const { return performance_name; };

protected:
  /// Set the name for writing performance information
  void setPerformanceName(std::string name) { performance_name = std::move(name); }
};

////////////////////////////////////////////
// Legacy interface
// These will be removed at some point

void laplace_tridag_coefs(int jx, int jy, int jz, dcomplex& a, dcomplex& b, dcomplex& c,
                          const Field2D* ccoef = nullptr, const Field2D* d = nullptr,
                          CELL_LOC loc = CELL_DEFAULT);

#endif // BOUT_LAPLACE_H
