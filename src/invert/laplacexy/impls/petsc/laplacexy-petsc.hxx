/**************************************************************************
 * Laplacian solver in 2D (X-Y)
 *
 * Equation solved is: 
 *
 * Div( A * Grad_perp(x) ) + B*x = b
 *
 * Intended for use in solving n = 0 component of potential
 * from inversion of vorticity equation
 * 
 **************************************************************************
 * Copyright 2015 B.Dudson
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

#ifndef BOUT_LAPLACE_XY_PETSC_H
#define BOUT_LAPLACE_XY_PETSC_H

#include "bout/build_defines.hxx"
#include "bout/invert/laplacexy.hxx"

#if !BOUT_HAS_PETSC

namespace {
RegisterUnavailableLaplaceXY
    registerlaplacexypetsc("petsc", "BOUT++ was not configured with PETSc");
}

#else // BOUT_HAS_PETSC

#include "bout/bout_types.hxx"
#include "bout/cyclic_reduction.hxx"
#include "bout/field2d.hxx"
#include "bout/mesh.hxx"
#include "bout/monitor.hxx"
#include "bout/petsclib.hxx"
#include "bout/utils.hxx"

#include <memory>
#include <string>

class Options;
class Solver;

class LaplaceXYpetsc : public LaplaceXY {
public:
  LaplaceXYpetsc(Mesh* m = nullptr, Options* opt = nullptr, CELL_LOC loc = CELL_CENTRE);
  LaplaceXYpetsc(const LaplaceXYpetsc&) = delete;
  LaplaceXYpetsc(LaplaceXYpetsc&&) = delete;
  LaplaceXYpetsc& operator=(const LaplaceXYpetsc&) = delete;
  LaplaceXYpetsc& operator=(LaplaceXYpetsc&&) = delete;
  ~LaplaceXYpetsc() override;

  /*!
   * Set coefficients (A, B) in equation:
   * Div( A * Grad_perp(x) ) + B*x = b
   */
  void setCoefs(const Field2D& A, const Field2D& B) override;

  /*!
   * Solve Laplacian in X-Y
   * 
   * Inputs
   * ======
   *
   * rhs  - The field to be inverted. This must be allocated
   *        and contain valid data.
   * x0   - Initial guess at the solution. If this is unallocated
   *        then an initial guess of zero will be used.
   * 
   * Returns
   * =======
   * 
   * The solution as a Field2D. On failure an exception will be raised
   * 
   */
  Field2D solve(const Field2D& rhs, const Field2D& x0) override;

  /*!
   * Preconditioner function
   * This is called by PETSc via a static function.
   * and should not be called by external users
   */
  int precon(Vec x, Vec y);

  /*!
   * If this method is called, save some performance monitoring information
   */
  void savePerformance(Solver& solver, const std::string& name = "");

private:
  PetscLib lib;       ///< Requires PETSc library
  Mat MatA = nullptr; ///< Matrix to be inverted
  Vec xs = nullptr;   ///< Solution vector
  Vec bs = nullptr;   ///< RHS vectors
  KSP ksp = nullptr;  ///< Krylov Subspace solver
  PC pc = nullptr;    ///< Preconditioner

  Mesh* localmesh; ///< The mesh this operates on, provides metrics and communication

  /// default prefix for writing performance logging variables
  std::string default_prefix;

  // Preconditioner
  int xstart, xend;
  int nloc, nsys;
  Matrix<BoutReal> acoef, bcoef, ccoef, xvals, bvals;
  std::unique_ptr<CyclicReduce<BoutReal>> cr; ///< Tridiagonal solver

  // Use finite volume or finite difference discretization
  bool finite_volume{true};

  // Y derivatives
  bool include_y_derivs; // Include Y derivative terms?

  // Boundary conditions
  bool x_inner_dirichlet;         // Dirichlet on inner X boundary?
  bool x_outer_dirichlet{false};  // Dirichlet on outer X boundary?
  std::string y_bndry{"neumann"}; // Boundary condition for y-boundary

  // Location of the rhs and solution
  CELL_LOC location;

  /*!
   * Number of grid points on this processor
   */
  int localSize();

  /*!
   * Return the communicator for XY
   */
  MPI_Comm communicator();

  /*!
   * Return the global index of a local (x,y) coordinate
   * including guard cells.
   * Boundary cells have a global index of -1
   *
   * To do this, a Field2D (indexXY) is used to store
   * the index as a floating point number which is then rounded
   * to an integer. Guard cells are filled by communication
   * so no additional logic is needed in Mesh.
   */
  int globalIndex(int x, int y);
  Field2D indexXY; ///< Global index (integer stored as BoutReal)

  // Save performance information?
  bool save_performance = false;

  // Running average of number of iterations taken for solve in each output timestep
  BoutReal average_iterations = 0.;

  // Variable to store the final result of average_iterations, since output is
  // written after all other monitors have been called, and average_iterations
  // must be reset in the monitor
  BoutReal output_average_iterations = 0.;

  // Running total of number of calls to the solver in each output timestep
  int n_calls = 0;

  // Utility methods
  void setPreallocationFiniteVolume(PetscInt* d_nnz, PetscInt* o_nnz);
  void setPreallocationFiniteDifference(PetscInt* d_nnz, PetscInt* o_nnz);
  void setMatrixElementsFiniteVolume(const Field2D& A, const Field2D& B);
  void setMatrixElementsFiniteDifference(const Field2D& A, const Field2D& B);
  void solveFiniteVolume(const Field2D& x0);
  void solveFiniteDifference(const Field2D& x0);

  // Monitor class used to reset performance-monitoring variables for a new
  // output timestep
  friend class LaplaceXYMonitor;
  class LaplaceXYMonitor : public Monitor {
  public:
    LaplaceXYMonitor(LaplaceXYpetsc& owner) : laplacexy(&owner) {}

    int call(Solver* /*solver*/, BoutReal /*time*/, int /*iter*/, int /*nout*/) override;
    void outputVars(Options& output_options, const std::string& time_dimension) override;

  private:
    // LaplaceXY object that this monitor belongs to
    LaplaceXYpetsc* laplacexy;
  };

  LaplaceXYMonitor monitor;
};

namespace {
const inline RegisterLaplaceXY<LaplaceXYpetsc> registerlaplacexypetsc{"petsc"};
}

#endif // BOUT_HAS_PETSC
#endif // BOUT_LAPLACE_XY_PETSC_H
