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

#ifndef BOUT_LAPLACE_XY_PETSC2_H
#define BOUT_LAPLACE_XY_PETSC2_H

#include "bout/build_defines.hxx"
#include "bout/invert/laplacexy.hxx"

#if !BOUT_HAS_PETSC
namespace {
RegisterUnavailableLaplaceXY
    registerlaplacexypetsc2("petsc2", "BOUT++ was not configured with PETSc");
}

#elif BOUT_USE_METRIC_3D
namespace {
RegisterUnavailableLaplaceXY
    registerlaplacexypetsc2("petsc2", "BOUT++ was configured with 3D metrics");
}

#else // BOUT_HAS_PETSC

#include "bout/bout_types.hxx"
#include "bout/field2d.hxx"
#include "bout/globalindexer.hxx"
#include "bout/mesh.hxx"
#include "bout/petsc_interface.hxx"
#include "bout/petsclib.hxx"

class LaplaceXYpetsc2 : public LaplaceXY {
public:
  LaplaceXYpetsc2(Mesh* m = nullptr, Options* opt = nullptr, CELL_LOC loc = CELL_CENTRE);
  LaplaceXYpetsc2(const LaplaceXYpetsc2&) = default;
  LaplaceXYpetsc2(LaplaceXYpetsc2&&) = delete;
  LaplaceXYpetsc2& operator=(const LaplaceXYpetsc2&) = default;
  LaplaceXYpetsc2& operator=(LaplaceXYpetsc2&&) = delete;
  ~LaplaceXYpetsc2() override;

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

private:
  PetscLib lib; ///< Requires PETSc library

  Mesh* localmesh; ///< The mesh this operates on, provides metrics and communication

  IndexerPtr<Field2D> indexConverter;
  PetscMatrix<Field2D> matrix; ///< Matrix to be inverted
  KSP ksp = nullptr;           ///< Krylov Subspace solver
  PC pc = nullptr;             ///< Preconditioner

  // Y derivatives
  bool include_y_derivs; // Include Y derivative terms?

  // Boundary conditions
  bool x_inner_dirichlet; // Dirichlet on inner X boundary?
  bool y_bndry_dirichlet; // Dirichlet on Y boundary?

  // Location of the rhs and solution
  CELL_LOC location;

  /*!
   * Return the communicator for XY
   */
  MPI_Comm communicator();
};

namespace {
const inline RegisterLaplaceXY<LaplaceXYpetsc2> registerlaplacexypetsc2{"petsc2"};
}

#endif // BOUT_HAS_PETSC
#endif // BOUT_LAPLACE_XY_PETSC2_H
