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

#ifndef LAPLACE_XY2_HYPRE_H
#define LAPLACE_XY2_HYPRE_H

#include "bout/build_defines.hxx"
#include "bout/invert/laplacexy.hxx"

#if !BOUT_HAS_HYPRE

namespace {
RegisterUnavailableLaplaceXY
    registerlaplacexyhypre("hypre", "BOUT++ was not configured with HYPRE");
}

#else // BOUT_HAS_HYPRE

class Mesh;

#include <bout/hypre_interface.hxx>

class LaplaceXY2Hypre : public LaplaceXY {
public:
  LaplaceXY2Hypre(Mesh* m = nullptr, Options* opt = nullptr, CELL_LOC loc = CELL_CENTRE);
  LaplaceXY2Hypre(const LaplaceXY2Hypre&) = delete;
  LaplaceXY2Hypre(LaplaceXY2Hypre&&) = delete;
  LaplaceXY2Hypre& operator=(const LaplaceXY2Hypre&) = delete;
  LaplaceXY2Hypre& operator=(LaplaceXY2Hypre&&) = delete;
  ~LaplaceXY2Hypre() override = default;

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
   * This is called by Hypre via a static function.
   * and should not be called by external users
   */
  // int precon(HYPRE_IJVector x, HYPRE_IJVector y);

private:
  Mesh* localmesh; ///< The mesh this operates on, provides metrics and communication
  IndexerPtr<Field2D> indexConverter;
  bout::HypreMatrix<Field2D> M;
  bout::HypreVector<Field2D> x;
  bout::HypreVector<Field2D> b;
  bout::HypreSystem<Field2D> linearSystem;

  // Y derivatives
  bool include_y_derivs; // Include Y derivative terms?

  // Boundary conditions
  bool x_inner_dirichlet; // Dirichlet on inner X boundary?
  bool y_bndry_dirichlet; // Dirichlet on Y boundary?

  bool print_timing;

  // Location of the rhs and solution
  CELL_LOC location;

  /*!
   * Return the communicator for XY
   */
  MPI_Comm communicator();
};

namespace {
const inline RegisterLaplaceXY<LaplaceXY2Hypre> registerlaplacexyhypre{"hypre"};
}

#endif // BOUT_HAS_HYPRE
#endif // LAPLACE_XY_HYPRE_H
