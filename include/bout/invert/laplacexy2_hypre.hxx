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

#if !BOUT_HAS_HYPRE
// If no Hypre

#warning LaplaceXY requires Hypre. No LaplaceXY available

#include <bout/mesh.hxx>
#include <boutexception.hxx>
#include <options.hxx>

/*!
 * Create a dummy class so that code will compile
 * without Hypre, but will throw an exception if
 * LaplaceXY is used.
 */
class LaplaceXY2Hypre {
public:
  LaplaceXY2Hypre(Mesh* UNUSED(m) = nullptr, Options* UNUSED(opt) = nullptr,
                  const CELL_LOC UNUSED(loc) = CELL_CENTRE) {
    throw BoutException("LaplaceXY requires Hypre. No LaplaceXY available");
  }
  void setCoefs(const Field2D& UNUSED(A), const Field2D& UNUSED(B)) {}
  const Field2D solve(const Field2D& UNUSED(rhs), const Field2D& UNUSED(x0)) {
    throw BoutException("LaplaceXY requires Hypre. No LaplaceXY available");
  }
};

#else // BOUT_HAS_HYPRE

#include "utils.hxx"
#include <bout/hypre_interface.hxx>
#include <bout/mesh.hxx>
#include <cyclic_reduction.hxx>

#include "HYPRE.h"
#include "HYPRE_parcsr_ls.h"

class LaplaceXY2Hypre {
public:
  /*!
   * Constructor
   */
  LaplaceXY2Hypre(Mesh* m = nullptr, Options* opt = nullptr,
                  const CELL_LOC loc = CELL_CENTRE);
  /*!
   * Destructor
   */
  ~LaplaceXY2Hypre();

  /*!
   * Set coefficients (A, B) in equation:
   * Div( A * Grad_perp(x) ) + B*x = b
   */
  void setCoefs(const Field2D& A, const Field2D& B);

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
  const Field2D solve(const Field2D& rhs, const Field2D& x0);

  /*!
   * Preconditioner function
   * This is called by Hypre via a static function.
   * and should not be called by external users
   */
  int precon(HYPRE_IJVector x, HYPRE_IJVector y);

  /// Set up a stencil describing the structure of the operator.
  static OperatorStencil<Ind2D> getStencil(Mesh* localmesh, RangeIterator lowerYBound,
                                           RangeIterator upperYBound);

private:
  Mesh* localmesh; ///< The mesh this operates on, provides metrics and communication
  RangeIterator lowerY, upperY; ///< Lower and upper boundary ranges
  IndexerPtr<Field2D> indexer;  ///< Converts local to global indices
  
  bout::HypreMatrix<Field2D> matrix; ///< Matrix to be inverted
  HYPRE_Solver solver;               ///< Solver

  // Y derivatives
  bool include_y_derivs; // Include Y derivative terms?

  // Boundary conditions
  bool x_inner_dirichlet; // Dirichlet on inner X boundary?
  bool x_outer_dirichlet; // Dirichlet on outer X boundary?
  bool y_bndry_dirichlet; // Dirichlet on Y boundary?

  // Location of the rhs and solution
  CELL_LOC location;

  /*!
   * Return the communicator for XY
   */
  MPI_Comm communicator();
};

#endif // BOUT_HAS_HYPRE
#endif // LAPLACE_XY_HYPRE_H
