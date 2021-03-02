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

#include <bout/build_defines.hxx>

#if not BOUT_HAS_HYPRE
// If no Hypre

#warning LaplaceXY requires Hypre. No LaplaceXY available

#include <bout/mesh.hxx>
#include <boutexception.hxx>
#include <options.hxx>
#include "bout/globalindexer.hxx"

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
   Field2D solve(Field2D& rhs, Field2D& x0);

  /*!
   * Preconditioner function
   * This is called by Hypre via a static function.
   * and should not be called by external users
   */
  //int precon(HYPRE_IJVector x, HYPRE_IJVector y);

private:
  Mesh* localmesh; ///< The mesh this operates on, provides metrics and communication
  Field2D f2dinit;                   ///< This is here just to initialise matrix
  IndexerPtr<Field2D> indexConverter;
  bout::HypreMatrix<Field2D> *M;
  bout::HypreVector<Field2D> *x;
  bout::HypreVector<Field2D> *b;
  bout::HypreSystem<Field2D> *linearSystem;

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


template <class T>
OperatorStencil<T> squareStencil(Mesh* localmesh) {
  OperatorStencil<T> stencil;
  IndexOffset<T> zero;
  std::set<IndexOffset<T>> offsets = {
      zero,
      zero.xp(),
      zero.xm(),
  };
  if (!std::is_same<T, IndPerp>::value) {
    offsets.insert(zero.yp());
    offsets.insert(zero.ym());
    offsets.insert(zero.xp().yp());
    offsets.insert(zero.xp().ym());
    offsets.insert(zero.xm().yp());
    offsets.insert(zero.xm().ym());
  }
  if (!std::is_same<T, Ind2D>::value) {
    offsets.insert(zero.zp());
    offsets.insert(zero.zm());
    offsets.insert(zero.xp().zp());
    offsets.insert(zero.xp().zm());
    offsets.insert(zero.xm().zp());
    offsets.insert(zero.xm().zm());
  }
  if (std::is_same<T, Ind3D>::value) {
    offsets.insert(zero.yp().zp());
    offsets.insert(zero.yp().zm());
    offsets.insert(zero.ym().zp());
    offsets.insert(zero.ym().zm());
  }
  std::vector<IndexOffset<T>> offsetsVec(offsets.begin(), offsets.end());
  stencil.add(
      [localmesh](T ind) -> bool {
        return (localmesh->xstart <= ind.x() && ind.x() <= localmesh->xend
                && (std::is_same<T, IndPerp>::value
                    || (localmesh->ystart <= ind.y() && ind.y() <= localmesh->yend))
                && (std::is_same<T, Ind2D>::value
                    || (localmesh->zstart <= ind.z() && ind.z() <= localmesh->zend)));
              },
      offsetsVec);
  stencil.add([](T UNUSED(ind)) -> bool { return true; }, {zero});
  return stencil;
}

#endif // BOUT_HAS_HYPRE
#endif // LAPLACE_XY_HYPRE_H
