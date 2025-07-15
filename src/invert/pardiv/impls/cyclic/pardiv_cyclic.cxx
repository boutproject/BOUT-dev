/************************************************************************
 * Inversion of parallel derivatives
 *
 * Inverts a matrix of the form
 *
 * A + Div_par( B Grad_par )
 *
 * Parallel algorithm, using Cyclic solver
 *
 * Known issues:
 * ------------
 *
 *
 *
 *
 **************************************************************************
 * Copyright 2010 - 2022 BOUT++ contributors
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
 ************************************************************************/

#include "pardiv_cyclic.hxx"
#include "bout/build_defines.hxx"

#if not BOUT_USE_METRIC_3D

#include <bout/boutexception.hxx>
#include <bout/constants.hxx>
#include <bout/cyclic_reduction.hxx>
#include <bout/derivs.hxx>
#include <bout/fft.hxx>
#include <bout/globals.hxx>
#include <bout/msg_stack.hxx>
#include <bout/surfaceiter.hxx>
#include <bout/utils.hxx>

#include <cmath>

InvertParDivCR::InvertParDivCR(Options* opt, CELL_LOC location, Mesh* mesh_in)
    : InvertParDiv(opt, location, mesh_in) {
  // Number of k equations to solve for each x location
  nsys = 1 + (localmesh->LocalNz) / 2;
}

Field3D InvertParDivCR::solve(const Field3D& f) {
  TRACE("InvertParDivCR::solve(Field3D)");
  ASSERT1(localmesh == f.getMesh());
  ASSERT1(location == f.getLocation());

  Field3D result = emptyFrom(f).setDirectionY(YDirectionType::Aligned);

  Coordinates* coord = f.getCoordinates();

  Field3D alignedField = toFieldAligned(f, "RGN_NOBNDRY");

  // Create cyclic reduction object
  auto cr = bout::utils::make_unique<CyclicReduce<dcomplex>>();

  // Find out if we are on a boundary
  int size = localmesh->LocalNy - 2 * localmesh->ystart;
  SurfaceIter surf(localmesh);
  for (surf.first(); !surf.isDone(); surf.next()) {
    int n = localmesh->LocalNy - 2 * localmesh->ystart;
    if (!surf.closed()) {
      // Open field line
      if (surf.firstY()) {
        if (location == CELL_YLOW) {
          // The 'boundary' includes the grid point at mesh->ystart
          n += localmesh->ystart - 1;
        } else {
          n += localmesh->ystart;
        }
      }
      if (surf.lastY()) {
        n += localmesh->ystart;
      }

      if (n > size) {
        size = n; // Maximum size
      }
    }
  }

  auto rhs = Matrix<dcomplex>(localmesh->LocalNy, nsys);
  auto rhsk = Matrix<dcomplex>(nsys, size);
  auto xk = Matrix<dcomplex>(nsys, size);
  auto a = Matrix<dcomplex>(nsys, size);
  auto b = Matrix<dcomplex>(nsys, size);
  auto c = Matrix<dcomplex>(nsys, size);

  const Field2D dy = coord->dy;
  const Field2D J = coord->J;
  const Field2D g_22 = coord->g_22;

  const auto zlength = getUniform(coord->zlength());
  // Loop over flux-surfaces
  for (surf.first(); !surf.isDone(); surf.next()) {
    int x = surf.xpos;

    // Test if open or closed field-lines
    BoutReal ts;
    bool closed = surf.closed(ts);

    // Number of rows
    int y0 = 0;
    int local_ystart = localmesh->ystart;
    size = localmesh->LocalNy - 2 * localmesh->ystart; // If no boundaries
    if (!closed) {
      if (surf.firstY()) {
        if (location == CELL_YLOW) {
          // The 'boundary' includes the grid point at mesh->ystart
          y0 += localmesh->ystart;
          size += localmesh->ystart - 1;
          local_ystart = localmesh->ystart + 1;
        } else {
          y0 += localmesh->ystart;
          size += localmesh->ystart;
        }
      }
      if (surf.lastY()) {
        size += localmesh->ystart;
      }
    }

    // Setup CyclicReduce object
    cr->setup(surf.communicator(), size);
    cr->setPeriodic(closed);

    // Take Fourier transform
    for (int y = 0; y < localmesh->LocalNy - localmesh->ystart - local_ystart; y++) {
      rfft(alignedField(x, y + local_ystart), localmesh->LocalNz, &rhs(y + y0, 0));
    }

    // Set up tridiagonal system
    for (int k = 0; k < nsys; k++) {
      for (int y = 0; y < localmesh->LocalNy - localmesh->ystart - local_ystart; y++) {

        // Interpolate Field2D onto left boundary
        auto Left = [&](const Field2D& field) {
          return 0.5 * (field(x, y + local_ystart) + field(x, y - 1 + local_ystart));
        };
        auto Right = [&](const Field2D& field) {
          return 0.5 * (field(x, y + local_ystart) + field(x, y + 1 + local_ystart));
        };

        BoutReal acoef = A(x, y + local_ystart); // Constant

        // Divergence form, at left and right boundaries
        BoutReal bcoef_L =
            Left(B) * Left(J)
            / (Left(dy) * Left(g_22) * J(x, y + local_ystart) * dy(x, y + local_ystart));
        BoutReal bcoef_R = Right(B) * Right(J)
                           / (Right(dy) * Right(g_22) * J(x, y + local_ystart)
                              * dy(x, y + local_ystart));

        //           const       div(grad)
        //           -----       ---------
        a(k, y + y0) = bcoef_L;
        b(k, y + y0) = acoef - (bcoef_L + bcoef_R);
        c(k, y + y0) = bcoef_R;

        rhsk(k, y + y0) = rhs(y + y0, k); // Transpose
      }
    }

    if (closed) {
      // Twist-shift
      int rank;
      int np;
      bout::globals::mpi->MPI_Comm_rank(surf.communicator(), &rank);
      bout::globals::mpi->MPI_Comm_size(surf.communicator(), &np);

      if (rank == 0) {
        for (int k = 0; k < nsys; k++) {
          BoutReal kwave = k * 2.0 * PI / zlength; // wave number is 1/[rad]
          dcomplex phase(cos(kwave * ts), -sin(kwave * ts));
          a(k, 0) *= phase;
        }
      }
      if (rank == np - 1) {
        for (int k = 0; k < nsys; k++) {
          BoutReal kwave = k * 2.0 * PI / zlength; // wave number is 1/[rad]
          dcomplex phase(cos(kwave * ts), sin(kwave * ts));
          c(k, localmesh->LocalNy - 2 * localmesh->ystart - 1) *= phase;
        }
      }
    } else {
      // Open surface, so may have boundaries
      if (surf.firstY()) {
        for (int k = 0; k < nsys; k++) {
          for (int y = 0; y < localmesh->ystart; y++) {
            a(k, y) = 0.;
            b(k, y) = 1.;
            c(k, y) = -1.;

            rhsk(k, y) = 0.;
          }
        }
      }
      if (surf.lastY()) {
        for (int k = 0; k < nsys; k++) {
          for (int y = size - localmesh->ystart; y < size; y++) {
            a(k, y) = -1.;
            b(k, y) = 1.;
            c(k, y) = 0.;

            rhsk(k, y) = 0.;
          }
        }
      }
    }

    // Solve cyclic tridiagonal system for each k
    cr->setCoefs(a, b, c);
    cr->solve(rhsk, xk);

    // Put back into rhs array
    for (int k = 0; k < nsys; k++) {
      for (int y = 0; y < size; y++) {
        rhs(y, k) = xk(k, y);
      }
    }

    // Inverse Fourier transform
    for (int y = 0; y < size; y++) {
      irfft(&rhs(y, 0), localmesh->LocalNz, result(x, y + local_ystart - y0));
    }
  }

  return fromFieldAligned(result, "RGN_NOBNDRY");
}

#endif // BOUT_USE_METRIC_3D
