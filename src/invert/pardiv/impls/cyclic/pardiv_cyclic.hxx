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
 * - For CELL_YLOW implementation, boundary conditions are only 1st order accurate.
 *   Should be OK for preconditioners, which are allowed to be less accurate.
 *   Only 1st-order accurate one-sided derivative is possible in a tri-diagonal matrix and
 *   staggered mesh requires one-sided derivative as boundary condition.
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

#ifndef __INV_PARDIV_CR_H__
#define __INV_PARDIV_CR_H__

#include "bout/build_config.hxx"
#include "bout/invert_pardiv.hxx"

#if BOUT_USE_METRIC_3D

namespace {
RegisterUnavailableInvertParDiv registerinvertpardivcyclic{
    PARDIVCYCLIC, "BOUT++ was configured with 3D metrics"};
}

#else

#include "dcomplex.hxx"
#include "utils.hxx"
#include <globals.hxx>

class InvertParDivCR : public InvertParDiv {
public:
  explicit InvertParDivCR(Options* opt, CELL_LOC location = CELL_CENTRE,
                          Mesh* mesh_in = bout::globals::mesh);

  using InvertParDiv::solve;
  const Field3D solve(const Field3D& f) override;

  using InvertParDiv::setCoefA;
  void setCoefA(const Field2D& f) override {
    ASSERT1(localmesh == f.getMesh());
    ASSERT1(location == f.getLocation());
    A = f;
  }
  using InvertParDiv::setCoefB;
  void setCoefB(const Field2D& f) override {
    ASSERT1(localmesh == f.getMesh());
    ASSERT1(location == f.getLocation());
    B = f;
  }

private:
  Field2D A{1.0}, B{0.0};

  int nsys;
};

namespace {
RegisterInvertParDiv<InvertParDivCR> registerinvertpardivcyclic{PARDIVCYCLIC};
}

#endif // BOUT_USE_METRIC_3D

#endif // __INV_PARDIV_CR_H__
