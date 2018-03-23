/**************************************************************************
 * Implementation of the Mesh class, handling input files compatible with
 * BOUT++.
 *
 * Changelog
 * ---------
 *
 * 2016-09 David Schwörer
 *           based on the BoutMesh
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

#include "aiolosmesh.hxx"

#include <boutexception.hxx>
#include <interpolation.hxx>
#include <output.hxx>
#include <utils.hxx>

#include <strings.h>

// Include the auto generated files
// Note: They MUST be included in this order. clang-format does break this
#include "generated_init.cxx"
#include "generated_stencils.cxx"
#include "generated_derivs.cxx"

AiolosMesh::AiolosMesh(GridDataSource *s, Options *options) : BoutMesh(s, options) {
  output_info.write("  Using Aiolos Mesh!\n");
  derivs_init(options);
}
AiolosMesh::~AiolosMesh() {}

BoutReal AiolosMesh::GlobalY(int jy) const {
  int gjy = BoutMesh::YGLOBAL(jy);
  output_debug.write(" %d->%d ", jy, gjy);

  return ((BoutReal)gjy) / ((BoutReal)(this->GlobalNy - 2 * this->ystart));
}

// Field3D AiolosMesh::
