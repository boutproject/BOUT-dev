/**************************************************************************
 * Implementation of the Mesh class. Based on the BoutMesh, but with
 * hopefuly faster derivatives.
 *
 * Changelog
 * ---------
 *
 * 2016..2018 David Schwörer
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

#pragma once

#include "../bout/boutmesh.hxx"

#include <stencils.hxx>

#include <cmath>
#include <list>
#include <vector>

using std::list;
using std::vector;

class AiolosMesh : public BoutMesh {
public:
  AiolosMesh(GridDataSource *s, Options *options = NULL);
  ~AiolosMesh();

  struct cart_diff_lookup_table {
    Mesh::deriv_func func; // Single-argument differencing function
    deriv_func norm;
    deriv_func on;
    deriv_func off;
  };

#include "aiolos_header.hxx"

  // virtual const Field3D interp_to(const Field3D &var, CELL_LOC loc) const;

  virtual const Field3D interp_to(const Field3D &f, CELL_LOC loc, REGION region) const override {
    ASSERT2(f.getMesh() == this);
    if (loc == f.getLocation() || loc == CELL_DEFAULT) {
      return f;
    } else {
      return interp_to_do(f, loc, region);
    }
  }
  virtual const Field2D interp_to(const Field2D &f, CELL_LOC loc, REGION region) const override {
    return f;
  }

  virtual void derivs_init(Options *option) override;

// to check in debugger we have the right mesh
#if CHECK > 1
  bool isAiolos = true;
#endif
private:
  const Field3D interp_to_do(const Field3D &f, CELL_LOC loc, REGION region) const;

#include "aiolos_derivs.hxx"

#include "aiolos_stencils.hxx"

#include "aiolos_interp_to.hxx"

};
