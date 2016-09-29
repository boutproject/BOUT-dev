/**************************************************************************
 * Implementation of the Mesh class, handling input files compatible with
 * BOUT / BOUT-06.
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

#include "cartesianmesh.hxx"

#include <boutexception.hxx>
#include <utils.hxx>
#include <output.hxx>
// Include the auto generated file
#include "generated_stencils.cxx"



CartesianMesh::CartesianMesh(GridDataSource *s, Options *options): BoutMesh(s,options){
}
CartesianMesh::~CartesianMesh(){
}

const Field3D CartesianMesh::applyYdiff(const Field3D &var, Mesh::deriv_func func, Mesh::inner_boundary_deriv_func func_in, Mesh::outer_boundary_deriv_func func_out, CELL_LOC loc ){
  cart_diff_lookup_table chosen;
  for (int i=0;;++i){
    output.write("maybe: %p == %p?\n",diff_lookup_y[i].func,func);
    if (diff_lookup_y[i].func==func){
      chosen=diff_lookup_y[i];
      break;
    }
    if (diff_lookup_y[i].func==NULL){
      throw BoutException("Diff method not found!");
    }
  }
  if ( loc == CELL_YLOW && var.getLocation() != CELL_YLOW ){
    return chosen.on(var);
  } else if ( var.getLocation() == CELL_YLOW  && ! (loc == CELL_CENTRE )){
    return chosen.off(var);
  } else {
    if (loc != var.getLocation() && loc != CELL_DEFAULT) throw BoutException("Didn't expect this - I wouldn't interpolate");
    return chosen.norm(var);
  }
}
