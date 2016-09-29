#pragma once

#include "../bout/boutmesh.hxx"

#include <stencils.hxx>

#include <list>
#include <vector>
#include <cmath>

using std::list;
using std::vector;

class CartesianMesh : public BoutMesh {
 public:
  CartesianMesh(GridDataSource *s, Options *options = NULL);
  ~CartesianMesh();

  
  typedef Field3D (*deriv_func)(const Field3D ); // f
  typedef Field3D (*upwind_func)(const Field3D , const Field3D); // v, f
  
  struct cart_diff_lookup_table {
    Mesh::deriv_func func;     // Single-argument differencing function
    deriv_func norm;
    deriv_func on;
    deriv_func off;
  };


  //virtual const Field2D applyXdiff(const Field2D &var, deriv_func func, inner_boundary_deriv_func func_in, outer_boundary_deriv_func func_out, CELL_LOC loc = CELL_DEFAULT);
  //virtual const Field3D applyXdiff(const Field3D &var, deriv_func func, inner_boundary_deriv_func func_in, outer_boundary_deriv_func func_out, CELL_LOC loc = CELL_DEFAULT);
  
  //virtual const Field2D applyYdiff(const Field2D &var, deriv_func func, inner_boundary_deriv_func func_in, outer_boundary_deriv_func func_out, CELL_LOC loc = CELL_DEFAULT);
  virtual const Field3D applyYdiff(const Field3D &var, Mesh::deriv_func func, Mesh::inner_boundary_deriv_func func_in, Mesh::outer_boundary_deriv_func func_out, CELL_LOC loc = CELL_DEFAULT);

  //virtual const Field3D applyZdiff(const Field3D &var, Mesh::deriv_func func, CELL_LOC loc = CELL_DEFAULT);


};
