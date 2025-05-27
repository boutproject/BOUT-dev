#pragma once
#ifndef BOUT_FIELDOPS_HXX
#define BOUT_FIELDOPS_HXX

#include "bout/bout_types.hxx"
#include "bout/field_accessor.hxx"

struct Add {
  __device__ inline BoutReal operator()(BoutReal a, BoutReal b) const { return a + b; }
};
struct Sub {
  __device__ inline BoutReal operator()(BoutReal a, BoutReal b) const { return a - b; }
};
struct Mul {
  __device__ inline BoutReal operator()(BoutReal a, BoutReal b) const { return a * b; }
};
struct Div {
  __device__ inline BoutReal operator()(BoutReal a, BoutReal b) const { return a / b; }
};

struct BinaryExpr {
  struct RegionIndices {
    int* data;
    int size;

    RegionIndices(int n) : size(n) {
      cudaMallocManaged(&data, n * sizeof(int));
      for (int i = 0; i < n; ++i)
        data[i] = 0;
    }
    ~RegionIndices() { cudaFree(data); }

    __device__ inline int operator()(int idx) const { return data[idx]; }
  };

  using FieldType = FieldAccessor<CELL_CENTRE, Field3D>;

  FieldType lhs;
  FieldType rhs;
  RegionIndices indices;
  Add op;

  Mesh* mesh;
  CELL_LOC location = CELL_CENTRE;
  DirectionTypes directions;

  template <typename IndType>
  BinaryExpr(FieldType lhs, FieldType rhs, Mesh* mesh, CELL_LOC location,
             DirectionTypes directions, const Region<IndType>& region)
      : lhs(lhs), rhs(rhs), mesh(mesh), location(location), directions(directions),
        indices(region.getIndices().size()) {
    // Copy the region indices into the managed array
    for (int i = 0; i < indices.size; ++i) {
      indices.data[i] = region.getIndices()[i].ind;
    }
  }

  __host__ __device__ inline int getSize() const { return indices.size; }
  __device__ inline int regionIdx(int idx) const { return indices(idx); }
  __device__ inline BoutReal operator()(int idx) const { return op(lhs(idx), rhs(idx)); }

  Mesh* getMesh() const { return mesh; }
  CELL_LOC getLocation() const { return location; }
  DirectionTypes getDirections() const { return directions; }
};

#endif // BOUT_EXPRESSION_HXX