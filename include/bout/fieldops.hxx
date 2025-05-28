#pragma once
#ifndef BOUT_FIELDOPS_HXX
#define BOUT_FIELDOPS_HXX

#include "bout/bout_types.hxx"

#include <cuda_runtime.h>
#include <optional>

class Mesh;
class Field3D;

#include <type_traits>

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

template <typename Expr>
__global__ static void evaluatorExpr(BoutReal* out, const Expr& expr) {
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  int stride = blockDim.x * gridDim.x;
  for (int i = tid; i < expr.getSize(); i += stride) {
    out[expr.regionIdx(i)] = expr(expr.regionIdx(i)); // single‐pass fusion
  }
}

template <typename L, typename R>
struct BinaryExpr {
  enum class Op { ADD, SUB, MUL, DIV };
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

  L lhs;
  R rhs;
  RegionIndices indices;
  Op op;

  Mesh* mesh;
  CELL_LOC location = CELL_CENTRE;
  DirectionTypes directions;
  std::optional<size_t> regionID;

  template <typename IndType>
  BinaryExpr(L lhs, R rhs, Op op, Mesh* mesh, CELL_LOC location,
             DirectionTypes directions, std::optional<size_t> regionID,
             const Region<IndType>& region)
      : lhs(lhs), rhs(rhs), op(op), mesh(mesh), location(location),
        directions(directions), regionID(regionID), indices(region.getIndices().size()) {
    // Copy the region indices into the managed array
    for (int i = 0; i < indices.size; ++i) {
      indices.data[i] = region.getIndices()[i].ind;
    }
  }

  __host__ inline int getSize() const { return indices.size; }

  struct View {
    L lhs;
    R rhs;
    int* indices;
    int size;
    Op op;

    __host__ __device__ inline int getSize() const { return size; }
    __device__ inline int regionIdx(int idx) const { return indices[idx]; }
    __device__ inline BoutReal operator()(int idx) const {
      switch (op) {
      case Op::ADD:
        return Add{}(lhs(idx), rhs(idx));
      case Op::SUB:
        return Sub{}(lhs(idx), rhs(idx));
      case Op::MUL:
        return Mul{}(lhs(idx), rhs(idx));
      case Op::DIV:
        return Div{}(lhs(idx), rhs(idx));
      }
    }
  };

  operator View() { return View{lhs, rhs, indices.data, indices.size, op}; }
  operator View() const { return View{lhs, rhs, indices.data, indices.size, op}; }

  void evaluate(BoutReal* data) const {}

  Mesh* getMesh() const { return mesh; }
  CELL_LOC getLocation() const { return location; }
  DirectionTypes getDirections() const { return directions; }
  std::optional<size_t> getRegionID() const { return regionID; };
};

//template <typename T>
//struct Expr {
//  using type = T;
//};
//
//template <>
//struct Expr<Field3D> {
//  using type = Field3D::View;
//};

// 1) detect our BinaryExpr<T,U> template
template <typename>
struct is_binary_expr : std::false_type {};
template <typename A, typename B>
struct is_binary_expr<BinaryExpr<A, B>> : std::true_type {};

// 2) detect “any subclass of Field”
//    assuming Field is your common base class
template <typename T>
constexpr bool is_field_v = std::is_base_of<Field3D, std::decay_t<T>>::value;

// 3) combine into “is one of our expression types”
template <typename T>
constexpr bool is_expr_v = is_field_v<T> || is_binary_expr<std::decay_t<T>>::value;

#endif // BOUT_EXPRESSION_HXX