#pragma once
#ifndef BOUT_FIELDOPS_HXX
#define BOUT_FIELDOPS_HXX

#include "bout/bout_types.hxx"

#include <cuda_runtime.h>
#include <optional>

class Mesh;
class Field3D;

#include <type_traits>

namespace bout {
namespace op {
  struct Add {
    template<typename LView, typename RView>
    __host__ __device__ inline BoutReal operator()(int idx, const LView &L, const RView &R) const {
      return L(idx) + R(idx);
    }
    __host__ __device__ inline BoutReal operator()(BoutReal a, BoutReal b) const { return a + b; }
  };
  struct Sub {
    template<typename LView, typename RView>
    __host__ __device__ inline BoutReal operator()(int idx, const LView &L, const RView &R) const { return L(idx) - R(idx); }
    __host__ __device__ inline BoutReal operator()(BoutReal a, BoutReal b) const { return a - b; }
  };                                                                                                
  struct Mul {                                                                                      
    template<typename LView, typename RView>
    __host__ __device__ inline BoutReal operator()(int idx, const LView &L, const RView &R) const { return L(idx) * R(idx); }
    __host__ __device__ inline BoutReal operator()(BoutReal a, BoutReal b) const { return a * b; }
  };                                                                                                
  struct Div {                                                                                      
    template<typename LView, typename RView>
    __host__ __device__ inline BoutReal operator()(int idx, const LView &L, const RView &R) const { return L(idx) / R(idx); }
    __host__ __device__ inline BoutReal operator()(BoutReal a, BoutReal b) const { return a / b; }
  };
};
};

template <typename Expr>
__global__ static void evaluatorExpr(BoutReal* out, const Expr& expr) {
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  int stride = blockDim.x * gridDim.x;
  for (int i = tid; i < expr.size(); i += stride) {
    out[expr.regionIdx(i)] = expr(expr.regionIdx(i)); // single‐pass fusion
  }
}

template <typename L, typename R, typename Func>
struct BinaryExpr {
  const L &LHS;
  const R &RHS;
  typename L::View lhs;
  typename R::View rhs;
  Array<int> indices;
  Func f;

  Mesh* mesh;
  CELL_LOC location = CELL_CENTRE;
  DirectionTypes directions;
  std::optional<size_t> regionID;

  template <typename IndType>
  BinaryExpr(const L &lhs, const R &rhs, Func f, Mesh* mesh, CELL_LOC location,
             DirectionTypes directions, std::optional<size_t> regionID,
             const Region<IndType>& region)
      : LHS(lhs), RHS(rhs), lhs(static_cast<typename L::View>(lhs)), rhs(static_cast<typename R::View>(rhs)),
        f(f), mesh(mesh), location(location), directions(directions), regionID(regionID),
        indices(region.getIndices().size()) {
    // Copy the region indices into the managed array
    for (int i = 0; i < indices.size(); ++i) {
      indices[i] = region.getIndices()[i].ind;
    }
  }

  inline int size() const { return indices.size(); }
  inline BoutReal operator()(int idx) const {
    return f(idx, lhs, rhs); // single‐pass fusion
  }
  inline int regionIdx(int idx) const { return indices[idx]; }

  struct View {
    typename L::View lhs;
    typename R::View rhs;
    const int* indices;
    int num_indices;
    Func f;

    __device__ inline int size() const { return num_indices; }
    __device__ inline int regionIdx(int idx) const { return indices[idx]; }
    __device__ inline BoutReal operator()(int idx) const {
      return f(idx, lhs, rhs); // single‐pass fusion
      //return f(lhs(idx), rhs(idx)); // single‐pass fusion
    }
  };

  operator View() { return View{lhs, rhs, &indices[0], indices.size(), f}; }
  operator View() const { return View{lhs, rhs, &indices[0], indices.size(), f}; }

  void evaluate(BoutReal* data) const {
    constexpr int THREADS = 256;
    int blocks = (size() + THREADS - 1) / THREADS;
    evaluatorExpr<<<blocks, THREADS>>>(&data[0], static_cast<View>(*this));
    cudaDeviceSynchronize();
    //for(int i=0; i<size(); ++i) {
    //  data[regionIdx(i)] = f(i, lhs, rhs); // single‐pass fusion
    //}
  }

  Mesh* getMesh() const { return mesh; }
  CELL_LOC getLocation() const { return location; }
  DirectionTypes getDirections() const { return directions; }
  std::optional<size_t> getRegionID() const { return regionID; };
};

#endif // BOUT_EXPRESSION_HXX