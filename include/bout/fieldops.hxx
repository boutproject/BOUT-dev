#pragma once
#include "bout/array.hxx"
#include <unordered_map>
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
struct Assign {
  int scale = 1;
  int offset = 0;
  template <typename Expr>
  __device__ void operator()(int idx, BoutReal* out, const Expr& expr) const {
    out[(idx * scale) + offset] = expr.lhs(idx) + expr.rhs(idx);
  }
};

struct Add {
  template <typename LView, typename RView>
  __device__ __forceinline__ BoutReal operator()(int idx, const LView& L,
                                                 const RView& R) const {
    return L(idx) + R(idx);
  }
  __device__ __forceinline__ BoutReal operator()(BoutReal a, BoutReal b) const {
    return a + b;
  }
};
  struct Sub {
    template <typename LView, typename RView>
    __device__ __forceinline__ BoutReal operator()(int idx, const LView& L,
                                                   const RView& R) const {
      return L(idx) - R(idx);
    }
    __device__ __forceinline__ BoutReal operator()(BoutReal a, BoutReal b) const {
      return a - b;
    }
  };
  struct Mul {
    template <typename LView, typename RView>
    __device__ __forceinline__ BoutReal operator()(int idx, const LView& L,
                                                   const RView& R) const {
      return L(idx) * R(idx);
    }
    __device__ __forceinline__ BoutReal operator()(BoutReal a, BoutReal b) const {
      return a * b;
    }
  };
  struct Div {
    template <typename LView, typename RView>
    __device__ __forceinline__ BoutReal operator()(int idx, const LView& L,
                                                   const RView& R) const {
      return L(idx) / R(idx);
    }
    __device__ __forceinline__ BoutReal operator()(BoutReal a, BoutReal b) const {
      return a / b;
    }
  };
};
};

template <typename Expr>
__global__ __launch_bounds__(256) static void evaluatorExpr(BoutReal* out,
                                                            const Expr expr) {
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  if (tid >= expr.size()) {
    return; // Out of bounds
  }
  int idx = expr.regionIdx(tid);
  out[idx] = expr(idx); // single‐pass fusion
  //int stride = blockDim.x * gridDim.x;
  //for (int i = tid, e = expr.size(); i < e; i += stride) {
  //  int idx = expr.regionIdx(i);
  //  out[idx] = expr(idx); // single‐pass fusion
  //}
}

template <typename Result, typename Expr>
__global__ __launch_bounds__(256) static void evaluatorExprWithResult(Result res,
                                                                      const Expr expr) {
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  if (tid >= expr.size()) {
    return; // Out of bounds
  }
  int idx = expr.regionIdx(tid);
  res[idx] = expr(idx); // single‐pass fusion
  //res(idx, expr(idx)); // single‐pass fusion
  //res(idx) = expr(idx); // single‐pass fusion
  //int stride = blockDim.x * gridDim.x;
  //for (int i = tid, e = expr.size(); i < e; i += stride) {
  //  int idx = expr.regionIdx(i);
  //  out[idx] = expr(idx); // single‐pass fusion
  //}
}

inline std::unordered_map<void*, Array<int>> regionIndicesCache;

template <typename L, typename R, typename Func>
struct BinaryExpr {
  typename L::View lhs;
  typename R::View rhs;
  Array<int> indices;
  Func f;

  Mesh* mesh;
  CELL_LOC location = CELL_CENTRE;
  DirectionTypes directions;
  std::optional<size_t> regionID;

  template <typename IndType>
  BinaryExpr(const typename L::View& lhs, const typename R::View& rhs, Func f, Mesh* mesh,
             CELL_LOC location, DirectionTypes directions, std::optional<size_t> regionID,
             const Region<IndType>& region)
      //: lhs(static_cast<typename L::View>(lhs)), rhs(static_cast<typename R::View>(rhs)),
      : lhs(lhs), rhs(rhs), f(f), mesh(mesh), location(location), directions(directions),
        regionID(regionID), indices(region.getIndices().size()) {
    // Copy the region indices into the managed array
    for (int i = 0; i < indices.size(); ++i) {
      indices[i] = region.getIndices()[i].ind;
    }
    //if (regionIndicesCache.find(static_cast<void*>(const_cast<Region<IndType>*>(&region)))
    //    != regionIndicesCache.end()) {
    //  // If we have already computed the indices for this region, use them
    //  indices =
    //      regionIndicesCache[static_cast<void*>(const_cast<Region<IndType>*>(&region))];
    //} else {
    //  // Otherwise, compute the indices and store them in the cache
    //  indices = Array<int>(region.getIndices().size());
    //  // Copy the region indices into the managed array
    //  for (int i = 0; i < indices.size(); ++i) {
    //    indices[i] = region.getIndices()[i].ind;
    //  }
    //  regionIndicesCache[static_cast<void*>(const_cast<Region<IndType>*>(&region))] =
    //      indices;
    //}
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
    int scale = 1;
    int offset = 0;

    View& setScale(int s) {
      scale = s;
      return *this;
    }
    View& setOffset(int o) {
      offset = o;
      return *this;
    }

    __device__ __forceinline__ int size() const { return num_indices; }
    __device__ __forceinline__ int regionIdx(int idx) const { return indices[idx]; }
    __device__ __forceinline__ BoutReal operator()(int idx) const {
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

  template <typename Result>
  void evaluateWithResult(const Result& res) const {
    constexpr int THREADS = 256;
    int blocks = (size() + THREADS - 1) / THREADS;
    evaluatorExprWithResult<<<blocks, THREADS>>>(res, static_cast<View>(*this));
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