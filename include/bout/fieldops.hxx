#pragma once
#ifndef BOUT_FIELDOPS_HXX
#define BOUT_FIELDOPS_HXX

#include "bout/array.hxx"
#include "bout/bout_types.hxx"

#include <cuda_runtime.h>
#include <optional>
#include <type_traits>
#include <unordered_map>

class Mesh;
class Field3D;
class Field2D;

template <typename T>
struct is_expr_field2d : std::false_type {};

template <typename T>
inline constexpr bool is_expr_field2d_v = is_expr_field2d<std::decay_t<T>>::value;

// Base template: nothing is an expression by default
template <typename T>
struct is_expr_field3d : std::false_type {};

template <typename T>
struct is_expr_fieldperp : std::false_type {};

template <typename T>
inline constexpr bool is_expr_fieldperp_v = is_expr_fieldperp<std::decay_t<T>>::value;

// Helper variable template
template <typename T>
inline constexpr bool is_expr_field3d_v = is_expr_field3d<std::decay_t<T>>::value;

template <typename T>
struct is_expr_constant : std::bool_constant<std::is_arithmetic_v<T>> {};

template <typename T>
inline constexpr bool is_expr_constant_v = is_expr_constant<std::decay_t<T>>::value;

template <typename T>
struct is_expr_constant<Constant<T>>
    : std::integral_constant<bool, is_expr_constant_v<std::decay_t<T>>> {};

constexpr int THREADS = 128;
namespace bout {
namespace op {
struct Assign {
  int scale = 1;
  int offset = 0;
  template <typename Expr>
  __host__ __device__ void operator()(int idx, BoutReal* out, const Expr& expr) const {
    out[(idx * scale) + offset] = expr.lhs(idx) + expr.rhs(idx);
  }
};

struct Add {
  template <typename LView, typename RView>
  __host__ __device__ __forceinline__ BoutReal operator()(int idx, const LView& L,
                                                          const RView& R) const {
    return L(idx) + R(idx);
  }
  __host__ __device__ __forceinline__ BoutReal operator()(BoutReal a, BoutReal b) const {
    return a + b;
  }
};
  struct Sub {
    template <typename LView, typename RView>
    __host__ __device__ __forceinline__ BoutReal operator()(int idx, const LView& L,
                                                            const RView& R) const {
      return L(idx) - R(idx);
    }
    __host__ __device__ __forceinline__ BoutReal operator()(BoutReal a,
                                                            BoutReal b) const {
      return a - b;
    }
  };
  struct Mul {
    template <typename LView, typename RView>
    __host__ __device__ __forceinline__ BoutReal operator()(int idx, const LView& L,
                                                            const RView& R) const {
      return L(idx) * R(idx);
    }
    __host__ __device__ __forceinline__ BoutReal operator()(BoutReal a,
                                                            BoutReal b) const {
      return a * b;
    }
  };
  struct Div {
    template <typename LView, typename RView>
    __host__ __device__ __forceinline__ BoutReal operator()(int idx, const LView& L,
                                                            const RView& R) const {
      return L(idx) / R(idx);
    }
    __host__ __device__ __forceinline__ BoutReal operator()(BoutReal a,
                                                            BoutReal b) const {
      return a / b;
    }
  };
};
};

template <typename Expr>
__global__ void __launch_bounds__(THREADS) evaluatorExpr(BoutReal* out, const Expr expr) {
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  int e = expr.size();

  // In-bounds version
  //if (tid < e) {
  //  int idx = expr.regionIdx(tid);
  //  out[idx] = expr(idx); // single‐pass fusion
  //}

  // Out-of-bounds version
  if (tid >= e) {
    return;
  }
  int idx = expr.regionIdx(tid);
  out[idx] = expr(idx); // single‐pass fusion

  // Grid-strided loop
  //int stride = blockDim.x * gridDim.x;
  //for (int i = tid; i < e; i += stride) {
  //  int idx = expr.regionIdx(i);
  //  out[idx] = expr(idx); // single‐pass fusion
  //}
}

inline std::unordered_map<void*, Array<int>> regionIndicesCache;

template <typename ResT, typename L, typename R, typename Func>
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
    //std::cout << "===PRE-sorting indices\n";
    //for (auto& ind : indices) {
    //  std::cout << ind << " ";
    //}
    //std::cout << "===end PRE\n";
    //std::sort(indices.begin(), indices.end(),
    //          [](const auto& a, const auto& b) { return a < b; });
    //std::cout << "===POST-sorting indices\n";
    //for (auto& ind : indices) {
    //  std::cout << ind << " ";
    //}
    //std::cout << "===end POST\n";
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

  BinaryExpr& operator=(BinaryExpr const&) = delete;
  BinaryExpr& operator=(BinaryExpr&&) = delete;

  inline int size() const { return indices.size(); }
  inline BoutReal operator()(int idx) const {
    return f(idx, lhs, rhs); // single‐pass fusion
  }
  inline int regionIdx(int idx) const { return indices[idx]; }

  //operator ResT() { return ResT{*this}; }
  struct View {
    typename L::View lhs;
    typename R::View rhs;
    const int* indices;
    int num_indices;
    Func f;
    int mul = 1;
    int div = 1;

    View& setScale(int mul, int div) {
      this->mul = mul;
      this->div = div;
      return *this;
    }
    __host__ __device__ __forceinline__ int size() const { return num_indices; }
    __host__ __device__ __forceinline__ int regionIdx(int idx) const {
      return indices[idx];
    }
    __host__ __device__ __forceinline__ BoutReal operator()(int idx) const {
      return f((idx * mul) / div, lhs, rhs); // single‐pass fusion
      //return f(lhs((idx * mul) / div), rhs((idx * mul) / div)); // single‐pass fusion
    }
  };

  operator View() { return View{lhs, rhs, &indices[0], indices.size(), f}; }
  operator View() const { return View{lhs, rhs, &indices[0], indices.size(), f}; }

  void evaluate(BoutReal* data) const {
    int blocks = (size() + THREADS - 1) / THREADS;
    evaluatorExpr<<<blocks, THREADS>>>(&data[0], static_cast<View>(*this));
    cudaDeviceSynchronize();
    // OpenMP impl.
    //int e = size();
    //#pragma omp parallel for
    //for (int i = 0; i < e; ++i) {
    //  int idx = regionIdx(i);
    //  data[idx] = operator()(idx); // single‐pass fusion
    //}
  }

  Mesh* getMesh() const { return mesh; }
  CELL_LOC getLocation() const { return location; }
  DirectionTypes getDirections() const { return directions; }
  std::optional<size_t> getRegionID() const { return regionID; };
};

#endif // BOUT_EXPRESSION_HX