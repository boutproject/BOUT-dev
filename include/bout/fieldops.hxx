#pragma once
#ifndef BOUT_FIELDOPS_HXX
#define BOUT_FIELDOPS_HXX

#include "bout/array.hxx"
#include "bout/assert.hxx"
#include "bout/bout_types.hxx"
#include "bout/build_config.hxx"
#include "bout/build_defines.hxx"
#include "bout/region.hxx"

#include <cstddef>
#include <limits>
#include <optional>
#include <type_traits>

#if BOUT_HAS_CUDA
#include <cuda_runtime.h>
#endif

class Mesh;
class Field3D;
class Field3DParallel;
class Field2D;
class FieldPerp;

namespace bout::detail {
const Region<Ind3D>& getField3DRegion(const Mesh* mesh, std::optional<size_t> regionID);
}

template <typename T>
struct is_expr_field2d : std::false_type {};

template <>
struct is_expr_field2d<Field2D> : std::true_type {};

template <typename T>
inline constexpr bool is_expr_field2d_v = is_expr_field2d<std::decay_t<T>>::value;

// Base template: nothing is an expression by default
template <typename T>
struct is_expr_field3d : std::false_type {};

template <>
struct is_expr_field3d<Field3D> : std::true_type {};

template <>
struct is_expr_field3d<Field3DParallel> : std::true_type {};

template <typename T>
struct is_expr_fieldperp : std::false_type {};

template <>
struct is_expr_fieldperp<FieldPerp> : std::true_type {};

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
  BOUT_HOST_DEVICE void operator()(int idx, BoutReal* out, const Expr& expr) const {
    out[(idx * scale) + offset] = expr.lhs(idx) + expr.rhs(idx);
  }
};

struct Add {
  template <typename LView, typename RView>
  BOUT_HOST_DEVICE BOUT_FORCEINLINE BoutReal operator()(int idx, const LView& L,
                                                        const RView& R) const {
    return L(idx) + R(idx);
  }
  BOUT_HOST_DEVICE BOUT_FORCEINLINE BoutReal operator()(BoutReal a, BoutReal b) const {
    return a + b;
  }
};
struct Sub {
  template <typename LView, typename RView>
  BOUT_HOST_DEVICE BOUT_FORCEINLINE BoutReal operator()(int idx, const LView& L,
                                                        const RView& R) const {
    return L(idx) - R(idx);
  }
  BOUT_HOST_DEVICE BOUT_FORCEINLINE BoutReal operator()(BoutReal a, BoutReal b) const {
    return a - b;
  }
};
struct Mul {
  template <typename LView, typename RView>
  BOUT_HOST_DEVICE BOUT_FORCEINLINE BoutReal operator()(int idx, const LView& L,
                                                        const RView& R) const {
    return L(idx) * R(idx);
  }
  BOUT_HOST_DEVICE BOUT_FORCEINLINE BoutReal operator()(BoutReal a, BoutReal b) const {
    return a * b;
  }
};
struct Div {
  template <typename LView, typename RView>
  BOUT_HOST_DEVICE BOUT_FORCEINLINE BoutReal operator()(int idx, const LView& L,
                                                        const RView& R) const {
    return L(idx) / R(idx);
  }
  BOUT_HOST_DEVICE BOUT_FORCEINLINE BoutReal operator()(BoutReal a, BoutReal b) const {
    return a / b;
  }
};
struct IfElse {
  bool condition;

  template <typename LView, typename RView>
  BOUT_HOST_DEVICE BOUT_FORCEINLINE BoutReal operator()(int idx, const LView& L,
                                                        const RView& R) const {
    return condition ? L(idx) : R(idx);
  }
  BOUT_HOST_DEVICE BOUT_FORCEINLINE BoutReal operator()(BoutReal a, BoutReal b) const {
    return condition ? a : b;
  }
};
}; // namespace op

namespace reduce {

struct Min {
  struct State {
    BoutReal value;
  };

  BOUT_HOST_DEVICE static State identity() {
    return {std::numeric_limits<BoutReal>::infinity()};
  }
  BOUT_HOST_DEVICE static void accumulate(State& state, BoutReal value) {
    state.value = value < state.value ? value : state.value;
  }
  BOUT_HOST_DEVICE static void combine(State& state, const State& other) {
    state.value = other.value < state.value ? other.value : state.value;
  }
  static BoutReal finalize(const State& state) { return state.value; }
};

struct Max {
  struct State {
    BoutReal value;
  };

  BOUT_HOST_DEVICE static State identity() {
    return {-std::numeric_limits<BoutReal>::infinity()};
  }
  BOUT_HOST_DEVICE static void accumulate(State& state, BoutReal value) {
    state.value = value > state.value ? value : state.value;
  }
  BOUT_HOST_DEVICE static void combine(State& state, const State& other) {
    state.value = other.value > state.value ? other.value : state.value;
  }
  static BoutReal finalize(const State& state) { return state.value; }
};

struct Mean {
  struct State {
    BoutReal sum;
    int count;
  };

  BOUT_HOST_DEVICE static State identity() { return {0.0, 0}; }
  BOUT_HOST_DEVICE static void accumulate(State& state, BoutReal value) {
    state.sum += value;
    state.count += 1;
  }
  BOUT_HOST_DEVICE static void combine(State& state, const State& other) {
    state.sum += other.sum;
    state.count += other.count;
  }
  static BoutReal finalize(const State& state) {
    return state.sum / static_cast<BoutReal>(state.count);
  }
};

} // namespace reduce
}; // namespace bout

template <typename ExprView>
struct ReductionView {
  ExprView expr;
  const int* indices;
  int num_indices;

  BOUT_HOST_DEVICE BOUT_FORCEINLINE int size() const { return num_indices; }
  BOUT_HOST_DEVICE BOUT_FORCEINLINE BoutReal valueAtRegionPos(int idx) const {
    return expr(indices[idx]);
  }
};

template <typename ExprView>
ReductionView<ExprView> makeReductionView(const ExprView& expr,
                                          const Array<int>& indices) {
  return ReductionView<ExprView>{expr, indices.size() > 0 ? &indices[0] : nullptr,
                                 indices.size()};
}

#if BOUT_HAS_CUDA && defined(__CUDACC__)
template <typename Expr>
__global__ void __launch_bounds__(THREADS) evaluatorExpr(BoutReal* out, const Expr expr) {
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  int e = expr.size();

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

template <typename Reducer, typename ExprView>
__global__ void __launch_bounds__(THREADS)
    reducerExpr(typename Reducer::State* partials, const ExprView expr) {
  using State = typename Reducer::State;

  __shared__ State shared[THREADS];

  const int tid = threadIdx.x;
  const int global = blockIdx.x * blockDim.x + tid;
  const int stride = blockDim.x * gridDim.x;

  State local = Reducer::identity();

  for (int i = global; i < expr.size(); i += stride) {
    Reducer::accumulate(local, expr.valueAtRegionPos(i));
  }

  shared[tid] = local;
  __syncthreads();

  for (int offset = blockDim.x / 2; offset > 0; offset /= 2) {
    if (tid < offset) {
      Reducer::combine(shared[tid], shared[tid + offset]);
    }
    __syncthreads();
  }

  if (tid == 0) {
    partials[blockIdx.x] = shared[0];
  }
}
#endif

#if BOUT_HAS_CUDA && defined(__CUDACC__)
struct StreamsRAII {
  std::vector<cudaStream_t> streams;

  cudaStream_t get() {
    cudaStream_t stream = 0;

    if (streams.empty()) {
      if (cudaStreamCreate(&stream) != cudaSuccess) {
        throw BoutException("Failed to create CUDA stream");
      }
    } else {
      stream = streams.back();
      streams.pop_back();
    }

    return stream;
  }

  void put(cudaStream_t stream) { streams.push_back(stream); }

  ~StreamsRAII() {
    for (auto& stream : streams) {
      cudaStreamDestroy(stream);
    }
  }

  StreamsRAII() = default;
  StreamsRAII(const StreamsRAII&) = delete;
  StreamsRAII(StreamsRAII&&) = delete;
  StreamsRAII& operator=(const StreamsRAII&) = delete;
  StreamsRAII& operator=(StreamsRAII&&) = delete;
};
inline struct StreamsRAII streams;
#endif

template <typename Reducer, typename ExprView>
auto reduceExpr(const ExprView& expr_view) -> typename Reducer::State {
  using State = typename Reducer::State;

  ASSERT1(expr_view.size() > 0);

#if BOUT_HAS_CUDA && defined(__CUDACC__)
  cudaStream_t stream = streams.get();
  int blocks = (expr_view.size() + THREADS - 1) / THREADS;
  blocks = blocks < 1024 ? blocks : 1024;
  Array<State> partials(blocks);

  reducerExpr<Reducer><<<blocks, THREADS, 0, stream>>>(&partials[0], expr_view);
  cudaStreamSynchronize(stream);
  streams.put(stream);

  State result = Reducer::identity();
  for (int i = 0; i < blocks; ++i) {
    Reducer::combine(result, partials[i]);
  }
  return result;
#else
  State result = Reducer::identity();
  for (int i = 0; i < expr_view.size(); ++i) {
    Reducer::accumulate(result, expr_view.valueAtRegionPos(i));
  }
  return result;
#endif
}

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
  std::optional<int> yindex;

  template <typename IndType>
  BinaryExpr(const typename L::View& lhs, const typename R::View& rhs, Func f, Mesh* mesh,
             CELL_LOC location, DirectionTypes directions, std::optional<size_t> regionID,
             const Region<IndType>& region, std::optional<int> yindex = std::nullopt)
      : lhs(lhs), rhs(rhs), indices(region.getLinearIndices()), f(f), mesh(mesh),
        location(location), directions(directions), regionID(regionID), yindex(yindex) {}

  BinaryExpr(const typename L::View& lhs, const typename R::View& rhs, Func f, Mesh* mesh,
             CELL_LOC location, DirectionTypes directions, std::optional<size_t> regionID,
             const Array<int>& indices, std::optional<int> yindex = std::nullopt)
      : lhs(lhs), rhs(rhs), indices(indices), f(f), mesh(mesh), location(location),
        directions(directions), regionID(regionID), yindex(yindex) {}

  BinaryExpr& operator=(const BinaryExpr&) = delete;
  BinaryExpr& operator=(BinaryExpr&&) = delete;

  BOUT_HOST_DEVICE BOUT_FORCEINLINE int size() const { return indices.size(); }
  BOUT_HOST_DEVICE BOUT_FORCEINLINE BoutReal operator()(int idx) const {
    return f(idx, lhs, rhs); // single‐pass fusion
  }
  template <typename IndType>
  BOUT_HOST_DEVICE BOUT_FORCEINLINE auto operator[](const IndType& d) const
      -> decltype(d.ind, BoutReal{}) {
    if constexpr (std::is_same_v<ResT, Field2D>) {
      return operator()(d.ind / d.nz);
    } else {
      return operator()(d.ind);
    }
  }
  BOUT_HOST_DEVICE BOUT_FORCEINLINE int regionIdx(int idx) const { return indices[idx]; }

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
    BOUT_HOST_DEVICE BOUT_FORCEINLINE int size() const { return num_indices; }
    BOUT_HOST_DEVICE BOUT_FORCEINLINE int regionIdx(int idx) const {
      return indices[idx];
    }
    BOUT_HOST_DEVICE BOUT_FORCEINLINE BoutReal operator()(int idx) const {
      return f((idx * mul) / div, lhs, rhs); // single‐pass fusion
      //return f(lhs((idx * mul) / div), rhs((idx * mul) / div)); // single‐pass fusion
    }
  };

  operator View() { return View{lhs, rhs, &indices[0], indices.size(), f}; }
  operator View() const { return View{lhs, rhs, &indices[0], indices.size(), f}; }

  void evaluate(BoutReal* data) const {
#if BOUT_HAS_CUDA && defined(__CUDACC__)
    cudaStream_t stream = streams.get();
    int blocks = (size() + THREADS - 1) / THREADS;
    evaluatorExpr<<<blocks, THREADS, 0, stream>>>(&data[0], static_cast<View>(*this));
    cudaStreamSynchronize(stream);
    streams.put(stream);
#else
    if constexpr (std::is_same_v<ResT, Field3D>) {
      const auto& region = bout::detail::getField3DRegion(mesh, regionID);
      for (auto block = region.getBlocks().cbegin(), end = region.getBlocks().cend();
           block < end; ++block) {
        const int block_begin = block->first.ind;
        const int block_end = block->second.ind;

        BOUT_ASSUME(block_begin >= 0);
        BOUT_ASSUME(block_begin <= block_end);

        for (int idx = block_begin; idx < block_end; ++idx) {
          data[idx] = operator()(idx);
        }
      }
    } else {
      int e = size();
      for (int i = 0; i < e; ++i) {
        int idx = regionIdx(i);
        data[idx] = operator()(idx); // single‐pass fusion
      }
    }
#endif
  }

  Mesh* getMesh() const { return mesh; }
  CELL_LOC getLocation() const { return location; }
  DirectionTypes getDirections() const { return directions; }
  std::optional<size_t> getRegionID() const { return regionID; };
  int getIndex() const { return yindex.value_or(-1); }
};

#endif // BOUT_FIELDSOPS_HXX
