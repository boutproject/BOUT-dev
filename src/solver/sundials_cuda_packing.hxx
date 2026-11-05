/**************************************************************************
 * CUDA packing for the SUNDIALS MPIPlusX backend.
 * Copyright 2026 BOUT++ contributors. Licensed under LGPL-3.0-or-later.
 **************************************************************************/
#ifndef BOUT_SUNDIALS_CUDA_PACKING_HXX
#define BOUT_SUNDIALS_CUDA_PACKING_HXX

#include "bout/build_defines.hxx"

#if BOUT_HAS_SUNDIALS_CUDA
#include "bout/bout_types.hxx"
#include "bout/boutexception.hxx"
#include <nvector/nvector_cuda.h>
#include <nvector/nvector_mpiplusx.h>
#include <type_traits>
#include <utility>
#include <vector>

namespace bout::detail {
// The map preserves Solver::loop_vars ordering, including per-field boundaries.
static __global__ void sundials_cuda_pack(BoutReal* field, BoutReal* vector,
                                          const int* map, int count, bool load) {
  const int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i < count) {
    const int index = map[2 * i];
    const int slot = map[2 * i + 1];
    if (load) {
      field[index] = vector[slot];
    } else {
      vector[slot] = field[index];
    }
  }
}
} // namespace bout::detail

struct SundialsCudaPacking {
  static_assert(std::is_same_v<BoutReal, sunrealtype>,
                "BOUT++ and SUNDIALS must use matching floating point types");

  static void check(cudaError_t error) {
    if (error != cudaSuccess) {
      throw BoutException("SUNDIALS CUDA packing: {}", cudaGetErrorString(error));
    }
  }

  struct FieldMap {
    std::vector<int> host;
    int* device = nullptr;
    int count = 0;

    FieldMap() = default;
    FieldMap(const FieldMap&) = delete;
    FieldMap& operator=(const FieldMap&) = delete;
    FieldMap(FieldMap&& other) noexcept
        : host(std::move(other.host)), device(std::exchange(other.device, nullptr)),
          count(other.count) {}
    ~FieldMap() {
      if (device) {
        cudaFree(device);
      }
    }

    void append(int index, int slot) {
      host.push_back(index);
      host.push_back(slot);
    }

    void upload() {
      count = host.size() / 2;
      if (count == 0) {
        return;
      }
      check(cudaMalloc(reinterpret_cast<void**>(&device), host.size() * sizeof(int)));
      check(cudaMemcpy(device, host.data(), host.size() * sizeof(int),
                       cudaMemcpyHostToDevice));
      std::vector<int>().swap(host);
    }

    void transfer(BoutReal* field, BoutReal* vector, bool load) const {
      if (count == 0) {
        return;
      }
      // BOUT++ CUDA fields use Umpire pinned memory. Resolve its device alias
      // rather than assuming the host and device virtual addresses are equal.
      cudaPointerAttributes attributes{};
      check(cudaPointerGetAttributes(&attributes, field));
      if (!attributes.devicePointer) {
        throw BoutException(
            "SUNDIALS CUDA packing requires GPU-accessible field storage");
      }
      auto* device_field = static_cast<BoutReal*>(attributes.devicePointer);
      bout::detail::sundials_cuda_pack<<<(count + 255) / 256, 256>>>(device_field, vector,
                                                                     device, count, load);
      check(cudaGetLastError());
    }
  };

  bool ready = false;
  std::vector<FieldMap> fields;
};
#else
// Keep the shared interface usable in CPU-only and older SUNDIALS builds.
struct SundialsCudaPacking {};
#endif

#endif
