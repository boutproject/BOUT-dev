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
#include "bout/field2d.hxx"
#include "bout/field3d.hxx"
#include "bout/profiling.hxx"
#include <memory>
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

  SundialsCudaPacking() = default;
  SundialsCudaPacking(const SundialsCudaPacking&) = delete;
  SundialsCudaPacking& operator=(const SundialsCudaPacking&) = delete;

  ~SundialsCudaPacking() {
    synchronize_noexcept();
    if (fields_read) {
      cudaEventDestroy(fields_read);
    }
    if (stream) {
      cudaStreamDestroy(stream);
    }
  }

  // All SUNDIALS streaming/reduction operations and packing use this stream.
  // N_VClone_Cuda clones both policies, so CVODE's internal vectors inherit it.
  void bind(N_Vector local) {
    if (!stream) {
      check(cudaStreamCreateWithFlags(&stream, cudaStreamNonBlocking));
    }
    if (!fields_read) {
      check(cudaEventCreateWithFlags(&fields_read, cudaEventDisableTiming));
    }
    auto* content = static_cast<N_VectorContent_Cuda>(local->content);
    std::unique_ptr<SUNCudaExecPolicy> streaming(
        content->stream_exec_policy->clone_new_stream(stream));
    std::unique_ptr<SUNCudaExecPolicy> reduction(
        content->reduce_exec_policy->clone_new_stream(stream));
    if (N_VSetKernelExecPolicy_Cuda(local, streaming.get(), reduction.get()) != 0) {
      throw BoutException("SUNDIALS CUDA execution policy setup failed");
    }
  }

  void check_stream(N_Vector local) const {
    const auto* content = static_cast<N_VectorContent_Cuda>(local->content);
    if (*content->stream_exec_policy->stream() != stream
        || *content->reduce_exec_policy->stream() != stream) {
      throw BoutException("SUNDIALS CUDA vector does not use the packing stream");
    }
  }

  // An event covers field reads only. Subsequent SUNDIALS vector operations may
  // keep running while the host prepares fields for the next callback.
  void record_field_reads() {
    check(cudaEventRecord(fields_read, stream));
    pending_field_reads = true;
  }

  void wait_for_fields() {
    if (pending_field_reads) {
      const bout::profiling::ProfileRange range("SundialsCudaPacking::wait_for_fields");
      check(cudaEventSynchronize(fields_read));
      pending_field_reads = false;
      retained_2d.clear();
      retained_3d.clear();
    }
  }

  void synchronize() {
    if (stream) {
      const bout::profiling::ProfileRange range("SundialsCudaPacking::wait_stream");
      check(cudaStreamSynchronize(stream));
      pending_field_reads = false;
      retained_2d.clear();
      retained_3d.clear();
    }
  }

  // Destruction and exception cleanup must finish queued accesses before field
  // storage, packing maps, or SUNDIALS vectors can be released.
  void synchronize_noexcept() noexcept {
    if (stream) {
      static_cast<void>(cudaStreamSynchronize(stream));
      pending_field_reads = false;
      retained_2d.clear();
      retained_3d.clear();
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

    void transfer(BoutReal* field, BoutReal* vector, bool load,
                  cudaStream_t stream) const {
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
      bout::detail::sundials_cuda_pack<<<(count + 255) / 256, 256, 0, stream>>>(
          device_field, vector, device, count, load);
      check(cudaGetLastError());
    }
  };

  bool ready = false;
  bool async_transfers = true;
  cudaStream_t stream = nullptr;
  cudaEvent_t fields_read = nullptr;
  bool pending_field_reads = false;
  // Shallow field copies keep allocations alive until the queued reads finish,
  // including when a model or solver is destroyed after a failed integration.
  std::vector<Field2D> retained_2d;
  std::vector<Field3D> retained_3d;
  std::vector<FieldMap> fields;
};
#else
// Keep the shared interface usable in CPU-only and older SUNDIALS builds.
struct SundialsCudaPacking {};
#endif

#endif
