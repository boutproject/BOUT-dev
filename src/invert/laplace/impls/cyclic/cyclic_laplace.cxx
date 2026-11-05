/**************************************************************************
 * Copyright 2013 B.D.Dudson
 *
 * Contact: Ben Dudson, benjamin.dudson@york.ac.uk
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
 */

#include "bout/build_defines.hxx"

#if not BOUT_USE_METRIC_3D

#include "cyclic_laplace.hxx"

#include "../../common_transform.hxx"

#include "bout/array.hxx"
#include "bout/dcomplex.hxx"
#include <bout/assert.hxx>
#include <bout/bout_types.hxx>
#include <bout/boutcomm.hxx>
#include <bout/boutexception.hxx>
#include <bout/constants.hxx>
#include <bout/fft.hxx>
#include <bout/globals.hxx>
#include <bout/mesh.hxx>
#include <bout/output.hxx>
#include <bout/profiling.hxx>
#include <bout/sys/timer.hxx>
#include <bout/utils.hxx>

#include "cyclic_gpu.hxx"

#include <algorithm>
#include <array>
#include <cstring>
#include <tuple>
#include <vector>

namespace {
#if BOUT_HAS_CUDA || BOUT_HAS_HIP
static_assert(sizeof(dcomplex) == sizeof(gpu::Complex),
              "dcomplex and gpu::Complex must have the same memory layout");
static_assert(sizeof(dcomplex) == sizeof(gpu::fftDoubleComplex),
              "dcomplex and gpu::fftDoubleComplex must have the same memory layout");

void checkGpu(gpu::runtimeError_t status, const char* call) {
  if (status != gpu::runtimeSuccess) {
    throw BoutException("{} runtime error in LaplaceCyclic {}: {}", gpu::backend, call,
                        gpu::runtimeGetErrorString(status));
  }
}

#if BOUT_HAS_CUDA && !BOUT_HAS_HIP
// Resolve mapped pinned storage without assuming identical host/device addresses.
template <class T>
T* fieldDeviceAlias(T* pointer) {
  cudaPointerAttributes attributes{};
  const auto status = cudaPointerGetAttributes(&attributes, pointer);
  if (status == cudaErrorInvalidValue) {
    // Ordinary host allocations in CUDA builds still use the staging path.
    cudaGetLastError();
    return nullptr;
  }
  checkGpu(status, "cudaPointerGetAttributes(field)");
  return static_cast<T*>(attributes.devicePointer);
}

__global__ void packFieldForFft(const BoutReal* rhs, const BoutReal* x0, BoutReal* packed,
                                int nx, int ny, int nz, int local_nx, int local_ny,
                                int local_nz, int xs, int ys, int zs, int inbndry,
                                int outbndry, bool inner_set, bool outer_set) {
  const int id = blockIdx.x * blockDim.x + threadIdx.x;
  if (id >= nx * ny * nz) {
    return;
  }
  const int ix = xs + id / (ny * nz);
  const int iy = ys + (id / nz) % ny;
  const int iz = zs + id % nz;
  const bool boundary =
      (ix < inbndry && inner_set) || (local_nx - ix - 1 < outbndry && outer_set);
  packed[id] = (boundary ? x0 : rhs)[(ix * local_ny + iy) * local_nz + iz];
}

__global__ void unpackFieldFromFft(const BoutReal* packed, BoutReal* field, int nx,
                                   int ny, int nz, int local_ny, int local_nz, int xs,
                                   int ys, int zs) {
  const int id = blockIdx.x * blockDim.x + threadIdx.x;
  if (id >= nx * ny * nz) {
    return;
  }
  const int ix = xs + id / (ny * nz);
  const int iy = ys + (id / nz) % ny;
  const int iz = zs + id % nz;
  field[(ix * local_ny + iy) * local_nz + iz] = packed[id];
}
#endif

void checkSparse(gpu::sparseStatus_t status, const char* call) {
  if (status != gpu::SPARSE_STATUS_SUCCESS) {
    throw BoutException("{} sparse error in LaplaceCyclic {}: status {}", gpu::backend,
                        call, static_cast<int>(status));
  }
}

void checkFft(gpu::fftResult status, const char* call) {
  if (status != gpu::FFT_SUCCESS) {
    throw BoutException("{} FFT error in LaplaceCyclic {}: status {}", gpu::backend, call,
                        static_cast<int>(status));
  }
}

template <class T>
class DeviceBuffer {
public:
  DeviceBuffer() = default;
  explicit DeviceBuffer(std::size_t count) { ensure(count); }
  ~DeviceBuffer() {
    if (data != nullptr) {
#if BOUT_HAS_UMPIRE
      umpire::ResourceManager::getInstance().deallocate(data);
#else
      gpu::runtimeFree(data);
#endif
    }
  }

  DeviceBuffer(const DeviceBuffer&) = delete;
  DeviceBuffer& operator=(const DeviceBuffer&) = delete;

  T* get() { return data; }
  const T* get() const { return data; }

  void ensure(std::size_t count) {
    if (count <= capacity) {
      return;
    }
    if (data != nullptr) {
#if BOUT_HAS_UMPIRE
      umpire::ResourceManager::getInstance().deallocate(data);
#else
      gpu::runtimeFree(data);
#endif
      data = nullptr;
    }
    allocate(count);
    capacity = count;
  }

private:
  void allocate(std::size_t count) {
    if (count == 0) {
      return;
    }
#if BOUT_HAS_UMPIRE
    auto allocator =
        umpire::ResourceManager::getInstance().getAllocator(umpire::resource::Device);
    data = static_cast<T*>(allocator.allocate(count * sizeof(T)));
#else
    checkGpu(gpu::runtimeMalloc(&data, count * sizeof(T)), "gpu::runtimeMalloc");
#endif
  }

  T* data{nullptr};
  std::size_t capacity{0};
};

class SparseHandle {
public:
  SparseHandle() { checkSparse(gpu::sparseCreate(&handle), "gpu::sparseCreate"); }
  ~SparseHandle() {
    if (handle != nullptr) {
      gpu::sparseDestroy(handle);
    }
  }

  SparseHandle(const SparseHandle&) = delete;
  SparseHandle& operator=(const SparseHandle&) = delete;

  operator gpu::sparseHandle_t() const { return handle; }

private:
  gpu::sparseHandle_t handle{nullptr};
};

class FftPlan {
public:
  FftPlan() = default;
  ~FftPlan() {
    if (plan != 0) {
      gpu::fftDestroy(plan);
    }
  }

  FftPlan(const FftPlan&) = delete;
  FftPlan& operator=(const FftPlan&) = delete;

  gpu::fftHandle get() const { return plan; }

  void ensure(int new_nz, int new_batch, gpu::fftType type) {
    if (plan != 0 && new_nz == nz && new_batch == batch && type == plan_type) {
      return;
    }
    if (plan != 0) {
      gpu::fftDestroy(plan);
      plan = 0;
    }

    int n[] = {new_nz};
    const int real_dist = new_nz;
    const int complex_dist = (new_nz / 2) + 1;
    if (type == gpu::FFT_D2Z) {
      checkFft(gpu::fftPlanMany(&plan, 1, n, nullptr, 1, real_dist, nullptr, 1,
                                complex_dist, type, new_batch),
               "gpu::fftPlanMany D2Z");
    } else {
      checkFft(gpu::fftPlanMany(&plan, 1, n, nullptr, 1, complex_dist, nullptr, 1,
                                real_dist, type, new_batch),
               "gpu::fftPlanMany Z2D");
    }
    nz = new_nz;
    batch = new_batch;
    plan_type = type;
  }

private:
  gpu::fftHandle plan{0};
  int nz{0};
  int batch{0};
  gpu::fftType plan_type{gpu::FFT_D2Z};
};

__device__ gpu::Complex cadd(gpu::Complex a, gpu::Complex b) {
  return gpu::makeComplex(gpu::real(a) + gpu::real(b), gpu::imag(a) + gpu::imag(b));
}

__device__ gpu::Complex csub(gpu::Complex a, gpu::Complex b) {
  return gpu::makeComplex(gpu::real(a) - gpu::real(b), gpu::imag(a) - gpu::imag(b));
}

__device__ gpu::Complex cmul(gpu::Complex a, gpu::Complex b) {
  return gpu::makeComplex(gpu::real(a) * gpu::real(b) - gpu::imag(a) * gpu::imag(b),
                          gpu::real(a) * gpu::imag(b) + gpu::imag(a) * gpu::real(b));
}

__device__ gpu::Complex cdiv(gpu::Complex a, gpu::Complex b) {
  const double denom = gpu::real(b) * gpu::real(b) + gpu::imag(b) * gpu::imag(b);
  return gpu::makeComplex(
      (gpu::real(a) * gpu::real(b) + gpu::imag(a) * gpu::imag(b)) / denom,
      (gpu::imag(a) * gpu::real(b) - gpu::real(a) * gpu::imag(b)) / denom);
}

__device__ gpu::Complex cneg(gpu::Complex a) {
  return gpu::makeComplex(-gpu::real(a), -gpu::imag(a));
}

__global__ void prepareTridiagonalBatch(const gpu::Complex* a, const gpu::Complex* b,
                                        const gpu::Complex* c, const gpu::Complex* rhs,
                                        gpu::Complex* dl, gpu::Complex* d,
                                        gpu::Complex* du, gpu::Complex* x, int nsys,
                                        int nx) {
  const int id = blockIdx.x * blockDim.x + threadIdx.x;
  const int total = nsys * nx;
  if (id >= total) {
    return;
  }

  const int i = id % nx;
  dl[id] = (i == 0) ? gpu::makeComplex(0.0, 0.0) : a[id];
  d[id] = b[id];
  du[id] = (i == nx - 1) ? gpu::makeComplex(0.0, 0.0) : c[id];
  x[id] = rhs[id];
}

__global__ void prepareCyclicBatch(const gpu::Complex* a, const gpu::Complex* b,
                                   const gpu::Complex* c, const gpu::Complex* rhs,
                                   gpu::Complex* dl, gpu::Complex* d, gpu::Complex* du,
                                   gpu::Complex* x, int nsys, int nx) {
  const int id = blockIdx.x * blockDim.x + threadIdx.x;
  const int total = 2 * nsys * nx;
  if (id >= total) {
    return;
  }

  const int batch = id / nx;
  const int system = batch % nsys;
  const int i = id % nx;
  const int src = system * nx + i;
  const bool correction_rhs = batch >= nsys;

  const gpu::Complex alpha = a[system * nx];
  const gpu::Complex beta = c[system * nx + nx - 1];
  const gpu::Complex gamma = cneg(b[system * nx]);

  dl[id] = (i == 0) ? gpu::makeComplex(0.0, 0.0) : a[src];
  du[id] = (i == nx - 1) ? gpu::makeComplex(0.0, 0.0) : c[src];

  gpu::Complex diag = b[src];
  if (i == 0) {
    diag = csub(diag, gamma);
  } else if (i == nx - 1) {
    diag = csub(diag, cdiv(cmul(alpha, beta), gamma));
  }
  d[id] = diag;

  if (correction_rhs) {
    if (i == 0) {
      x[id] = gamma;
    } else if (i == nx - 1) {
      x[id] = beta;
    } else {
      x[id] = gpu::makeComplex(0.0, 0.0);
    }
  } else {
    x[id] = rhs[src];
  }
}

__global__ void finishCyclicBatch(const gpu::Complex* a, const gpu::Complex* b,
                                  const gpu::Complex* c, gpu::Complex* x, int nsys,
                                  int nx) {
  const int id = blockIdx.x * blockDim.x + threadIdx.x;
  const int total = nsys * nx;
  if (id >= total) {
    return;
  }

  const int system = id / nx;
  const int i = id % nx;
  const int correction = (system + nsys) * nx + i;

  const gpu::Complex alpha = a[system * nx];
  const gpu::Complex gamma = cneg(b[system * nx]);

  const gpu::Complex x0 = x[system * nx];
  const gpu::Complex xn = x[system * nx + nx - 1];
  const gpu::Complex z0 = x[(system + nsys) * nx];
  const gpu::Complex zn = x[(system + nsys) * nx + nx - 1];

  const gpu::Complex numerator = cadd(x0, cdiv(cmul(alpha, xn), gamma));
  const gpu::Complex denominator =
      cadd(gpu::makeComplex(1.0, 0.0), cadd(z0, cdiv(cmul(alpha, zn), gamma)));
  const gpu::Complex factor = cdiv(numerator, denominator);

  x[id] = csub(x[id], cmul(factor, x[correction]));
}

__global__ void fftOutputToSparseRhs(const gpu::fftDoubleComplex* fft_out,
                                     gpu::Complex* rhs, int nx, int ny, int nz, int nmode,
                                     int nmodes) {
  const int id = blockIdx.x * blockDim.x + threadIdx.x;
  const int total = nx * ny * nmode;
  if (id >= total) {
    return;
  }

  const int ix = id % nx;
  const int tmp = id / nx;
  const int kz = tmp % nmode;
  const int iy = tmp / nmode;
  const int fft_index = (ix * ny + iy) * nmodes + kz;
  const int rhs_index = (iy * nmode + kz) * nx + ix;
  const double scale = 1.0 / static_cast<double>(nz);

  rhs[rhs_index] = gpu::makeComplex(gpu::real(fft_out[fft_index]) * scale,
                                    gpu::imag(fft_out[fft_index]) * scale);
}

__global__ void sparseXToIfftInput(const gpu::Complex* x, gpu::fftDoubleComplex* ifft_in,
                                   int nx, int ny, int nmode, int nmodes, bool zero_dc) {
  const int id = blockIdx.x * blockDim.x + threadIdx.x;
  const int total = nx * ny * nmodes;
  if (id >= total) {
    return;
  }

  const int kz = id % nmodes;
  const int tmp = id / nmodes;
  const int iy = tmp % ny;
  const int ix = tmp / ny;

  gpu::fftDoubleComplex value = gpu::makeComplex(0.0, 0.0);
  if (kz < nmode && !(zero_dc && kz == 0)) {
    const int x_index = (iy * nmode + kz) * nx + ix;
    value = x[x_index];
  }
  ifft_in[id] = value;
}

__global__ void subtractPeriodicXAverage(gpu::Complex* x, int nx, int ny, int nmode) {
  extern __shared__ double partial[];
  const int iy = blockIdx.x;
  const int thread = threadIdx.x;

  double sum = 0.0;
  for (int ix = thread; ix < nx; ix += blockDim.x) {
    sum += gpu::real(x[(iy * nmode) * nx + ix]);
  }
  partial[thread] = sum;
  __syncthreads();

  for (int stride = blockDim.x / 2; stride > 0; stride /= 2) {
    if (thread < stride) {
      partial[thread] += partial[thread + stride];
    }
    __syncthreads();
  }

  const double avg = partial[0] / static_cast<double>(nx);
  for (int ix = thread; ix < nx; ix += blockDim.x) {
    gpu::Complex& value = x[(iy * nmode) * nx + ix];
    value = gpu::makeComplex(gpu::real(value) - avg, gpu::imag(value));
  }
}

__global__ void buildTridagMatrices(
    gpu::Complex* a, gpu::Complex* b, gpu::Complex* c, gpu::Complex* rhs,
    const double* acoef, const double* c1coef, const double* c2coef, const double* dcoef,
    const double* g11, const double* g33, const double* g13, const double* G1,
    const double* G3, const double* dx, const double* int_shift_torsion, int nx, int ny,
    int nmode, int local_nx, int local_ny, int xs, int ys, double zlength, bool all_terms,
    bool nonuniform, bool inc_int_shear, bool pin_zero_mode, bool periodic, int inbndry,
    int outbndry, bool inner_rhs, bool outer_rhs) {
  const int id = blockIdx.x * blockDim.x + threadIdx.x;
  const int total = nx * ny * nmode;
  if (id >= total) {
    return;
  }

  const int ix_local = id % nx;
  const int tmp = id / nx;
  const int kz = tmp % nmode;
  const int iy_local = tmp / nmode;
  const int jx = xs + ix_local;
  const int jy = ys + iy_local;
  const int field_index = jx * local_ny + jy;
  const double kwave = kz * 2.0 * PI / zlength;

  // Match tridagMatrix's half-cell Dirichlet rows, including all boundary
  // guard cells. No metric/coefficient evaluation is needed for these rows.
  if (!periodic && (ix_local < inbndry || ix_local >= nx - outbndry)) {
    const bool inner = ix_local < inbndry;
    a[id] = gpu::makeComplex(inner ? 0.0 : 0.5, 0.0);
    b[id] = gpu::makeComplex(0.5, 0.0);
    c[id] = gpu::makeComplex(inner ? 0.5 : 0.0, 0.0);
    if (!(inner ? inner_rhs : outer_rhs)) {
      rhs[id] = gpu::makeComplex(0.0, 0.0);
    }
    return;
  }

  if (pin_zero_mode && kz == 0 && ix_local == 0) {
    a[id] = gpu::makeComplex(0.0, 0.0);
    b[id] = gpu::makeComplex(1.0, 0.0);
    c[id] = gpu::makeComplex(0.0, 0.0);
    rhs[id] = gpu::makeComplex(0.0, 0.0);
    return;
  }

  double coef1 = g11[field_index];
  double coef2 = g33[field_index];
  double coef3 = 2.0 * g13[field_index];
  double coef4 = all_terms ? G1[field_index] : 0.0;
  double coef5 = all_terms ? G3[field_index] : 0.0;

  const double d = dcoef[field_index];
  coef1 *= d;
  coef2 *= d;
  coef3 *= d;
  coef4 *= d;
  coef5 *= d;

  if (nonuniform && jx != 0 && jx != local_nx - 1) {
    const double dx_center = dx[field_index];
    const double dx_plus = dx[(jx + 1) * local_ny + jy];
    const double dx_minus = dx[(jx - 1) * local_ny + jy];
    coef4 -= 0.5 * ((dx_plus - dx_minus) / (dx_center * dx_center)) * coef1;
  }

  if (jx > 0 && jx < local_nx - 1) {
    const double dc2dx_over_c1 =
        (c2coef[(jx + 1) * local_ny + jy] - c2coef[(jx - 1) * local_ny + jy])
        / (2.0 * dx[field_index] * c1coef[field_index]);
    coef4 += g11[field_index] * dc2dx_over_c1;
    coef5 += g13[field_index] * dc2dx_over_c1;
  }

  if (inc_int_shear) {
    const double shift = int_shift_torsion[field_index];
    coef2 += g11[field_index] * shift * shift;
    coef3 = 0.0;
  }

  const double dx_center = dx[field_index];
  coef1 /= dx_center * dx_center;
  coef3 /= 2.0 * dx_center;
  coef4 /= 2.0 * dx_center;

  b[id] = gpu::makeComplex(-2.0 * coef1 - kwave * kwave * coef2 + acoef[field_index],
                           kwave * coef5);
  a[id] = gpu::makeComplex(coef1 - coef4, -kwave * coef3);
  c[id] = gpu::makeComplex(coef1 + coef4, kwave * coef3);
}

} // namespace

class LaplaceCyclicGpuScratch {
public:
  SparseHandle handle;
  Array<BoutReal> real_host;
  DeviceBuffer<gpu::Complex> a;
  DeviceBuffer<gpu::Complex> b;
  DeviceBuffer<gpu::Complex> c;
  DeviceBuffer<gpu::Complex> rhs;
  DeviceBuffer<gpu::Complex> dl;
  DeviceBuffer<gpu::Complex> d;
  DeviceBuffer<gpu::Complex> du;
  DeviceBuffer<gpu::Complex> x;
  DeviceBuffer<gpu::fftDoubleReal> real_device;
  DeviceBuffer<gpu::fftDoubleComplex> spectral_device;
  DeviceBuffer<double> acoef;
  DeviceBuffer<double> c1coef;
  DeviceBuffer<double> c2coef;
  DeviceBuffer<double> dcoef;
  DeviceBuffer<double> g11;
  DeviceBuffer<double> g33;
  DeviceBuffer<double> g13;
  DeviceBuffer<double> G1;
  DeviceBuffer<double> G3;
  DeviceBuffer<double> dx;
  DeviceBuffer<double> int_shift_torsion;
  DeviceBuffer<char> buffer;
  FftPlan forward_plan;
  FftPlan backward_plan;
  int real_host_size{0};
  std::size_t buffer_size{0};
  // Content snapshots detect in-place edits as well as coefficient replacement.
  std::array<std::vector<double>, 11> cached_field_values;

  void ensureFft(int nz, int batch) {
    const auto real_size = batch * nz;
    const auto spectral_size = batch * ((nz / 2) + 1);
    real_device.ensure(static_cast<std::size_t>(real_size));
    spectral_device.ensure(static_cast<std::size_t>(spectral_size));
    forward_plan.ensure(nz, batch, gpu::FFT_D2Z);
    backward_plan.ensure(nz, batch, gpu::FFT_Z2D);
  }
};

namespace {

const double* field2DData(const Field2D& field) {
  const auto view = static_cast<Field2D::View>(field);
  return view.data;
}

// Keep repeatedly accessed coefficients in device memory. Compare host contents
// rather than addresses: fields and metrics can be modified in place.
const double* field2DDeviceData(DeviceBuffer<double>& destination, const Field2D& field,
                                std::size_t size, const char* name,
                                std::vector<double>& cached_values) {
  const auto* host = field2DData(field);
  const auto bytes = size * sizeof(double);
  if (cached_values.size() != size
      || std::memcmp(cached_values.data(), host, bytes) != 0) {
    destination.ensure(size);
    checkGpu(gpu::runtimeMemcpy(destination.get(), host, bytes,
                                gpu::runtimeMemcpyHostToDevice),
             name);
    // Update only after a successful upload so failures cannot validate stale data.
    cached_values.assign(host, host + size);
  }
  return destination.get();
}

void solveWithGpuSparseDevice(int nsys, int nx, bool periodic,
                              LaplaceCyclicGpuScratch& scratch) {
  const bout::profiling::ProfileRange solve_range("LaplaceCyclic::gpu_sparse");
  const int solve_batches = periodic ? 2 * nsys : nsys;
  const std::size_t matrix_values = static_cast<std::size_t>(nsys) * nx;
  const std::size_t solve_values = static_cast<std::size_t>(solve_batches) * nx;

  {
    const bout::profiling::ProfileRange range("LaplaceCyclic::sparse_buffers");
    scratch.a.ensure(matrix_values);
    scratch.b.ensure(matrix_values);
    scratch.c.ensure(matrix_values);
    scratch.dl.ensure(solve_values);
    scratch.d.ensure(solve_values);
    scratch.du.ensure(solve_values);
    scratch.x.ensure(solve_values);
  }

  {
    const bout::profiling::ProfileRange range("LaplaceCyclic::sparse_prepare");
    const int block_size = 256;
    const int prepare_blocks =
        (static_cast<int>(solve_values) + block_size - 1) / block_size;
    if (periodic) {
      prepareCyclicBatch<<<prepare_blocks, block_size>>>(
          scratch.a.get(), scratch.b.get(), scratch.c.get(), scratch.rhs.get(),
          scratch.dl.get(), scratch.d.get(), scratch.du.get(), scratch.x.get(), nsys, nx);
    } else {
      prepareTridiagonalBatch<<<prepare_blocks, block_size>>>(
          scratch.a.get(), scratch.b.get(), scratch.c.get(), scratch.rhs.get(),
          scratch.dl.get(), scratch.d.get(), scratch.du.get(), scratch.x.get(), nsys, nx);
    }
    checkGpu(gpu::runtimeGetLastError(), "prepare sparse batch");
  }

  size_t buffer_size = 0;
  {
    const bout::profiling::ProfileRange range("LaplaceCyclic::sparse_workspace");
    checkSparse(gpu::sparseZgtsv2StridedBatch_bufferSizeExt(
                    scratch.handle, nx, scratch.dl.get(), scratch.d.get(),
                    scratch.du.get(), scratch.x.get(), solve_batches, nx, &buffer_size),
                "gpu::sparseZgtsv2StridedBatch_bufferSizeExt");

    if (buffer_size > scratch.buffer_size) {
      scratch.buffer.ensure(buffer_size);
      scratch.buffer_size = buffer_size;
    }
  }
  {
    const bout::profiling::ProfileRange range("LaplaceCyclic::sparse_execute");
    checkSparse(gpu::sparseZgtsv2StridedBatch(scratch.handle, nx, scratch.dl.get(),
                                              scratch.d.get(), scratch.du.get(),
                                              scratch.x.get(), solve_batches, nx,
                                              scratch.buffer.get()),
                "gpu::sparseZgtsv2StridedBatch");
  }

  if (periodic) {
    const bout::profiling::ProfileRange range("LaplaceCyclic::sparse_finish_cyclic");
    const int block_size = 256;
    const int finish_blocks =
        (static_cast<int>(matrix_values) + block_size - 1) / block_size;
    finishCyclicBatch<<<finish_blocks, block_size>>>(
        scratch.a.get(), scratch.b.get(), scratch.c.get(), scratch.x.get(), nsys, nx);
    checkGpu(gpu::runtimeGetLastError(), "finish cyclic sparse batch");
  }
}

void solveWithGpuSparse(const Matrix<dcomplex>& a, const Matrix<dcomplex>& b,
                        const Matrix<dcomplex>& c, const Matrix<dcomplex>& rhs,
                        Matrix<dcomplex>& x, bool periodic,
                        LaplaceCyclicGpuScratch& scratch) {
  const bout::profiling::ProfileRange solve_range("LaplaceCyclic::gpu_sparse");
  const int nsys = std::get<0>(a.shape());
  const int nx = std::get<1>(a.shape());
  ASSERT2(std::get<0>(b.shape()) == nsys);
  ASSERT2(std::get<1>(b.shape()) == nx);
  ASSERT2(std::get<0>(c.shape()) == nsys);
  ASSERT2(std::get<1>(c.shape()) == nx);
  ASSERT2(std::get<0>(rhs.shape()) == nsys);
  ASSERT2(std::get<1>(rhs.shape()) == nx);
  ASSERT2(std::get<0>(x.shape()) == nsys);
  ASSERT2(std::get<1>(x.shape()) == nx);

  const int solve_batches = periodic ? 2 * nsys : nsys;
  const std::size_t matrix_values = static_cast<std::size_t>(nsys) * nx;
  const std::size_t solve_values = static_cast<std::size_t>(solve_batches) * nx;

  {
    const bout::profiling::ProfileRange range("LaplaceCyclic::sparse_buffers");
    scratch.a.ensure(matrix_values);
    scratch.b.ensure(matrix_values);
    scratch.c.ensure(matrix_values);
    scratch.rhs.ensure(matrix_values);
    scratch.dl.ensure(solve_values);
    scratch.d.ensure(solve_values);
    scratch.du.ensure(solve_values);
    scratch.x.ensure(solve_values);
  }

  {
    const bout::profiling::ProfileRange range("LaplaceCyclic::sparse_upload");
    checkGpu(gpu::runtimeMemcpy(
                 scratch.a.get(), reinterpret_cast<const gpu::Complex*>(a.begin()),
                 matrix_values * sizeof(gpu::Complex), gpu::runtimeMemcpyHostToDevice),
             "copy a to device");
    checkGpu(gpu::runtimeMemcpy(
                 scratch.b.get(), reinterpret_cast<const gpu::Complex*>(b.begin()),
                 matrix_values * sizeof(gpu::Complex), gpu::runtimeMemcpyHostToDevice),
             "copy b to device");
    checkGpu(gpu::runtimeMemcpy(
                 scratch.c.get(), reinterpret_cast<const gpu::Complex*>(c.begin()),
                 matrix_values * sizeof(gpu::Complex), gpu::runtimeMemcpyHostToDevice),
             "copy c to device");
    checkGpu(gpu::runtimeMemcpy(
                 scratch.rhs.get(), reinterpret_cast<const gpu::Complex*>(rhs.begin()),
                 matrix_values * sizeof(gpu::Complex), gpu::runtimeMemcpyHostToDevice),
             "copy rhs to device");
  }

  {
    const bout::profiling::ProfileRange range("LaplaceCyclic::sparse_prepare");
    const int block_size = 256;
    const int prepare_blocks =
        (static_cast<int>(solve_values) + block_size - 1) / block_size;
    if (periodic) {
      prepareCyclicBatch<<<prepare_blocks, block_size>>>(
          scratch.a.get(), scratch.b.get(), scratch.c.get(), scratch.rhs.get(),
          scratch.dl.get(), scratch.d.get(), scratch.du.get(), scratch.x.get(), nsys, nx);
    } else {
      prepareTridiagonalBatch<<<prepare_blocks, block_size>>>(
          scratch.a.get(), scratch.b.get(), scratch.c.get(), scratch.rhs.get(),
          scratch.dl.get(), scratch.d.get(), scratch.du.get(), scratch.x.get(), nsys, nx);
    }
    checkGpu(gpu::runtimeGetLastError(), "prepare sparse batch");
  }

  size_t buffer_size = 0;
  {
    const bout::profiling::ProfileRange range("LaplaceCyclic::sparse_workspace");
    checkSparse(gpu::sparseZgtsv2StridedBatch_bufferSizeExt(
                    scratch.handle, nx, scratch.dl.get(), scratch.d.get(),
                    scratch.du.get(), scratch.x.get(), solve_batches, nx, &buffer_size),
                "gpu::sparseZgtsv2StridedBatch_bufferSizeExt");

    if (buffer_size > scratch.buffer_size) {
      scratch.buffer.ensure(buffer_size);
      scratch.buffer_size = buffer_size;
    }
  }
  {
    const bout::profiling::ProfileRange range("LaplaceCyclic::sparse_execute");
    checkSparse(gpu::sparseZgtsv2StridedBatch(scratch.handle, nx, scratch.dl.get(),
                                              scratch.d.get(), scratch.du.get(),
                                              scratch.x.get(), solve_batches, nx,
                                              scratch.buffer.get()),
                "gpu::sparseZgtsv2StridedBatch");
  }

  if (periodic) {
    const bout::profiling::ProfileRange range("LaplaceCyclic::sparse_finish_cyclic");
    const int block_size = 256;
    const int finish_blocks =
        (static_cast<int>(matrix_values) + block_size - 1) / block_size;
    finishCyclicBatch<<<finish_blocks, block_size>>>(
        scratch.a.get(), scratch.b.get(), scratch.c.get(), scratch.x.get(), nsys, nx);
    checkGpu(gpu::runtimeGetLastError(), "finish cyclic sparse batch");
  }

  {
    const bout::profiling::ProfileRange range("LaplaceCyclic::sparse_download");
    checkGpu(gpu::runtimeMemcpy(reinterpret_cast<gpu::Complex*>(x.begin()),
                                scratch.x.get(), matrix_values * sizeof(gpu::Complex),
                                gpu::runtimeMemcpyDeviceToHost),
             "copy solution to host");
  }
}

bool canUseGpuSolve(const Mesh& mesh) { return mesh.getNXPE() == 1; }

bool canUseResidentFftSolve(const Mesh& mesh, int inner_flags, int outer_flags) {
  // The device boundary rows below implement second-order Dirichlet conditions.
  // Other boundary operators retain the general CPU construction path.
  constexpr int supported = INVERT_BNDRY_ONE | INVERT_SET | INVERT_RHS;
  return mesh.getNXPE() == 1
         && (mesh.periodicX
             || ((inner_flags & ~supported) == 0 && (outer_flags & ~supported) == 0));
}
#endif
} // namespace

LaplaceCyclic::LaplaceCyclic(Options* opt, const CELL_LOC loc, Mesh* mesh_in,
                             Solver* UNUSED(solver))
    : Laplacian(opt, loc, mesh_in), Acoef(0.0), C1coef(1.0), C2coef(1.0), Dcoef(1.0) {

  bout::fft::assertZSerial(*localmesh, "`cyclic` inversion");

  Acoef.setLocation(location);
  C1coef.setLocation(location);
  C2coef.setLocation(location);
  Dcoef.setLocation(location);

  // Get options

  dst = (*opt)["dst"]
            .doc("Use Discrete Sine Transform in Z to enforce Dirichlet boundaries in Z")
            .withDefault<bool>(false);
#if BOUT_HAS_HIP
  use_gpu_sparse =
      (*opt)["use_hipsparse"]
#else
  use_gpu_sparse =
      (*opt)["use_cusparse"]
#endif
          .doc("Use GPU sparse batched tridiagonal solves for the cyclic "
               "Laplacian when CUDA or HIP is enabled and X is not distributed")
          .withDefault<bool>(true);
  compare_device_tridag =
      (*opt)["compare_device_tridag"]
          .doc("Compare GPU-built cyclic Laplacian tridiagonal matrices against the "
               "CPU tridagMatrix implementation once, when using the resident GPU path")
          .withDefault<bool>(true);

  if (dst) {
    nmode = localmesh->LocalNz - 2;
  } else {
    // Number of Z modes. maxmode set in invert_laplace.cxx from options
    nmode = maxmode + 1;
  }

  // Note nmode == nsys of cyclic_reduction

  // Allocate arrays

  xs = localmesh->xstart; // Starting X index
  if (localmesh->firstX()
      && !localmesh->periodicX) { // Only want to include guard cells at boundaries
                                  // (unless periodic in x)
    xs = 0;
  }
  xe = localmesh->xend; // Last X index
  if (localmesh->lastX()
      && !localmesh->periodicX) { // Only want to include guard cells at boundaries
                                  // (unless periodic in x)
    xe = localmesh->LocalNx - 1;
  }
  int n = xe - xs + 1; // Number of X points on this processor,
                       // including boundaries but not guard cells

  a.reallocate(nmode, n);
  b.reallocate(nmode, n);
  c.reallocate(nmode, n);
  xcmplx.reallocate(nmode, n);
  bcmplx.reallocate(nmode, n);

  // Create a cyclic reduction object, operating on dcomplex values
  cr = new CyclicReduce<dcomplex>(localmesh->getXcomm(), n);
  cr->setPeriodic(localmesh->periodicX);
}

LaplaceCyclic::~LaplaceCyclic() {
  // Delete tridiagonal solver
  delete cr;
}

FieldPerp LaplaceCyclic::solve(const FieldPerp& rhs, const FieldPerp& x0) {
  ASSERT1(localmesh == rhs.getMesh() && localmesh == x0.getMesh());
  ASSERT1(rhs.getLocation() == location);
  ASSERT1(x0.getLocation() == location);

  // Get the width of the boundary

  // If the flags to assign that only one guard cell should be used is set
  int inbndry = localmesh->xstart;
  int outbndry = localmesh->xstart;
  if (isGlobalFlagSet(INVERT_BOTH_BNDRY_ONE) || (localmesh->xstart < 2)) {
    inbndry = outbndry = 1;
  }
  if (isInnerBoundaryFlagSet(INVERT_BNDRY_ONE)) {
    inbndry = 1;
  }
  if (isOuterBoundaryFlagSet(INVERT_BNDRY_ONE)) {
    outbndry = 1;
  }

  if (dst) {
    const DSTTransform transform(
        *localmesh, nmode, xs, xe, 0, 0, localmesh->zstart, localmesh->zend, inbndry,
        outbndry, isInnerBoundaryFlagSetOnFirstX(INVERT_SET),
        isOuterBoundaryFlagSetOnLastX(INVERT_SET), isGlobalFlagSet(INVERT_ZERO_DC));

    auto matrices = transform.forward(*this, rhs, x0, Acoef, C1coef, C2coef, Dcoef);

    // Solve tridiagonal systems
#if BOUT_HAS_CUDA || BOUT_HAS_HIP
    if (use_gpu_sparse && canUseGpuSolve(*localmesh)) {
      if (!gpu_scratch) {
        gpu_scratch = std::make_unique<LaplaceCyclicGpuScratch>();
      }
      solveWithGpuSparse(matrices.a, matrices.b, matrices.c, matrices.bcmplx, xcmplx,
                         localmesh->periodicX, *gpu_scratch);
    } else
#endif
    {
      cr->setCoefs(a, b, c);
      cr->solve(bcmplx, xcmplx);
    }

    return transform.backward(rhs, xcmplx);
  }

  const FFTTransform transform(
      *localmesh, nmode, xs, xe, 0, 0, localmesh->zstart, localmesh->zend, inbndry,
      outbndry, isInnerBoundaryFlagSetOnFirstX(INVERT_SET),
      isOuterBoundaryFlagSetOnLastX(INVERT_SET), isGlobalFlagSet(INVERT_ZERO_DC));

  auto matrices = transform.forward(*this, rhs, x0, Acoef, C1coef, C2coef, Dcoef);

  // Solve tridiagonal systems
#if BOUT_HAS_CUDA || BOUT_HAS_HIP
  if (use_gpu_sparse && canUseGpuSolve(*localmesh)) {
    if (!gpu_scratch) {
      gpu_scratch = std::make_unique<LaplaceCyclicGpuScratch>();
    }
    solveWithGpuSparse(matrices.a, matrices.b, matrices.c, matrices.bcmplx, xcmplx,
                       localmesh->periodicX, *gpu_scratch);
  } else
#endif
  {
    cr->setCoefs(matrices.a, matrices.b, matrices.c);
    cr->solve(matrices.bcmplx, xcmplx);
  }

  if (localmesh->periodicX) {
    // Subtract X average of kz=0 mode
    std::array<BoutReal, 2> local = {
        0.0,                               // index 0 = sum of coefficients
        static_cast<BoutReal>(xe - xs + 1) // number of grid cells
    };
    for (int ix = xs; ix <= xe; ix++) {
      local[0] += xcmplx(0, ix - xs).real();
    }
    std::array<BoutReal, 2> global{};
    MPI_Allreduce(local.data(), global.data(), 2, MPI_DOUBLE, MPI_SUM,
                  localmesh->getXcomm());
    const BoutReal avg = global[0] / global[1];
    for (int ix = xs; ix <= xe; ix++) {
      xcmplx(0, ix - xs) -= avg;
    }
  }

  return transform.backward(rhs, xcmplx);
}

Field3D LaplaceCyclic::solve(const Field3D& rhs, const Field3D& x0) {

  const bout::profiling::ProfileRange solve_range("LaplaceCyclic::solve3D");

  ASSERT1(rhs.getLocation() == location);
  ASSERT1(x0.getLocation() == location);
  ASSERT1(localmesh == rhs.getMesh() && localmesh == x0.getMesh());

  Timer timer("invert");

  Field3D x{emptyFrom(rhs)}; // Result

  // Get the width of the boundary

  // If the flags to assign that only one guard cell should be used is set
  int inbndry = localmesh->xstart;
  int outbndry = localmesh->xstart;
  if (isGlobalFlagSet(INVERT_BOTH_BNDRY_ONE) || (localmesh->xstart < 2)) {
    inbndry = outbndry = 1;
  }
  if (isInnerBoundaryFlagSet(INVERT_BNDRY_ONE)) {
    inbndry = 1;
  }
  if (isOuterBoundaryFlagSet(INVERT_BNDRY_ONE)) {
    outbndry = 1;
  }

  int nx = xe - xs + 1; // Number of X points on this processor

  // Get range of Y indices
  int ys = localmesh->ystart;
  int ye = localmesh->yend;

  if (localmesh->hasBndryLowerY()) {
    if (include_yguards) {
      ys = 0; // Mesh contains a lower boundary and we are solving in the guard cells
    }

    ys += extra_yguards_lower;
  }
  if (localmesh->hasBndryUpperY()) {
    if (include_yguards) {
      // Contains upper boundary and we are solving in the guard cells
      ye = localmesh->LocalNy - 1;
    }
    ye -= extra_yguards_upper;
  }

  const int ny = (ye - ys + 1); // Number of Y points
  const int nsys = nmode * ny;  // Number of systems of equations to solve

  // This is just to silence static analysis
  ASSERT0(ny > 0);

  auto xcmplx3D = [&]() {
    const bout::profiling::ProfileRange range("LaplaceCyclic::host_matrix_allocate");
    return Matrix<dcomplex>(nsys, nx);
  }();

  if (dst) {
    const DSTTransform transform(
        *localmesh, nmode, xs, xe, ys, ye, localmesh->zstart, localmesh->zend, inbndry,
        outbndry, isInnerBoundaryFlagSetOnFirstX(INVERT_SET),
        isOuterBoundaryFlagSetOnLastX(INVERT_SET), isGlobalFlagSet(INVERT_ZERO_DC));

    auto matrices = [&]() {
      const bout::profiling::ProfileRange range("LaplaceCyclic::transform_forward");
      return transform.forward(*this, rhs, x0, Acoef, C1coef, C2coef, Dcoef);
    }();

    // Solve tridiagonal systems
#if BOUT_HAS_CUDA || BOUT_HAS_HIP
    if (use_gpu_sparse && canUseGpuSolve(*localmesh)) {
      if (!gpu_scratch) {
        gpu_scratch = std::make_unique<LaplaceCyclicGpuScratch>();
      }
      solveWithGpuSparse(matrices.a, matrices.b, matrices.c, matrices.bcmplx, xcmplx3D,
                         localmesh->periodicX, *gpu_scratch);
    } else
#endif
    {
      {
        const bout::profiling::ProfileRange range("LaplaceCyclic::cpu_set_coefficients");
        cr->setCoefs(matrices.a, matrices.b, matrices.c);
      }
      {
        const bout::profiling::ProfileRange range("LaplaceCyclic::cpu_solve");
        cr->solve(matrices.bcmplx, xcmplx3D);
      }
    }

    {
      const bout::profiling::ProfileRange range("LaplaceCyclic::transform_backward");
      return transform.backward(rhs, xcmplx3D);
    }
  }
#if BOUT_HAS_CUDA || BOUT_HAS_HIP
  if (use_gpu_sparse
      && canUseResidentFftSolve(*localmesh, getInnerBoundaryFlags(),
                                getOuterBoundaryFlags())) {
    const bout::profiling::ProfileRange resident_range("LaplaceCyclic::resident_fft");
    if (!gpu_scratch) {
      gpu_scratch = std::make_unique<LaplaceCyclicGpuScratch>();
    }
    auto& scratch = *gpu_scratch;

    {
      const int nz = localmesh->zend - localmesh->zstart + 1;
      const int nxny = nx * ny;
      const int nmodes = (nz / 2) + 1;
      {
        const bout::profiling::ProfileRange range("LaplaceCyclic::fft_setup");
        scratch.ensureFft(nz, nxny);
      }

#if BOUT_HAS_CUDA && !BOUT_HAS_HIP
      const auto* rhs_device = fieldDeviceAlias(&rhs(0, 0, 0));
      const bool inner_set = isInnerBoundaryFlagSetOnFirstX(INVERT_SET);
      const bool outer_set = isOuterBoundaryFlagSetOnLastX(INVERT_SET);
      const auto* x0_device =
          (inner_set || outer_set) ? fieldDeviceAlias(&x0(0, 0, 0)) : rhs_device;
      if (rhs_device && x0_device) {
        const bout::profiling::ProfileRange range("LaplaceCyclic::pack_rhs_device");
        // Field producers may use other streams; preserve host-path ordering.
        checkGpu(cudaDeviceSynchronize(), "wait for FFT input fields");
        const int total = nxny * nz;
        packFieldForFft<<<(total + 255) / 256, 256>>>(
            rhs_device, x0_device, scratch.real_device.get(), nx, ny, nz,
            localmesh->LocalNx, localmesh->LocalNy, localmesh->LocalNz, xs, ys,
            localmesh->zstart, inbndry, outbndry, inner_set, outer_set);
        checkGpu(cudaGetLastError(), "pack FFT input from mapped fields");
      } else
#endif
      {
        if (nxny * nz > scratch.real_host_size) {
          scratch.real_host.reallocate(nxny * nz);
          scratch.real_host_size = nxny * nz;
        }
        {
          const bout::profiling::ProfileRange range("LaplaceCyclic::pack_rhs_host");
          for (int ind = 0; ind < nxny; ++ind) {
            const int ix = xs + (ind / ny);
            const int iy = ys + (ind % ny);
            const BoutReal* input =
                (((ix < inbndry) and isInnerBoundaryFlagSetOnFirstX(INVERT_SET))
                 || ((localmesh->LocalNx - ix - 1 < outbndry)
                     and isOuterBoundaryFlagSetOnLastX(INVERT_SET)))
                    ? &(x0(ix, iy, localmesh->zstart))
                    : &(rhs(ix, iy, localmesh->zstart));
            std::copy(input, input + nz, std::begin(scratch.real_host) + ind * nz);
          }
        }

        {
          const bout::profiling::ProfileRange range("LaplaceCyclic::fft_upload");
          checkGpu(
              gpu::runtimeMemcpy(scratch.real_device.get(), std::begin(scratch.real_host),
                                 static_cast<std::size_t>(nxny) * nz * sizeof(BoutReal),
                                 gpu::runtimeMemcpyHostToDevice),
              "copy resident rfft input to device");
        }
      } // Host staging fallback
      {
        const bout::profiling::ProfileRange range("LaplaceCyclic::fft_forward");
        checkFft(gpu::fftExecD2Z(scratch.forward_plan.get(), scratch.real_device.get(),
                                 scratch.spectral_device.get()),
                 "resident gpu::fftExecD2Z");
      }

      {
        const bout::profiling::ProfileRange range("LaplaceCyclic::fft_to_sparse");
        scratch.rhs.ensure(static_cast<std::size_t>(nsys) * nx);
        const int block_size = 256;
        const int total = nx * ny * nmode;
        const int blocks = (total + block_size - 1) / block_size;
        fftOutputToSparseRhs<<<blocks, block_size>>>(
            scratch.spectral_device.get(), scratch.rhs.get(), nx, ny, nz, nmode, nmodes);
        checkGpu(gpu::runtimeGetLastError(), "transpose rfft output to GPU sparse rhs");
      }
    }

    {
      const Coordinates* localcoords = localmesh->getCoordinates(location);
      std::array<const double*, 11> coefficient_data{};
      {
        const bout::profiling::ProfileRange range("LaplaceCyclic::coefficients_access");
        const std::size_t field_size =
            static_cast<std::size_t>(localmesh->LocalNx) * localmesh->LocalNy;
        coefficient_data[0] =
            field2DDeviceData(scratch.acoef, Acoef, field_size, "copy acoef to device",
                              scratch.cached_field_values[0]);
        coefficient_data[1] =
            field2DDeviceData(scratch.c1coef, C1coef, field_size, "copy c1coef to device",
                              scratch.cached_field_values[1]);
        coefficient_data[2] =
            field2DDeviceData(scratch.c2coef, C2coef, field_size, "copy c2coef to device",
                              scratch.cached_field_values[2]);
        coefficient_data[3] =
            field2DDeviceData(scratch.dcoef, Dcoef, field_size, "copy dcoef to device",
                              scratch.cached_field_values[3]);
        coefficient_data[4] =
            field2DDeviceData(scratch.g11, localcoords->g11(), field_size,
                              "copy g11 to device", scratch.cached_field_values[4]);
        coefficient_data[5] =
            field2DDeviceData(scratch.g33, localcoords->g33(), field_size,
                              "copy g33 to device", scratch.cached_field_values[5]);
        coefficient_data[6] =
            field2DDeviceData(scratch.g13, localcoords->g13(), field_size,
                              "copy g13 to device", scratch.cached_field_values[6]);
        coefficient_data[7] =
            field2DDeviceData(scratch.G1, localcoords->G1(), field_size,
                              "copy G1 to device", scratch.cached_field_values[7]);
        coefficient_data[8] =
            field2DDeviceData(scratch.G3, localcoords->G3(), field_size,
                              "copy G3 to device", scratch.cached_field_values[8]);
        coefficient_data[9] =
            field2DDeviceData(scratch.dx, localcoords->dx(), field_size,
                              "copy dx to device", scratch.cached_field_values[9]);
        coefficient_data[10] = field2DDeviceData(
            scratch.int_shift_torsion, localcoords->IntShiftTorsion(), field_size,
            "copy IntShiftTorsion to device", scratch.cached_field_values[10]);

        scratch.a.ensure(static_cast<std::size_t>(nsys) * nx);
        scratch.b.ensure(static_cast<std::size_t>(nsys) * nx);
        scratch.c.ensure(static_cast<std::size_t>(nsys) * nx);
      }

      const BoutReal zlength = getUniform(localcoords->zlength());
      {
        const bout::profiling::ProfileRange range("LaplaceCyclic::build_matrices");
        const int block_size = 256;
        const int total = nx * ny * nmode;
        const int blocks = (total + block_size - 1) / block_size;
        const bool pin_zero_mode = localmesh->periodicX && localmesh->firstX();
        buildTridagMatrices<<<blocks, block_size>>>(
            scratch.a.get(), scratch.b.get(), scratch.c.get(), scratch.rhs.get(),
            coefficient_data[0], coefficient_data[1], coefficient_data[2],
            coefficient_data[3], coefficient_data[4], coefficient_data[5],
            coefficient_data[6], coefficient_data[7], coefficient_data[8],
            coefficient_data[9], coefficient_data[10], nx, ny, nmode, localmesh->LocalNx,
            localmesh->LocalNy, xs, ys, zlength, all_terms, nonuniform,
            localmesh->IncIntShear, pin_zero_mode, localmesh->periodicX, inbndry,
            outbndry, isInnerBoundaryFlagSet(INVERT_RHS | INVERT_SET),
            isOuterBoundaryFlagSet(INVERT_RHS | INVERT_SET));
        checkGpu(gpu::runtimeGetLastError(), "build tridag matrices on device");
      }

      if (compare_device_tridag && !compared_device_tridag) {
        const bout::profiling::ProfileRange range(
            "LaplaceCyclic::compare_device_tridag_once");
        Matrix<dcomplex> cpu_a(nsys, nx);
        Matrix<dcomplex> cpu_b(nsys, nx);
        Matrix<dcomplex> cpu_c(nsys, nx);
        Matrix<dcomplex> cpu_rhs(nsys, nx);
        Matrix<dcomplex> gpu_a(nsys, nx);
        Matrix<dcomplex> gpu_b(nsys, nx);
        Matrix<dcomplex> gpu_c(nsys, nx);
        Matrix<dcomplex> gpu_rhs(nsys, nx);

        const auto matrix_bytes =
            static_cast<std::size_t>(nsys) * nx * sizeof(gpu::Complex);
        checkGpu(gpu::runtimeMemcpy(reinterpret_cast<gpu::Complex*>(gpu_a.begin()),
                                    scratch.a.get(), matrix_bytes,
                                    gpu::runtimeMemcpyDeviceToHost),
                 "copy GPU tridag a to host for comparison");
        checkGpu(gpu::runtimeMemcpy(reinterpret_cast<gpu::Complex*>(gpu_b.begin()),
                                    scratch.b.get(), matrix_bytes,
                                    gpu::runtimeMemcpyDeviceToHost),
                 "copy GPU tridag b to host for comparison");
        checkGpu(gpu::runtimeMemcpy(reinterpret_cast<gpu::Complex*>(gpu_c.begin()),
                                    scratch.c.get(), matrix_bytes,
                                    gpu::runtimeMemcpyDeviceToHost),
                 "copy GPU tridag c to host for comparison");
        checkGpu(gpu::runtimeMemcpy(reinterpret_cast<gpu::Complex*>(gpu_rhs.begin()),
                                    scratch.rhs.get(), matrix_bytes,
                                    gpu::runtimeMemcpyDeviceToHost),
                 "copy GPU tridag rhs to host for comparison");
        if (localmesh->periodicX) {
          std::copy(gpu_rhs.begin(), gpu_rhs.end(), cpu_rhs.begin());
        } else {
          const FFTTransform reference_transform(
              *localmesh, nmode, xs, xe, ys, ye, localmesh->zstart, localmesh->zend,
              inbndry, outbndry, isInnerBoundaryFlagSetOnFirstX(INVERT_SET),
              isOuterBoundaryFlagSetOnLastX(INVERT_SET), isGlobalFlagSet(INVERT_ZERO_DC));
          auto reference =
              reference_transform.forward(*this, rhs, x0, Acoef, C1coef, C2coef, Dcoef);
          std::copy(reference.bcmplx.begin(), reference.bcmplx.end(), cpu_rhs.begin());
        }

        for (int ind = 0; ind < nsys; ind++) {
          const int iy = ys + (ind / nmode);
          const int kz = ind % nmode;
          const BoutReal kwave = kz * 2.0 * PI / zlength;
          tridagMatrix(&cpu_a(ind, 0), &cpu_b(ind, 0), &cpu_c(ind, 0), &cpu_rhs(ind, 0),
                       iy, kz, kwave, &Acoef, &C1coef, &C2coef, &Dcoef, false);
        }

        BoutReal max_a = 0.0;
        BoutReal max_b = 0.0;
        BoutReal max_c = 0.0;
        BoutReal max_rhs = 0.0;
        int max_a_index = 0;
        int max_b_index = 0;
        int max_c_index = 0;
        int max_rhs_index = 0;
        const int values = nsys * nx;
        for (int ind = 0; ind < values; ind++) {
          const BoutReal diff_a = std::abs(cpu_a.begin()[ind] - gpu_a.begin()[ind]);
          const BoutReal diff_b = std::abs(cpu_b.begin()[ind] - gpu_b.begin()[ind]);
          const BoutReal diff_c = std::abs(cpu_c.begin()[ind] - gpu_c.begin()[ind]);
          const BoutReal diff_rhs = std::abs(cpu_rhs.begin()[ind] - gpu_rhs.begin()[ind]);
          if (diff_a > max_a) {
            max_a = diff_a;
            max_a_index = ind;
          }
          if (diff_b > max_b) {
            max_b = diff_b;
            max_b_index = ind;
          }
          if (diff_c > max_c) {
            max_c = diff_c;
            max_c_index = ind;
          }
          if (diff_rhs > max_rhs) {
            max_rhs = diff_rhs;
            max_rhs_index = ind;
          }
        }

        const int compare_nmode = nmode;
        const int compare_xs = xs;
        const int compare_ys = ys;
        auto describe_index = [nx, compare_xs, compare_ys,
                               compare_nmode](int flat_index) {
          const int ix_local = flat_index % nx;
          const int system = flat_index / nx;
          const int kz = system % compare_nmode;
          const int iy_local = system / compare_nmode;
          return std::tuple<int, int, int>{compare_ys + iy_local, kz,
                                           compare_xs + ix_local};
        };
        const auto [a_iy, a_kz, a_ix] = describe_index(max_a_index);
        const auto [b_iy, b_kz, b_ix] = describe_index(max_b_index);
        const auto [c_iy, c_kz, c_ix] = describe_index(max_c_index);
        const auto [rhs_iy, rhs_kz, rhs_ix] = describe_index(max_rhs_index);

        output.write("\n\tLaplaceCyclic GPU tridag comparison:"
                     "\n\t  max |a_cpu-a_gpu| = {:.16e} at iy={}, kz={}, ix={}"
                     "\n\t  max |b_cpu-b_gpu| = {:.16e} at iy={}, kz={}, ix={}"
                     "\n\t  max |c_cpu-c_gpu| = {:.16e} at iy={}, kz={}, ix={}"
                     "\n\t  max |rhs_cpu-rhs_gpu| = {:.16e} at iy={}, kz={}, ix={}\n",
                     max_a, a_iy, a_kz, a_ix, max_b, b_iy, b_kz, b_ix, max_c, c_iy, c_kz,
                     c_ix, max_rhs, rhs_iy, rhs_kz, rhs_ix);
        compared_device_tridag = true;
      }
    }

    {
      solveWithGpuSparseDevice(nsys, nx, localmesh->periodicX, scratch);
    }

    if (localmesh->periodicX) {
      const bout::profiling::ProfileRange range("LaplaceCyclic::subtract_x_average");
      const int block_size = 256;
      subtractPeriodicXAverage<<<ny, block_size, block_size * sizeof(double)>>>(
          scratch.x.get(), nx, ny, nmode);
      checkGpu(gpu::runtimeGetLastError(), "subtract periodic X average");
    }

    {
      const int nz = localmesh->zend - localmesh->zstart + 1;
      const int nxny = nx * ny;
      const int nmodes = (nz / 2) + 1;
      {
        const bout::profiling::ProfileRange range("LaplaceCyclic::sparse_to_fft");
        const int block_size = 256;
        const int total = nxny * nmodes;
        const int blocks = (total + block_size - 1) / block_size;
        sparseXToIfftInput<<<blocks, block_size>>>(
            scratch.x.get(), scratch.spectral_device.get(), nx, ny, nmode, nmodes,
            isGlobalFlagSet(INVERT_ZERO_DC));
        checkGpu(gpu::runtimeGetLastError(),
                 "transpose GPU sparse solution to irfft input");
      }

      {
        const bout::profiling::ProfileRange range("LaplaceCyclic::fft_backward");
        checkFft(gpu::fftExecZ2D(scratch.backward_plan.get(),
                                 scratch.spectral_device.get(),
                                 scratch.real_device.get()),
                 "resident gpu::fftExecZ2D");
      }
#if BOUT_HAS_CUDA && !BOUT_HAS_HIP
      if (auto* x_device = fieldDeviceAlias(&x(0, 0, 0))) {
        const bout::profiling::ProfileRange range(
            "LaplaceCyclic::unpack_solution_device");
        const int total = nxny * nz;
        unpackFieldFromFft<<<(total + 255) / 256, 256>>>(
            scratch.real_device.get(), x_device, nx, ny, nz, localmesh->LocalNy,
            localmesh->LocalNz, xs, ys, localmesh->zstart);
        checkGpu(cudaGetLastError(), "unpack FFT output to mapped field");
        // The returned field is immediately accessible to host boundary/comm code.
        checkGpu(cudaStreamSynchronize(0), "wait for mapped FFT output");
      } else
#endif
      {
        if (nxny * nz > scratch.real_host_size) {
          scratch.real_host.reallocate(nxny * nz);
          scratch.real_host_size = nxny * nz;
        }
        {
          const bout::profiling::ProfileRange range("LaplaceCyclic::fft_download");
          checkGpu(
              gpu::runtimeMemcpy(std::begin(scratch.real_host), scratch.real_device.get(),
                                 static_cast<std::size_t>(nxny) * nz * sizeof(BoutReal),
                                 gpu::runtimeMemcpyDeviceToHost),
              "copy resident irfft output to host");
        }

        {
          const bout::profiling::ProfileRange range(
              "LaplaceCyclic::unpack_solution_host");
          for (int ind = 0; ind < nxny; ++ind) {
            const int ix = xs + (ind / ny);
            const int iy = ys + (ind % ny);
            std::copy(std::begin(scratch.real_host) + ind * nz,
                      std::begin(scratch.real_host) + (ind + 1) * nz,
                      &(x(ix, iy, localmesh->zstart)));
          }
        }
      } // Host staging fallback
    }

    return x;
  }
#endif
  const FFTTransform transform(
      *localmesh, nmode, xs, xe, ys, ye, localmesh->zstart, localmesh->zend, inbndry,
      outbndry, isInnerBoundaryFlagSetOnFirstX(INVERT_SET),
      isOuterBoundaryFlagSetOnLastX(INVERT_SET), isGlobalFlagSet(INVERT_ZERO_DC));

  auto matrices = [&]() {
    const bout::profiling::ProfileRange range("LaplaceCyclic::transform_forward");
    return transform.forward(*this, rhs, x0, Acoef, C1coef, C2coef, Dcoef);
  }();

  // Solve tridiagonal systems
#if BOUT_HAS_CUDA || BOUT_HAS_HIP
  if (use_gpu_sparse && canUseGpuSolve(*localmesh)) {
    if (!gpu_scratch) {
      gpu_scratch = std::make_unique<LaplaceCyclicGpuScratch>();
    }
    solveWithGpuSparse(matrices.a, matrices.b, matrices.c, matrices.bcmplx, xcmplx3D,
                       localmesh->periodicX, *gpu_scratch);
  } else
#endif
  {
    {
      const bout::profiling::ProfileRange range("LaplaceCyclic::cpu_set_coefficients");
      cr->setCoefs(matrices.a, matrices.b, matrices.c);
    }
    {
      const bout::profiling::ProfileRange range("LaplaceCyclic::cpu_solve");
      cr->solve(matrices.bcmplx, xcmplx3D);
    }
  }

  if (localmesh->periodicX) {
    // Subtract X average of kz=0 mode
    std::vector<BoutReal> local(ny + 1, 0.0);
    for (int y = 0; y < ny; y++) {
      for (int ix = xs; ix <= xe; ix++) {
        local[y] += xcmplx3D(y * nmode, ix - xs).real();
      }
    }
    local[ny] = static_cast<BoutReal>(xe - xs + 1);

    // Global reduce
    std::vector<BoutReal> global(ny + 1, 0.0);
    MPI_Allreduce(local.data(), global.data(), ny + 1, MPI_DOUBLE, MPI_SUM,
                  localmesh->getXcomm());
    // Subtract average from kz=0 modes
    for (int y = 0; y < ny; y++) {
      BoutReal avg = global[y] / global[ny];
      for (int ix = xs; ix <= xe; ix++) {
        xcmplx3D(y * nmode, ix - xs) -= avg;
      }
    }
  }

  {
    const bout::profiling::ProfileRange range("LaplaceCyclic::transform_backward");
    return transform.backward(rhs, xcmplx3D);
  }
}

void LaplaceCyclic ::verify_solution(const Matrix<dcomplex>& a_ver,
                                     const Matrix<dcomplex>& b_ver,
                                     const Matrix<dcomplex>& c_ver,
                                     const Matrix<dcomplex>& r_ver,
                                     const Matrix<dcomplex>& x_sol, const int nsys) {
  output.write("Verify solution\n");
  const int nx = xe - xs + 1; // Number of X points on this processor,
                              // including boundaries but not guard cells
  const int xproc = localmesh->getXProcIndex();
  const int yproc = localmesh->getYProcIndex();
  const int nprocs = localmesh->getNXPE();
  const int myrank = yproc * nprocs + xproc;
  Matrix<dcomplex> y_ver(nsys, nx + 2);
  Matrix<dcomplex> error(nsys, nx + 2);

  MPI_Status status;
  Array<MPI_Request> request(4);
  Array<dcomplex> sbufup(nsys);
  Array<dcomplex> sbufdown(nsys);
  Array<dcomplex> rbufup(nsys);
  Array<dcomplex> rbufdown(nsys);

  // nsys = nmode * ny;  // Number of systems of equations to solve
  Matrix<dcomplex> x_ver(nsys, nx + 2);

  for (int kz = 0; kz < nsys; kz++) {
    for (int ix = 0; ix < nx; ix++) {
      x_ver(kz, ix + 1) = x_sol(kz, ix);
    }
  }

  if (xproc > 0) {
    MPI_Irecv(&rbufdown[0], nsys, MPI_DOUBLE_COMPLEX, myrank - 1, 901, BoutComm::get(),
              &request[1]);
    for (int kz = 0; kz < nsys; kz++) {
      sbufdown[kz] = x_ver(kz, 1);
    }
    MPI_Isend(&sbufdown[0], nsys, MPI_DOUBLE_COMPLEX, myrank - 1, 900, BoutComm::get(),
              &request[0]);
  }
  if (xproc < nprocs - 1) {
    MPI_Irecv(&rbufup[0], nsys, MPI_DOUBLE_COMPLEX, myrank + 1, 900, BoutComm::get(),
              &request[3]);
    for (int kz = 0; kz < nsys; kz++) {
      sbufup[kz] = x_ver(kz, nx);
    }
    MPI_Isend(&sbufup[0], nsys, MPI_DOUBLE_COMPLEX, myrank + 1, 901, BoutComm::get(),
              &request[2]);
  }

  if (xproc > 0) {
    MPI_Wait(&request[0], &status);
    MPI_Wait(&request[1], &status);
    for (int kz = 0; kz < nsys; kz++) {
      x_ver(kz, 0) = rbufdown[kz];
    }
  }
  if (xproc < nprocs - 1) {
    MPI_Wait(&request[2], &status);
    MPI_Wait(&request[3], &status);
    for (int kz = 0; kz < nsys; kz++) {
      x_ver(kz, nx + 1) = rbufup[kz];
    }
  }

  BoutReal max_error = 0.0;
  for (int kz = 0; kz < nsys; kz++) {
    for (int i = 0; i < nx; i++) {
      y_ver(kz, i) = a_ver(kz, i) * x_ver(kz, i) + b_ver(kz, i) * x_ver(kz, i + 1)
                     + c_ver(kz, i) * x_ver(kz, i + 2);
      error(kz, i) = y_ver(kz, i) - r_ver(kz, i);
      max_error = std::max(max_error, std::abs(error(kz, i)));
      output.write("abs error {}, r={}, y={}, kz {}, i {},  a={}, b={}, c={}, x-= {}, "
                   "x={}, x+ = {}\n",
                   error(kz, i).real(), r_ver(kz, i).real(), y_ver(kz, i).real(), kz, i,
                   a_ver(kz, i).real(), b_ver(kz, i).real(), c_ver(kz, i).real(),
                   x_ver(kz, i).real(), x_ver(kz, i + 1).real(), x_ver(kz, i + 2).real());
    }
  }
  output.write("max abs error {}\n", max_error);
}

#endif // BOUT_USE_METRIC_3D
