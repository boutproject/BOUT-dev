#include "common_transform.hxx"

#include "bout/array.hxx"
#include "bout/bout_types.hxx"
#include "bout/constants.hxx"
#include "bout/dcomplex.hxx"
#include "bout/fft.hxx"
#include "bout/field2d.hxx"
#include "bout/field3d.hxx"
#include "bout/fieldperp.hxx"
#include "bout/invert_laplace.hxx"
#include "bout/mesh.hxx"
#include "bout/openmpwrap.hxx"
#include "bout/utils.hxx"

#if BOUT_HAS_CUDA
#include <cuComplex.h>
#include <cuda_runtime.h>
#include <cufft.h>
#endif

#include <cstdint>
#include <algorithm>
#include <iterator>
#include <memory>

namespace {
#if BOUT_HAS_CUDA
static_assert(sizeof(dcomplex) == sizeof(cufftDoubleComplex),
              "dcomplex and cufftDoubleComplex must have the same memory layout");

void checkCuda(cudaError_t status, const char* call) {
  if (status != cudaSuccess) {
    throw BoutException("CUDA error in FFTTransform {}: {}", call,
                        cudaGetErrorString(status));
  }
}

void checkCufft(cufftResult status, const char* call) {
  if (status != CUFFT_SUCCESS) {
    throw BoutException("cuFFT error in FFTTransform {}: status {}", call,
                        static_cast<int>(status));
  }
}

template <class T>
class CudaBuffer {
public:
  CudaBuffer() = default;
  ~CudaBuffer() {
    if (data != nullptr) {
      cudaFree(data);
    }
  }

  CudaBuffer(const CudaBuffer&) = delete;
  CudaBuffer& operator=(const CudaBuffer&) = delete;

  T* get() { return data; }
  const T* get() const { return data; }

  void ensure(std::size_t count) {
    if (count <= capacity) {
      return;
    }
    if (data != nullptr) {
      cudaFree(data);
      data = nullptr;
    }
    checkCuda(cudaMalloc(&data, count * sizeof(T)), "cudaMalloc");
    capacity = count;
  }

private:
  T* data{nullptr};
  std::size_t capacity{0};
};

class CufftPlan {
public:
  CufftPlan() = default;
  ~CufftPlan() {
    if (plan != 0) {
      cufftDestroy(plan);
    }
  }

  CufftPlan(const CufftPlan&) = delete;
  CufftPlan& operator=(const CufftPlan&) = delete;

  cufftHandle get() const { return plan; }

  void ensure(int new_nz, int new_batch, cufftType type) {
    if (plan != 0 && new_nz == nz && new_batch == batch && type == plan_type) {
      return;
    }
    if (plan != 0) {
      cufftDestroy(plan);
      plan = 0;
    }

    int n[] = {new_nz};
    const int real_dist = new_nz;
    const int complex_dist = (new_nz / 2) + 1;
    if (type == CUFFT_D2Z) {
      checkCufft(cufftPlanMany(&plan, 1, n, nullptr, 1, real_dist, nullptr, 1,
                               complex_dist, type, new_batch),
                 "cufftPlanMany D2Z");
    } else {
      checkCufft(cufftPlanMany(&plan, 1, n, nullptr, 1, complex_dist, nullptr, 1,
                               real_dist, type, new_batch),
                 "cufftPlanMany Z2D");
    }
    nz = new_nz;
    batch = new_batch;
    plan_type = type;
  }

private:
  cufftHandle plan{0};
  int nz{0};
  int batch{0};
  cufftType plan_type{CUFFT_D2Z};
};

class CufftScratch {
public:
  void ensure(int nz, int batch) {
    const int nmodes = (nz / 2) + 1;
    const auto real_size = batch * nz;
    const auto complex_size = batch * nmodes;
    if (real_size > real_host_size) {
      real_host.reallocate(real_size);
      real_host_size = real_size;
    }
    if (complex_size > complex_host_size) {
      complex_host.reallocate(complex_size);
      complex_host_size = complex_size;
    }
    real_device.ensure(static_cast<std::size_t>(real_size));
    complex_device.ensure(static_cast<std::size_t>(complex_size));
    forward_plan.ensure(nz, batch, CUFFT_D2Z);
    backward_plan.ensure(nz, batch, CUFFT_Z2D);
  }

  Array<BoutReal> real_host;
  Array<dcomplex> complex_host;
  CudaBuffer<cufftDoubleReal> real_device;
  CudaBuffer<cufftDoubleComplex> complex_device;
  CufftPlan forward_plan;
  CufftPlan backward_plan;
  int real_host_size{0};
  int complex_host_size{0};
};

CufftScratch& cufftScratch() {
  static CufftScratch scratch;
  return scratch;
}
#endif
} // namespace

FFTTransform::FFTTransform(const Mesh& mesh, int nmode, int xs, int xe, int ys, int ye,
                           int zs, int ze, int inbndry, int outbndry,
                           bool inner_boundary_set_on_first_x,
                           bool outer_boundary_set_on_last_x, bool zero_DC)
    : zlength(getUniform(mesh.getCoordinatesConst()->zlength())), nmode(nmode), xs(xs),
      xe(xe), ys(ys), ye(ye), zs(zs), ze(ze), nx(xe - xs + 1), ny(ye - ys + 1),
      nz(ze - zs + 1), nxny(nx * ny), nsys(nmode * ny), inbndry(inbndry),
      outbndry(outbndry), inner_boundary_set_on_first_x(inner_boundary_set_on_first_x),
      outer_boundary_set_on_last_x(outer_boundary_set_on_last_x), zero_DC(zero_DC),
      localmesh(&mesh) {}

auto FFTTransform::forward(const Laplacian& laplacian, const Field3D& rhs,
                           const Field3D& x0, const Field2D& Acoef, const Field2D& C1coef,
                           const Field2D& C2coef,
                           const Field2D& Dcoef) const -> Matrices {

  Matrices result(nsys, nx);

#if BOUT_HAS_CUDA
  {
    auto& scratch = cufftScratch();
    const int nmodes = (nz / 2) + 1;
    scratch.ensure(nz, nxny);

    for (int ind = 0; ind < nxny; ++ind) {
      const int ix = xs + (ind / ny);
      const int iy = ys + (ind % ny);
      const BoutReal* input =
          (((ix < inbndry) and inner_boundary_set_on_first_x)
           || ((localmesh->LocalNx - ix - 1 < outbndry)
               and outer_boundary_set_on_last_x))
              ? &(x0(ix, iy, zs))
              : &(rhs(ix, iy, zs));
      std::copy(input, input + nz, std::begin(scratch.real_host) + ind * nz);
    }

    checkCuda(cudaMemcpy(scratch.real_device.get(), std::begin(scratch.real_host),
                         static_cast<std::size_t>(nxny) * nz * sizeof(BoutReal),
                         cudaMemcpyHostToDevice),
              "copy rfft input to device");
    checkCufft(cufftExecD2Z(scratch.forward_plan.get(), scratch.real_device.get(),
                            scratch.complex_device.get()),
               "cufftExecD2Z");
    checkCuda(cudaMemcpy(reinterpret_cast<cufftDoubleComplex*>(
                             std::begin(scratch.complex_host)),
                         scratch.complex_device.get(),
                         static_cast<std::size_t>(nxny) * nmodes
                             * sizeof(cufftDoubleComplex),
                         cudaMemcpyDeviceToHost),
              "copy rfft output to host");

    const BoutReal fac = 1.0 / nz;
    for (int ind = 0; ind < nxny; ++ind) {
      const int ix = xs + (ind / ny);
      const int iy = ys + (ind % ny);
      for (int kz = 0; kz < nmode; kz++) {
        result.bcmplx(((iy - ys) * nmode) + kz, ix - xs) =
            scratch.complex_host[ind * nmodes + kz] * fac;
      }
    }
  }
#else
  BOUT_OMP_PERF(parallel)
  {
    /// Create a local thread-scope working array
    // ZFFT routine expects input of this length
    auto k1d = Array<dcomplex>((nz / 2) + 1);

    // Loop over X and Y indices, including boundaries but not guard cells
    // (unless periodic in x)
    {
      BOUT_OMP_PERF(for)
      for (int ind = 0; ind < nxny; ++ind) {
        const int ix = xs + (ind / ny);
        const int iy = ys + (ind % ny);

        // Take FFT in Z direction, apply shift, and put result in k1d

        if (((ix < inbndry) and inner_boundary_set_on_first_x)
            || ((localmesh->LocalNx - ix - 1 < outbndry)
                and outer_boundary_set_on_last_x)) {
          // Use the values in x0 in the boundary
          rfft(&(x0(ix, iy, zs)), nz, std::begin(k1d));
        } else {
          rfft(&(rhs(ix, iy, zs)), nz, std::begin(k1d));
        }

        // Copy into array, transposing so kz is first index
        for (int kz = 0; kz < nmode; kz++) {
          result.bcmplx(((iy - ys) * nmode) + kz, ix - xs) = k1d[kz];
        }
      }
    }
  }
#endif

    // Get elements of the tridiagonal matrix
    // including boundary conditions
  BOUT_OMP_PERF(parallel)
  {
    {
      BOUT_OMP_PERF(for nowait)
      for (int ind = 0; ind < nsys; ind++) {
        const int iy = ys + (ind / nmode);
        const int kz = ind % nmode;

        const BoutReal kwave = kz * 2.0 * PI / zlength; // wave number is 1/[rad]
        laplacian.tridagMatrix(&result.a(ind, 0), &result.b(ind, 0),
                               &result.c(ind, 0), &result.bcmplx(ind, 0), iy, kz,
                               kwave, &Acoef, &C1coef, &C2coef, &Dcoef, false);
      }
    }
  }
  return result;
}

auto FFTTransform::backward(const Field3D& rhs,
                            const Matrix<dcomplex>& xcmplx3D) const -> Field3D {
  Field3D x{emptyFrom(rhs)};

#if BOUT_HAS_CUDA
  {
    auto& scratch = cufftScratch();
    const int nmodes = (nz / 2) + 1;
    scratch.ensure(nz, nxny);

    std::fill(std::begin(scratch.complex_host), std::end(scratch.complex_host),
              dcomplex{0.0, 0.0});
    for (int ind = 0; ind < nxny; ++ind) {
      const int ix = xs + (ind / ny);
      const int iy = ys + (ind % ny);

      if (zero_DC) {
        scratch.complex_host[ind * nmodes] = 0.0;
      }
      for (int kz = static_cast<int>(zero_DC); kz < nmode; kz++) {
        scratch.complex_host[ind * nmodes + kz] =
            xcmplx3D(((iy - ys) * nmode) + kz, ix - xs);
      }
    }

    checkCuda(cudaMemcpy(scratch.complex_device.get(),
                         reinterpret_cast<cufftDoubleComplex*>(
                             std::begin(scratch.complex_host)),
                         static_cast<std::size_t>(nxny) * nmodes
                             * sizeof(cufftDoubleComplex),
                         cudaMemcpyHostToDevice),
              "copy irfft input to device");
    checkCufft(cufftExecZ2D(scratch.backward_plan.get(), scratch.complex_device.get(),
                            scratch.real_device.get()),
               "cufftExecZ2D");
    checkCuda(cudaMemcpy(std::begin(scratch.real_host), scratch.real_device.get(),
                         static_cast<std::size_t>(nxny) * nz * sizeof(BoutReal),
                         cudaMemcpyDeviceToHost),
              "copy irfft output to host");

    for (int ind = 0; ind < nxny; ++ind) {
      const int ix = xs + (ind / ny);
      const int iy = ys + (ind % ny);
      std::copy(std::begin(scratch.real_host) + ind * nz,
                std::begin(scratch.real_host) + (ind + 1) * nz, &(x(ix, iy, zs)));
    }
  }
#else
  // FFT back to real space
  BOUT_OMP_PERF(parallel)
  {
    /// Create a local thread-scope working array
    // ZFFT routine expects input of this length
    auto k1d = Array<dcomplex>((nz / 2) + 1);

    {
      BOUT_OMP_PERF(for nowait)
      for (int ind = 0; ind < nxny; ++ind) { // Loop over X and Y
        const int ix = xs + (ind / ny);
        const int iy = ys + (ind % ny);

        if (zero_DC) {
          k1d[0] = 0.;
        }

        for (int kz = static_cast<int>(zero_DC); kz < nmode; kz++) {
          k1d[kz] = xcmplx3D(((iy - ys) * nmode) + kz, ix - xs);
        }

        for (int kz = nmode; kz < (nz / 2) + 1; kz++) {
          k1d[kz] = 0.0; // Filtering out all higher harmonics
        }

        irfft(std::begin(k1d), nz, &(x(ix, iy, zs)));
      }
    }
  }
#endif

  return x;
}

auto FFTTransform::forward(const Laplacian& laplacian, const FieldPerp& rhs,
                           const FieldPerp& x0, const Field2D& Acoef,
                           const Field2D& C1coef, const Field2D& C2coef,
                           const Field2D& Dcoef) const -> Matrices {

  Matrices result(nmode, nx);
  const int jy = rhs.getIndex();

  BOUT_OMP_PERF(parallel)
  {
    /// Create a local thread-scope working array
    // ZFFT routine expects input of this length
    auto k1d = Array<dcomplex>((nz / 2) + 1);

    // Loop over X indices, including boundaries but not guard
    // cells (unless periodic in x)
    BOUT_OMP_PERF(for)
    for (int ix = xs; ix <= xe; ix++) {
      // Take FFT in Z direction, apply shift, and put result in k1d

      if (((ix < inbndry) and inner_boundary_set_on_first_x)
          || ((localmesh->LocalNx - ix - 1 < outbndry)
              and outer_boundary_set_on_last_x)) {
        // Use the values in x0 in the boundary
        rfft(&(x0(ix, zs)), nz, std::begin(k1d));
      } else {
        rfft(&(rhs(ix, zs)), nz, std::begin(k1d));
      }

      // Copy into array, transposing so kz is first index
      for (int kz = 0; kz < nmode; kz++) {
        result.bcmplx(kz, ix - xs) = k1d[kz];
      }
    }

    // Get elements of the tridiagonal matrix
    // including boundary conditions
    BOUT_OMP_PERF(for nowait)
    for (int kz = 0; kz < nmode; kz++) {
      const BoutReal kwave = kz * 2.0 * PI / zlength; // wave number is 1/[rad]
      laplacian.tridagMatrix(&result.a(kz, 0), &result.b(kz, 0), &result.c(kz, 0),
                             &result.bcmplx(kz, 0), jy, kz, kwave, &Acoef, &C1coef,
                             &C2coef, &Dcoef, false);
    }
  }
  return result;
}

auto FFTTransform::backward(const FieldPerp& rhs,
                            const Matrix<dcomplex>& xcmplx) const -> FieldPerp {
  FieldPerp x{emptyFrom(rhs)};

  // FFT back to real space
  BOUT_OMP_PERF(parallel)
  {
    /// Create a local thread-scope working array
    // ZFFT routine expects input of this length
    auto k1d = Array<dcomplex>((nz / 2) + 1);

    BOUT_OMP_PERF(for nowait)
    for (int ix = xs; ix <= xe; ++ix) {
      if (zero_DC) {
        k1d[0] = 0.;
      }

      for (int kz = static_cast<int>(zero_DC); kz < nmode; kz++) {
        k1d[kz] = xcmplx(kz, ix - xs);
      }

      for (int kz = nmode; kz < (nz / 2) + 1; kz++) {
        k1d[kz] = 0.0; // Filtering out all higher harmonics
      }

      irfft(std::begin(k1d), nz, &(x(ix, zs)));
    }
  }

  checkData(x);

  return x;
}

DSTTransform::DSTTransform(const Mesh& mesh, int nmode, int xs, int xe, int ys, int ye,
                           int zs, int ze, int inbndry, int outbndry,
                           bool inner_boundary_set_on_first_x,
                           bool outer_boundary_set_on_last_x, bool zero_DC)
    : zlength(getUniform(mesh.getCoordinatesConst()->zlength())), nmode(nmode), xs(xs),
      xe(xe), ys(ys), ye(ye), zs(zs), ze(ze), nx(xe - xs + 1), ny(ye - ys + 1),
      nz(ze - zs + 1), nxny(nx * ny), nsys(nmode * ny), inbndry(inbndry),
      outbndry(outbndry), inner_boundary_set_on_first_x(inner_boundary_set_on_first_x),
      outer_boundary_set_on_last_x(outer_boundary_set_on_last_x), zero_DC(zero_DC),
      localmesh(&mesh) {}

auto DSTTransform::forward(const Laplacian& laplacian, const Field3D& rhs,
                           const Field3D& x0, const Field2D& Acoef, const Field2D& C1coef,
                           const Field2D& C2coef,
                           const Field2D& Dcoef) const -> Matrices {

  Matrices result(nsys, nx);

  BOUT_OMP_PERF(parallel)
  {
    /// Create a local thread-scope working array
    // ZDST routine expects input of this length
    auto k1d = Array<dcomplex>(nz);

    // Loop over X and Y indices, including boundaries but not guard cells
    // (unless periodic in x)
    BOUT_OMP_PERF(for)
    for (int ind = 0; ind < nxny; ++ind) {
      const int ix = xs + (ind / ny);
      const int iy = ys + (ind % ny);

      // Take DST in Z direction, apply shift, and put result in k1d

      if (((ix < inbndry) and inner_boundary_set_on_first_x)
          || ((localmesh->LocalNx - ix - 1 < outbndry)
              and outer_boundary_set_on_last_x)) {
        // Use the values in x0 in the boundary
        bout::fft::DST(&(x0(ix, iy, zs + 1)), nz - 2, std::begin(k1d));
      } else {
        bout::fft::DST(&(rhs(ix, iy, zs + 1)), nz - 2, std::begin(k1d));
      }

      // Copy into array, transposing so kz is first index
      for (int kz = 0; kz < nmode; kz++) {
        result.bcmplx(((iy - ys) * nmode) + kz, ix - xs) = k1d[kz];
      }
    }

    // Get elements of the tridiagonal matrix
    // including boundary conditions
    BOUT_OMP_PERF(for nowait)
    for (int ind = 0; ind < nsys; ind++) {
      const int iy = ys + (ind / nmode);
      const int kz = ind % nmode;

      const BoutReal kwave = kz * 2.0 * PI / zlength; // wave number is 1/[rad]
      laplacian.tridagMatrix(&result.a(ind, 0), &result.b(ind, 0), &result.c(ind, 0),
                             &result.bcmplx(ind, 0), iy, kz, kwave, &Acoef, &C1coef,
                             &C2coef, &Dcoef, false);
    }
  }
  return result;
}

auto DSTTransform::backward(const Field3D& rhs,
                            const Matrix<dcomplex>& xcmplx3D) const -> Field3D {
  Field3D x{emptyFrom(rhs)};

  // DST back to real space
  BOUT_OMP_PERF(parallel)
  {
    /// Create a local thread-scope working array
    // ZDST routine expects input of this length
    auto k1d = Array<dcomplex>((nz / 2) + 1);

    BOUT_OMP_PERF(for nowait)
    for (int ind = 0; ind < nxny; ++ind) { // Loop over X and Y
      const int ix = xs + (ind / ny);
      const int iy = ys + (ind % ny);

      if (zero_DC) {
        k1d[0] = 0.;
      }

      for (int kz = static_cast<int>(zero_DC); kz < nmode; kz++) {
        k1d[kz] = xcmplx3D(((iy - ys) * nmode) + kz, ix - xs);
      }

      for (int kz = nmode; kz < (nz / 2) + 1; kz++) {
        k1d[kz] = 0.0; // Filtering out all higher harmonics
      }

      bout::fft::DST_rev(std::begin(k1d), nz - 2, &(x(ix, iy, zs + 1)));

      x(ix, iy, zs) = -x(ix, iy, zs + 2);
      x(ix, iy, ze) = -x(ix, iy, ze - 2);
    }
  }

  return x;
}

auto DSTTransform::forward(const Laplacian& laplacian, const FieldPerp& rhs,
                           const FieldPerp& x0, const Field2D& Acoef,
                           const Field2D& C1coef, const Field2D& C2coef,
                           const Field2D& Dcoef) const -> Matrices {

  Matrices result(nmode, nx);
  const int jy = rhs.getIndex();

  BOUT_OMP_PERF(parallel)
  {
    /// Create a local thread-scope working array
    // ZDST routine expects input of this length
    auto k1d = Array<dcomplex>((nz / 2) + 1);

    // Loop over X indices, including boundaries but not guard
    // cells (unless periodic in x)
    BOUT_OMP_PERF(for)
    for (int ix = xs; ix <= xe; ix++) {
      // Take DST in Z direction, apply shift, and put result in k1d

      if (((ix < inbndry) and inner_boundary_set_on_first_x)
          || ((localmesh->LocalNx - ix - 1 < outbndry)
              and outer_boundary_set_on_last_x)) {
        // Use the values in x0 in the boundary
        rfft(&(x0(ix, zs)), nz, std::begin(k1d));
      } else {
        rfft(&(rhs(ix, zs)), nz, std::begin(k1d));
      }

      // Copy into array, transposing so kz is first index
      for (int kz = 0; kz < nmode; kz++) {
        result.bcmplx(kz, ix - xs) = k1d[kz];
      }
    }

    // Get elements of the tridiagonal matrix
    // including boundary conditions
    BOUT_OMP_PERF(for nowait)
    for (int kz = 0; kz < nmode; kz++) {
      const BoutReal kwave = kz * 2.0 * PI / zlength; // wave number is 1/[rad]
      laplacian.tridagMatrix(&result.a(kz, 0), &result.b(kz, 0), &result.c(kz, 0),
                             &result.bcmplx(kz, 0), jy, kz, kwave, &Acoef, &C1coef,
                             &C2coef, &Dcoef, false);
    }
  }
  return result;
}

auto DSTTransform::backward(const FieldPerp& rhs,
                            const Matrix<dcomplex>& xcmplx) const -> FieldPerp {
  FieldPerp x{emptyFrom(rhs)};

  // DST back to real space
  BOUT_OMP_PERF(parallel)
  {
    /// Create a local thread-scope working array
    // ZDST routine expects input of this length
    auto k1d = Array<dcomplex>((nz / 2) + 1);

    BOUT_OMP_PERF(for nowait)
    for (int ix = xs; ix < xe; ++ix) {
      if (zero_DC) {
        k1d[0] = 0.;
      }

      for (int kz = static_cast<int>(zero_DC); kz < nmode; kz++) {
        k1d[kz] = xcmplx(kz, ix - xs);
      }

      for (int kz = nmode; kz < (nz / 2) + 1; kz++) {
        k1d[kz] = 0.0; // Filtering out all higher harmonics
      }

      irfft(std::begin(k1d), nz, &(x(ix, zs)));
    }
  }

  checkData(x);

  return x;
}
