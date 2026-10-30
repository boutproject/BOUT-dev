// Private CUDA/HIP adapter for the cyclic Laplacian solver.
#pragma once

#include <bout/build_defines.hxx>
#include <cstddef>

#if BOUT_HAS_CUDA || BOUT_HAS_HIP
#if BOUT_HAS_HIP
#include <hip/hip_complex.h>
#include <hip/hip_runtime.h>
#include <hipfft/hipfft.h>
#include <hipsparse/hipsparse.h>
#else
#include <cuComplex.h>
#include <cuda_runtime.h>
#include <cufft.h>
#include <cusparse.h>
#endif

namespace gpu {
#if BOUT_HAS_HIP
inline constexpr const char* backend = "HIP";
using runtimeError_t = hipError_t;
using Complex = hipDoubleComplex;
using sparseStatus_t = hipsparseStatus_t;
using sparseHandle_t = hipsparseHandle_t;
using fftResult = hipfftResult;
using fftHandle = hipfftHandle;
using fftType = hipfftType;
using fftDoubleComplex = hipfftDoubleComplex;
using fftDoubleReal = hipfftDoubleReal;
inline constexpr auto runtimeSuccess = hipSuccess;
inline constexpr auto runtimeMemcpyHostToDevice = hipMemcpyHostToDevice;
inline constexpr auto runtimeMemcpyDeviceToHost = hipMemcpyDeviceToHost;
inline constexpr auto SPARSE_STATUS_SUCCESS = HIPSPARSE_STATUS_SUCCESS;
inline constexpr auto FFT_SUCCESS = HIPFFT_SUCCESS;
inline constexpr auto FFT_D2Z = HIPFFT_D2Z;
inline constexpr auto FFT_Z2D = HIPFFT_Z2D;
inline constexpr auto runtimeGetErrorString = hipGetErrorString;
inline constexpr auto runtimeFree = hipFree;
inline constexpr auto runtimeMemcpy = hipMemcpy;
inline constexpr auto runtimeGetLastError = hipGetLastError;
inline constexpr auto sparseCreate = hipsparseCreate;
inline constexpr auto sparseDestroy = hipsparseDestroy;
inline constexpr auto sparseZgtsv2StridedBatch_bufferSizeExt =
    hipsparseZgtsv2StridedBatch_bufferSizeExt;
inline constexpr auto sparseZgtsv2StridedBatch = hipsparseZgtsv2StridedBatch;
inline constexpr auto fftDestroy = hipfftDestroy;
inline constexpr auto fftPlanMany = hipfftPlanMany;
inline constexpr auto fftExecD2Z = hipfftExecD2Z;
inline constexpr auto fftExecZ2D = hipfftExecZ2D;
template <class T>
inline auto runtimeMalloc(T** ptr, std::size_t bytes) {
  return hipMalloc(reinterpret_cast<void**>(ptr), bytes);
}
__host__ __device__ inline Complex makeComplex(double re, double im) {
  return make_hipDoubleComplex(re, im);
}
__host__ __device__ inline double real(Complex value) { return hipCreal(value); }
__host__ __device__ inline double imag(Complex value) { return hipCimag(value); }
#else
inline constexpr const char* backend = "CUDA";
using runtimeError_t = cudaError_t;
using Complex = cuDoubleComplex;
using sparseStatus_t = cusparseStatus_t;
using sparseHandle_t = cusparseHandle_t;
using fftResult = cufftResult;
using fftHandle = cufftHandle;
using fftType = cufftType;
using fftDoubleComplex = cufftDoubleComplex;
using fftDoubleReal = cufftDoubleReal;
inline constexpr auto runtimeSuccess = cudaSuccess;
inline constexpr auto runtimeMemcpyHostToDevice = cudaMemcpyHostToDevice;
inline constexpr auto runtimeMemcpyDeviceToHost = cudaMemcpyDeviceToHost;
inline constexpr auto SPARSE_STATUS_SUCCESS = CUSPARSE_STATUS_SUCCESS;
inline constexpr auto FFT_SUCCESS = CUFFT_SUCCESS;
inline constexpr auto FFT_D2Z = CUFFT_D2Z;
inline constexpr auto FFT_Z2D = CUFFT_Z2D;
inline constexpr auto runtimeGetErrorString = cudaGetErrorString;
inline constexpr auto runtimeFree = cudaFree;
inline constexpr auto runtimeMemcpy = cudaMemcpy;
inline constexpr auto runtimeGetLastError = cudaGetLastError;
inline constexpr auto sparseCreate = cusparseCreate;
inline constexpr auto sparseDestroy = cusparseDestroy;
inline constexpr auto sparseZgtsv2StridedBatch_bufferSizeExt =
    cusparseZgtsv2StridedBatch_bufferSizeExt;
inline constexpr auto sparseZgtsv2StridedBatch = cusparseZgtsv2StridedBatch;
inline constexpr auto fftDestroy = cufftDestroy;
inline constexpr auto fftPlanMany = cufftPlanMany;
inline constexpr auto fftExecD2Z = cufftExecD2Z;
inline constexpr auto fftExecZ2D = cufftExecZ2D;
template <class T>
inline auto runtimeMalloc(T** ptr, std::size_t bytes) {
  return cudaMalloc(reinterpret_cast<void**>(ptr), bytes);
}
__host__ __device__ inline Complex makeComplex(double re, double im) {
  return make_cuDoubleComplex(re, im);
}
__host__ __device__ inline double real(Complex value) { return cuCreal(value); }
__host__ __device__ inline double imag(Complex value) { return cuCimag(value); }
#endif
} // namespace gpu
#endif
