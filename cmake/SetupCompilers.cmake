# Set specific options for CUDA if enabled
if (ENABLE_RAJA AND ENABLE_CUDA)
  # RAJA requires some experimental features
  set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -arch ${CUDA_ARCH} --expt-extended-lambda --expt-relaxed-constexpr")
endif ()

# TODO Ensure openmp flags are not enabled twice!
if (ENABLE_OPENMP AND ENABLE_CUDA)
  # CMAKE_CUDA_FLAGS does not pass OpenMP_CXX_FLAGS to the host compiler by default
  set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -Xcompiler ${OpenMP_CXX_FLAGS}")
endif ()
