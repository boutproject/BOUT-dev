# This might not be entirely sensible, but helps CMake to find the
# correct MPI, workaround for https://gitlab.kitware.com/cmake/cmake/issues/18895
find_program(MPIEXEC_EXECUTABLE NAMES mpiexec mpirun)
find_package(MPI REQUIRED)

# Set specific options for CUDA if enabled
if(ENABLE_CUDA)
   enable_language(CUDA)
   set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -arch ${CUDA_ARCH} -ccbin ${CMAKE_CXX_COMPILER}")
   if (ENABLE_RAJA)
   # RAJA uses lambda expressions
      set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} --expt-extended-lambda --expt-relaxed-constexpr")
   endif ()

# TODO Ensure openmp flags are not enabled twice!
   if (ENABLE_OPENMP)
   # CMAKE_CUDA_FLAGS does not pass OpenMP_CXX_FLAGS to the host compiler by default
      set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -Xcompiler ${OpenMP_CXX_FLAGS}")
   endif ()
   set(BOUT_USE_CUDA ON)
endif()

