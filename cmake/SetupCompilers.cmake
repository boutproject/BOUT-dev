# Note: Currently BOUT++ always needs MPI. This option just determines
# whether the find_* routines are used
option(BOUT_ENABLE_MPI "Enable MPI support" ON)
if(BOUT_ENABLE_MPI)
  # This might not be entirely sensible, but helps CMake to find the
  # correct MPI, workaround for https://gitlab.kitware.com/cmake/cmake/issues/18895
  find_program(MPIEXEC_EXECUTABLE NAMES mpiexec mpirun)
  find_package(MPI REQUIRED)
endif()
set(BOUT_USE_MPI ${BOUT_ENABLE_MPI})

option(BOUT_ENABLE_OPENMP "Enable OpenMP support" OFF)
set(BOUT_OPENMP_SCHEDULE
    static
    CACHE STRING "Set OpenMP schedule"
)
set_property(
  CACHE BOUT_OPENMP_SCHEDULE PROPERTY STRINGS static dynamic guided auto
)
if(BOUT_ENABLE_OPENMP)
  find_package(OpenMP REQUIRED)
  set(possible_openmp_schedules static dynamic guided auto)
  if(NOT BOUT_OPENMP_SCHEDULE IN_LIST possible_openmp_schedules)
    message(
      FATAL_ERROR
        "BOUT_OPENMP_SCHEDULE must be one of ${possible_openmp_schedules}; got ${BOUT_OPENMP_SCHEDULE}"
    )
  endif()
  message(STATUS "OpenMP schedule: ${BOUT_OPENMP_SCHEDULE}")
endif()
set(BOUT_USE_OPENMP ${BOUT_ENABLE_OPENMP})
message(STATUS "Enable OpenMP: ${BOUT_ENABLE_OPENMP}")

option(BOUT_ENABLE_CUDA "Enable CUDA support" OFF)
option(BOUT_ENABLE_HIP "Enable AMD HIP support" OFF)
if(BOUT_ENABLE_CUDA AND BOUT_ENABLE_HIP)
  message(
    FATAL_ERROR "BOUT_ENABLE_CUDA and BOUT_ENABLE_HIP are mutually exclusive"
  )
endif()
if(BOUT_ENABLE_HIP)
  if(CMAKE_VERSION VERSION_LESS 3.21)
    message(FATAL_ERROR "BOUT_ENABLE_HIP requires CMake 3.21 or newer")
  endif()
  if(NOT BOUT_ENABLE_RAJA OR NOT BOUT_ENABLE_UMPIRE)
    message(
      FATAL_ERROR
        "BOUT_ENABLE_HIP requires BOUT_ENABLE_RAJA and BOUT_ENABLE_UMPIRE"
    )
  endif()
  # Keep compiler detection and compilation on the same ROCm installation.
  # Otherwise a Spack view can be mistaken for an implicit include directory.
  if(CMAKE_HIP_COMPILER_ROCM_ROOT)
    string(APPEND CMAKE_HIP_FLAGS_INIT
           " --rocm-path=${CMAKE_HIP_COMPILER_ROCM_ROOT}"
    )
  endif()
  enable_language(HIP)
  set(CMAKE_HIP_STANDARD 20)
  set(CMAKE_HIP_STANDARD_REQUIRED ON)
  set(CMAKE_HIP_EXTENSIONS OFF)
endif()
set(BOUT_HAS_HIP ${BOUT_ENABLE_HIP})
set(CUDA_ARCH
    "compute_70,code=sm_70"
    CACHE STRING "CUDA architecture"
)
if(BOUT_ENABLE_CUDA)
  # Set specific options for CUDA if enabled
  enable_language(CUDA)
  set(CMAKE_CUDA_FLAGS
      "${CMAKE_CUDA_FLAGS} -gencode arch=${CUDA_ARCH} -ccbin ${CMAKE_CXX_COMPILER}"
  )
  if(BOUT_ENABLE_RAJA)
    # RAJA uses lambda expressions
    set(CMAKE_CUDA_FLAGS
        "${CMAKE_CUDA_FLAGS} --expt-extended-lambda --expt-relaxed-constexpr"
    )
  endif()

  # TODO Ensure openmp flags are not enabled twice!
  if(BOUT_ENABLE_OPENMP)
    # CMAKE_CUDA_FLAGS does not pass OpenMP_CXX_FLAGS to the host compiler by default
    set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -Xcompiler ${OpenMP_CXX_FLAGS}")
  endif()
endif()
set(BOUT_HAS_CUDA ${BOUT_ENABLE_CUDA})
