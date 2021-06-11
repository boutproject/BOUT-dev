set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}/cmake/")


set(BOUT_DEPENDS "")

# determined in SetupCompilers.cmake
if (HAVE_MPI)
   list(APPEND BOUT_DEPENDS mpi)
   target_link_libraries(bout++ PUBLIC MPI::MPI_CXX mpark_variant)
endif ()

# determined in SetupCompilers.cmake
if (HAVE_OPENMP)
   list(APPEND BOUT_DEPENDS openmp)
   target_link_libraries(bout++ PUBLIC OpenMP::OpenMP_CXX)
   set(possible_openmp_schedules static dynamic guided auto)
   if (NOT BOUT_OPENMP_SCHEDULE IN_LIST possible_openmp_schedules)
    message(FATAL_ERROR "BOUT_OPENMP_SCHEDULE must be one of ${possible_openmp_schedules}; got ${BOUT_OPENMP_SCHEDULE}")
   endif()
   message(STATUS "OpenMP schedule: ${BOUT_OPENMP_SCHEDULE}")
endif()
message(STATUS "Enable OpenMP: ${BOUT_ENABLE_OPENMP}")
set(BOUT_USE_OPENMP ${BOUT_ENABLE_OPENMP})


# determined in SetupCompilers.cmake
if (BOUT_USE_CUDA)
   list(APPEND BOUT_DEPENDS cuda)
   enable_language(CUDA)
   message(STATUS "BOUT_USE_CUDA ${CMAKE_CUDA_COMPILER}")
   set_source_files_properties(${BOUT_SOURCES} PROPERTIES LANGUAGE CUDA )
   # CMAKE 3.14 if we don't use deprecated FindCUDA, then need to compute CUDA_TOOLKIT_ROOT_DIR
   #cmake 3.17 has FindCUDAToolkit
   get_filename_component(cuda_bin_dir ${CMAKE_CUDA_COMPILER} DIRECTORY)
   set(BOUT_CUDA_LIB_DIR ${cuda_bin_dir}/../lib64 CACHE STRING "CUDA Library DIR")
   message(STATUS "CUDA LIBRARY DIR: ${BOUT_CUDA_LIB_DIR}")
   set_target_properties(bout++ PROPERTIES CUDA_STANDARD 14)
   set_target_properties(bout++ PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
   set_target_properties(bout++ PROPERTIES POSITION_INDEPENDENT_CODE ON)
   set_target_properties(bout++ PROPERTIES LINKER_LANGUAGE CUDA)
   #target_compile_definitions(bout++ PUBLIC "BOUT_USE_CUDA") 
endif ()   

# Caliper
set(BOUT_HAS_CALIPER OFF)
if (ENABLE_CALIPER)
  find_package(caliper REQUIRED)
  list(APPEND BOUT_DEPENDS caliper)
  set (BOUT_HAS_CALIPER ON)
  target_compile_definitions(bout++ PUBLIC "BOUT_HAS_CALIPER")
  target_include_directories(bout++ PUBLIC ${caliper_INCLUDE_DIR})
  target_link_libraries(bout++ PUBLIC caliper)
endif ()

# UMPIRE
set(BOUT_HAS_UMPIRE OFF)
if (ENABLE_UMPIRE)
  find_package(UMPIRE REQUIRED)
  list(APPEND BOUT_DEPENDS umpire)
  set (BOUT_HAS_UMPIRE ON)
  target_compile_definitions(bout++ PUBLIC "BOUT_HAS_UMPIRE")
  target_include_directories(bout++ PUBLIC ${UMPIRE_INCLUDE_DIRS}/include)
  target_link_libraries(bout++ PUBLIC umpire)
endif ()

# RAJA
set(BOUT_HAS_RAJA OFF)   
if (ENABLE_RAJA)
  find_package(RAJA REQUIRED)
  list(APPEND BOUT_DEPENDS RAJA)
  message(STATUS "RAJA_CONFIG:" ${RAJA_CONFIG}) 
  string(FIND ${RAJA_CONFIG} "raja" loc)
  math(EXPR value "${loc} + 5" OUTPUT_FORMAT DECIMAL)
  string(SUBSTRING ${RAJA_CONFIG} 0 ${value}  RAJA_PATH)
  message(STATUS "RAJA_PATH" ${RAJA_PATH})
  target_include_directories(bout++ PUBLIC ${RAJA_PATH}/include)
  target_link_libraries(bout++ PUBLIC RAJA)
  set(BOUT_HAS_RAJA ON)
  target_compile_definitions(bout++ PUBLIC "BOUT_HAS_RAJA")
endif ()


#HAVE_HYPRE
set(BOUT_HAS_HYPRE OFF)
set(HYPRE_DIR "" CACHE STRING "Point to HYPRE Install")
if (ENABLE_HYPRE OR HYPRE_DIR)
  enable_language(C)
  find_package(HYPRE REQUIRED)
  list(APPEND BOUT_DEPENDS HYPRE)
  set (BOUT_HAS_HYPRE ON)
  set (ENABLE_HYPRE ON) # may just have HYPRE_DIR
  target_compile_definitions(bout++ PUBLIC "BOUT_HAS_HYPRE")
  #target_include_directories(bout++ PUBLIC ${HYPRE_DIR}/include)
  target_link_libraries(bout++ PUBLIC HYPRE::HYPRE)
  if (HYPRE_CUDA)
     target_compile_definitions(bout++ PUBLIC "HYPRE_USING_CUDA;HYPRE_USING_UNIFIED_MEMORY")
     target_link_libraries(bout++ PUBLIC "${BOUT_CUDA_LIB_DIR}/libcusparse_static.a;${BOUT_CUDA_LIB_DIR}/libcurand_static.a;${BOUT_CUDA_LIB_DIR}/libculibos.a;${BOUT_CUDA_LIB_DIR}/libcublas_static.a;${BOUT_CUDA_LIB_DIR}/libcublasLt_static.a")
  endif ()
endif ()


#PETSc
option(BOUT_USE_PETSC "Enable support for PETSc time solvers and inversions" On)
if (BOUT_USE_PETSC)
  if (NOT CMAKE_SYSTEM_NAME STREQUAL "CrayLinuxEnvironment")
    # Cray wrappers sort this out for us
    find_package(PETSc REQUIRED)
    list(APPEND BOUT_DEPENDS PETSc)
    target_link_libraries(bout++ PUBLIC PETSc::PETSc)
  endif()
endif()
message(STATUS "PETSc support: ${BOUT_USE_PETSC}")
set(BOUT_HAS_PETSC ${BOUT_USE_PETSC})

cmake_dependent_option(BOUT_USE_SYSTEM_MPARK_VARIANT "Use external installation of mpark.variant" OFF
   "BOUT_UPDATE_GIT_SUBMODULE OR EXISTS ${PROJECT_SOURCE_DIR}/externalpackages/mpark.variant/CMakeLists.txt" ON)
message(STATUS "BOUT_USE_SYSTEM_MPARK_VARIANT : ${BOUT_USE_SYSTEM_MPARK_VARIANT}") 
if(BOUT_USE_SYSTEM_MPARK_VARIANT)
  message(STATUS "Using external mpark.variant")
  find_package(mpark_variant REQUIRED)
else()
  message(STATUS "Using mpark.variant submodule")
  bout_update_submodules()
  add_subdirectory(externalpackages/mpark.variant)
  if(NOT TARGET mpark_variant)
    message(FATAL_ERROR "mpark_variant not found! Have you disabled the git submodules (BOUT_UPDATE_GIT_SUBMODULE)?")
  endif()
endif()
target_link_libraries(bout++ PUBLIC mpark_variant)

cmake_dependent_option(BOUT_USE_SYSTEM_FMT "Use external installation of fmt" OFF
  "BOUT_UPDATE_GIT_SUBMODULE OR EXISTS ${PROJECT_SOURCE_DIR}/externalpackages/fmt/CMakeLists.txt" ON)
if(BOUT_USE_SYSTEM_FMT)
  message(STATUS "Using external fmt")
  find_package(fmt 6 REQUIRED)
else()
  message(STATUS "Using fmt submodule")
  bout_update_submodules()
  # Need to install fmt alongside BOUT++
  set(FMT_INSTALL ON CACHE BOOL "")
  add_subdirectory(externalpackages/fmt)
  if(NOT TARGET fmt::fmt)
    message(FATAL_ERROR "fmt not found! Have you disabled the git submodules (BOUT_UPDATE_GIT_SUBMODULE)?")
  endif()
endif()
target_link_libraries(bout++ PUBLIC fmt::fmt)


option(BOUT_USE_PVODE "Enable support for bundled PVODE" ON)
if (BOUT_USE_PVODE)
  add_subdirectory(externalpackages/PVODE)
  target_link_libraries(bout++ PUBLIC pvode pvpre)
endif()
message(STATUS "PVODE support: ${BOUT_USE_PVODE}")
set(BOUT_HAS_PVODE ${BOUT_USE_PVODE})

#option(BOUT_USE_NETCDF "Enable support for NetCDF output" ON)
#if (BOUT_USE_NETCDF)
#  find_package(NetCDF REQUIRED)
#  target_link_libraries(bout++ PUBLIC NetCDF::NetCDF_CXX)
#endif()
#message(STATUS "NetCDF support: ${BOUT_USE_NETCDF}")
#set(BOUT_HAS_NETCDF ${BOUT_USE_NETCDF})

option(BOUT_USE_HDF5 "Enable support for HDF5 output" OFF)
if (BOUT_USE_HDF5)
  find_package(HDF5 REQUIRED COMPONENTS CXX)
  target_link_libraries(bout++ PUBLIC "${HDF5_CXX_LIBRARIES}")
  target_include_directories(bout++ PUBLIC "${HDF5_CXX_INCLUDE_DIRS}")
endif()
message(STATUS "HDF5 support: ${BOUT_USE_HDF5}")
set(BOUT_HAS_HDF5 ${BOUT_USE_HDF5})

option(BOUT_USE_FFTW "Enable support for FFTW" ON)
if (BOUT_USE_FFTW)
  find_package(FFTW REQUIRED)
  target_link_libraries(bout++ PUBLIC FFTW::FFTW)
endif()
message(STATUS "FFTW support: ${BOUT_USE_FFTW}")
set(BOUT_HAS_FFTW ${BOUT_USE_FFTW})

option(BOUT_USE_LAPACK "Enable support for LAPACK" ON)
if (BOUT_USE_LAPACK)
  if (NOT CMAKE_SYSTEM_NAME STREQUAL "CrayLinuxEnvironment")
    # Cray wrappers sort this out for us
    find_package(LAPACK REQUIRED)
    target_link_libraries(bout++ PUBLIC "${LAPACK_LIBRARIES}")
  endif()
endif()
message(STATUS "LAPACK support: ${BOUT_USE_LAPACK}")
set(BOUT_HAS_LAPACK ${BOUT_USE_LAPACK})


option(BOUT_USE_SLEPC "Enable support for SLEPc eigen solver" OFF)
if (BOUT_USE_SLEPC)
  find_package(SLEPc REQUIRED)
  target_link_libraries(bout++ PUBLIC SLEPc::SLEPc)
endif()
message(STATUS "SLEPc support: ${BOUT_USE_SLEPC}")
set(BOUT_HAS_SLEPC ${BOUT_USE_SLEPC})

option(BOUT_USE_SUNDIALS "Enable support for SUNDIALS time solvers" On)
if (BOUT_USE_SUNDIALS)
  find_package(SUNDIALS REQUIRED)
  target_link_libraries(bout++ PUBLIC SUNDIALS::cvode)
  target_link_libraries(bout++ PUBLIC SUNDIALS::ida)
  target_link_libraries(bout++ PUBLIC SUNDIALS::arkode)
endif()
message(STATUS "SUNDIALS support: ${BOUT_USE_SUNDIALS}")
set(BOUT_HAS_SUNDIALS ${BOUT_USE_SUNDIALS})
set(BOUT_HAS_ARKODE ${BOUT_USE_SUNDIALS})
set(BOUT_HAS_CVODE ${BOUT_USE_SUNDIALS})
set(BOUT_HAS_IDA ${BOUT_USE_SUNDIALS})

option(BOUT_USE_NLS "Enable Native Language Support" ON)
if (BOUT_USE_NLS)
  find_package(Gettext)
  if (GETTEXT_FOUND)
	find_package(Intl)
	if (Intl_FOUND)
	  target_link_libraries(bout++
		PUBLIC ${Intl_LIBRARIES})
	  target_include_directories(bout++
		PUBLIC ${Intl_INCLUDE_DIRS})
	endif()
  endif()
endif()
set(BOUT_HAS_GETTEXT ${BOUT_USE_NLS})

option(BOUT_USE_SCOREP "Enable support for Score-P based instrumentation" OFF)
if (BOUT_USE_SCOREP)
  message(STATUS "Score-P support enabled. Please make sure you are calling CMake like so:

  SCOREP_WRAPPER=off cmake -DCMAKE_C_COMPILER=scorep-mpicc -DCMAKE_CXX_COMPILER=scorep-mpicxx <other CMake options>
")
endif()
set(BOUT_HAS_SCOREP ${BOUT_USE_SCOREP})

message(STATUS "BOUT_DEPENDS: ${BOUT_DEPENDS}")
