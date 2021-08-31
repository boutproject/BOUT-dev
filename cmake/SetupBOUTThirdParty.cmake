set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}/cmake/")


set(BOUT_DEPENDS "")

# determined in SetupCompilers.cmake
if (BOUT_USE_MPI)
  list(APPEND BOUT_DEPENDS mpi)
  target_link_libraries(bout++ PUBLIC MPI::MPI_CXX)
endif ()

# determined in SetupCompilers.cmake
if (BOUT_USE_OPENMP)
  list(APPEND BOUT_DEPENDS openmp)
  target_link_libraries(bout++ PUBLIC OpenMP::OpenMP_CXX)
endif()

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
endif ()   

# Caliper
option(BOUT_ENABLE_CALIPER "Enable Caliper" OFF)
if (BOUT_ENABLE_CALIPER)
  find_package(caliper REQUIRED)
  list(APPEND BOUT_DEPENDS caliper)
  target_include_directories(bout++ PUBLIC ${caliper_INCLUDE_DIR})
  target_link_libraries(bout++ PUBLIC caliper)
  target_compile_definitions(bout++ PUBLIC "BOUT_HAS_CALIPER")
endif ()
set(BOUT_HAS_CALIPER ${BOUT_ENABLE_CALIPER})

# UMPIRE
option(BOUT_ENABLE_UMPIRE "Enable UMPIRE memory management" OFF)
if (BOUT_ENABLE_UMPIRE)
  find_package(UMPIRE REQUIRED)
  list(APPEND BOUT_DEPENDS umpire)
  target_include_directories(bout++ PUBLIC ${UMPIRE_INCLUDE_DIRS}/include)
  target_link_libraries(bout++ PUBLIC umpire)
  target_compile_definitions(bout++ PUBLIC "BOUT_HAS_UMPIRE")
endif ()
set(BOUT_HAS_UMPIRE ${BOUT_ENABLE_UMPIRE})

# RAJA
option(BOUT_ENABLE_RAJA "Enable RAJA" OFF)
if (BOUT_ENABLE_RAJA)
  find_package(RAJA REQUIRED)
  list(APPEND BOUT_DEPENDS RAJA)
  message(STATUS "RAJA_CONFIG:" ${RAJA_CONFIG}) 
  string(FIND ${RAJA_CONFIG} "raja" loc)
  math(EXPR value "${loc} + 5" OUTPUT_FORMAT DECIMAL)
  string(SUBSTRING ${RAJA_CONFIG} 0 ${value}  RAJA_PATH)
  message(STATUS "RAJA_PATH" ${RAJA_PATH})
  target_include_directories(bout++ PUBLIC ${RAJA_PATH}/include)
  target_link_libraries(bout++ PUBLIC RAJA)
  target_compile_definitions(bout++ PUBLIC "BOUT_HAS_RAJA")
endif ()
set(BOUT_HAS_RAJA ${BOUT_ENABLE_RAJA})

# Hypre
option(BOUT_USE_HYPRE "Enable support for Hypre solvers" OFF)
set(HYPRE_DIR "" CACHE STRING "Point to HYPRE Install")
if (BOUT_USE_HYPRE)
  enable_language(C)
  find_package(HYPRE REQUIRED)
  list(APPEND BOUT_DEPENDS HYPRE)
  target_compile_definitions(bout++ PUBLIC "BOUT_HAS_HYPRE")
  target_link_libraries(bout++ PUBLIC HYPRE::HYPRE)
  if (HYPRE_CUDA)
     target_compile_definitions(bout++ PUBLIC "HYPRE_USING_CUDA;HYPRE_USING_UNIFIED_MEMORY")
     target_link_libraries(bout++ PUBLIC "${BOUT_CUDA_LIB_DIR}/libcusparse_static.a;${BOUT_CUDA_LIB_DIR}/libcurand_static.a;${BOUT_CUDA_LIB_DIR}/libculibos.a;${BOUT_CUDA_LIB_DIR}/libcublas_static.a;${BOUT_CUDA_LIB_DIR}/libcublasLt_static.a")
  endif ()
endif ()
set(BOUT_HAS_HYPRE ${BOUT_USE_HYPRE})

# PETSc
option(BOUT_USE_PETSC "Enable support for PETSc time solvers and inversions" OFF)
if (BOUT_USE_PETSC)
  if (NOT CMAKE_SYSTEM_NAME STREQUAL "CrayLinuxEnvironment")
    # Cray wrappers sort this out for us
    find_package(PETSc REQUIRED)
    list(APPEND BOUT_DEPENDS PETSc)
    target_link_libraries(bout++ PUBLIC PETSc::PETSc)
    string(JOIN " " CONFIG_PETSC_LIBRARIES ${PETSC_LIBRARIES})
    set(CONFIG_LDFLAGS "${CONFIG_PETSC_LIBRARIES} ${CONFIG_LDFLAGS}")
  endif()
endif()
message(STATUS "PETSc support: ${BOUT_USE_PETSC}")
set(BOUT_HAS_PETSC ${BOUT_USE_PETSC})


cmake_dependent_option(BOUT_USE_SYSTEM_MPARK_VARIANT "Use external installation of mpark.variant" OFF
  "BOUT_UPDATE_GIT_SUBMODULE OR EXISTS externalpackages/mpark.variant/CMakeLists.txt" ON)

if(BOUT_USE_SYSTEM_MPARK_VARIANT)
  message(STATUS "Using external mpark.variant")
  find_package(mpark_variant REQUIRED)
  get_target_property(MPARK_VARIANT_INCLUDE_PATH mpark_variant INTERFACE_INCLUDE_DIRECTORIES)
else()
  message(STATUS "Using mpark.variant submodule")
  bout_update_submodules()
  add_subdirectory(externalpackages/mpark.variant)
  if(NOT TARGET mpark_variant)
    message(FATAL_ERROR "mpark_variant not found! Have you disabled the git submodules (BOUT_UPDATE_GIT_SUBMODULE)?")
  endif()
  set(MPARK_VARIANT_INCLUDE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/externalpackages/mpark.variant/include")
  set(CONFIG_CFLAGS "${CONFIG_CFLAGS} -I\${MPARK_VARIANT_INCLUDE_PATH}")
endif()
target_link_libraries(bout++ PUBLIC mpark_variant)

cmake_dependent_option(BOUT_USE_SYSTEM_FMT "Use external installation of fmt" OFF
  "BOUT_UPDATE_GIT_SUBMODULE OR EXISTS externalpackages/fmt/CMakeLists.txt" ON)

if(BOUT_USE_SYSTEM_FMT)
  message(STATUS "Using external fmt")
  find_package(fmt 6 REQUIRED)
  get_target_property(FMT_INCLUDE_PATH fmt::fmt INTERFACE_INCLUDE_DIRECTORIES)
else()
  message(STATUS "Using fmt submodule")
  bout_update_submodules()
  # Need to install fmt alongside BOUT++
  set(FMT_INSTALL ON CACHE BOOL "")
  set(FMT_DEBUG_POSTFIX "" CACHE STRING "")
  add_subdirectory(externalpackages/fmt)
  if(NOT TARGET fmt::fmt)
    message(FATAL_ERROR "fmt not found! Have you disabled the git submodules (BOUT_UPDATE_GIT_SUBMODULE)?")
  endif()
  # Build the library in <build dir>/lib: this makes updating the path
  # for bout-config much easier
  set_target_properties(fmt PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/lib"
    ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/lib")
  set(FMT_INCLUDE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/externalpackages/fmt/include")
  set(CONFIG_CFLAGS "${CONFIG_CFLAGS} -I\${FMT_INCLUDE_PATH}")
  set(CONFIG_LDFLAGS "${CONFIG_LDFLAGS} -lfmt")
endif()
target_link_libraries(bout++ PUBLIC fmt::fmt)

option(BOUT_USE_PVODE "Enable support for bundled PVODE" ON)
if (BOUT_USE_PVODE)
  add_subdirectory(externalpackages/PVODE)
  target_link_libraries(bout++ PUBLIC pvode pvpre)
  # Build the libraries in <build dir>/lib: this makes updating the
  # path for bout-config much easier
  set_target_properties(pvode pvpre PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/lib"
    ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/lib")
  set(CONFIG_LDFLAGS "${CONFIG_LDFLAGS} -lpvode -lpvpre")
endif()
message(STATUS "PVODE support: ${BOUT_USE_PVODE}")
set(BOUT_HAS_PVODE ${BOUT_USE_PVODE})

option(BOUT_USE_NETCDF "Enable support for NetCDF output" ON)
option(BOUT_DOWNLOAD_NETCDF_CXX4 "Download and build netCDF-cxx4" OFF)
if (BOUT_USE_NETCDF)
  if (BOUT_DOWNLOAD_NETCDF_CXX4)
    include(FetchContent)
    FetchContent_Declare(
      netcdf-cxx4
      GIT_REPOSITORY https://github.com/ZedThree/netcdf-cxx4
      GIT_TAG        "ad3e50953190615cb69dcc8a4652f9a88a8499cf"
      )
    # Don't build the netcdf tests, they have lots of warnings
    set(NCXX_ENABLE_TESTS OFF CACHE BOOL "" FORCE)
    # Use our own FindnetCDF module which uses nc-config
    find_package(netCDF REQUIRED)
    FetchContent_MakeAvailable(netcdf-cxx4)
  else()
    find_package(netCDFCxx REQUIRED)
    set(CONFIG_LDFLAGS "${CONFIG_LDFLAGS} ${netCDF_CXX_LIBRARY} ${netCDF_LIBRARIES}")
  endif()
  target_link_libraries(bout++ PUBLIC netCDF::netcdf-cxx4)
endif()
message(STATUS "NetCDF support: ${BOUT_USE_NETCDF}")
set(BOUT_HAS_NETCDF ${BOUT_USE_NETCDF})

option(BOUT_USE_FFTW "Enable support for FFTW" ON)
if (BOUT_USE_FFTW)
  find_package(FFTW REQUIRED)
  target_link_libraries(bout++ PUBLIC FFTW::FFTW)
  set(CONFIG_LDFLAGS "${CONFIG_LDFLAGS} ${FFTW_LIBRARIES}")
endif()
message(STATUS "FFTW support: ${BOUT_USE_FFTW}")
set(BOUT_HAS_FFTW ${BOUT_USE_FFTW})

option(BOUT_USE_LAPACK "Enable support for LAPACK" ON)
if (BOUT_USE_LAPACK)
  if (NOT CMAKE_SYSTEM_NAME STREQUAL "CrayLinuxEnvironment")
    # Cray wrappers sort this out for us
    find_package(LAPACK REQUIRED)
    target_link_libraries(bout++ PUBLIC "${LAPACK_LIBRARIES}")
    string(JOIN " " CONFIG_LAPACK_LIBRARIES ${LAPACK_LIBRARIES})
    set(CONFIG_LDFLAGS "${CONFIG_LDFLAGS} ${CONFIG_LAPACK_LIBRARIES}")
  endif()
endif()
message(STATUS "LAPACK support: ${BOUT_USE_LAPACK}")
set(BOUT_HAS_LAPACK ${BOUT_USE_LAPACK})

option(BOUT_USE_SLEPC "Enable support for SLEPc eigen solver" OFF)
if (BOUT_USE_SLEPC)
  find_package(SLEPc REQUIRED)
  target_link_libraries(bout++ PUBLIC SLEPc::SLEPc)
  string(JOIN " " CONFIG_SLEPC_LIBRARIES ${SLEPC_LIBRARIES})
  set(CONFIG_LDFLAGS "${CONFIG_LDFLAGS} ${CONFIG_SLEPC_LIBRARIES}")
endif()
message(STATUS "SLEPc support: ${BOUT_USE_SLEPC}")
set(BOUT_HAS_SLEPC ${BOUT_USE_SLEPC})

option(BOUT_DOWNLOAD_SUNDIALS "Download and build SUNDIALS" OFF)
# Force BOUT_USE_SUNDIALS if we're downloading it!
cmake_dependent_option(BOUT_USE_SUNDIALS "Enable support for SUNDIALS time solvers" OFF
  "NOT BOUT_DOWNLOAD_SUNDIALS" ON)
if (BOUT_USE_SUNDIALS)
  if (BOUT_DOWNLOAD_SUNDIALS)
    message(STATUS "Downloading and configuring SUNDIALS")
    include(FetchContent)
    FetchContent_Declare(
      sundials
      GIT_REPOSITORY https://github.com/ZedThree/sundials
      GIT_TAG        "4f3bb8281c7b27343bcb95386ebbb665fb6196a5"
      )
    set(EXAMPLES_ENABLE_C OFF CACHE BOOL "" FORCE)
    set(EXAMPLES_INSTALL OFF CACHE BOOL "" FORCE)
    set(BOUT_ENABLE_MPI ON CACHE BOOL "" FORCE)
    set(BOUT_ENABLE_OPENMP OFF CACHE BOOL "" FORCE)
    set(BUILD_STATIC_LIBS OFF CACHE BOOL "" FORCE)
    FetchContent_MakeAvailable(sundials)
    message(STATUS "SUNDIALS done configuring")
  else()
    find_package(SUNDIALS REQUIRED)
  endif()
  target_link_libraries(bout++ PUBLIC SUNDIALS::nvecparallel)
  target_link_libraries(bout++ PUBLIC SUNDIALS::cvode)
  target_link_libraries(bout++ PUBLIC SUNDIALS::ida)
  target_link_libraries(bout++ PUBLIC SUNDIALS::arkode)
  set(CONFIG_LDFLAGS "${CONFIG_LDFLAGS} ${SUNDIALS_cvode_LIBRARY} ${SUNDIALS_ida_LIBRARY} ${SUNDIALS_arkode_LIBRARY} ${SUNDIALS_nvecparallel_LIBRARY}")
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

option(BOUT_USE_UUID_SYSTEM_GENERATOR "Enable support for using a system UUID generator" ON)
if (BOUT_USE_UUID_SYSTEM_GENERATOR)
  find_package(Libuuid QUIET)
  if (Libuuid_FOUND)
    target_link_libraries(bout++
      PUBLIC Libuuid::libuuid)
    set(CONFIG_LDFLAGS "${CONFIG_LDFLAGS} ${Libuuid_LIBRARIES}")
  else()
    message(STATUS "libuuid not found, using fallback UUID generator")
    set(BOUT_USE_UUID_SYSTEM_GENERATOR FALSE)
  endif()
endif()
message(STATUS "UUID_SYSTEM_GENERATOR: ${BOUT_USE_UUID_SYSTEM_GENERATOR}")
set(BOUT_HAS_UUID_SYSTEM_GENERATOR ${BOUT_USE_UUID_SYSTEM_GENERATOR})

message(STATUS "BOUT_DEPENDS: ${BOUT_DEPENDS}")
