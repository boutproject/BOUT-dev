cmake_minimum_required(VERSION 3.9...3.12)

if(${CMAKE_VERSION} VERSION_LESS 3.12)
  cmake_policy(VERSION ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION})
else()
  cmake_policy(VERSION 3.12)
endif()

# CMake currently doesn't support proper semver
# Set the version here, strip any extra tags to use in `project`
set(BOUT_FULL_VERSION 5.0.0-alpha)
string(REGEX REPLACE "^([0-9]+\.[0-9]+\.[0-9]+)-.*" "\\1" BOUT_CMAKE_ACCEPTABLE_VERSION ${BOUT_FULL_VERSION})
string(REGEX REPLACE "^[0-9]+\.[0-9]+\.[0-9]+-(.*)" "\\1" BOUT_VERSION_TAG ${BOUT_FULL_VERSION})

project(BOUT++
  DESCRIPTION "Fluid PDE solver framework"
  VERSION ${BOUT_CMAKE_ACCEPTABLE_VERSION}
  LANGUAGES CXX)

# This might not be entirely sensible, but helps CMake to find the
# correct MPI, workaround for https://gitlab.kitware.com/cmake/cmake/issues/18895
find_program(MPIEXEC_EXECUTABLE NAMES mpiexec mpirun)
find_package(MPI REQUIRED)

include(CMakeDependentOption)

# Override default
option(INSTALL_GTEST "Enable installation of googletest. (Projects embedding googletest may want to turn this OFF.)" OFF)

set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

option(BOUT_UPDATE_GIT_SUBMODULE "Check submodules are up-to-date during build" ON)
# Adapted from https://cliutils.gitlab.io/modern-cmake/chapters/projects/submodule.html
# Update submodules as needed
function(bout_update_submodules)
  if(NOT BOUT_UPDATE_GIT_SUBMODULE)
    return()
  endif()
  find_package(Git QUIET)
  if(GIT_FOUND AND EXISTS "${PROJECT_SOURCE_DIR}/.git")
    message(STATUS "Submodule update")
    execute_process(COMMAND ${GIT_EXECUTABLE} -c submodule.recurse=false submodule update --init --recursive
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
      RESULT_VARIABLE GIT_SUBMOD_RESULT)
    if(NOT GIT_SUBMOD_RESULT EQUAL "0")
      message(FATAL_ERROR "git submodule update --init failed with ${GIT_SUBMOD_RESULT}, please checkout submodules")
    endif()
  endif()
endfunction()

set(BOUT_SOURCES
  ./include/boundary_factory.hxx
  ./include/boundary_op.hxx
  ./include/boundary_region.hxx
  ./include/boundary_standard.hxx
  ./include/bout.hxx
  ./include/bout/array.hxx
  ./include/bout/assert.hxx
  ./include/bout/constants.hxx
  ./include/bout/coordinates.hxx
  ./include/bout/deprecated.hxx
  ./include/bout/deriv_store.hxx
  ./include/bout/expr.hxx
  ./include/bout/field_visitor.hxx
  ./include/bout/fieldgroup.hxx
  ./include/bout/format.hxx
  ./include/bout/fv_ops.hxx
  ./include/bout/generic_factory.hxx
  ./include/bout/globalfield.hxx
  ./include/bout/griddata.hxx
  ./include/bout/index_derivs.hxx
  ./include/bout/index_derivs_interface.hxx
  ./include/bout/invert/laplacexy.hxx
  ./include/bout/invert/laplacexz.hxx
  ./include/bout/invertable_operator.hxx
  ./include/bout/macro_for_each.hxx
  ./include/bout/mesh.hxx
  ./include/bout/monitor.hxx
  ./include/bout/mpi_wrapper.hxx
  ./include/bout/openmpwrap.hxx
  ./include/bout/paralleltransform.hxx
  ./include/bout/petsclib.hxx
  ./include/bout/petsc_interface.hxx
  ./include/bout/physicsmodel.hxx
  ./include/bout/region.hxx
  ./include/bout/rkscheme.hxx
  ./include/bout/rvec.hxx
  ./include/bout/scorepwrapper.hxx
  ./include/bout/slepclib.hxx
  ./include/bout/snb.hxx
  ./include/bout/solver.hxx
  ./include/bout/surfaceiter.hxx
  ./include/bout/operatorstencil.hxx
  ./include/bout/sys/expressionparser.hxx
  ./include/bout/sys/generator_context.hxx
  ./include/bout/sys/gettext.hxx
  ./include/bout/sys/range.hxx
  ./include/bout/sys/timer.hxx
  ./include/bout/sys/type_name.hxx
  ./include/bout/sys/uncopyable.hxx
  ./include/bout/sys/variant.hxx
  ./include/bout/template_combinations.hxx
  ./include/bout_types.hxx
  ./include/boutcomm.hxx
  ./include/boutexception.hxx
  ./include/boutmain.hxx
  ./include/cyclic_reduction.hxx
  ./include/datafile.hxx
  ./include/dataformat.hxx
  ./include/dcomplex.hxx
  ./include/derivs.hxx
  ./include/difops.hxx
  ./include/fft.hxx
  ./include/field.hxx
  ./include/field2d.hxx
  ./include/field3d.hxx
  ./include/field_data.hxx
  ./include/field_factory.hxx
  ./include/fieldperp.hxx
  ./include/globals.hxx
  ./include/gyro_average.hxx
  ./include/initialprofiles.hxx
  ./include/interpolation.hxx
  ./include/interpolation_xz.hxx
  ./include/invert_laplace.hxx
  ./include/invert_parderiv.hxx
  ./include/lapack_routines.hxx
  ./include/mask.hxx
  ./include/msg_stack.hxx
  ./include/multiostream.hxx
  ./include/options.hxx
  ./include/options_netcdf.hxx
  ./include/optionsreader.hxx
  ./include/output.hxx
  ./include/parallel_boundary_op.hxx
  ./include/parallel_boundary_region.hxx
  ./include/smoothing.hxx
  ./include/sourcex.hxx
  ./include/stencils.hxx
  ./include/unused.hxx
  ./include/utils.hxx
  ./include/vecops.hxx
  ./include/vector2d.hxx
  ./include/vector3d.hxx
  ./include/where.hxx
  ./src/bout++.cxx
  ./src/field/field.cxx
  ./src/field/field2d.cxx
  ./src/field/field3d.cxx
  ./src/field/field_data.cxx
  ./src/field/field_factory.cxx
  ./src/field/fieldgenerators.cxx
  ./src/field/fieldgenerators.hxx
  ./src/field/fieldgroup.cxx
  ./src/field/fieldperp.cxx
  ./src/field/generated_fieldops.cxx
  ./src/field/globalfield.cxx
  ./src/field/initialprofiles.cxx
  ./src/field/vecops.cxx
  ./src/field/vector2d.cxx
  ./src/field/vector3d.cxx
  ./src/field/where.cxx
  ./src/fileio/datafile.cxx
  ./src/fileio/dataformat.cxx
  ./src/fileio/formatfactory.cxx
  ./src/fileio/formatfactory.hxx
  ./src/fileio/impls/emptyformat.hxx
  ./src/fileio/impls/hdf5/h5_format.cxx
  ./src/fileio/impls/hdf5/h5_format.hxx
  ./src/fileio/impls/netcdf/nc_format.cxx
  ./src/fileio/impls/netcdf/nc_format.hxx
  ./src/fileio/impls/netcdf4/ncxx4.cxx
  ./src/fileio/impls/netcdf4/ncxx4.hxx
  ./src/fileio/impls/pnetcdf/pnetcdf.cxx
  ./src/fileio/impls/pnetcdf/pnetcdf.hxx
  ./src/invert/fft_fftw.cxx
  ./src/invert/lapack_routines.cxx
  ./src/invert/laplace/impls/cyclic/cyclic_laplace.cxx
  ./src/invert/laplace/impls/cyclic/cyclic_laplace.hxx
  ./src/invert/laplace/impls/multigrid/multigrid_alg.cxx
  ./src/invert/laplace/impls/multigrid/multigrid_laplace.cxx
  ./src/invert/laplace/impls/multigrid/multigrid_laplace.hxx
  ./src/invert/laplace/impls/multigrid/multigrid_solver.cxx
  ./src/invert/laplace/impls/naulin/naulin_laplace.cxx
  ./src/invert/laplace/impls/naulin/naulin_laplace.hxx
  ./src/invert/laplace/impls/pdd/pdd.cxx
  ./src/invert/laplace/impls/pdd/pdd.hxx
  ./src/invert/laplace/impls/petsc/petsc_laplace.cxx
  ./src/invert/laplace/impls/petsc/petsc_laplace.hxx
  ./src/invert/laplace/impls/petsc3damg/petsc3damg.cxx
  ./src/invert/laplace/impls/petsc3damg/petsc3damg.hxx
 ./src/invert/laplace/impls/serial_band/serial_band.cxx
  ./src/invert/laplace/impls/serial_band/serial_band.hxx
  ./src/invert/laplace/impls/serial_tri/serial_tri.cxx
  ./src/invert/laplace/impls/serial_tri/serial_tri.hxx
  ./src/invert/laplace/impls/shoot/shoot_laplace.cxx
  ./src/invert/laplace/impls/shoot/shoot_laplace.hxx
  ./src/invert/laplace/impls/spt/spt.cxx
  ./src/invert/laplace/impls/spt/spt.hxx
  ./src/invert/laplace/invert_laplace.cxx
  ./src/invert/laplacexy/laplacexy.cxx
  ./src/invert/laplacexz/impls/cyclic/laplacexz-cyclic.cxx
  ./src/invert/laplacexz/impls/cyclic/laplacexz-cyclic.hxx
  ./src/invert/laplacexz/impls/petsc/laplacexz-petsc.cxx
  ./src/invert/laplacexz/impls/petsc/laplacexz-petsc.hxx
  ./src/invert/laplacexz/laplacexz.cxx
  ./src/invert/parderiv/impls/cyclic/cyclic.cxx
  ./src/invert/parderiv/impls/cyclic/cyclic.hxx
  ./src/invert/parderiv/invert_parderiv.cxx
  ./src/mesh/boundary_factory.cxx
  ./src/mesh/boundary_region.cxx
  ./src/mesh/boundary_standard.cxx
  ./src/mesh/coordinates.cxx
  ./src/mesh/data/gridfromfile.cxx
  ./src/mesh/data/gridfromoptions.cxx
  ./src/mesh/difops.cxx
  ./src/mesh/fv_ops.cxx
  ./src/mesh/impls/bout/boutmesh.cxx
  ./src/mesh/impls/bout/boutmesh.hxx
  ./src/mesh/index_derivs.cxx
  ./src/mesh/interpolation_xz.cxx
  ./src/mesh/interpolation/bilinear_xz.cxx
  ./src/mesh/interpolation/hermite_spline_xz.cxx
  ./src/mesh/interpolation/hermite_spline_z.cxx
  ./src/mesh/interpolation/interpolation_z.cxx
  ./src/mesh/interpolation/lagrange_4pt_xz.cxx
  ./src/mesh/interpolation/monotonic_hermite_spline_xz.cxx
  ./src/mesh/mesh.cxx
  ./src/mesh/parallel/fci.cxx
  ./src/mesh/parallel/fci.hxx
  ./src/mesh/parallel/identity.cxx
  ./src/mesh/parallel/shiftedmetric.cxx
  ./src/mesh/parallel/shiftedmetricinterp.cxx
  ./src/mesh/parallel/shiftedmetricinterp.hxx
  ./src/mesh/parallel_boundary_op.cxx
  ./src/mesh/parallel_boundary_region.cxx
  ./src/mesh/surfaceiter.cxx
  ./src/physics/gyro_average.cxx
  ./src/physics/physicsmodel.cxx
  ./src/physics/smoothing.cxx
  ./src/physics/snb.cxx
  ./src/physics/sourcex.cxx
  ./src/solver/impls/adams_bashforth/adams_bashforth.cxx
  ./src/solver/impls/adams_bashforth/adams_bashforth.hxx
  ./src/solver/impls/arkode/arkode.cxx
  ./src/solver/impls/arkode/arkode.hxx
  ./src/solver/impls/cvode/cvode.cxx
  ./src/solver/impls/cvode/cvode.hxx
  ./src/solver/impls/euler/euler.cxx
  ./src/solver/impls/euler/euler.hxx
  ./src/solver/impls/ida/ida.cxx
  ./src/solver/impls/ida/ida.hxx
  ./src/solver/impls/imex-bdf2/imex-bdf2.cxx
  ./src/solver/impls/imex-bdf2/imex-bdf2.hxx
  ./src/solver/impls/karniadakis/karniadakis.cxx
  ./src/solver/impls/karniadakis/karniadakis.hxx
  ./src/solver/impls/petsc/petsc.cxx
  ./src/solver/impls/petsc/petsc.hxx
  ./src/solver/impls/power/power.cxx
  ./src/solver/impls/power/power.hxx
  ./src/solver/impls/pvode/pvode.cxx
  ./src/solver/impls/pvode/pvode.hxx
  ./src/solver/impls/rk3-ssp/rk3-ssp.cxx
  ./src/solver/impls/rk3-ssp/rk3-ssp.hxx
  ./src/solver/impls/rk4/rk4.cxx
  ./src/solver/impls/rk4/rk4.hxx
  ./src/solver/impls/rkgeneric/impls/cashkarp/cashkarp.cxx
  ./src/solver/impls/rkgeneric/impls/cashkarp/cashkarp.hxx
  ./src/solver/impls/rkgeneric/impls/rk4simple/rk4simple.cxx
  ./src/solver/impls/rkgeneric/impls/rk4simple/rk4simple.hxx
  ./src/solver/impls/rkgeneric/impls/rkf34/rkf34.cxx
  ./src/solver/impls/rkgeneric/impls/rkf34/rkf34.hxx
  ./src/solver/impls/rkgeneric/impls/rkf45/rkf45.cxx
  ./src/solver/impls/rkgeneric/impls/rkf45/rkf45.hxx
  ./src/solver/impls/rkgeneric/rkgeneric.cxx
  ./src/solver/impls/rkgeneric/rkgeneric.hxx
  ./src/solver/impls/rkgeneric/rkscheme.cxx
  ./src/solver/impls/slepc/slepc.cxx
  ./src/solver/impls/slepc/slepc.hxx
  ./src/solver/impls/snes/snes.cxx
  ./src/solver/impls/snes/snes.hxx
  ./src/solver/impls/split-rk/split-rk.cxx
  ./src/solver/impls/split-rk/split-rk.hxx
  ./src/solver/solver.cxx
  ./src/sys/bout_types.cxx
  ./src/sys/boutcomm.cxx
  ./src/sys/boutexception.cxx
  ./src/sys/derivs.cxx
  ./src/sys/expressionparser.cxx
  ./src/sys/generator_context.cxx
  ./src/sys/msg_stack.cxx
  ./src/sys/options.cxx
  ./src/sys/options/optionparser.hxx
  ./src/sys/options/options_ini.cxx
  ./src/sys/options/options_ini.hxx
  ./src/sys/options/options_netcdf.cxx
  ./src/sys/optionsreader.cxx
  ./src/sys/output.cxx
  ./src/sys/petsclib.cxx
  ./src/sys/range.cxx
  ./src/sys/slepclib.cxx
  ./src/sys/timer.cxx
  ./src/sys/type_name.cxx
  ./src/sys/utils.cxx
  ${CMAKE_CURRENT_BINARY_DIR}/include/bout/revision.hxx
  ${CMAKE_CURRENT_BINARY_DIR}/include/bout/version.hxx
  )

add_library(bout++
  ${BOUT_SOURCES}
  )
add_library(bout++::bout++ ALIAS bout++)
target_link_libraries(bout++ PUBLIC MPI::MPI_CXX)
target_include_directories(bout++ PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/include>
  $<INSTALL_INTERFACE:include>
  )

target_compile_features(bout++ PUBLIC cxx_std_14)
set_target_properties(bout++ PROPERTIES CXX_EXTENSIONS OFF)

##################################################
# Components of the version number
# Pre-release identifier (BOUT_VERSION_TAG) set above
set(BOUT_VERSION ${BOUT_FULL_VERSION})
set(BOUT_VERSION_MAJOR ${PROJECT_VERSION_MAJOR})
set(BOUT_VERSION_MINOR ${PROJECT_VERSION_MINOR})
set(BOUT_VERSION_PATCH ${PROJECT_VERSION_PATCH})

# Get the git commit
include(GetGitRevisionDescription)
get_git_head_revision(GIT_REFSPEC BOUT_REVISION)
if(BOUT_REVISION STREQUAL "GITDIR-NOTFOUND")
  set(BOUT_REVISION "Unknown")
endif()
message(STATUS "Git revision: ${BOUT_REVISION}")

# Build the file containing the version information
configure_file(
  "${PROJECT_SOURCE_DIR}/include/bout/version.hxx.in"
  "${PROJECT_BINARY_DIR}/include/bout/version.hxx")
# Build the file containing just the commit hash
# This will be rebuilt on every commit!
configure_file(
  "${PROJECT_SOURCE_DIR}/include/bout/revision.hxx.in"
  "${PROJECT_BINARY_DIR}/include/bout/revision.hxx")

##################################################

cmake_dependent_option(BOUT_USE_SYSTEM_MPARK_VARIANT "Use external installation of mpark.variant" OFF
  "BOUT_UPDATE_GIT_SUBMODULE OR EXISTS externalpackages/mpark.variant/CMakeLists.txt" ON)

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
  "BOUT_UPDATE_GIT_SUBMODULE OR EXISTS externalpackages/fmt/CMakeLists.txt" ON)

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

option(BOUT_ENABLE_WARNINGS "Enable compiler warnings" ON)
if (BOUT_ENABLE_WARNINGS)
  target_compile_options(bout++ PRIVATE
    $<$<OR:$<CXX_COMPILER_ID:GNU>,$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:AppleClang>>:
      -Wall -Wextra >
    $<$<CXX_COMPILER_ID:MSVC>:
      /W4 >
   )

 include(EnableCXXWarningIfSupport)
 # Note we explicitly turn off -Wcast-function-type as PETSc *requires*
 # we cast a function to the wrong type in MatFDColoringSetFunction
 target_enable_cxx_warning_if_supported(bout++
   FLAGS -Wnull-dereference -Wno-cast-function-type
   )

endif()

# Compile time features

set(CHECK_LEVELS 0 1 2 3 4)
set(CHECK 2 CACHE STRING "Set run-time checking level")
set_property(CACHE CHECK PROPERTY STRINGS ${CHECK_LEVELS})
if (NOT CHECK IN_LIST CHECK_LEVELS)
  message(FATAL_ERROR "CHECK must be one of ${CHECK_LEVELS}")
endif()
message(STATUS "Runtime checking level: CHECK=${CHECK}")
target_compile_definitions(bout++ PUBLIC "CHECK=${CHECK}")
set(BOUT_CHECK_LEVEL ${CHECK})

cmake_dependent_option(BOUT_ENABLE_OUTPUT_DEBUG "Enable extra debug output" OFF
  "CHECK LESS 3" ON)
message(STATUS "Extra debug output: BOUT_USE_OUTPUT_DEBUG=${BOUT_ENABLE_OUTPUT_DEBUG}")
set(BOUT_USE_OUTPUT_DEBUG ${BOUT_ENABLE_OUTPUT_DEBUG})

option(BOUT_ENABLE_SIGNAL "SegFault handling" ON)
message(STATUS "Signal handling: BOUT_USE_SIGNAL=${BOUT_ENABLE_SIGNAL}")
set(BOUT_USE_SIGNAL ${BOUT_ENABLE_SIGNAL})

option(BOUT_ENABLE_COLOR "Output coloring" ON)
message(STATUS "Output coloring: BOUT_USE_COLOR=${BOUT_ENABLE_COLOR}")
set(BOUT_USE_COLOR ${BOUT_ENABLE_COLOR})

option(BOUT_ENABLE_TRACK "Field name tracking" ON)
message(STATUS "Field name tracking: BOUT_USE_TRACK=${BOUT_ENABLE_TRACK}")
set(BOUT_USE_TRACK ${BOUT_ENABLE_TRACK})

option(BOUT_ENABLE_SIGFPE "Signalling floating point exceptions" OFF)
message(STATUS "Signalling floating point exceptions: BOUT_USE_SIGFPE=${BOUT_ENABLE_SIGFPE}")
set(BOUT_USE_SIGFPE ${BOUT_ENABLE_SIGFPE})

option(BOUT_ENABLE_BACKTRACE "Enable backtrace" ON)
if (BOUT_ENABLE_BACKTRACE)
  find_program(ADDR2LINE_FOUND addr2line)
  if (NOT ADDR2LINE_FOUND)
    message(FATAL_ERROR "addr2line not found")
  endif()
  target_link_libraries(bout++ PUBLIC ${CMAKE_DL_LIBS})
endif()
message(STATUS "Enable backtrace: BOUT_USE_BACKTRACE=${BOUT_ENABLE_BACKTRACE}")
set(BOUT_USE_BACKTRACE ${BOUT_ENABLE_BACKTRACE})

option(BOUT_ENABLE_OPENMP "Enable OpenMP support" OFF)
set(BOUT_OPENMP_SCHEDULE static CACHE STRING "Set OpenMP schedule")
set_property(CACHE BOUT_OPENMP_SCHEDULE PROPERTY STRINGS static dynamic guided auto)
if (BOUT_ENABLE_OPENMP)
  find_package(OpenMP REQUIRED)
  target_link_libraries(bout++ PUBLIC OpenMP::OpenMP_CXX)
  set(possible_openmp_schedules static dynamic guided auto)
  if (NOT BOUT_OPENMP_SCHEDULE IN_LIST possible_openmp_schedules)
    message(FATAL_ERROR "BOUT_OPENMP_SCHEDULE must be one of ${possible_openmp_schedules}; got ${BOUT_OPENMP_SCHEDULE}")
  endif()
  message(STATUS "OpenMP schedule: ${BOUT_OPENMP_SCHEDULE}")
endif()
message(STATUS "Enable OpenMP: ${BOUT_ENABLE_OPENMP}")
set(BOUT_USE_OPENMP ${BOUT_ENABLE_OPENMP})

# Optional dependencies

option(BOUT_USE_PVODE "Enable support for bundled PVODE" ON)
if (BOUT_USE_PVODE)
  add_subdirectory(externalpackages/PVODE)
  target_link_libraries(bout++ PUBLIC pvode pvpre)
endif()
message(STATUS "PVODE support: ${BOUT_USE_PVODE}")
set(BOUT_HAS_PVODE ${BOUT_USE_PVODE})

option(BOUT_USE_NETCDF "Enable support for NetCDF output" ON)
if (BOUT_USE_NETCDF)
  find_package(NetCDF REQUIRED)
  target_link_libraries(bout++ PUBLIC NetCDF::NetCDF_CXX)
endif()
message(STATUS "NetCDF support: ${BOUT_USE_NETCDF}")
set(BOUT_HAS_NETCDF ${BOUT_USE_NETCDF})

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

option(BOUT_USE_PETSC "Enable support for PETSc time solvers and inversions" OFF)
if (BOUT_USE_PETSC)
  if (NOT CMAKE_SYSTEM_NAME STREQUAL "CrayLinuxEnvironment")
    # Cray wrappers sort this out for us
    find_package(PETSc REQUIRED)
    target_link_libraries(bout++ PUBLIC PETSc::PETSc)
  endif()
endif()
message(STATUS "PETSc support: ${BOUT_USE_PETSC}")
set(BOUT_HAS_PETSC ${BOUT_USE_PETSC})

option(BOUT_USE_SLEPC "Enable support for SLEPc eigen solver" OFF)
if (BOUT_USE_SLEPC)
  find_package(SLEPc REQUIRED)
  target_link_libraries(bout++ PUBLIC SLEPc::SLEPc)
endif()
message(STATUS "SLEPc support: ${BOUT_USE_SLEPC}")
set(BOUT_HAS_SLEPC ${BOUT_USE_SLEPC})

option(BOUT_USE_SUNDIALS "Enable support for SUNDIALS time solvers" OFF)
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

include(CheckCXXSourceCompiles)
check_cxx_source_compiles("int main() { const char* name = __PRETTY_FUNCTION__; }"
  HAS_PRETTY_FUNCTION)
set(BOUT_HAS_PRETTY_FUNCTION ${HAS_PRETTY_FUNCTION})

# We want to compile the actual flags used into the library so we can
# see them at runtime. This needs a few steps:

# 1. Get the macro definitions. They come as a ;-separated list and
#    without the -D. We also need to also stick a -D on the front of
#    the first item
get_property(BOUT_COMPILE_DEFINITIONS
  TARGET bout++
  PROPERTY COMPILE_DEFINITIONS)
string(REPLACE ";" " -D" BOUT_COMPILE_DEFINITIONS "${BOUT_COMPILE_DEFINITIONS}")
string(CONCAT BOUT_COMPILE_DEFINITIONS " -D" "${BOUT_COMPILE_DEFINITIONS}")

# 2. Get the compiler options. Again, they come as a ;-separated
#    list. Note that they don't include optimisation or debug flags:
#    they're in the CMAKE_CXX_FLAGS* variables
get_property(BOUT_COMPILE_OPTIONS
  TARGET bout++
  PROPERTY COMPILE_OPTIONS)
string(REPLACE ";" " " BOUT_COMPILE_OPTIONS "${BOUT_COMPILE_OPTIONS}")

# 3. The optimisation and/or debug flags are in the CMAKE_CXX_FLAGS*
#    variables. We need both the common flags as well as those for the
#    build type actually being used. Note: this might behave weirdly
#    on Windows. Might need to expand CMAKE_CONFIGURATION_TYPES
#    instead?
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE_UPPER "DEBUG")
else()
  string(TOUPPER "${CMAKE_BUILD_TYPE}" CMAKE_BUILD_TYPE_UPPER)
endif()
string(CONCAT BOUT_COMPILE_BUILD_FLAGS
  " "
  "${CMAKE_CXX_FLAGS}"
  "${CMAKE_CXX_FLAGS_${CMAKE_BUILD_TYPE_UPPER}}")

# 4. Now we join all the flags from the first three steps together
string(CONCAT BOUT_FLAGS_STRING
  "${BOUT_COMPILE_OPTIONS}"
  "${BOUT_COMPILE_DEFINITIONS}"
  "${BOUT_COMPILE_BUILD_FLAGS}")

# 5. Finally actually add the flags as a define
target_compile_definitions(bout++
  PRIVATE BOUT_FLAGS_STRING=${BOUT_FLAGS_STRING})

##################################################
# Tests

# Copy FILENAME from source directory to build directory
macro(bout_test_copy_file FILENAME)
  configure_file(
      ${CMAKE_CURRENT_SOURCE_DIR}/${FILENAME}
      ${CMAKE_CURRENT_BINARY_DIR}/${FILENAME}
      COPYONLY)
endmacro()

# Add a new integrated test. By default, the executable is named like
# the first source, stripped of its file extension.
#
# Required arguments:
#
# - TESTNAME: name of the test
#
# - SOURCES: list of source files
#
# Optional arguments:
#
# - USE_RUNTEST: if given, the test uses `./runtest` as the test
#   command, otherwise it uses the executable
#
# - USE_DATA_BOUT_INP: if given, copy `data/BOUT.inp`
#
# - EXTRA_FILES: any extra files that are required to run the test
#
# - REQUIRES: list of variables that must be truthy to enable test
#   (note: use `CONFLICTS` to negate the variable, rather than `NOT
#   VARIABLE`)
#
# - CONFLICTS: list of variables that must be falsey to enable test
#
# - EXECUTABLE_NAME: name of the executable, if different from the
#   first source name

function(bout_add_integrated_test TESTNAME)
  set(options USE_RUNTEST USE_DATA_BOUT_INP)
  set(oneValueArgs EXECUTABLE_NAME)
  set(multiValueArgs SOURCES EXTRA_FILES REQUIRES CONFLICTS)
  cmake_parse_arguments(BOUT_TEST_OPTIONS "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

  foreach (REQUIREMENT IN LISTS BOUT_TEST_OPTIONS_REQUIRES)
    if (NOT ${REQUIREMENT})
      message(STATUS "Not building test ${TESTNAME}, requirement not met: ${REQUIREMENT}")
      return()
    endif()
  endforeach()

  foreach (CONFLICT IN LISTS BOUT_TEST_OPTIONS_CONFLICTS)
    if (${CONFLICT})
      message(STATUS "Not building test ${TESTNAME}, conflicts with: ${CONFLICT}")
      return()
    endif()
  endforeach()

  add_executable(${TESTNAME} ${BOUT_TEST_OPTIONS_SOURCES})
  target_link_libraries(${TESTNAME} bout++)
  target_include_directories(${TESTNAME} PRIVATE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)

  # Set the name of the executable. We either take it as an option,
  # or use the first source file, stripping the file suffix
  if (BOUT_TEST_OPTIONS_EXECUTABLE_NAME)
    set_target_properties(${TESTNAME} PROPERTIES OUTPUT_NAME ${BOUT_TEST_OPTIONS_EXECUTABLE_NAME})
  else()
    # If more than one source file, just get the first one
    list(LENGTH ${BOUT_TEST_OPTIONS_SOURCES} BOUT_SOURCES_LENGTH)
    if (BOUT_SOURCES_LENGTH GREATER 0)
      list(GET ${BOUT_TEST_OPTIONS_SOURCES} 0 BOUT_TEST_FIRST_SOURCE)
    else()
      set(BOUT_TEST_FIRST_SOURCE ${BOUT_TEST_OPTIONS_SOURCES})
    endif()
    # Strip the directory and file extension from the source file
    get_filename_component(BOUT_TEST_EXECUTABLE_NAME ${BOUT_TEST_FIRST_SOURCE} NAME_WE)
    set_target_properties(${TESTNAME} PROPERTIES OUTPUT_NAME ${BOUT_TEST_EXECUTABLE_NAME})
  endif()

  # Set the actual test command
  if (BOUT_TEST_OPTIONS_USE_RUNTEST)
    add_test(NAME ${TESTNAME} COMMAND ./runtest)
    bout_test_copy_file(runtest)
  else()
    add_test(NAME ${TESTNAME} COMMAND ${TESTNAME})
  endif()

  # Copy the input file if needed
  if (BOUT_TEST_OPTIONS_USE_DATA_BOUT_INP)
    bout_test_copy_file(data/BOUT.inp)
  endif()

  # Copy any other needed files
  if (BOUT_TEST_OPTIONS_EXTRA_FILES)
    foreach (FILE ${BOUT_TEST_OPTIONS_EXTRA_FILES})
      bout_test_copy_file("${FILE}")
    endforeach()
  endif()

  set_target_properties(${TESTNAME} PROPERTIES FOLDER tests/integrated)

  # Add the test to the build-check-integrated-tests target
  add_dependencies(build-check-integrated-tests ${TESTNAME})
endfunction()

option(PACKAGE_TESTS "Build the tests" ON)
if(PACKAGE_TESTS)
  enable_testing()

  # Targets for just building the tests
  # Tests need to add themselves as dependencies to these targets
  add_custom_target(build-check-unit-tests)
  add_custom_target(build-check-integrated-tests)

  # Build all the tests
  add_custom_target(build-check)
  add_dependencies(build-check build-check-unit-tests build-check-integrated-tests)

  add_subdirectory(tests/unit EXCLUDE_FROM_ALL)
  add_subdirectory(tests/integrated EXCLUDE_FROM_ALL)

  # Targets for running the tests
  add_custom_target(check-unit-tests
    COMMAND ctest -R serial_tests --output-on-failure)
  add_dependencies(check-unit-tests build-check-unit-tests)

  add_custom_target(check-integrated-tests
    COMMAND ctest -R "test-" --output-on-failure)
  add_dependencies(check-integrated-tests build-check-integrated-tests)

  # Run all the tests
  add_custom_target(check)
  add_dependencies(check check-unit-tests check-integrated-tests)
endif()

option(BOUT_ENABLE_METRIC_3D "Enable 3D metric support" OFF)
if(BOUT_ENABLE_METRIC_3D)
  set(BOUT_METRIC_TYPE "3D")
else()
  set(BOUT_METRIC_TYPE "2D")
endif()
set(BOUT_USE_METRIC_3D ${BOUT_ENABLE_METRIC_3D})
##################################################
# Generate the build config header

if (EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/include/bout/build_defines.hxx")
  message(FATAL_ERROR "Generated build_defines.hxx header already exists; please remove '${CMAKE_CURRENT_SOURCE_DIR}/include/bout/build_defines.hxx' before continuing")
endif()

configure_file(cmake_build_defines.hxx.in include/bout/build_defines.hxx)

##################################################
# Installation

include(GNUInstallDirs)
install(TARGETS bout++
  EXPORT bout++Targets
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
  ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}"
  RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}"
  INCLUDES DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
  )
# Repo files
install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
  FILES_MATCHING PATTERN "*.hxx")
# Generated headers
install(DIRECTORY "${PROJECT_BINARY_DIR}/include/"
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

include(CMakePackageConfigHelpers)
write_basic_package_version_file(
  bout++ConfigVersion.cmake
  VERSION ${PACKAGE_VERSION}
  COMPATIBILITY SameMajorVersion
  )

install(EXPORT bout++Targets
  FILE bout++Targets.cmake
  NAMESPACE bout++::
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/bout++"
  )

configure_package_config_file(bout++Config.cmake.in bout++Config.cmake
  INSTALL_DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/bout++Config.cmake"
  )

configure_package_config_file(tools/pylib/boutconfig/__init__.py.cin tools/pylib/boutconfig/__init__.py
  INSTALL_DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/tools/pylib/boutconfig/__init__.py"
  )
install(
  FILES
      "${CMAKE_CURRENT_BINARY_DIR}/bout++Config.cmake"
      "${CMAKE_CURRENT_BINARY_DIR}/bout++ConfigVersion.cmake"
      "${CMAKE_CURRENT_SOURCE_DIR}/cmake/CorrectWindowsPaths.cmake"
      "${CMAKE_CURRENT_SOURCE_DIR}/cmake/FindFFTW.cmake"
      "${CMAKE_CURRENT_SOURCE_DIR}/cmake/FindNetCDF.cmake"
      "${CMAKE_CURRENT_SOURCE_DIR}/cmake/FindPackageMultipass.cmake"
      "${CMAKE_CURRENT_SOURCE_DIR}/cmake/FindPETSc.cmake"
      "${CMAKE_CURRENT_SOURCE_DIR}/cmake/FindScoreP.cmake"
      "${CMAKE_CURRENT_SOURCE_DIR}/cmake/FindSLEPc.cmake"
      "${CMAKE_CURRENT_SOURCE_DIR}/cmake/FindSUNDIALS.cmake"
      "${CMAKE_CURRENT_SOURCE_DIR}/cmake/ResolveCompilerPaths.cmake"
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/bout++"
  )

export(EXPORT bout++Targets
  FILE "${CMAKE_CURRENT_BINARY_DIR}/bout++Targets.cmake"
  NAMESPACE bout++::
  )

export(PACKAGE bout)

##################################################
# Configure summary

message("
   --------------------------------
     BOUT++ Configuration Summary
   --------------------------------

   Bundled PVODE support    : ${BOUT_HAS_PVODE}
   PETSc support            : ${BOUT_HAS_PETSC}
   SLEPc support            : ${BOUT_HAS_SLEPC}
   SUNDIALS support         : ${BOUT_HAS_SUNDIALS}
   NetCDF support           : ${BOUT_HAS_NETCDF}
   HDF5 support             : ${BOUT_HAS_HDF5}
   FFTW support             : ${BOUT_HAS_FFTW}
   LAPACK support           : ${BOUT_HAS_LAPACK}
   OpenMP support           : ${BOUT_USE_OPENMP}
   Natural language support : ${BOUT_HAS_GETTEXT}
   ScoreP support           : ${BOUT_HAS_SCOREP}
   Extra debug output       : ${BOUT_USE_OUTPUT_DEBUG}
   Runtime error check level: ${BOUT_CHECK_LEVEL}
   Signal handling          : ${BOUT_USE_SIGNAL}
   Output coloring          : ${BOUT_USE_COLOR}
   Field name tracking      : ${BOUT_USE_TRACK}
   Floating point exceptions: ${BOUT_USE_SIGFPE}
   Backtrace enabled        : ${BOUT_USE_BACKTRACE}
   Metric type              : ${BOUT_METRIC_TYPE}

   === Python ===

   Make sure that the tools/pylib directory is in your PYTHONPATH
   e.g. by adding to your ~/.bashrc file

       export PYTHONPATH=$PWD/tools/pylib/:\$PYTHONPATH

*** Now run `cmake --build .` to compile BOUT++ ***
")
