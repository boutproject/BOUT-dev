cmake_minimum_required(VERSION 3.1...3.12)

if(${CMAKE_VERSION} VERSION_LESS 3.11)
  cmake_policy(VERSION ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION})
else()
  cmake_policy(VERSION 3.12)
endif()

project(BOUT++
  DESCRIPTION "Fluid PDE solver framework"
  VERSION 4.2.2
  LANGUAGES CXX)

# This might not be entirely sensible, but helps CMake to find the
# correct MPI, workaround for https://gitlab.kitware.com/cmake/cmake/issues/18895
find_program(MPIEXEC_EXECUTABLE NAMES mpiexec mpirun)
find_package(MPI REQUIRED)

set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

add_subdirectory(externalpackages/mpark.variant)
add_subdirectory(externalpackages/PVODE)

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
  ./include/bout/openmpwrap.hxx
  ./include/bout/paralleltransform.hxx
  ./include/bout/petsclib.hxx
  ./include/bout/physicsmodel.hxx
  ./include/bout/region.hxx
  ./include/bout/rkscheme.hxx
  ./include/bout/rvec.hxx
  ./include/bout/scorepwrapper.hxx
  ./include/bout/slepclib.hxx
  ./include/bout/solver.hxx
  ./include/bout/solverfactory.hxx
  ./include/bout/surfaceiter.hxx
  ./include/bout/sys/expressionparser.hxx
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
  ./include/interpolation_factory.hxx
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
  ./src/invert/laplace/impls/mumps/mumps_laplace.cxx
  ./src/invert/laplace/impls/mumps/mumps_laplace.hxx
  ./src/invert/laplace/impls/naulin/naulin_laplace.cxx
  ./src/invert/laplace/impls/naulin/naulin_laplace.hxx
  ./src/invert/laplace/impls/pdd/pdd.cxx
  ./src/invert/laplace/impls/pdd/pdd.hxx
  ./src/invert/laplace/impls/petsc/petsc_laplace.cxx
  ./src/invert/laplace/impls/petsc/petsc_laplace.hxx
  ./src/invert/laplace/impls/serial_band/serial_band.cxx
  ./src/invert/laplace/impls/serial_band/serial_band.hxx
  ./src/invert/laplace/impls/serial_tri/serial_tri.cxx
  ./src/invert/laplace/impls/serial_tri/serial_tri.hxx
  ./src/invert/laplace/impls/shoot/shoot_laplace.cxx
  ./src/invert/laplace/impls/shoot/shoot_laplace.hxx
  ./src/invert/laplace/impls/spt/spt.cxx
  ./src/invert/laplace/impls/spt/spt.hxx
  ./src/invert/laplace/invert_laplace.cxx
  ./src/invert/laplace/laplacefactory.cxx
  ./src/invert/laplace/laplacefactory.hxx
  ./src/invert/laplacexy/laplacexy.cxx
  ./src/invert/laplacexz/impls/cyclic/laplacexz-cyclic.cxx
  ./src/invert/laplacexz/impls/cyclic/laplacexz-cyclic.hxx
  ./src/invert/laplacexz/impls/petsc/laplacexz-petsc.cxx
  ./src/invert/laplacexz/impls/petsc/laplacexz-petsc.hxx
  ./src/invert/laplacexz/laplacexz.cxx
  ./src/invert/parderiv/impls/cyclic/cyclic.cxx
  ./src/invert/parderiv/impls/cyclic/cyclic.hxx
  ./src/invert/parderiv/invert_parderiv.cxx
  ./src/invert/parderiv/parderiv_factory.cxx
  ./src/invert/parderiv/parderiv_factory.hxx
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
  ./src/mesh/interpolation.cxx
  ./src/mesh/interpolation/bilinear.cxx
  ./src/mesh/interpolation/hermite_spline.cxx
  ./src/mesh/interpolation/interpolation_factory.cxx
  ./src/mesh/interpolation/lagrange_4pt.cxx
  ./src/mesh/interpolation/monotonic_hermite_spline.cxx
  ./src/mesh/mesh.cxx
  ./src/mesh/meshfactory.cxx
  ./src/mesh/meshfactory.hxx
  ./src/mesh/parallel/fci.cxx
  ./src/mesh/parallel/fci.hxx
  ./src/mesh/parallel/identity.cxx
  ./src/mesh/parallel/shiftedmetric.cxx
  ./src/mesh/parallel_boundary_op.cxx
  ./src/mesh/parallel_boundary_region.cxx
  ./src/mesh/surfaceiter.cxx
  ./src/physics/gyro_average.cxx
  ./src/physics/physicsmodel.cxx
  ./src/physics/smoothing.cxx
  ./src/physics/sourcex.cxx
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
  ./src/solver/impls/rkgeneric/rkschemefactory.cxx
  ./src/solver/impls/rkgeneric/rkschemefactory.hxx
  ./src/solver/impls/slepc/slepc.cxx
  ./src/solver/impls/slepc/slepc.hxx
  ./src/solver/impls/snes/snes.cxx
  ./src/solver/impls/snes/snes.hxx
  ./src/solver/impls/split-rk/split-rk.cxx
  ./src/solver/impls/split-rk/split-rk.hxx
  ./src/solver/solver.cxx
  ./src/solver/solverfactory.cxx
  ./src/sys/bout_types.cxx
  ./src/sys/boutcomm.cxx
  ./src/sys/boutexception.cxx
  ./src/sys/derivs.cxx
  ./src/sys/expressionparser.cxx
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
  )

add_library(bout++
  ${BOUT_SOURCES}
  )
add_library(bout++::bout++ ALIAS bout++)
target_link_libraries(bout++ PUBLIC MPI::MPI_CXX pvode pvpre mpark_variant)
target_include_directories(bout++ PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>
  )

target_compile_definitions(bout++
  PUBLIC "BOUT_VERSION_STRING=\"4.2.2\""
  PUBLIC "BOUT_VERSION_DOUBLE=4.22"
  PUBLIC "BOUT_HAS_PVODE"
  )
set(BOUT_HAS_PVODE TRUE)

target_compile_features(bout++ PUBLIC cxx_std_11)
set_target_properties(bout++ PROPERTIES CXX_EXTENSIONS OFF)

# Compile time features

set(CHECK_LEVELS 0 1 2 3 4)
set(CHECK 3 CACHE STRINGS "Set run-time checking level")
set_property(CACHE CHECK PROPERTY STRINGS ${CHECK_LEVELS})
if (NOT CHECK IN_LIST CHECK_LEVELS)
  message(FATAL_ERROR "CHECK must be one of ${CHECK_LEVELS}")
endif()
message(STATUS "Runtime checking level: CHECK=${CHECK}")
target_compile_definitions(bout++
  PUBLIC "CHECK=${CHECK}"
  PUBLIC "BOUT_CHECK=${CHECK}")
set(BOUT_CHECK_LEVEL ${CHECK})

option(DEBUG_ENABLED "Enable extra debug output" OFF)
if (DEBUG_ENABLED)
  target_compile_definitions(bout++
    PUBLIC "DEBUG_ENABLED"
    PUBLIC "BOUT_DEBUG_ENABLED")
endif()
message(STATUS "Extra debug output: DEBUG_ENABLED=${DEBUG_ENABLED}")
set(BOUT_HAS_OUTPUT_DEBUG ${DEBUG_ENABLED})

option(ENABLE_SIGNAL "SegFault handling" ON)
if (ENABLE_SIGNAL)
  target_compile_definitions(bout++
    PUBLIC "SIGHANDLE"
    PUBLIC "BOUT_SIGHANDLE")
endif()
message(STATUS "Signal handling: SIGHANDLE=${ENABLE_SIGNAL}")
set(BOUT_USE_SIGNAL ${ENABLE_SIGNAL})

option(ENABLE_COLOR "Output coloring" ON)
if (ENABLE_COLOR)
  target_compile_definitions(bout++
    PUBLIC "LOGCOLOR"
    PUBLIC "BOUT_LOGCOLOR")
endif()
message(STATUS "Output coloring: LOGCOLOR=${ENABLE_COLOR}")
set(BOUT_USE_COLOR ${ENABLE_COLOR})

option(ENABLE_TRACK "Field name tracking" ON)
if (ENABLE_TRACK)
  target_compile_definitions(bout++
    PUBLIC "TRACK"
    PUBLIC "BOUT_TRACK")
endif()
message(STATUS "Field name tracking: TRACK=${ENABLE_TRACK}")
set(BOUT_USE_TRACK ${ENABLE_TRACK})

option(ENABLE_SIGFPE "Signalling floating point exceptions" ON)
if (ENABLE_SIGFPE)
  target_compile_definitions(bout++
    PUBLIC "BOUT_FPE")
endif()
message(STATUS "Signalling floating point exceptions: BOUT_FPE=${ENABLE_SIGFPE}")
set(BOUT_USE_SIGFPE ${ENABLE_SIGFPE})

option(ENABLE_BACKTRACE "Enable backtrace" ON)
if (ENABLE_BACKTRACE)
  find_program(ADDR2LINE_FOUND addr2line)
  if (NOT ADDR2LINE_FOUND)
    message(FATAL_ERROR "addr2line not found")
  endif()
  target_compile_definitions(bout++
    PUBLIC "BACKTRACE"
    PUBLIC "BOUT_BACKTRACE")
  target_link_libraries(bout++ PUBLIC ${CMAKE_DL_LIBS})
endif()
message(STATUS "Enable backtrace: BACKTRACE=${ENABLE_BACKTRACE}")
set(BOUT_USE_BACKTRACE ${ENABLE_BACKTRACE})

option(ENABLE_OPENMP "Enable OpenMP support" OFF)
if (ENABLE_OPENMP)
  find_package(OpenMP REQUIRED)
  target_link_libraries(bout++ PUBLIC OpenMP::OpenMP_CXX)
  set(possible_openmp_schedules static dynamic guided auto)
  set(OPENMP_SCHEDULE static CACHE STRINGS "Set OpenMP schedule")
  if (NOT OPENMP_SCHEDULE IN_LIST possible_openmp_schedules)
    message(FATAL_ERROR "OPENMP_SCHEDULE must be one of ${possible_openmp_schedules}")
  endif()
  target_compile_definitions(bout++
    PUBLIC "OPENMP_SCHEDULE=${OPENMP_SCHEDULE}"
    PUBLIC "BOUT_OPENMP_SCHEDULE=${OPENMP_SCHEDULE}")
  message(STATUS "OpenMP schedule: ${OPENMP_SCHEDULE}")
endif()
message(STATUS "Enable OpenMP: ${ENABLE_OPENMP}")
set(BOUT_USE_OPENMP ${ENABLE_OPENMP})

include(GetGitRevisionDescription)
get_git_head_revision(GIT_REFSPEC GIT_SHA1)
message(STATUS "Git revision: ${GIT_SHA1}")
target_compile_definitions(bout++
  PUBLIC "GIT_REVISION=${GIT_SHA1}")
set(BOUT_GIT_REVISION ${GIT_SHA1})

# Optional dependencies

option(USE_NETCDF "Enable support for NetCDF output" ON)
if (USE_NETCDF)
  find_package(NetCDF REQUIRED)
  target_compile_definitions(bout++
    PUBLIC "NCDF4"
    PUBLIC "BOUT_HAS_NETCDF")
  target_link_libraries(bout++ PUBLIC NetCDF::NetCDF_CXX)
endif()
message(STATUS "NetCDF support: ${USE_NETCDF}")
set(BOUT_HAS_NETCDF ${USE_NETCDF})

option(USE_FFTW "Enable support for FFTW" ON)
if (USE_FFTW)
  find_package(FFTW REQUIRED)
  target_compile_definitions(bout++
    PUBLIC "BOUT_HAS_FFTW")
  target_link_libraries(bout++ PUBLIC FFTW::FFTW)
endif()
message(STATUS "FFTW support: ${USE_FFTW}")
set(BOUT_HAS_FFTW ${USE_FFTW})

option(USE_LAPACK "Enable support for LAPACK" ON)
if (USE_LAPACK)
  if (NOT CMAKE_SYSTEM_NAME STREQUAL "CrayLinuxEnvironment")
    # Cray wrappers sort this out for us
    find_package(LAPACK REQUIRED)
    target_link_libraries(bout++ PUBLIC "${LAPACK_LIBRARIES}")
  endif()
  target_compile_definitions(bout++
    PUBLIC "LAPACK"
    PUBLIC "BOUT_HAS_LAPACK")
endif()
message(STATUS "LAPACK support: ${USE_LAPACK}")
set(BOUT_HAS_LAPACK ${USE_LAPACK})

option(USE_PETSC "Enable support for PETSc time solvers and inversions" OFF)
if (USE_PETSC)
  if (NOT CMAKE_SYSTEM_NAME STREQUAL "CrayLinuxEnvironment")
    # Cray wrappers sort this out for us
    find_package(PETSc REQUIRED)
    target_link_libraries(bout++ PUBLIC PETSc::PETSc)
  endif()
  target_compile_definitions(bout++
    PUBLIC "BOUT_HAS_PETSC")
endif()
message(STATUS "PETSc support: ${USE_PETSC}")
set(BOUT_HAS_PETSC ${USE_PETSC})

option(USE_SLEPC "Enable support for SLEPc eigen solver" OFF)
if (USE_SLEPC)
  find_package(SLEPc REQUIRED)
  target_compile_definitions(bout++
    PUBLIC "BOUT_HAS_SLEPC")
  target_link_libraries(bout++ PUBLIC SLEPc::SLEPc)
endif()
message(STATUS "SLEPc support: ${USE_SLEPC}")
set(BOUT_HAS_SLEPC ${USE_SLEPC})

option(USE_SUNDIALS "Enable support for SUNDIALS time solvers" OFF)
if (USE_SUNDIALS)
  find_package(SUNDIALS)
  target_compile_definitions(bout++
    PUBLIC "BOUT_HAS_CVODE"
    PUBLIC "BOUT_HAS_IDA"
    PUBLIC "BOUT_HAS_ARKODE"
    PUBLIC "BOUT_HAS_SUNDIALS")
  target_link_libraries(bout++ PUBLIC SUNDIALS_cvode::cvode)
  target_link_libraries(bout++ PUBLIC SUNDIALS_ida::ida)
  target_link_libraries(bout++ PUBLIC SUNDIALS_arkode::arkode)
endif()
message(STATUS "SUNDIALS support: ${USE_SUNDIALS}")
set(BOUT_HAS_SUNDIALS ${USE_SUNDIALS})

option(USE_NLS "Enable Native Language Support" ON)
find_package(Gettext)
if (GETTEXT_FOUND)
  target_compile_definitions(bout++
    PUBLIC "BOUT_HAS_GETTEXT")
  find_package(Intl)
  if (Intl_FOUND)
    target_link_libraries(bout++
      PUBLIC ${Intl_LIBRARIES})
    target_include_directories(bout++
      PUBLIC ${Intl_INCLUDE_DIRS})
  endif()
endif()
set(BOUT_HAS_GETTEXT ${GETTEXT_FOUND})

option(USE_SCOREP "Enable support for Score-P based instrumentation" OFF)
if (USE_SCOREP)
  find_package(ScoreP REQUIRED)
  if (ScoreP_FOUND)
    target_link_libraries(bout++
      PUBLIC ScoreP::ScoreP)
    target_compile_definitions(bout++
      PUBLIC "BOUT_HAS_SCOREP")

    if (NOT DEFINED CMAKE_CXX_COMPILER_LAUNCHER)
      set_target_properties(bout++ PROPERTIES
        CXX_COMPILER_LAUNCHER "${ScoreP_COMPILER_LAUNCHER}")
    else()
      if (NOT "${CMAKE_CXX_COMPILER_LAUNCHER}" STREQUAL "${ScoreP_COMPILER_LAUNCHER}")
        message(FATAL_ERROR "Score-P conflicts with an existing compiler wrapper: ${CMAKE_CXX_COMPILER_LAUNCHER}")
      endif()
    endif()
  endif()
endif()
set(BOUT_HAS_SCOREP ${USE_SCOREP})

include(CheckCXXSourceCompiles)
check_cxx_source_compiles("int main() { const char* name = __PRETTY_FUNCTION__; }"
  HAS_PRETTY_FUNCTION)
set(BOUT_HAS_PRETTY_FUNCTION ${HAS_PRETTY_FUNCTION})
if (HAS_PRETTY_FUNCTION)
  target_compile_definitions(bout++
    PUBLIC "HAS_PRETTY_FUNCTION"
    PUBLIC "BOUT_HAS_PRETTY_FUNCTION")
endif()

# Copy FILENAME from source directory to build directory
# and add dependency on TARGET
macro(bout_test_copy_file TARGET FILENAME)
  add_custom_command(TARGET "${TARGET}" POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy
    ${CMAKE_CURRENT_SOURCE_DIR}/"${FILENAME}"
    ${CMAKE_CURRENT_BINARY_DIR}/"${FILENAME}")
endmacro()

# Add a new integrated test
#
# Required arguments:
# - TESTNAME: name of the test
# - SOURCES: list of source files
# Optional arguments:
# - USE_RUNTEST: if given, the test uses `./runtest` as the test
#   command, otherwise it uses the executable
# - USE_DATA_BOUT_INP: if given, copy `data/BOUT.inp`
# - REQUIRES: list of required dependencies
function(bout_add_integrated_test TESTNAME)
  set(options USE_RUNTEST USE_DATA_BOUT_INP)
  set(multiValueArgs SOURCES EXTRA_FILES REQUIRES)
  cmake_parse_arguments(BOUT_TEST_OPTIONS "${options}" "" "${multiValueArgs}" ${ARGN})

  foreach (REQUIREMENT "${BOUT_TEST_OPTIONS_REQUIRES}")
    if (REQUIREMENT AND NOT BOUT_HAS_${REQUIREMENT})
      message(STATUS "Not building test ${TESTNAME}, requirement not met: ${REQUIREMENT}")
      return()
    endif()
  endforeach()

  add_executable(${TESTNAME} ${BOUT_TEST_OPTIONS_SOURCES})
  target_link_libraries(${TESTNAME} bout++)
  target_include_directories(${TESTNAME} PRIVATE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)

  if (BOUT_TEST_OPTIONS_USE_RUNTEST)
    add_test(NAME ${TESTNAME} COMMAND ./runtest)
    bout_test_copy_file("${TESTNAME}" runtest)
  else()
    add_test(NAME ${TESTNAME} COMMAND ${TESTNAME})
  endif()
  if (BOUT_TEST_OPTIONS_USE_DATA_BOUT_INP)
    bout_test_copy_file("${TESTNAME}" data/BOUT.inp)
  endif()
  if (BOUT_TEST_OPTIONS_EXTRA_FILES)
    foreach (FILE ${BOUT_TEST_OPTIONS_EXTRA_FILES})
      bout_test_copy_file("${TESTNAME}" "${FILE}")
    endforeach()
  endif()
  set_target_properties(${TESTNAME} PROPERTIES FOLDER tests/integrated)
endfunction()

option(PACKAGE_TESTS "Build the tests" ON)
if(PACKAGE_TESTS)
  enable_testing()
  add_subdirectory(tests/unit)
  add_subdirectory(tests/integrated)
endif()

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
install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

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
install(
  FILES
      "${CMAKE_CURRENT_BINARY_DIR}/bout++Config.cmake"
      "${CMAKE_CURRENT_BINARY_DIR}/bout++ConfigVersion.cmake"
      "${CMAKE_CURRENT_SOURCE_DIR}/cmake/CorrectWindowsPaths.cmake"
      "${CMAKE_CURRENT_SOURCE_DIR}/cmake/FindFFTW.cmake"
      "${CMAKE_CURRENT_SOURCE_DIR}/cmake/FindNetCDF.cmake"
      "${CMAKE_CURRENT_SOURCE_DIR}/cmake/FindPackageMultipass.cmake"
      "${CMAKE_CURRENT_SOURCE_DIR}/cmake/FindPETSc.cmake"
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

   PETSc support            : ${BOUT_HAS_PETSC}
   SLEPc support            : ${BOUT_HAS_SLEPC}
   SUNDIALS support         : ${BOUT_HAS_SUNDIALS}
   NetCDF support           : ${BOUT_HAS_NETCDF}
   FFTW support             : ${BOUT_HAS_FFTW}
   LAPACK support           : ${BOUT_HAS_LAPACK}
   OpenMP support           : ${BOUT_USE_OPENMP}
   Natural language support : ${BOUT_HAS_GETTEXT}
   ScoreP support           : ${BOUT_HAS_SCOREP}
   Extra debug output       : ${BOUT_HAS_OUTPUT_DEBUG}
   CHECK level              : ${BOUT_CHECK_LEVEL}
   Signal handling          : ${BOUT_USE_SIGNAL}
   Output coloring          : ${BOUT_USE_COLOR}
   Field name tracking      : ${BOUT_USE_TRACK}
   Floating point exceptions: ${BOUT_USE_SIGFPE}
   Backtrace enabled        : ${BOUT_USE_BACKTRACE}

   === Python ===

   Make sure that the tools/pylib directory is in your PYTHONPATH
   e.g. by adding to your ~/.bashrc file

       export PYTHONPATH=$PWD/tools/pylib/:\$PYTHONPATH

*** Now run `cmake --build .` to compile BOUT++ ***
")
