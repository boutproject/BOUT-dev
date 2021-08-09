cmake_minimum_required(VERSION 3.9...3.12)

if(${CMAKE_VERSION} VERSION_LESS 3.12)
  cmake_policy(VERSION ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION})
else()
  cmake_policy(VERSION 3.12)
endif()
cmake_policy(SET CMP0025 NEW)
cmake_policy(SET CMP0048 NEW)
cmake_policy(SET CMP0057 NEW)
cmake_policy(SET CMP0074 NEW)

set(CXX_STANDARD "14")

# CMake currently doesn't support proper semver
# Set the version here, strip any extra tags to use in `project`
set(BOUT_FULL_VERSION 5.0.0-alpha)
string(REGEX REPLACE "^([0-9]+\.[0-9]+\.[0-9]+)-.*" "\\1" BOUT_CMAKE_ACCEPTABLE_VERSION ${BOUT_FULL_VERSION})
string(REGEX REPLACE "^[0-9]+\.[0-9]+\.[0-9]+-(.*)" "\\1" BOUT_VERSION_TAG ${BOUT_FULL_VERSION})

project(BOUT++
  DESCRIPTION "Fluid PDE solver framework"
  VERSION ${BOUT_CMAKE_ACCEPTABLE_VERSION}
  LANGUAGES CXX)

set(ENABLE_OPENMP On CACHE STRING "Enable OpenMP")
set(ENABLE_CUDA Off CACHE STRING "Enable CUDA")
option(ENABLE_UMPIRE "Enable UMPIRE memory management" Off)
option(ENABLE_RAJA "Enable RAJA" Off)
option(ENABLE_CALIPER "Enable Caliper" Off)

include_directories(${CMAKE_BINARY_DIR}/include)
include(CMakeDependentOption)

set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})
<<<<<<< HEAD
include(cmake/SetupCompilers.cmake)

# Override default
option(INSTALL_GTEST "Enable installation of googletest. (Projects embedding googletest may want to turn this OFF.)" ON)
option(BOUT_UPDATE_GIT_SUBMODULE "Check submodules are up-to-date during build" OFF)
option(BOUT_USE_SYSTEM_MPARK_VARIANT "Use external installation of mpark.variant" OFF)
option(BOUT_USE_SYSTEM_FMT "Use external installation of fmt" OFF)
=======
include(BOUT++functions)
>>>>>>> next

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

set(BOUT_HEADERS
  ./include/boundary_factory.hxx
  ./include/boundary_op.hxx
  ./include/boundary_region.hxx
  ./include/boundary_standard.hxx
  ./include/bout.hxx
  ./include/bout/array.hxx
  ./include/bout/assert.hxx
  ./include/bout/bout_enum_class.hxx
  ./include/bout/constants.hxx
  ./include/bout/coordinates.hxx
  ./include/bout/deprecated.hxx
  ./include/bout/deriv_store.hxx
  ./include/bout/expr.hxx
  ./include/bout/fieldgroup.hxx
  ./include/bout/format.hxx
  ./include/bout/fv_ops.hxx
  ./include/bout/generic_factory.hxx
  ./include/bout/globalfield.hxx
  ./include/bout/globalindexer.hxx
  ./include/bout/griddata.hxx
  ./include/bout/hypre_interface.hxx
  ./include/bout/hyprelib.hxx
  ./include/bout/index_derivs.hxx
  ./include/bout/index_derivs_interface.hxx
  ./include/bout/invert/laplacexy2_hypre.hxx
  #./include/bout/invert/laplacexy2.hxx
  ./include/bout/invert/laplacexy.hxx
  ./include/bout/invert/laplacexz.hxx
  ./include/bout/invertable_operator.hxx
  ./include/bout/macro_for_each.hxx
  ./include/bout/mesh.hxx
  ./include/bout/monitor.hxx
  ./include/bout/mpi_wrapper.hxx
  ./include/bout/openmpwrap.hxx
  ./include/bout/operatorstencil.hxx
  ./include/bout/paralleltransform.hxx
  ./include/bout/petsclib.hxx
  ./include/bout/petsc_interface.hxx
  ./include/bout/physicsmodel.hxx
  ./include/bout/region.hxx
  ./include/bout/rkscheme.hxx
  ./include/bout/rvec.hxx
  ./include/bout/scorepwrapper.hxx
  ./include/bout/single_index_ops.hxx
  ./include/bout/field2D_accessor.hxx
  ./include/bout/slepclib.hxx
  ./include/bout/snb.hxx
  ./include/bout/solver.hxx
  ./include/bout/solverfactory.hxx
  ./include/bout/surfaceiter.hxx
  ./include/bout/output_bout_types.hxx
  ./include/bout/sys/expressionparser.hxx
  ./include/bout/sys/generator_context.hxx
  ./include/bout/sys/gettext.hxx
  ./include/bout/sys/range.hxx
  ./include/bout/sys/timer.hxx
  ./include/bout/sys/type_name.hxx
  ./include/bout/sys/uncopyable.hxx
  ./include/bout/sys/variant.hxx
  ./include/bout/sys/Collectives.h
  ./include/bout/sys/ExecutionPolicy.h
  ./include/bout/sys/ForAll.h
  ./include/bout/sys/NVTXUtilities.h
  ./include/bout/template_combinations.hxx
  ./include/bout/traits.hxx
  ./include/bout_types.hxx
  ./include/boutcomm.hxx
  ./include/boutexception.hxx
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
  ./include/interpolation_z.hxx
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
  ./include/StackTrace.hxx
  ./include/stencils.hxx
  ./include/unused.hxx
  ./include/utils.hxx
  ./include/vecops.hxx
  ./include/vector2d.hxx
  ./include/vector3d.hxx
  ./include/where.hxx
  ./src/field/fieldgenerators.hxx
  ./src/fileio/formatfactory.hxx
  ./src/fileio/impls/emptyformat.hxx
  ./src/fileio/impls/hdf5/h5_format.hxx
  ./src/fileio/impls/netcdf/nc_format.hxx
  ./src/fileio/impls/netcdf4/ncxx4.hxx
  ./src/fileio/impls/pnetcdf/pnetcdf.hxx
  ./src/invert/laplace/impls/cyclic/cyclic_laplace.hxx
  ./src/invert/laplace/impls/multigrid/multigrid_laplace.hxx
  ./src/invert/laplace/impls/naulin/naulin_laplace.hxx
  ./src/invert/laplace/impls/pdd/pdd.hxx
  ./src/invert/laplace/impls/petsc/petsc_laplace.hxx
  ./src/invert/laplace/impls/petsc3damg/petsc3damg.hxx
  ./src/invert/laplace/impls/serial_band/serial_band.hxx
  ./src/invert/laplace/impls/serial_tri/serial_tri.hxx
  ./src/invert/laplace/impls/shoot/shoot_laplace.hxx
  ./src/invert/laplace/impls/spt/spt.hxx
  ./src/invert/laplace/impls/hypre3d/hypre3d_laplace.hxx
  ./src/invert/laplacexz/impls/cyclic/laplacexz-cyclic.hxx
  ./src/invert/laplacexz/impls/petsc/laplacexz-petsc.hxx
  ./src/invert/parderiv/impls/cyclic/cyclic.hxx
  ./src/mesh/impls/bout/boutmesh.hxx
  ./src/mesh/parallel/fci.hxx
  ./src/mesh/parallel/shiftedmetricinterp.hxx
  ./src/solver/impls/adams_bashforth/adams_bashforth.hxx
  ./src/solver/impls/arkode/arkode.hxx
  ./src/solver/impls/cvode/cvode.hxx
  ./src/solver/impls/euler/euler.hxx
  ./src/solver/impls/ida/ida.hxx
  ./src/solver/impls/imex-bdf2/imex-bdf2.hxx
  ./src/solver/impls/karniadakis/karniadakis.hxx
  ./src/solver/impls/petsc/petsc.hxx
  ./src/solver/impls/power/power.hxx
  ./src/solver/impls/pvode/pvode.hxx
  ./src/solver/impls/rk3-ssp/rk3-ssp.hxx
  ./src/solver/impls/rk4/rk4.hxx
  ./src/solver/impls/rkgeneric/impls/cashkarp/cashkarp.hxx
  ./src/solver/impls/rkgeneric/impls/rk4simple/rk4simple.hxx
  ./src/solver/impls/rkgeneric/impls/rkf34/rkf34.hxx
  ./src/solver/impls/rkgeneric/impls/rkf45/rkf45.hxx
  ./src/solver/impls/rkgeneric/rkgeneric.hxx
  ./src/solver/impls/slepc/slepc.hxx
  ./src/solver/impls/snes/snes.hxx
  ./src/solver/impls/split-rk/split-rk.hxx
  ./src/sys/options/optionparser.hxx
  ./src/sys/options/options_ini.hxx
  ${CMAKE_CURRENT_BINARY_DIR}/include/bout/revision.hxx
  ${CMAKE_CURRENT_BINARY_DIR}/include/bout/version.hxx
)


set(BOUT_SOURCES
  ./src/bout++.cxx
  ./src/bout++-time.hxx
  ./src/field/field.cxx
  ./src/field/field2d.cxx
  ./src/field/field3d.cxx
  ./src/field/field_data.cxx
  ./src/field/field_factory.cxx
  ./src/field/fieldgenerators.cxx
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
<<<<<<< HEAD
  ./src/fileio/impls/hdf5/h5_format.cxx
=======
  ./src/fileio/formatfactory.hxx
  ./src/fileio/impls/emptyformat.hxx
>>>>>>> next
  ./src/fileio/impls/netcdf/nc_format.cxx
  ./src/fileio/impls/netcdf4/ncxx4.cxx
  ./src/fileio/impls/pnetcdf/pnetcdf.cxx
  ./src/invert/fft_fftw.cxx
  ./src/invert/lapack_routines.cxx
  ./src/invert/laplace/impls/cyclic/cyclic_laplace.cxx
<<<<<<< HEAD
=======
  ./src/invert/laplace/impls/cyclic/cyclic_laplace.hxx
  ./src/invert/laplace/impls/iterative_parallel_tri/iterative_parallel_tri.cxx
  ./src/invert/laplace/impls/iterative_parallel_tri/iterative_parallel_tri.hxx
>>>>>>> next
  ./src/invert/laplace/impls/multigrid/multigrid_alg.cxx
  ./src/invert/laplace/impls/multigrid/multigrid_laplace.cxx
  ./src/invert/laplace/impls/multigrid/multigrid_solver.cxx
  ./src/invert/laplace/impls/naulin/naulin_laplace.cxx
<<<<<<< HEAD
=======
  ./src/invert/laplace/impls/naulin/naulin_laplace.hxx
  ./src/invert/laplace/impls/pcr/pcr.cxx
  ./src/invert/laplace/impls/pcr/pcr.hxx
>>>>>>> next
  ./src/invert/laplace/impls/pdd/pdd.cxx
  ./src/invert/laplace/impls/petsc/petsc_laplace.cxx
  ./src/invert/laplace/impls/petsc3damg/petsc3damg.cxx
<<<<<<< HEAD
  ./src/invert/laplace/impls/serial_band/serial_band.cxx
  ./src/invert/laplace/impls/serial_tri/serial_tri.cxx
  ./src/invert/laplace/impls/shoot/shoot_laplace.cxx
=======
  ./src/invert/laplace/impls/petsc3damg/petsc3damg.hxx
  ./src/invert/laplace/impls/serial_band/serial_band.cxx
  ./src/invert/laplace/impls/serial_band/serial_band.hxx
  ./src/invert/laplace/impls/serial_tri/serial_tri.cxx
  ./src/invert/laplace/impls/serial_tri/serial_tri.hxx
>>>>>>> next
  ./src/invert/laplace/impls/spt/spt.cxx
  ./src/invert/laplace/impls/hypre3d/hypre3d_laplace.cxx
  ./src/invert/laplace/invert_laplace.cxx
  ./src/invert/laplacexy/laplacexy.cxx
  #./src/invert/laplacexy2/laplacexy2.cxx
  ./src/invert/laplacexy2/laplacexy2_hypre.cxx 
  ./src/invert/laplacexz/impls/cyclic/laplacexz-cyclic.cxx
  ./src/invert/laplacexz/impls/petsc/laplacexz-petsc.cxx
  ./src/invert/laplacexz/laplacexz.cxx
  ./src/invert/parderiv/impls/cyclic/cyclic.cxx
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
  ./src/mesh/parallel/identity.cxx
  ./src/mesh/parallel/shiftedmetric.cxx
  ./src/mesh/parallel/shiftedmetricinterp.cxx
  ./src/mesh/parallel_boundary_op.cxx
  ./src/mesh/parallel_boundary_region.cxx
  ./src/mesh/surfaceiter.cxx
  ./src/physics/gyro_average.cxx
  ./src/physics/physicsmodel.cxx
  ./src/physics/smoothing.cxx
  ./src/physics/snb.cxx
  ./src/physics/sourcex.cxx
  ./src/solver/impls/adams_bashforth/adams_bashforth.cxx
  ./src/solver/impls/arkode/arkode.cxx
  ./src/solver/impls/cvode/cvode.cxx
  ./src/solver/impls/euler/euler.cxx
  ./src/solver/impls/ida/ida.cxx
  ./src/solver/impls/imex-bdf2/imex-bdf2.cxx
<<<<<<< HEAD
  ./src/solver/impls/karniadakis/karniadakis.cxx
=======
  ./src/solver/impls/imex-bdf2/imex-bdf2.hxx
>>>>>>> next
  ./src/solver/impls/petsc/petsc.cxx
  ./src/solver/impls/power/power.cxx
  ./src/solver/impls/pvode/pvode.cxx
  ./src/solver/impls/rk3-ssp/rk3-ssp.cxx
  ./src/solver/impls/rk4/rk4.cxx
  ./src/solver/impls/rkgeneric/impls/cashkarp/cashkarp.cxx
  ./src/solver/impls/rkgeneric/impls/rk4simple/rk4simple.cxx
  ./src/solver/impls/rkgeneric/impls/rkf34/rkf34.cxx
  ./src/solver/impls/rkgeneric/impls/rkf45/rkf45.cxx
  ./src/solver/impls/rkgeneric/rkgeneric.cxx
  ./src/solver/impls/rkgeneric/rkscheme.cxx
  ./src/solver/impls/slepc/slepc.cxx
  ./src/solver/impls/snes/snes.cxx
  ./src/solver/impls/split-rk/split-rk.cxx
  ./src/solver/solver.cxx
  ./src/sys/bout_types.cxx
  ./src/sys/boutcomm.cxx
  ./src/sys/boutexception.cxx
  ./src/sys/derivs.cxx
  ./src/sys/expressionparser.cxx
  ./src/sys/generator_context.cxx
  ./src/sys/hyprelib.cxx
  ./src/sys/msg_stack.cxx
  ./src/sys/options.cxx
  ./src/sys/options/options_ini.cxx
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

<<<<<<< HEAD
if(BUILD_SHARED_LIBS)
   add_library(bout++ SHARED ${BOUT_SOURCES} ${BOUT_HEADERS})
else ()
   add_library(bout++ STATIC ${BOUT_SOURCES} ${BOUT_HEADERS})
endif ()

set_target_properties(bout++ PROPERTIES LINKER_LANGUAGE CXX)
=======
option(BOUT_GENERATE_FIELDOPS "Automatically re-generate the Field arithmetic operators from the Python templates. \
Requires Python, clang-format, and Jinja2. Turn this OFF to skip generating them if, for example, \
you are unable to install the Jinja2 Python module. This is only important for BOUT++ developers." ON)

find_package(Python)
find_package(ClangFormat)
if (BOUT_GENERATE_FIELDOPS AND Python_FOUND AND ClangFormat_FOUND)
  add_custom_command( OUTPUT ${CMAKE_CURRENT_SOURCE_DIR}/src/field/generated_fieldops.cxx
    COMMAND ${Python_EXECUTABLE} gen_fieldops.py --filename generated_fieldops.cxx.tmp
    COMMAND ${ClangFormat_BIN} generated_fieldops.cxx.tmp -i
    COMMAND ${CMAKE_COMMAND} -E rename generated_fieldops.cxx.tmp generated_fieldops.cxx
    DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/src/field/gen_fieldops.jinja  ${CMAKE_CURRENT_SOURCE_DIR}/src/field/gen_fieldops.py
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/src/field/
    COMMENT "Generating source code" )
else()
  message(AUTHOR_WARNING "'src/field/generated_fieldops.cxx' will not be \
regenerated when you make changes to either \
'src/field/gen_fieldops.py' or 'src/field/gen_fieldops.jinja'. \
This is because either Python or clang-format is missing \
(see above messages for more information) \
or BOUT_GENERATE_FIELDOPS is OFF (current value: ${BOUT_GENERATE_FIELDOPS}). \
This warning is only important for BOUT++ developers and can otherwise be \
safely ignored.")
endif()

# Ensure that the compile date/time is up-to-date when any of the sources change
add_custom_command(
  OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/bout++-time.cxx
  COMMAND ${CMAKE_COMMAND} -P "${CMAKE_CURRENT_LIST_DIR}/cmake/GenerateDateTimeFile.cmake"
  DEPENDS ${BOUT_SOURCES}
  MAIN_DEPENDENCY "${CMAKE_CURRENT_LIST_DIR}/cmake/GenerateDateTimeFile.cmake"
  )

add_library(bout++
  ${BOUT_SOURCES}
  ${CMAKE_CURRENT_BINARY_DIR}/bout++-time.cxx
  )
>>>>>>> next
add_library(bout++::bout++ ALIAS bout++)

target_include_directories(bout++ PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/include>
  $<INSTALL_INTERFACE:include>
  )
set(BOUT_LIB_PATH "${CMAKE_CURRENT_BINARY_DIR}/lib")
set_target_properties(bout++ PROPERTIES
  LIBRARY_OUTPUT_DIRECTORY "${BOUT_LIB_PATH}"
  ARCHIVE_OUTPUT_DIRECTORY "${BOUT_LIB_PATH}")

# Set some variables for the bout-config script
set(CONFIG_LDFLAGS "${CONFIG_LDFLAGS} -L\$BOUT_LIB_PATH -lbout++")
set(BOUT_INCLUDE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/include")
set(CONFIG_CFLAGS "${CONFIG_CFLAGS} -I\${BOUT_INCLUDE_PATH} -I${CMAKE_CURRENT_BINARY_DIR}/include")


target_compile_features(bout++ PUBLIC cxx_std_14)
set_target_properties(bout++ PROPERTIES CXX_EXTENSIONS OFF)

<<<<<<< HEAD
# Optional dependencies
include(cmake/SetupBOUTThirdParty.cmake)

#if(BOUT_HAS_HYPRE)
#   target_compile_definitions(bout++ PUBLIC "BOUT_HAS_HYPRE")
#endif ()
=======
# Various sanitizers, including coverage and address sanitizer
include(cmake/Sanitizers.cmake)
enable_sanitizers(bout++)
>>>>>>> next

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

<<<<<<< HEAD
# temporarily turn this off - need to reconcile with nvcc
option(ENABLE_WARNINGS "Enable compiler warnings" Off)
#option(BOUT_ENABLE_WARNINGS "Enable compiler warnings" ON)
=======
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

option(BOUT_ENABLE_WARNINGS "Enable compiler warnings" ON)
>>>>>>> next
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

if (CHECK GREATER 1)
  set(bout_use_msgstack_default ON)
else()
  set(bout_use_msgstack_default OFF)
endif()
set(BOUT_ENABLE_MSGSTACK ${bout_use_msgstack_default} CACHE BOOL "Enable debug message stack")
message(STATUS "Message stack: BOUT_USE_MSGSTACK=${BOUT_ENABLE_MSGSTACK}")
set(BOUT_USE_MSGSTACK ${BOUT_ENABLE_MSGSTACK})

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
  set(CONFIG_LDFLAGS "${CONFIG_LDFLAGS} -l${CMAKE_DL_LIBS}")
endif()
message(STATUS "Enable backtrace: BOUT_USE_BACKTRACE=${BOUT_ENABLE_BACKTRACE}")
set(BOUT_USE_BACKTRACE ${BOUT_ENABLE_BACKTRACE})

<<<<<<< HEAD
=======
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

option(BOUT_USE_PETSC "Enable support for PETSc time solvers and inversions" OFF)
if (BOUT_USE_PETSC)
  if (NOT CMAKE_SYSTEM_NAME STREQUAL "CrayLinuxEnvironment")
    # Cray wrappers sort this out for us
    find_package(PETSc REQUIRED)
    target_link_libraries(bout++ PUBLIC PETSc::PETSc)
    string(JOIN " " CONFIG_PETSC_LIBRARIES ${PETSC_LIBRARIES})
    set(CONFIG_LDFLAGS "${CONFIG_PETSC_LIBRARIES} ${CONFIG_LDFLAGS}")
  endif()
endif()
message(STATUS "PETSc support: ${BOUT_USE_PETSC}")
set(BOUT_HAS_PETSC ${BOUT_USE_PETSC})

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
    set(ENABLE_MPI ON CACHE BOOL "" FORCE)
    set(ENABLE_OPENMP OFF CACHE BOOL "" FORCE)
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

>>>>>>> next
include(CheckCXXSourceCompiles)
check_cxx_source_compiles("int main() { const char* name = __PRETTY_FUNCTION__; }"
  HAS_PRETTY_FUNCTION)
set(BOUT_HAS_PRETTY_FUNCTION ${HAS_PRETTY_FUNCTION})

# Locations of the various Python modules, including the generated boutconfig module
set(BOUT_PYTHONPATH "${CMAKE_CURRENT_BINARY_DIR}/tools/pylib:${CMAKE_CURRENT_SOURCE_DIR}/tools/pylib")
# Variables for boutconfig module -- note that these will contain
# generator expressions and CMake targets, and not generally be very
# useful
get_target_property(BOUT_LIBS bout++ INTERFACE_LINK_LIBRARIES)
get_target_property(BOUT_CFLAGS bout++ INTERFACE_INCLUDE_DIRECTORIES)

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

<<<<<<< HEAD
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
#
# - EXECUTABLE_NAME: name of the executable, if different from the
#   first source name

function(bout_add_integrated_test TESTNAME)
  set(options USE_RUNTEST USE_DATA_BOUT_INP)
  set(oneValueArgs EXECUTABLE_NAME)
  set(multiValueArgs SOURCES EXTRA_FILES REQUIRES)
  cmake_parse_arguments(BOUT_TEST_OPTIONS "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

  foreach (REQUIREMENT IN LISTS BOUT_TEST_OPTIONS_REQUIRES)
    if (NOT ${REQUIREMENT})
      message(STATUS "Not building test ${TESTNAME}, requirement not met: ${REQUIREMENT}")
      return()
    endif()
  endforeach()

  add_executable(${TESTNAME} ${BOUT_TEST_OPTIONS_SOURCES})
  target_link_libraries(${TESTNAME} bout++)
  target_include_directories(${TESTNAME} PRIVATE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)

  if(BOUT_USE_CUDA)
     set_source_files_properties(${BOUT_TEST_OPTIONS_SOURCES} PROPERTIES LANGUAGE CUDA )
     target_include_directories(${TESTNAME} PUBLIC ${BOUT_CUDA_LIB_DIR}/../include ) 
    if(BOUT_HAS_HYPRE)
     target_compile_definitions(${TESTNAME} PUBLIC "HYPRE_USING_CUDA;HYPRE_USING_UNIFIED_MEMORY")
    endif ()
  endif ()

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

option(PACKAGE_TESTS "Build the tests" OFF)
=======
option(PACKAGE_TESTS "Build the tests" ON)
option(BOUT_RUN_ALL_TESTS "Run all of the tests (this can be slow!)" OFF)
>>>>>>> next
if(PACKAGE_TESTS)
  enable_testing()

  # Targets for just building the tests
  # Tests need to add themselves as dependencies to these targets
  add_custom_target(build-check-unit-tests)
  add_custom_target(build-check-integrated-tests)
  add_custom_target(build-check-mms-tests)

  # Build all the tests
  add_custom_target(build-check)
  add_dependencies(build-check build-check-unit-tests build-check-integrated-tests build-check-mms-tests)

<<<<<<< HEAD
  #add_subdirectory(tests/unit EXCLUDE_FROM_ALL)
  #add_subdirectory(tests/integrated EXCLUDE_FROM_ALL)
  add_subdirectory(tests/unit)
  #add_subdirectory(tests/integrated)
=======
  add_subdirectory(tests/unit EXCLUDE_FROM_ALL)
  add_subdirectory(tests/integrated EXCLUDE_FROM_ALL)
  add_subdirectory(tests/MMS EXCLUDE_FROM_ALL)
>>>>>>> next

  # Targets for running the tests
  add_custom_target(check-unit-tests
    COMMAND ctest -R serial_tests --output-on-failure)
  add_dependencies(check-unit-tests build-check-unit-tests)

  add_custom_target(check-integrated-tests
    COMMAND ctest -R "test-" --output-on-failure)
  add_dependencies(check-integrated-tests build-check-integrated-tests)

  add_custom_target(check-mms-tests
    COMMAND ctest -R "MMS-" --output-on-failure)
  add_dependencies(check-mms-tests build-check-mms-tests)

  # Run all the tests
<<<<<<< HEAD
  add_custom_target(check_bout)
  add_dependencies(check_bout check-unit-tests check-integrated-tests)
=======
  add_custom_target(check)
  add_dependencies(check check-unit-tests check-integrated-tests check-mms-tests)
endif()

option(BOUT_BUILD_EXAMPLES "Build the examples" OFF)
if(BOUT_BUILD_EXAMPLES)
  add_subdirectory(examples)
endif()

##################################################
# Documentation

option(BOUT_BUILD_DOCS "Build the documentation" OFF)
if (BOUT_BUILD_DOCS)
  add_subdirectory(manual EXCLUDE_FROM_ALL)
>>>>>>> next
endif()

if(BOUT_HAS_HYPRE)
   add_subdirectory(examples/laplacexy/simple-hypre)
endif ()

add_subdirectory(examples/laplacexy/simple)
add_subdirectory(examples/hasegawa-wakatani-3d)
add_subdirectory(examples/elm-pb)
add_subdirectory(tests/integrated/test-laplace-petsc3d)
add_subdirectory(tests/integrated/test-laplace-hypre3d)
add_subdirectory(tests/integrated/test-laplacexy2-hypre)

##################################################
# Generate the build config header

if (EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/include/bout/build_defines.hxx")
  message(FATAL_ERROR "Generated build_defines.hxx header already exists; please remove '${CMAKE_CURRENT_SOURCE_DIR}/include/bout/build_defines.hxx' before continuing")
endif()

configure_file(cmake_build_defines.hxx.in include/bout/build_defines.hxx)


##################################################
# Generate the bout-config script

# Set some variables to match autotools so we can use the same input file
set(CXX "${MPI_CXX_COMPILER}")
set(PYTHONCONFIGPATH "${BOUT_PYTHONPATH}")
set(BOUT_HAS_LEGACY_NETCDF OFF)
set(BOUT_HAS_PNETCDF OFF)

# For shared libraries we only need to know how to link against BOUT++,
# while for static builds we need the dependencies too
if (BUILD_SHARED_LIBS)
  # Include rpath linker flag so user doesn't need to set LD_LIBRARY_PATH
  set(CONFIG_LDFLAGS "${CMAKE_SHARED_LIBRARY_RUNTIME_CXX_FLAG}\$BOUT_LIB_PATH -L\$BOUT_LIB_PATH -lbout++")
else()
  set(CONFIG_LDFLAGS "${CONFIG_LDFLAGS}")
endif()

# This version of the file allows the build directory to be used directly
configure_file(bin/bout-config.in bin/bout-config @ONLY)

# We need to generate a separate version for installation, with the
# correct install paths. So first we need to replace the build
# directory library path with the installation path
string(REPLACE
  "${CMAKE_BINARY_DIR}/lib" "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}"
  CONFIG_LDFLAGS "${CONFIG_LDFLAGS}")
# Update mpark.variant and fmt include paths if we're building them
if (NOT BOUT_USE_SYSTEM_MPARK_VARIANT)
  set(MPARK_VARIANT_INCLUDE_PATH "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_INCLUDEDIR}")
endif()
if (NOT BOUT_USE_SYSTEM_FMT)
  set(FMT_INCLUDE_PATH "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_INCLUDEDIR}")
endif()
set(BOUT_INCLUDE_PATH "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_INCLUDEDIR}")
# We don't need the build include path any more
string(REPLACE "-I${CMAKE_CURRENT_BINARY_DIR}/include" "" CONFIG_CFLAGS "${CONFIG_CFLAGS}")

# This version now has the correct paths to use the final installation
configure_file(bin/bout-config.in bin/bout-config-install @ONLY)

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

# The various helper scripts
install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/bin/"
  USE_SOURCE_PERMISSIONS
  DESTINATION "${CMAKE_INSTALL_BINDIR}")

# The installed version of bout-config needs renaming when we install
# it. Note this MUST be done after the installation of bin/, to make
# sure we clobber any versions of bout-config hanging around from an
# autotools build
install(PROGRAMS "${CMAKE_CURRENT_BINARY_DIR}/bin/bout-config-install"
  DESTINATION "${CMAKE_INSTALL_BINDIR}"
  RENAME "bout-config"
  )

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

# CMake configuration files
install(
  FILES
      "${CMAKE_CURRENT_BINARY_DIR}/bout++Config.cmake"
      "${CMAKE_CURRENT_BINARY_DIR}/bout++ConfigVersion.cmake"
      "${CMAKE_CURRENT_SOURCE_DIR}/cmake/BOUT++functions.cmake"
      "${CMAKE_CURRENT_SOURCE_DIR}/cmake/CorrectWindowsPaths.cmake"
      "${CMAKE_CURRENT_SOURCE_DIR}/cmake/FindClangFormat.cmake"
      "${CMAKE_CURRENT_SOURCE_DIR}/cmake/FindFFTW.cmake"
<<<<<<< HEAD
      "${CMAKE_CURRENT_SOURCE_DIR}/cmake/FindHYPRE.cmake"
      "${CMAKE_CURRENT_SOURCE_DIR}/cmake/FindNetCDF.cmake"
=======
      "${CMAKE_CURRENT_SOURCE_DIR}/cmake/FindnetCDF.cmake"
      "${CMAKE_CURRENT_SOURCE_DIR}/cmake/FindnetCDFCxx.cmake"
>>>>>>> next
      "${CMAKE_CURRENT_SOURCE_DIR}/cmake/FindPackageMultipass.cmake"
      "${CMAKE_CURRENT_SOURCE_DIR}/cmake/FindLibuuid.cmake"
      "${CMAKE_CURRENT_SOURCE_DIR}/cmake/FindPETSc.cmake"
      "${CMAKE_CURRENT_SOURCE_DIR}/cmake/FindScoreP.cmake"
      "${CMAKE_CURRENT_SOURCE_DIR}/cmake/FindSLEPc.cmake"
      "${CMAKE_CURRENT_SOURCE_DIR}/cmake/FindSUNDIALS.cmake"
      "${CMAKE_CURRENT_SOURCE_DIR}/cmake/FindSphinx.cmake"
      "${CMAKE_CURRENT_SOURCE_DIR}/cmake/ResolveCompilerPaths.cmake"
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/bout++"
  )

configure_package_config_file(tools/pylib/boutconfig/__init__.py.cin tools/pylib/boutconfig/__init__.py
  INSTALL_DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/tools/pylib/boutconfig/__init__.py"
  )

install(
  FILES "${CMAKE_CURRENT_BINARY_DIR}/tools/pylib/boutconfig/__init__.py"
  DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/boutconfig"
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
   HYPRE support            : ${BOUT_HAS_HYPRE}
   NetCDF support           : ${BOUT_HAS_NETCDF}
   FFTW support             : ${BOUT_HAS_FFTW}
   LAPACK support           : ${BOUT_HAS_LAPACK}
   OpenMP support           : ${BOUT_USE_OPENMP}
   Natural language support : ${BOUT_HAS_GETTEXT}
   ScoreP support           : ${BOUT_HAS_SCOREP}
   System UUID generator    : ${BOUT_HAS_UUID_SYSTEM_GENERATOR}
   Extra debug output       : ${BOUT_USE_OUTPUT_DEBUG}
   Runtime error check level: ${BOUT_CHECK_LEVEL}
   Signal handling          : ${BOUT_USE_SIGNAL}
   Output coloring          : ${BOUT_USE_COLOR}
   Field name tracking      : ${BOUT_USE_TRACK}
   Floating point exceptions: ${BOUT_USE_SIGFPE}
   Backtrace enabled        : ${BOUT_USE_BACKTRACE}
<<<<<<< HEAD
   RAJA enabled             : ${BOUT_HAS_RAJA}
   Umpire enabled           : ${BOUT_HAS_UMPIRE}
   CUDA enabled             : ${BOUT_USE_CUDA}
=======
   Sanitizers enabled       : ${BOUT_USE_SANITIZERS}
>>>>>>> next

   === Python ===

   Make sure that the tools/pylib directory is in your PYTHONPATH
   e.g. by adding to your ~/.bashrc file

       export PYTHONPATH=${BOUT_PYTHONPATH}:\$PYTHONPATH

*** Now run `cmake --build .` to compile BOUT++ ***
")
