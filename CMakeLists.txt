cmake_minimum_required(VERSION 3.17...3.20)
if(${CMAKE_VERSION} VERSION_LESS 3.20)
  cmake_policy(VERSION ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION})
else()
  cmake_policy(VERSION 3.17)
endif()

# Find python based on version
if(${CMAKE_VERSION} VERSION_GREATER_EQUAL 3.15)
  cmake_policy(SET CMP0094 NEW)
endif()
cmake_policy(SET CMP0095 NEW)
# We support the new syntax, so avoid deprecation warning
if(${CMAKE_VERSION} VERSION_GREATER_EQUAL 3.22)
  cmake_policy(SET CMP0127 NEW)
endif()

if ("${CMAKE_CURRENT_BINARY_DIR}" STREQUAL "${CMAKE_CURRENT_SOURCE_DIR}")
  option(BOUT_ALLOW_INSOURCE_BUILD "Whether BOUT++ should really allow to build in source." OFF)
  if (NOT ${BOUT_ALLOW_INSOURCE_BUILD})
    message(FATAL_ERROR "BOUT++ does not recommend in source builds. Try building out of source, e.g. with `cmake -S . -B build` or set -DBOUT_ALLOW_INSOURCE_BUILD=ON - but things may break!")
  endif()
endif()

# CMake currently doesn't support proper semver
# Set the version here, strip any extra tags to use in `project`
# We try to use git to get a full description, inspired by setuptools_scm
set(_bout_previous_version "5.1.1")
set(_bout_next_version "5.1.2")
execute_process(
  COMMAND "git" describe --tags --match=v${_bout_previous_version}
  COMMAND sed -e s/${_bout_previous_version}-/${_bout_next_version}.dev/ -e s/-/+/
  RESULTS_VARIABLE error_codes
  OUTPUT_VARIABLE BOUT_FULL_VERSION
  OUTPUT_STRIP_TRAILING_WHITESPACE
  )

foreach(error_code ${error_codes})
  if (NOT ${error_code} STREQUAL 0)
    set(BOUT_FULL_VERSION ${_bout_next_version})
  endif()
endforeach()
# Remove leading "v"
string(REGEX REPLACE "^v(.*)" "\\1" BOUT_FULL_VERSION ${BOUT_FULL_VERSION})
# Remove trailing tag
string(REGEX REPLACE "^([0-9]+\.[0-9]+\.[0-9]+)\..*" "\\1" BOUT_CMAKE_ACCEPTABLE_VERSION ${BOUT_FULL_VERSION})
# Get the trailing tag
string(REGEX REPLACE "^[0-9]+\.[0-9]+\.[0-9]+\.(.*)" "\\1" BOUT_VERSION_TAG ${BOUT_FULL_VERSION})

message(STATUS "Configuring BOUT++ version ${BOUT_FULL_VERSION}")
project(BOUT++
  DESCRIPTION "Fluid PDE solver framework"
  VERSION ${BOUT_CMAKE_ACCEPTABLE_VERSION}
  LANGUAGES CXX)

include(CMakeDependentOption)

option(BUILD_SHARED_LIBS "Build shared libs" ON)

# Override default
option(INSTALL_GTEST "Enable installation of googletest. (Projects embedding googletest may want to turn this OFF.)" OFF)

set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})
include(BOUT++functions)

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
      WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
      RESULT_VARIABLE GIT_SUBMOD_RESULT)
    if(NOT GIT_SUBMOD_RESULT EQUAL "0")
      message(FATAL_ERROR "git submodule update --init failed with ${GIT_SUBMOD_RESULT}, please checkout submodules")
    endif()
  endif()
endfunction()

set(BOUT_SOURCES
  ./include/bout/array.hxx
  ./include/bout/assert.hxx
  ./include/bout/boundary_factory.hxx
  ./include/bout/boundary_op.hxx
  ./include/bout/boundary_region.hxx
  ./include/bout/boundary_standard.hxx
  ./include/bout/bout.hxx
  ./include/bout/bout_enum_class.hxx
  ./include/bout/bout_types.hxx
  ./include/bout/build_config.hxx
  ./include/bout/boutcomm.hxx
  ./include/bout/boutexception.hxx
  ./include/bout/caliper_wrapper.hxx
  ./include/bout/constants.hxx
  ./include/bout/coordinates.hxx
  ./include/bout/coordinates_accessor.hxx
  ./include/bout/cyclic_reduction.hxx
  ./include/bout/dcomplex.hxx
  ./include/bout/deriv_store.hxx
  ./include/bout/derivs.hxx
  ./include/bout/difops.hxx
  ./include/bout/expr.hxx
  ./include/bout/fft.hxx
  ./include/bout/field.hxx
  ./include/bout/field2d.hxx
  ./include/bout/field3d.hxx
  ./include/bout/field_accessor.hxx
  ./include/bout/field_data.hxx
  ./include/bout/field_factory.hxx
  ./include/bout/fieldgroup.hxx
  ./include/bout/fieldperp.hxx
  ./include/bout/format.hxx
  ./include/bout/fv_ops.hxx
  ./include/bout/generic_factory.hxx
  ./include/bout/globalfield.hxx
  ./include/bout/globalindexer.hxx
  ./include/bout/globals.hxx
  ./include/bout/griddata.hxx
  ./include/bout/gyro_average.hxx
  ./include/bout/hypre_interface.hxx
  ./include/bout/index_derivs.hxx
  ./include/bout/index_derivs_interface.hxx
  ./include/bout/initialprofiles.hxx
  ./include/bout/interpolation.hxx
  ./include/bout/interpolation_xz.hxx
  ./include/bout/interpolation_z.hxx
  ./include/bout/invert/laplacexy.hxx
  ./include/bout/invert/laplacexz.hxx
  ./include/bout/invert_laplace.hxx
  ./include/bout/invert_parderiv.hxx
  ./include/bout/invert_pardiv.hxx
  ./include/bout/invertable_operator.hxx
  ./include/bout/lapack_routines.hxx
  ./include/bout/macro_for_each.hxx
  ./include/bout/mask.hxx
  ./include/bout/mesh.hxx
  ./include/bout/monitor.hxx
  ./include/bout/mpi_wrapper.hxx
  ./include/bout/msg_stack.hxx
  ./include/bout/multiostream.hxx
  ./include/bout/openmpwrap.hxx
  ./include/bout/operatorstencil.hxx
  ./include/bout/options.hxx
  ./include/bout/options_netcdf.hxx
  ./include/bout/optionsreader.hxx
  ./include/bout/output.hxx
  ./include/bout/output_bout_types.hxx
  ./include/bout/parallel_boundary_op.hxx
  ./include/bout/parallel_boundary_region.hxx
  ./include/bout/paralleltransform.hxx
  ./include/bout/petsc_interface.hxx
  ./include/bout/petsclib.hxx
  ./include/bout/physicsmodel.hxx
  ./include/bout/rajalib.hxx
  ./include/bout/region.hxx
  ./include/bout/rkscheme.hxx
  ./include/bout/rvec.hxx
  ./include/bout/scorepwrapper.hxx
  ./include/bout/single_index_ops.hxx
  ./include/bout/slepclib.hxx
  ./include/bout/smoothing.hxx
  ./include/bout/snb.hxx
  ./include/bout/solver.hxx
  ./include/bout/solverfactory.hxx
  ./include/bout/sourcex.hxx
  ./include/bout/stencils.hxx
  ./include/bout/sundials_backports.hxx
  ./include/bout/surfaceiter.hxx
  ./include/bout/sys/expressionparser.hxx
  ./include/bout/sys/generator_context.hxx
  ./include/bout/sys/gettext.hxx
  ./include/bout/sys/range.hxx
  ./include/bout/sys/timer.hxx
  ./include/bout/sys/type_name.hxx
  ./include/bout/sys/uncopyable.hxx
  ./include/bout/sys/uuid.h
  ./include/bout/sys/variant.hxx
  ./include/bout/template_combinations.hxx
  ./include/bout/traits.hxx
  ./include/bout/unused.hxx
  ./include/bout/utils.hxx
  ./include/bout/vecops.hxx
  ./include/bout/vector2d.hxx
  ./include/bout/vector3d.hxx
  ./include/bout/where.hxx
  ./src/bout++.cxx
  ./src/bout++-time.hxx
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
  ./src/invert/fft_fftw.cxx
  ./src/invert/lapack_routines.cxx
  ./src/invert/laplace/impls/cyclic/cyclic_laplace.cxx
  ./src/invert/laplace/impls/cyclic/cyclic_laplace.hxx
  ./src/invert/laplace/impls/iterative_parallel_tri/iterative_parallel_tri.cxx
  ./src/invert/laplace/impls/iterative_parallel_tri/iterative_parallel_tri.hxx
  ./src/invert/laplace/impls/multigrid/multigrid_alg.cxx
  ./src/invert/laplace/impls/multigrid/multigrid_laplace.cxx
  ./src/invert/laplace/impls/multigrid/multigrid_laplace.hxx
  ./src/invert/laplace/impls/multigrid/multigrid_solver.cxx
  ./src/invert/laplace/impls/naulin/naulin_laplace.cxx
  ./src/invert/laplace/impls/naulin/naulin_laplace.hxx
  ./src/invert/laplace/impls/pcr/pcr.cxx
  ./src/invert/laplace/impls/pcr/pcr.hxx
  ./src/invert/laplace/impls/pcr_thomas/pcr_thomas.cxx
  ./src/invert/laplace/impls/pcr_thomas/pcr_thomas.hxx
  ./src/invert/laplace/impls/petsc/petsc_laplace.cxx
  ./src/invert/laplace/impls/petsc/petsc_laplace.hxx
  ./src/invert/laplace/impls/petsc3damg/petsc3damg.cxx
  ./src/invert/laplace/impls/petsc3damg/petsc3damg.hxx
  ./src/invert/laplace/impls/serial_band/serial_band.cxx
  ./src/invert/laplace/impls/serial_band/serial_band.hxx
  ./src/invert/laplace/impls/serial_tri/serial_tri.cxx
  ./src/invert/laplace/impls/serial_tri/serial_tri.hxx
  ./src/invert/laplace/impls/spt/spt.cxx
  ./src/invert/laplace/impls/spt/spt.hxx
  ./src/invert/laplace/impls/hypre3d/hypre3d_laplace.cxx
  ./src/invert/laplace/impls/hypre3d/hypre3d_laplace.hxx
  ./src/invert/laplace/invert_laplace.cxx
  ./src/invert/laplacexy/laplacexy.cxx
  ./include/bout/invert/laplacexy2.hxx
  ./src/invert/laplacexy2/laplacexy2.cxx
  ./include/bout/invert/laplacexy2_hypre.hxx
  ./src/invert/laplacexy2/laplacexy2_hypre.cxx 
  ./src/invert/laplacexz/impls/cyclic/laplacexz-cyclic.cxx
  ./src/invert/laplacexz/impls/cyclic/laplacexz-cyclic.hxx
  ./src/invert/laplacexz/impls/petsc/laplacexz-petsc.cxx
  ./src/invert/laplacexz/impls/petsc/laplacexz-petsc.hxx
  ./src/invert/laplacexz/laplacexz.cxx
  ./src/invert/parderiv/impls/cyclic/cyclic.cxx
  ./src/invert/parderiv/impls/cyclic/cyclic.hxx
  ./src/invert/parderiv/invert_parderiv.cxx
  ./src/invert/pardiv/impls/cyclic/pardiv_cyclic.cxx
  ./src/invert/pardiv/impls/cyclic/pardiv_cyclic.hxx
  ./src/invert/pardiv/invert_pardiv.cxx
  ./src/mesh/boundary_factory.cxx
  ./src/mesh/boundary_region.cxx
  ./src/mesh/boundary_standard.cxx
  ./src/mesh/coordinates.cxx
  ./src/mesh/coordinates_accessor.cxx
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
  ./include/bout/hyprelib.hxx
  ./src/sys/hyprelib.cxx
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


find_package(Python3)
find_package(ClangFormat)

if (Python3_FOUND AND ClangFormat_FOUND)
  set(BOUT_GENERATE_FIELDOPS_DEFAULT ON)
else()
  set(BOUT_GENERATE_FIELDOPS_DEFAULT OFF)
endif()

execute_process(COMMAND ${Python3_EXECUTABLE} -c "import zoidberg"
  RESULT_VARIABLE zoidberg_FOUND)
if (zoidberg_FOUND EQUAL 0)
  set(zoidberg_FOUND ON)
else()
  set(zoidberg_FOUND OFF)
endif()

option(BOUT_GENERATE_FIELDOPS "Automatically re-generate the Field arithmetic operators from the Python templates. \
Requires Python3, clang-format, and Jinja2. Turn this OFF to skip generating them if, for example, \
you are unable to install the Jinja2 Python module. This is only important for BOUT++ developers." ${BOUT_GENERATE_FIELDOPS_DEFAULT})

if (BOUT_GENERATE_FIELDOPS)
  if (NOT Python3_FOUND)
    message(FATAL_ERROR "python not found, but you have requested to generate code!")
  endif()
  if (NOT ClangFormat_FOUND)
    message(FATAL_ERROR "clang-format not found, but you have requested to generate code!")
  endif()
  add_custom_command( OUTPUT ${CMAKE_CURRENT_SOURCE_DIR}/src/field/generated_fieldops.cxx
    COMMAND ${Python3_EXECUTABLE} gen_fieldops.py --filename generated_fieldops.cxx.tmp
    COMMAND ${ClangFormat_BIN} generated_fieldops.cxx.tmp -i
    COMMAND ${CMAKE_COMMAND} -E rename generated_fieldops.cxx.tmp generated_fieldops.cxx
    DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/src/field/gen_fieldops.jinja  ${CMAKE_CURRENT_SOURCE_DIR}/src/field/gen_fieldops.py
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/src/field/
    COMMENT "Generating source code" )
else()
  message(AUTHOR_WARNING "'src/field/generated_fieldops.cxx' will not be \
regenerated when you make changes to either \
'src/field/gen_fieldops.py' or 'src/field/gen_fieldops.jinja'. \
This is because either Python3 or clang-format is missing \
(see above messages for more information) \
or BOUT_GENERATE_FIELDOPS is OFF (current value: ${BOUT_GENERATE_FIELDOPS}). \
This warning is only important for BOUT++ developers and can otherwise be \
safely ignored.")
endif()

include(GNUInstallDirs)

# use, i.e. don't skip the full RPATH for the build tree
set(CMAKE_SKIP_BUILD_RPATH FALSE)

# when building, don't use the install RPATH already
# (but later on when installing)
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}")

# add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)


execute_process(COMMAND ${Python3_EXECUTABLE} -c "import site ; print('/'.join(site.getusersitepackages().split('/')[-2:]))"
  RESULT_VARIABLE PYTHON_WORKING
  OUTPUT_VARIABLE PYTHON_SITEPATH_SUFFIX
  OUTPUT_STRIP_TRAILING_WHITESPACE
)
set(CMAKE_INSTALL_PYTHON_SITEARCH lib/${PYTHON_SITEPATH_SUFFIX} CACHE STRING "Location to install python arch-specific modules")

set(ON_OFF_AUTO ON OFF AUTO)
set(BOUT_ENABLE_PYTHON AUTO CACHE STRING "Build the Python interface")
set_property(CACHE BOUT_ENABLE_PYTHON PROPERTY STRINGS ${ON_OFF_AUTO})
if (NOT BOUT_ENABLE_PYTHON IN_LIST ON_OFF_AUTO)
  message(FATAL_ERROR "BOUT_ENABLE_PYTHON must be one of ${ON_OFF_AUTO}")
endif()
if (BOUT_ENABLE_PYTHON OR BOUT_ENABLE_PYTHON STREQUAL "AUTO")
  add_subdirectory(tools/pylib/_boutpp_build)
else()
  set(BOUT_ENABLE_PYTHON OFF)
endif()
set(BOUT_USE_PYTHON ${BOUT_ENABLE_PYTHON})

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
add_library(bout++::bout++ ALIAS bout++)
target_link_libraries(bout++ PUBLIC MPI::MPI_CXX)
target_include_directories(bout++ PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/include>
  $<INSTALL_INTERFACE:include>
  )
set(BOUT_LIB_PATH "${CMAKE_CURRENT_BINARY_DIR}/lib")
set_target_properties(bout++ PROPERTIES
  LIBRARY_OUTPUT_DIRECTORY "${BOUT_LIB_PATH}"
  ARCHIVE_OUTPUT_DIRECTORY "${BOUT_LIB_PATH}"
  SOVERSION 5.1.0)

# Set some variables for the bout-config script
set(CONFIG_LDFLAGS "${CONFIG_LDFLAGS} -L\$BOUT_LIB_PATH -lbout++")
set(BOUT_INCLUDE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/include")
set(CONFIG_CFLAGS "${CONFIG_CFLAGS} -I\${BOUT_INCLUDE_PATH} -I${CMAKE_CURRENT_BINARY_DIR}/include ${CMAKE_CXX_FLAGS}")

target_compile_features(bout++ PUBLIC cxx_std_14)
set_target_properties(bout++ PROPERTIES CXX_EXTENSIONS OFF)

# Optional compiler features
include(cmake/SetupCompilers.cmake)

# Optional dependencies
include(cmake/SetupBOUTThirdParty.cmake)

# Various sanitizers, including coverage and address sanitizer
include(cmake/Sanitizers.cmake)
enable_sanitizers(bout++)

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

option(BOUT_ENABLE_WARNINGS "Enable compiler warnings" ON)
if (BOUT_ENABLE_WARNINGS)
  target_compile_options(bout++ PRIVATE
    $<$<NOT:$<COMPILE_LANGUAGE:CUDA>>:
    $<$<OR:$<CXX_COMPILER_ID:GNU>,$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:AppleClang>>:
      -Wall -Wextra > >
    $<$<CXX_COMPILER_ID:MSVC>:
      /W4 >  
    $<$<COMPILE_LANGUAGE:CUDA>:-Xcompiler=-Wall -Xcompiler=-Wextra >
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
    message(FATAL_ERROR "addr2line not found. Disable backtrace by setting BOUT_ENABLE_BACKTRACE=Off")
  endif()
  target_link_libraries(bout++ PUBLIC ${CMAKE_DL_LIBS})
  set(CONFIG_LDFLAGS "${CONFIG_LDFLAGS} -l${CMAKE_DL_LIBS}")
endif()
message(STATUS "Enable backtrace: BOUT_USE_BACKTRACE=${BOUT_ENABLE_BACKTRACE}")
set(BOUT_USE_BACKTRACE ${BOUT_ENABLE_BACKTRACE})

option(BOUT_ENABLE_METRIC_3D "Enable 3D metric support" OFF)
if(BOUT_ENABLE_METRIC_3D)
  set(BOUT_METRIC_TYPE "3D")
else()
  set(BOUT_METRIC_TYPE "2D")
endif()
set(BOUT_USE_METRIC_3D ${BOUT_ENABLE_METRIC_3D})

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

include(BuildType)
# Here CMAKE_BUILD_TYPE is always set
string(TOUPPER "${CMAKE_BUILD_TYPE}" CMAKE_BUILD_TYPE_UPPER)
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

# Are we building BOUT++ directly, or as part of another project
string(COMPARE EQUAL
  "${PROJECT_NAME}" "${CMAKE_PROJECT_NAME}"
  PROJECT_IS_TOP_LEVEL
)
option(BOUT_TESTS "Build the tests" ${PROJECT_IS_TOP_LEVEL})
option(BOUT_ENABLE_ALL_TESTS "Enable running all of the tests, rather then the standard selection of fast tests" OFF)
if(BOUT_TESTS)
  enable_testing()
  # Targets for just building the tests
  # Tests need to add themselves as dependencies to these targets
  add_custom_target(build-check-unit-tests)
  add_custom_target(build-check-integrated-tests)
  add_custom_target(build-check-mms-tests)

  # Build all the tests
  add_custom_target(build-check)
  add_dependencies(build-check build-check-unit-tests build-check-integrated-tests build-check-mms-tests)

  add_subdirectory(tests/unit EXCLUDE_FROM_ALL)
  add_subdirectory(tests/integrated EXCLUDE_FROM_ALL)
  add_subdirectory(tests/MMS EXCLUDE_FROM_ALL)

  # Targets for running the tests
  if (BOUT_ENABLE_UNIT_TESTS)
    add_custom_target(check-unit-tests
      COMMAND ctest -R serial_tests --output-on-failure)
    add_dependencies(check-unit-tests build-check-unit-tests)
  endif()

  add_custom_target(check-integrated-tests
    COMMAND ctest -R "test-" --output-on-failure)
  add_dependencies(check-integrated-tests build-check-integrated-tests)

  add_custom_target(check-mms-tests
    COMMAND ctest -R "MMS-" --output-on-failure)
  add_dependencies(check-mms-tests build-check-mms-tests)

  # Run all the tests
  add_custom_target(check)
  add_dependencies(check check-integrated-tests check-mms-tests)
  if (BOUT_ENABLE_UNIT_TESTS)
    add_dependencies(check check-unit-tests)
  endif()

endif()

option(BOUT_BUILD_EXAMPLES "Build the examples" ON)
if(BOUT_BUILD_EXAMPLES)
  add_custom_target(build-all-examples)
  add_subdirectory(examples EXCLUDE_FROM_ALL)
endif()


##################################################
# L10N: localisation - include translations

find_package(Gettext)

if (GETTEXT_FOUND)
  #add_custom_target(mofiles ALL)
  set(bout_langs es de fr zh_CN zh_TW)

  foreach(_lang IN LISTS bout_langs)
    set(_gmoFile ${CMAKE_CURRENT_BINARY_DIR}/locale/${_lang}/libbout.gmo)
    set(_poFile ${CMAKE_CURRENT_SOURCE_DIR}/locale/${_lang}/libbout.po)
    add_custom_command(OUTPUT ${_gmoFile} _mo_file_${_lang}
      COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_CURRENT_BINARY_DIR}/locale/${_lang}/
      COMMAND ${GETTEXT_MSGFMT_EXECUTABLE} -o ${_gmoFile} ${_poFile}
      WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
      DEPENDS ${_poFile}
    )

    list(APPEND _gmoFiles ${_gmoFile})

    install(FILES ${_gmoFile} DESTINATION ${CMAKE_INSTALL_LOCALEDIR}/${_lang}/LC_MESSAGES/ RENAME libbout.mo)
  endforeach()
  add_custom_target(mofiles ALL
    DEPENDS ${_gmoFiles})
endif()


##################################################
# Documentation

option(BOUT_BUILD_DOCS "Build the documentation" OFF)
if (BOUT_BUILD_DOCS)
  add_subdirectory(manual EXCLUDE_FROM_ALL)
endif()


add_custom_target(dist
  COMMAND ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/tools/pylib/_boutpp_build/backend.py dist
  # there is no cmake equivalent to `mv` - so only works on systems that are not inentionally non-POSIX complient
  COMMAND mv BOUT++-v${BOUT_FULL_VERSION}.tar.gz ${CMAKE_CURRENT_BINARY_DIR}/
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  )
##################################################
# Generate the build config header

if (EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/include/bout/build_defines.hxx")
  # If we do in source builds, this is fine
  if (NOT ${CMAKE_CURRENT_SOURCE_DIR} STREQUAL ${CMAKE_CURRENT_BINARY_DIR})
    message(FATAL_ERROR "Generated build_defines.hxx header already exists; please remove '${CMAKE_CURRENT_SOURCE_DIR}/include/bout/build_defines.hxx' before continuing")
  endif()
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
  set(CONFIG_LDFLAGS "${CMAKE_SHARED_LIBRARY_RUNTIME_CXX_FLAG}\$BOUT_LIB_PATH -L\$BOUT_LIB_PATH -lbout++ -lfmt")
else()
  set(CONFIG_LDFLAGS "${CONFIG_LDFLAGS}")
endif()

set(ISINSTALLED "FALSE")

set(_CONFIG_LDFLAGS)
string(REPLACE " " ";" CONFIG_LDFLAGS_LIST ${CONFIG_LDFLAGS})
foreach (flag ${CONFIG_LDFLAGS_LIST})
  string(REGEX MATCH "^-.*$" isopt "${flag}")
  if (isopt)
    # message("${flag} is an option")
    set(_CONFIG_LDFLAGS "${_CONFIG_LDFLAGS} ${flag}")
    # All good
  else()
    if(EXISTS ${flag})
      # message("${flag} is a file")
      set(_CONFIG_LDFLAGS "${_CONFIG_LDFLAGS} ${flag}")
    else()
      string(FIND "${flag}" "::" hascolcol)
      if (${hascolcol} EQUAL -1)
	message("Fixing ${flag} to -l${flag}")
	set(_CONFIG_LDFLAGS "${_CONFIG_LDFLAGS} -l${flag}")
      else()
	string(REGEX MATCH "[^:]*$" flag2 "${flag}")
	message("Fixing ${flag} to -l${flag2}")
	set(_CONFIG_LDFLAGS "${_CONFIG_LDFLAGS} -l${flag2}")
      endif()
    endif()
  endif()
endforeach()
set( CONFIG_LDFLAGS ${_CONFIG_LDFLAGS})

# This version of the file allows the build directory to be used directly
configure_file(bin/bout-config.in bin/bout-config @ONLY)
configure_file(tools/pylib/boutconfig/__init__.py.cin tools/pylib/boutconfig/__init__.py @ONLY)
configure_file(bout++Config.cmake.in bout++Config.cmake @ONLY)

# We need to generate a separate version for installation, with the
# correct install paths. So first we need to replace the build
# directory library path with the installation path
string(REPLACE
  "${CMAKE_BINARY_DIR}/lib" "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}"
  CONFIG_LDFLAGS "${CONFIG_LDFLAGS}")
set(BOUT_LIB_PATH "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}")
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

set(PYTHONCONFIGPATH "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_PYTHON_SITEARCH}")
set(ISINSTALLED "TRUE")

# This version now has the correct paths to use the final installation
configure_file(bin/bout-config.in bin/bout-config-install @ONLY)
configure_file(tools/pylib/boutconfig/__init__.py.cin tools/pylib/boutconfig/__init__.py-install @ONLY)
configure_file(bout++Config.cmake.in bout++Config.cmake-install @ONLY)

##################################################
# Installation

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
  DESTINATION "${CMAKE_INSTALL_BINDIR}"
  REGEX "bout-squashoutput" EXCLUDE
  REGEX "bout-config\.in" EXCLUDE
  REGEX "bout-pylib-cmd-to-bin" EXCLUDE
)

# The installed version of bout-config needs renaming when we install
# it. Note this MUST be done after the installation of bin/, to make
# sure we clobber any versions of bout-config hanging around from an
# autotools build
install(PROGRAMS "${CMAKE_CURRENT_BINARY_DIR}/bin/bout-config-install"
  DESTINATION "${CMAKE_INSTALL_BINDIR}"
  RENAME "bout-config"
  )

install(
  FILES "${CMAKE_CURRENT_BINARY_DIR}/tools/pylib/boutconfig/__init__.py-install"
  DESTINATION "${CMAKE_INSTALL_PYTHON_SITEARCH}/boutconfig"
  RENAME "__init__.py"
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

# CMake configuration files
install(
  FILES
      "${CMAKE_CURRENT_BINARY_DIR}/bout++ConfigVersion.cmake"
      "${CMAKE_CURRENT_SOURCE_DIR}/cmake/BOUT++functions.cmake"
      "${CMAKE_CURRENT_SOURCE_DIR}/cmake/CorrectWindowsPaths.cmake"
      "${CMAKE_CURRENT_SOURCE_DIR}/cmake/FindClangFormat.cmake"
      "${CMAKE_CURRENT_SOURCE_DIR}/cmake/FindFFTW.cmake"
      "${CMAKE_CURRENT_SOURCE_DIR}/cmake/FindHYPRE.cmake"
      "${CMAKE_CURRENT_SOURCE_DIR}/cmake/FindnetCDF.cmake"
      "${CMAKE_CURRENT_SOURCE_DIR}/cmake/FindnetCDFCxx.cmake"
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
install(
  FILES
      "${CMAKE_CURRENT_BINARY_DIR}/bout++Config.cmake-install"
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/bout++"
  RENAME "bout++Config.cmake"
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
   RAJA enabled             : ${BOUT_HAS_RAJA}
   Umpire enabled           : ${BOUT_HAS_UMPIRE}
   Caliper enabled          : ${BOUT_HAS_CALIPER}
   CUDA enabled             : ${BOUT_HAS_CUDA}
   Metric type              : ${BOUT_METRIC_TYPE}
   Python API support       : ${BOUT_USE_PYTHON}
   Sanitizers enabled       : ${BOUT_USE_SANITIZERS}

   === Python ===

   Make sure that the tools/pylib directory is in your PYTHONPATH
   e.g. by adding to your ~/.bashrc file

       export PYTHONPATH=${BOUT_PYTHONPATH}:\$PYTHONPATH

*** Now run `cmake --build ${CMAKE_BINARY_DIR}` to compile BOUT++ ***
")
