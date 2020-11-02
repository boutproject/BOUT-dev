<<<<<<< HEAD
#add_subdirectory("${PROJECT_SOURCE_DIR}/externalpackages/googletest"
#  "externalpackages/googletest")
=======
#bout_update_submodules()
#add_subdirectory("${PROJECT_SOURCE_DIR}/externalpackages/googletest"
#  "externalpackages/googletest")

if(NOT TARGET gtest)
  message(FATAL_ERROR "googletest not found! Have you disabled the git submodules (GIT_SUBMODULE)?")
endif()
>>>>>>> next-outerloop-GPU-umpire

mark_as_advanced(
  BUILD_GMOCK BUILD_GTEST BUILD_SHARED_LIBS
  gmock_build_tests gtest_build_samples gtest_build_tests
  gtest_disable_pthreads gtest_force_shared_crt gtest_hide_internal_symbols
  )

set(serial_tests_sources
  ./bout_test_main.cxx
  ./field/test_field.cxx
  ./field/test_field2d.cxx
  ./field/test_field3d.cxx
  ./field/test_field_factory.cxx
  ./field/test_fieldgroup.cxx
  ./field/test_fieldperp.cxx
  ./field/test_initialprofiles.cxx
  ./field/test_vector2d.cxx
  ./field/test_vector3d.cxx
  ./field/test_where.cxx
  ./include/bout/test_array.cxx
  ./include/bout/test_assert.cxx
  ./include/bout/test_deriv_store.cxx
  ./include/bout/test_generic_factory.cxx
  ./include/bout/test_hypre_interface.cxx
  ./include/bout/test_macro_for_each.cxx
  ./include/bout/test_monitor.cxx
  # ./include/bout/test_petsc_indexer.cxx
  ./include/bout/test_petsc_matrix.cxx
  ./include/bout/test_petsc_setters.cxx
  ./include/bout/test_petsc_vector.cxx
  ./include/bout/test_region.cxx
  ./include/bout/test_template_combinations.cxx
  ./include/test_cyclic_reduction.cxx
  ./include/test_derivs.cxx
  ./include/test_mask.cxx
  ./invert/test_fft.cxx
  ./mesh/data/test_gridfromoptions.cxx
  ./mesh/parallel/test_shiftedmetric.cxx
  ./mesh/test_boundary_factory.cxx
  ./mesh/test_boutmesh.cxx
  ./mesh/test_coordinates.cxx
  ./mesh/test_interpolation.cxx
  ./mesh/test_mesh.cxx
  ./mesh/test_paralleltransform.cxx
  ./solver/test_fakesolver.cxx
  ./solver/test_fakesolver.hxx
  ./solver/test_solver.cxx
  ./solver/test_solverfactory.cxx
  ./sys/test_boutexception.cxx
  ./sys/test_expressionparser.cxx
  ./sys/test_msg_stack.cxx
  ./sys/test_options.cxx
  ./sys/test_options_fields.cxx
  ./sys/test_options_netcdf.cxx
  ./sys/test_optionsreader.cxx
  ./sys/test_output.cxx
  ./sys/test_range.cxx
  ./sys/test_timer.cxx
  ./sys/test_type_name.cxx
  ./sys/test_utils.cxx
  ./sys/test_variant.cxx
  ./sys/test_raja.cxx
  ./test_extras.cxx
  ./test_extras.hxx
  ./src/test_bout++.cxx)

set(serial_tests_depends bout++ gtest gmock)
if (ENABLE_OPENMP)
  set(serial_tests_depends ${serial_tests_depends} openmp)
  set(OPENMP_SCHEDULE static CACHE STRING "Set OpenMP schedule")
endif ()

if (ENABLE_CUDA)
  set(serial_tests_depends ${serial_tests_depends} cuda)
endif ()

message(STATUS "serial_tests depends ${serial_tests_depends}")
blt_add_executable(
  NAME serial_tests
  SOURCES ${serial_tests_sources}
  DEPENDS_ON ${serial_tests_depends})


target_include_directories(serial_tests PUBLIC .
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../../include>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../../externalpackages/fmt/include>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../../externalpackages/mpark.variant/include>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/../../include>
  $<INSTALL_INTERFACE:include>
  )

<<<<<<< HEAD
set(serial_tests_depends bout++ gtest gmock PETSc HYPRE)
set(HYPRE_USING_CUDA_DEFINE "")
set(HYPRE_USING_CUBLAS_DEFINE "")
set(MEMORY_DEFINE "")
if (ENABLE_OPENMP)
  set(serial_tests_depends ${serial_tests_depends} openmp)
  set(OPENMP_SCHEDULE static CACHE STRING "Set OpenMP schedule")
endif ()

if (ENABLE_CUDA)
  set(serial_tests_depends ${serial_tests_depends} cuda)
  set(HYPRE_USING_CUDA_DEFINE "HYPRE_USING_CUDA")
  set(MEMORY_DEFINE "HYPRE_USING_UNIFIED_MEMORY")
endif ()

message(STATUS "serial_tests depends ${serial_tests_depends}")
blt_add_executable(
  NAME serial_tests
  SOURCES ${serial_tests_sources}
   INCLUDES 
     ${PROJECT_SOURCE_DIR}/include
     ${PROJECT_SOURCE_DIR}/externalpackages/fmt/include
     ${PROJECT_SOURCE_DIR}/externalpackages/mpark.variant/include
     ${CUDA_TOOLKIT_ROOT_DIR}/include
   DEFINES 
      "BOUT_VERSION_STRING=\"${BOUT_FULL_VERSION}\""
      "BOUT_VERSION_DOUBLE=${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}${PROJECT_VERSION_PATCH}"
      BOUT_HAS_HYPRE
      ${HYPRE_USING_CUDA_DEFINE}
      ${HYPRE_USING_CUBLAS_DEFINE}
      ${MEMORY_DEFINE}
  DEPENDS_ON ${serial_tests_depends})


target_include_directories(serial_tests PUBLIC .
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../../include>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../../externalpackages/fmt/include>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../../externalpackages/mpark.variant/include>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/../../include>
  $<INSTALL_INTERFACE:include>
  )

=======
>>>>>>> next-outerloop-GPU-umpire
target_compile_definitions(bout++
    PUBLIC "OPENMP_SCHEDULE=${OPENMP_SCHEDULE}"
    PUBLIC "BOUT_OPENMP_SCHEDULE=${OPENMP_SCHEDULE}")
set_target_properties(serial_tests PROPERTIES LINKER_LANGUAGE CXX)
#target_link_libraries(serial_tests gtest gmock bout++::bout++)
#add_test(NAME serial_tests COMMAND serial_tests)
set_target_properties(serial_tests PROPERTIES FOLDER tests/unit)
add_dependencies(build-check-unit-tests serial_tests)
<<<<<<< HEAD
target_link_libraries(
   serial_tests PUBLIC 
   ${CUDA_TOOLKIT_ROOT_DIR}/lib64/libcusparse_static.a
   ${CUDA_TOOLKIT_ROOT_DIR}/lib64/libcurand_static.a
   ${CUDA_TOOLKIT_ROOT_DIR}/lib64/libculibos.a
   ${CUDA_TOOLKIT_ROOT_DIR}/lib64/libcublas_static.a
   ${CUDA_TOOLKIT_ROOT_DIR}/lib64/libcublasLt_static.a
)



blt_add_test(
  NAME serial_tests
  COMMAND serial_tests)

#target_include_directories(serial_tests PUBLIC .)
#target_link_libraries(serial_tests gtest gmock bout++::bout++)
#add_test(NAME serial_tests COMMAND serial_tests)
#set_target_properties(serial_tests PROPERTIES FOLDER tests/unit)
#add_dependencies(build-check-unit-tests serial_tests)
=======
blt_add_test(
  NAME serial_tests
  COMMAND serial_tests)
>>>>>>> next-outerloop-GPU-umpire
