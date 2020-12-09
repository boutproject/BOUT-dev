cmake_minimum_required(VERSION 3.13)

project(invertable_operator LANGUAGES CXX)

if (NOT TARGET bout++::bout++)
  find_package(bout++ REQUIRED)
endif()

if(NOT BOUT_HAS_PETSC)
  if ("${PROJECT_SOURCE_DIR}" STREQUAL "${CMAKE_CURRENT_SOURCE_DIR}")
    message(STATUS "The invertable operator example requires PETSc. Please compile BOUT++ with PETSc to use it")
    return()
  else()
    message(FATAL_ERROR "The invertable operator example requires PETSc. Please compile BOUT++ with PETSc to use it")
  endif()
endif()

bout_add_example(invertable_operator SOURCES invertable_operator.cxx)
