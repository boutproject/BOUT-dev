cmake_minimum_required(VERSION 3.13)

project(invertable_operator LANGUAGES CXX)

if (NOT TARGET bout++::bout++)
  find_package(bout++ REQUIRED)
endif()

if(NOT BOUT_HAS_PETSC)
  message(FATAL_ERROR "This example requires PETSc. Please compile BOUT++ with PETSc")
endif()
add_executable(invertable_operator invertable_operator.cxx)
target_link_libraries(invertable_operator PRIVATE bout++::bout++)
