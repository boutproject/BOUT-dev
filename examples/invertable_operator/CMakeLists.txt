cmake_minimum_required(VERSION 3.13)

project(invertable_operator LANGUAGES CXX)

if (NOT TARGET bout++::bout++)
  find_package(bout++ REQUIRED)
endif()

bout_add_example(invertable_operator
  SOURCES invertable_operator.cxx
  REQUIRES BOUT_HAS_PETSC)
