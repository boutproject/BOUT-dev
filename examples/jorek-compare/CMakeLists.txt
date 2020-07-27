cmake_minimum_required(VERSION 3.13)

project(jorek-compare LANGUAGES CXX)

if (NOT TARGET bout++::bout++)
  find_package(bout++ REQUIRED)
endif()

bout_add_example(jorek-compare
  SOURCES jorek_compare.cxx
  EXTRA_FILES d3d_119919.nc)
