cmake_minimum_required(VERSION 3.13)

project(gravity_reduced LANGUAGES CXX)

if (NOT TARGET bout++::bout++)
  find_package(bout++ REQUIRED)
endif()

bout_add_example(gravity_reduced
  SOURCES gravity_reduced.cxx
  EXTRA_FILES slab_grid.nc)
