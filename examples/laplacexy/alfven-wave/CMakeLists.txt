cmake_minimum_required(VERSION 3.13)

project(laplacexy-alfven-wave LANGUAGES CXX)

if (NOT TARGET bout++::bout++)
  find_package(bout++ REQUIRED)
endif()

bout_add_example(laplacexy-alfven-wave
  SOURCES alfven.cxx
  DATA_DIRS cbm18 data
  EXTRA_FILES cbm18_dens8.grid_nx68ny64.nc d3d_119919.nc)
