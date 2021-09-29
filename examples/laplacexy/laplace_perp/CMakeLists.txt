cmake_minimum_required(VERSION 3.13)

project(laplacexy-laplace_perp LANGUAGES CXX)

if (NOT TARGET bout++::bout++)
  find_package(bout++ REQUIRED)
endif()

bout_add_example(laplacexy-laplace_perp
  SOURCES test.cxx
  EXTRA_FILES cbm18_dens8.grid_nx68ny64.nc 
  DATA_DIRS square torus)
