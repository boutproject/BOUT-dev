cmake_minimum_required(VERSION 3.13)

project(advdiff LANGUAGES CXX)

if (NOT TARGET bout++::bout++)
  find_package(bout++ REQUIRED)
endif()

bout_add_example(advdiff
  SOURCES advdiff.cxx
  EXTRA_FILES slab.grd.nc)
