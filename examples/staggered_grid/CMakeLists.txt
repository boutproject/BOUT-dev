cmake_minimum_required(VERSION 3.13)

project(staggered_grid LANGUAGES CXX)

if (NOT TARGET bout++::bout++)
  find_package(bout++ REQUIRED)
endif()

bout_add_example(staggered_grid
  SOURCES test_staggered.cxx
  EXTRA_FILES generate.py runandplot test-staggered.nc
  DATA_DIRS data test)
