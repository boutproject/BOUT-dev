cmake_minimum_required(VERSION 3.13)

project(conduction LANGUAGES CXX)

if (NOT TARGET bout++::bout++)
  find_package(bout++ REQUIRED)
endif()

bout_add_example(conduction
  SOURCES conduction.cxx
  DATA_DIRS data fromfile
  EXTRA_FILES generate.py)

bout_add_example(dump_coords
  SOURCES dump_coords.cxx)

bout_add_example(old_coord_dump
  SOURCES old_coord_dump.cxx)
