cmake_minimum_required(VERSION 3.13)

project(conduction LANGUAGES CXX)

if (NOT TARGET bout++::bout++)
  find_package(bout++ REQUIRED)
endif()

bout_add_example(conduction
  SOURCES conduction.cxx
  DATA_DIRS data fromfile
  EXTRA_FILES generate.py)
