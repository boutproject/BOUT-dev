cmake_minimum_required(VERSION 3.13)

project(finite-volume-fluid LANGUAGES CXX)

if (NOT TARGET bout++::bout++)
  find_package(bout++ REQUIRED)
endif()

bout_add_example(finite-volume-fluid
  SOURCES fluid.cxx
  DATA_DIRS data mms
  EXTRA_FILES mms.py)
