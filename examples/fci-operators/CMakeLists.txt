cmake_minimum_required(VERSION 3.13)

project(fci-operators LANGUAGES CXX)

if(NOT TARGET bout++::bout++)
  find_package(bout++ REQUIRED)
endif()

bout_add_example(
  fci-operators
  SOURCES fci_operators_example.cxx
  DATA_DIRS data
  EXTRA_FILES rotating-ellipse-20x8x20.fci.nc makeplots.py
)
