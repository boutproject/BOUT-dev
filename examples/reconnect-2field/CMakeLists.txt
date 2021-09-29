cmake_minimum_required(VERSION 3.13)

project(reconnect-2field LANGUAGES CXX)

if (NOT TARGET bout++::bout++)
  find_package(bout++ REQUIRED)
endif()

bout_add_example(reconnect-2field
  SOURCES 2field.cxx
  EXTRA_FILES slab_68x32.nc)
