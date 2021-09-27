cmake_minimum_required(VERSION 3.13)

project(wave-slab LANGUAGES CXX)

if (NOT TARGET bout++::bout++)
  find_package(bout++ REQUIRED)
endif()

bout_add_example(wave-slab
  SOURCES wave_slab.cxx
  EXTRA_FILES generate.py)
