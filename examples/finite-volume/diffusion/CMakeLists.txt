cmake_minimum_required(VERSION 3.13)

project(finite-volume-diffusion LANGUAGES CXX)

if (NOT TARGET bout++::bout++)
  find_package(bout++ REQUIRED)
endif()

bout_add_example(finite-volume-diffusion
  SOURCES diffusion.cxx
  EXTRA_FILES mms.py)
