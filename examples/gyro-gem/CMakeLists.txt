cmake_minimum_required(VERSION 3.13)

project(gyro-gem LANGUAGES CXX)

if (NOT TARGET bout++::bout++)
  find_package(bout++ REQUIRED)
endif()

bout_add_example(gyro-gem
  SOURCES gem.cxx
  EXTRA_FILES cyclone_68x32.nc)
