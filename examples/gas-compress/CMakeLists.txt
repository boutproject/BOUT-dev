cmake_minimum_required(VERSION 3.13)

project(gas-compress LANGUAGES CXX)

if (NOT TARGET bout++::bout++)
  find_package(bout++ REQUIRED)
endif()

bout_add_example(gas-compress
  SOURCES gas_compress.cxx gas_compress.hxx
  DATA_DIRS rayleigh-taylor sod-shock
  EXTRA_FILES rt.grd.nc sod.grd.nc)
