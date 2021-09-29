cmake_minimum_required(VERSION 3.13)

project(tokamak-2fluid LANGUAGES CXX)

if (NOT TARGET bout++::bout++)
  find_package(bout++ REQUIRED)
endif()

bout_add_example(tokamak-2fluid
  SOURCES 2fluid.cxx
  DATA_DIRS d3d-119919
            d3d-129131
            data
  EXTRA_FILES uedge.grd.nc
              uedge.grd_129131_newpproc.nc
              data/uedge.grd.nc)
