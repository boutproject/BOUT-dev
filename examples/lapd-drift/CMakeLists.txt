cmake_minimum_required(VERSION 3.13)

project(lapd-drift LANGUAGES CXX)

if (NOT TARGET bout++::bout++)
  find_package(bout++ REQUIRED)
endif()

bout_add_example(lapd-drift
  SOURCES lapd_drift.cxx
  DATA_DIRS lapd
            data
            pisces
  EXTRA_FILES BOUT.inp
              BOUT.inp.nn
              BOUT.inp.nn_zem
              uedge.grd.nc)
