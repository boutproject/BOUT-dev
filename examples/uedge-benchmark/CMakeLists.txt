cmake_minimum_required(VERSION 3.13)

project(uedge-benchmark LANGUAGES CXX)

if (NOT TARGET bout++::bout++)
  find_package(bout++ REQUIRED)
endif()

bout_add_example(uedge-benchmark
  SOURCES ue_bmark.cxx
  EXTRA_FILES uedge.grd_Up_Ni_Tei_2d.nc)
