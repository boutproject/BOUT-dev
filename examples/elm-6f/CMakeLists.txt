cmake_minimum_required(VERSION 3.13)

project(elm-6f LANGUAGES CXX)

if (NOT TARGET bout++::bout++)
  find_package(bout++ REQUIRED)
endif()

bout_add_example(elm_6f_v0
  SOURCES elm_6f_v0.cxx
  EXTRA_FILES cbm18_dens8.grid_nx68ny64.nc
              cbm18_dens6-0.5BS.grid.nc
	      job_script.sh
	      data/BOUT.inp)
