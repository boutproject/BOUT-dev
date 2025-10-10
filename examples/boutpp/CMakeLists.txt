cmake_minimum_required(VERSION 3.13)

if (NOT TARGET bout++::bout++)
  find_package(bout++ REQUIRED)
endif()

bout_copy_file(runexample)
bout_copy_file(blob2d.py)
bout_copy_file(simulation.py)
bout_copy_file(data/BOUT.inp)
