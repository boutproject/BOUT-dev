cmake_minimum_required(VERSION 3.13)

project(conduction-snb LANGUAGES CXX)

if (NOT TARGET bout++::bout++)
  find_package(bout++ REQUIRED)
endif()

bout_add_example(conduction-snb
  SOURCES conduction-snb.cxx
  EXTRA_FILES fit_temperature.py sinusoid.py snb.csv spitzer-harm.csv step.py temperature.csv vfp.csv
  DATA_DIRS data step)
