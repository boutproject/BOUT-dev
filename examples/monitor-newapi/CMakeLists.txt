cmake_minimum_required(VERSION 3.13)

project(monitor-newapi LANGUAGES CXX)

if (NOT TARGET bout++::bout++)
  find_package(bout++ REQUIRED)
endif()

bout_add_example(monitor-newapi SOURCES monitor.cxx)
