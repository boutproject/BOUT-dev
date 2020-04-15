cmake_minimum_required(VERSION 3.13)

project(advdiff2 LANGUAGES CXX)

if (NOT TARGET bout++::bout++)
  find_package(bout++ REQUIRED)
endif()

add_executable(advdiff2 
  globals.cxx
  globals.hxx
  header.hxx
  init.cxx
  run.cxx)
target_link_libraries(advdiff2 PRIVATE bout++::bout++)
