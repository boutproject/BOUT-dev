cmake_minimum_required(VERSION 3.13)

project(2Dturbulence_multigrid LANGUAGES CXX)

if (NOT TARGET bout++::bout++)
  find_package(bout++ REQUIRED)
endif()

add_executable(esel esel.cxx)
target_link_libraries(esel PRIVATE bout++::bout++)
