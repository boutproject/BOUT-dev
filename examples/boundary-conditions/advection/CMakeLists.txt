cmake_minimum_required(VERSION 3.13)

project(advection LANGUAGES CXX)

if (NOT TARGET bout++::bout++)
  find_package(bout++ REQUIRED)
endif()

add_executable(advection advection.cxx)
target_link_libraries(advection PRIVATE bout++::bout++)
