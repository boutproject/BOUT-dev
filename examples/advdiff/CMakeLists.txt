cmake_minimum_required(VERSION 3.13)

project(advdiff LANGUAGES CXX)

if (NOT TARGET bout++::bout++)
  find_package(bout++ REQUIRED)
endif()

add_executable(advdiff advdiff.cxx)
target_link_libraries(advdiff PRIVATE bout++::bout++)
