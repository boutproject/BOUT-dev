cmake_minimum_required(VERSION 3.13)

project(jorek-compare LANGUAGES CXX)

if (NOT TARGET bout++::bout++)
  find_package(bout++ REQUIRED)
endif()

add_executable(jorek-compare jorek_compare.cxx)
target_link_libraries(jorek-compare PRIVATE bout++::bout++)
