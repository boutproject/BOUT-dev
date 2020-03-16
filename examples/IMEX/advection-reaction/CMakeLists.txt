cmake_minimum_required(VERSION 3.13)

project(advection-reaction LANGUAGES CXX)

if (NOT TARGET bout++::bout++)
  find_package(bout++ REQUIRED)
endif()

add_executable(split_operator split_operator.cxx)
target_link_libraries(split_operator PRIVATE bout++::bout++)
