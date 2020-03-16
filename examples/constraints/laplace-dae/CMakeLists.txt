cmake_minimum_required(VERSION 3.13)

project(constraints-laplace-dae LANGUAGES CXX)

if (NOT TARGET bout++::bout++)
  find_package(bout++ REQUIRED)
endif()

add_executable(constraints-laplace-dae laplace_dae.cxx)
target_link_libraries(constraints-laplace-dae PRIVATE bout++::bout++)
