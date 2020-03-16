cmake_minimum_required(VERSION 3.13)

project(advection-diffusion LANGUAGES CXX)

if (NOT TARGET bout++::bout++)
  find_package(bout++ REQUIRED)
endif()

add_executable(imex imex.cxx)
target_link_libraries(imex PRIVATE bout++::bout++)
