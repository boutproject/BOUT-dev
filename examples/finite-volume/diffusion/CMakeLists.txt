cmake_minimum_required(VERSION 3.13)

project(finite-volume-diffusion LANGUAGES CXX)

if (NOT TARGET bout++::bout++)
  find_package(bout++ REQUIRED)
endif()

add_executable(finite-volume-diffusion diffusion.cxx)
target_link_libraries(finite-volume-diffusion PRIVATE bout++::bout++)
