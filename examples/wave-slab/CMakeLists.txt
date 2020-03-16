cmake_minimum_required(VERSION 3.13)

project(wave-slab LANGUAGES CXX)

if (NOT TARGET bout++::bout++)
  find_package(bout++ REQUIRED)
endif()

add_executable(wave-slab wave_slab.cxx)
target_link_libraries(wave-slab PRIVATE bout++::bout++)
