cmake_minimum_required(VERSION 3.13)

project(diffusion-nl LANGUAGES CXX)

if (NOT TARGET bout++::bout++)
  find_package(bout++ REQUIRED)
endif()

add_executable(diffusion-nl diffusion-nl.cxx)
target_link_libraries(diffusion-nl PRIVATE bout++::bout++)
