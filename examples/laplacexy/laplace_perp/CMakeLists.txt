cmake_minimum_required(VERSION 3.13)

project(laplacexy-laplace_perp LANGUAGES CXX)

if (NOT TARGET bout++::bout++)
  find_package(bout++ REQUIRED)
endif()

add_executable(laplacexy-laplace_perp test.cxx)
target_link_libraries(laplacexy-laplace_perp PRIVATE bout++::bout++)
