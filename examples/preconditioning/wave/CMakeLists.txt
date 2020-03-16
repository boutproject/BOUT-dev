cmake_minimum_required(VERSION 3.13)

project(preconditioning-wave LANGUAGES CXX)

if (NOT TARGET bout++::bout++)
  find_package(bout++ REQUIRED)
endif()

add_executable(preconditioning-wave test_precon.cxx)
target_link_libraries(preconditioning-wave PRIVATE bout++::bout++)
