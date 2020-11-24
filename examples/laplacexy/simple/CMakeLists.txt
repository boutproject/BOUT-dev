cmake_minimum_required(VERSION 3.13)

cmake_minimum_required(VERSION 3.13)

project(test_laplacexy LANGUAGES CXX C)

set(test_laplacexy_sources
   test-laplacexy.cxx
)

set(test_laplacexy_depends bout++ PETSc)

add_executable(test_laplacexy ${test_laplacexy_sources})
target_link_libraries(test_laplacexy PRIVATE bout++::bout++)
target_include_directories(test_laplacexy PUBLIC 
      ${PROJECT_SOURCE_DIR}/include
      ${PROJECT_SOURCE_DIR}/externalpackages/fmt/include
      ${PROJECT_SOURCE_DIR}/externalpackages/mpark.variant/include
      )

message(STATUS "test_laplacexy_depends: ${test_laplacexy_depends}")


