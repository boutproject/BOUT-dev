cmake_minimum_required(VERSION 3.13)

project(mhd LANGUAGES CXX C)

set(mhd_sources
   mhd.cxx
)

set(mhd_depends bout++ PETSc)

if(BOUT_HAS_RAJA)
   list(APPEND mhd_depends RAJA)
   find_package(RAJA REQUIRED)
endif ()

if(BOUT_HAS_UMPIRE)
   list(APPEND mhd_depends umpire)
   find_package(umpire REQUIRED)
endif ()

add_executable(mhd ${mhd_sources})
target_link_libraries(mhd PRIVATE bout++::bout++)
target_include_directories(mhd PUBLIC
      ${PROJECT_SOURCE_DIR}/include
      ${PROJECT_SOURCE_DIR}/externalpackages/fmt/include
      ${PROJECT_SOURCE_DIR}/externalpackages/mpark.variant/include
      )
if(BOUT_USE_CUDA)
   list(APPEND mhd_depends cuda)
   set_source_files_properties(${mhd_sources} PROPERTIES LANGUAGE CUDA )
   message(STATUS "mhd needs CUDA LIBRARY DIR: ${BOUT_CUDA_LIB_DIR}")
   target_include_directories(mhd PUBLIC ${BOUT_CUDA_LIB_DIR}/../include )
   target_compile_definitions(mhd PUBLIC "BOUT_USE_CUDA")
endif ()

if(BOUT_HAS_RAJA)
   target_compile_definitions(mhd PUBLIC "BOUT_HAS_RAJA")
endif ()

message(STATUS "mhd_depends: ${mhd_depends}")
