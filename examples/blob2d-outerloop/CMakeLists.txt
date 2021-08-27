cmake_minimum_required(VERSION 3.13)


project(blob2d-outerloop LANGUAGES CXX C)

set(bd_sources
   blob2d.cxx
)

set(bd_depends bout++ PETSc)
if(BOUT_HAS_RAJA)
   list(APPEND bd_depends RAJA)
   find_package(RAJA REQUIRED)
endif ()

if(BOUT_HAS_UMPIRE)
   list(APPEND bd_depends umpire)
   find_package(umpire REQUIRED)
endif ()

add_executable(blob2d-outerloop ${bd_sources})
target_link_libraries(blob2d-outerloop PRIVATE bout++::bout++)
target_include_directories(blob2d-outerloop PUBLIC
      ${PROJECT_SOURCE_DIR}/include
      ${PROJECT_SOURCE_DIR}/externalpackages/fmt/include
      ${PROJECT_SOURCE_DIR}/externalpackages/mpark.variant/include
      )
if(BOUT_USE_CUDA)
   list(APPEND bd_depends cuda)
   set_source_files_properties(${bd_sources} PROPERTIES LANGUAGE CUDA )
   message(STATUS "blob2d-outerloop needs CUDA LIBRARY DIR: ${BOUT_CUDA_LIB_DIR}")
   target_include_directories(bd PUBLIC ${BOUT_CUDA_LIB_DIR}/../include )

   target_compile_definitions(bd PUBLIC "BOUT_USE_CUDA")
endif ()

if(BOUT_HAS_RAJA)
   target_compile_definitions(bd PUBLIC "BOUT_HAS_RAJA")
endif ()

message(STATUS "bd_depends: ${bd_depends}")
