cmake_minimum_required(VERSION 3.13)

set(BLT_CXX_STD "c++14" CACHE STRING "Version of C++ standard")

set(mdh_sources
   mhd.cxx
)

set(mdh_depends bout++ PETSc)

if(HAVE_CUDA)
   list(APPEND mdh_depends cuda nvToolsExt)
endif ()

if(HAVE_RAJA)
   list(APPEND mdh_depends RAJA)
endif ()

message(STATUS "CUDA_TOOLKIT_ROOT_DIR : ${CUDA_TOOLKIT_ROOT_DIR}")
message(STATUS "mdh_depends: ${mdh_depends}")
blt_add_executable(
   NAME mdh
   SOURCES ${mdh_sources}
   INCLUDES 
     ${PROJECT_SOURCE_DIR}/include
     ${PROJECT_SOURCE_DIR}/externalpackages/fmt/include
     ${PROJECT_SOURCE_DIR}/externalpackages/mpark.variant/include
     ${CUDA_TOOLKIT_ROOT_DIR}/include
   DEFINES 
      "BOUT_VERSION_STRING=\"${BOUT_FULL_VERSION}\""
      "BOUT_VERSION_DOUBLE=${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}${PROJECT_VERSION_PATCH}"
   DEPENDS_ON
     ${mdh_depends})

