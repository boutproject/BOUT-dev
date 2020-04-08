include(CheckIncludeFiles)
include(CheckTypeSize)
include(CheckFunctionExists)
include(CheckCXXSourceCompiles)

if (CMAKE_BUILD_TYPE STREQUAL "Debug")
  set (DEBUG_INITIALIZE_UNDEFINED On)
  set (DEBUG_CHECK_ASSERTIONS On)
  set (DEBUG_CHECK_DIM_ASSERTIONS On)
endif()

#HAVE_CMATH
check_include_files("math.h" HAVE_CMATH)

#HAVE_CMATH_ISNAN
check_function_exists(std::isnan "cmath" HAVE_CMATH_ISNAN)

#HAVE_CTIME
check_include_files("time.h" HAVE_CTIME)

#HAVE_EXCEPTION_HANDLING

#HAVE_INLINE_ISNAND
check_function_exists(__inline_isnand "math.h" HAVE_INLINE_ISNAND)

#HAVE_ISNAN
check_cxx_source_compiles("int test = std::isnan(0.0)" HAVE_ISNAN)

#HAVE_ISNAND

#HAVE_MALLINFO
check_function_exists(mallinfo HAVE_MALLINFO)

#HAVE_MALLOC_H
check_include_files(malloc.h HAVE_MALLOC_H)

#HAVE_SYS_TIMES_H
check_include_files(sys/times.h HAVE_SYS_TIMES_H)

#HAVE_UNISTD_H
check_include_files(unistd.h HAVE_UNISTD_H)

#IOMANIP_HEADER_FILE
set(IOSTREAM_HEADER_FILE "<iomanip>")

#IOSTREAM_HEADER_FILE
set(IOSTREAM_HEADER_FILE "<iostream>")

#LACKS_SSTREAM
check_include_files(sstream HAVE_SSTREAM_H)

#LACKS_TEMPLATE_COMPLEX

#OPT_BUILD


#OSTRINGSTREAM_TYPE_IS_BROKEN

#OSTRSTREAM_TYPE_IS_BROKEN

#STL_SSTREAM_HEADER_FILE

configure_file(${PROJECT_SOURCE_DIR}/config/BOUT_config.h.cmake.in ${CMAKE_BINARY_DIR}/include/BOUT_config.h)

install(FILES ${CMAKE_BINARY_DIR}/include/BOUT_config.h
  DESTINATION ${CMAKE_INSTALL_PREFIX}/include)

