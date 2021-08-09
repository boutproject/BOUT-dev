# Find Clang format
#
# Taken from https://github.com/ttroy50/cmake-examples commit 64bd54a
# This file is under MIT Licence

if (NOT ClangFormat_BIN_NAME)
  set(ClangFormat_BIN_NAME clang-format)
endif()

# if custom path check there first
if (ClangFormat_ROOT_DIR)
  find_program(ClangFormat_BIN
    NAMES
    ${ClangFormat_BIN_NAME}
    PATHS
    "${ClangFormat_ROOT_DIR}"
    NO_DEFAULT_PATH)
endif()

find_program(ClangFormat_BIN NAMES ${ClangFormat_BIN_NAME})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(
  ClangFormat
  DEFAULT_MSG
  ClangFormat_BIN)

mark_as_advanced(
  ClangFormat_BIN)
