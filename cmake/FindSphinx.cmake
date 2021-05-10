# FindSphinx
# ----------
#
# Find the Sphinx documentation generator
#
# This module will define the following variables:
#
# ::
#
#   Sphinx_FOUND - true if Sphinx was found
#   Sphinx_EXECUTABLE - Path to the ``sphinx-build`` executable

# Taken from
# https://devblogs.microsoft.com/cppblog/clear-functional-c-documentation-with-sphinx-breathe-doxygen-cmake/

#Look for an executable called sphinx-build
find_program(SPHINX_EXECUTABLE
             NAMES sphinx-build sphinx-build-3
             DOC "Path to sphinx-build executable")

include(FindPackageHandleStandardArgs)

#Handle standard arguments to find_package like REQUIRED and QUIET
find_package_handle_standard_args(Sphinx
                                  "Failed to find sphinx-build executable"
                                  SPHINX_EXECUTABLE)
