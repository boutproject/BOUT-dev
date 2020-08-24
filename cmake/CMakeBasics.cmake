include(cmake/thirdparty/SetupBOUTThirdParty.cmake)

include(cmake/SetupCompilers.cmake)

# This needs to happen last, once all we have defined options and checked for
# third-party packages.
include(cmake/CMakeConfigureFile.cmake)
