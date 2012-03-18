#
# $Id: ConfigCMake.cmake 9131 2011-09-27 01:24:06Z remko $
#
# Useful CMake variables.
#
# There are three configuration files:
#   1) "ConfigDefault.cmake" - is version controlled and used to add new default
#      variables and set defaults for everyone.
#   2) "ConfigUser.cmake" in the source tree - is not version controlled
#      (currently listed in svn:ignore property) and used to override defaults on
#      a per-user basis.
#   3) "ConfigUser.cmake" in the build tree - is used to override
#      "ConfigUser.cmake" in the source tree.
#
# NOTE: If you want to change CMake behaviour just for yourself then copy
#      "ConfigUserTemplate.cmake" to "ConfigUser.cmake" and then edit
#      "ConfigUser.cmake" (not "ConfigDefault.cmake" or "ConfigUserTemplate.cmake").
#
include ("${CMAKE_SOURCE_DIR}/cmake/ConfigDefault.cmake")

# If "ConfigUser.cmake" doesn't exist then create one for convenience.
if (EXISTS "${CMAKE_SOURCE_DIR}/cmake/ConfigUser.cmake")
	include ("${CMAKE_SOURCE_DIR}/cmake/ConfigUser.cmake")
endif (EXISTS "${CMAKE_SOURCE_DIR}/cmake/ConfigUser.cmake")

# If you've got a 'ConfigUser.cmake' in the build tree then that overrides the
# one in the source tree.
if (EXISTS "${CMAKE_BINARY_DIR}/cmake/ConfigUser.cmake")
	include ("${CMAKE_BINARY_DIR}/cmake/ConfigUser.cmake")
endif (EXISTS "${CMAKE_BINARY_DIR}/cmake/ConfigUser.cmake")

###########################################################
# Do any needed processing of the configuration variables #
###########################################################

# Build type
if (NOT CMAKE_BUILD_TYPE)
	set (CMAKE_BUILD_TYPE Release)
endif (NOT CMAKE_BUILD_TYPE)

if (CMAKE_BUILD_TYPE MATCHES "Debug|RelWithDebInfo")
	set (DEBUG_BUILD TRUE)
endif (CMAKE_BUILD_TYPE MATCHES "Debug|RelWithDebInfo")
