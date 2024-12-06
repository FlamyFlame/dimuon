#
# This is a template for a CMakeLists.txt file that can be used in a client
# project (work area) to set up building ATLAS packages against the configured
# release.
#

# Set the minimum required CMake version:
cmake_minimum_required( VERSION 3.15 FATAL_ERROR )

project( WorkDir VERSION 1.0.0 ) # could add: LANGUAGES CXX Fortran

#register the package_filter file as a dependent
set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS $ENV{ACM_PACKAGE_FILTER_FILE} )
set( ATLAS_PACKAGE_FILTER_FILE $ENV{ACM_PACKAGE_FILTER_FILE} )

# If there's a directory called AtlasCMake in the project,
# and the user didn't specify AtlasCMake_DIR yet, then let's
# give it a default value.
if( IS_DIRECTORY ${CMAKE_SOURCE_DIR}/Build/AtlasCMake AND
      NOT AtlasCMake_DIR AND NOT ENV{AtlasCMake_DIR} )
   set( AtlasCMake_DIR ${CMAKE_SOURCE_DIR}/Build/AtlasCMake )
endif()

# If there's a directory called AtlasLCG in the project,
# and the user didn't specify LCG_DIR yet, then let's
# give it a default value.
if( IS_DIRECTORY ${CMAKE_SOURCE_DIR}/Build/AtlasLCG AND
      NOT LCG_DIR AND NOT ENV{LCG_DIR} )
   set( LCG_DIR ${CMAKE_SOURCE_DIR}/Build/AtlasLCG )
endif()

# Pick up a local version of the AtlasCMake code if it exists:
find_package( AtlasCMake QUIET )

# Find the project that we depend on:
find_package( $ENV{AtlasProject} )

# Set up CTest:
atlas_ctest_setup()

# These next lines are a temporary fix for issue: ATLINFR-2388
if( DEFINED ENV{WorkDir_PLATFORM} )
    set( ATLAS_PLATFORM $ENV{WorkDir_PLATFORM} )
endif()

# Set up a work directory project - hack for special case of 21.2 and 21.0 needing old signature:
if ( $ENV{AtlasBuildBranch} STREQUAL "21.2" OR $ENV{AtlasBuildBranch} STREQUAL "21.0" )
   atlas_project(WorkDir 1.0.0 USE $ENV{AtlasProject} $ENV{AtlasVersion})
else()
   atlas_project(USE $ENV{AtlasProject} $ENV{AtlasVersion})
endif()

# Set up the runtime environment setup script(s):
lcg_generate_env( SH_FILE ${CMAKE_BINARY_DIR}/${ATLAS_PLATFORM}/env_setup.sh )
install( FILES ${CMAKE_BINARY_DIR}/${ATLAS_PLATFORM}/env_setup.sh
   DESTINATION . )

# Set up CPack:
atlas_cpack_setup()
