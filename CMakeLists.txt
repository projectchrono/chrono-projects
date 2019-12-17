#=============================================================================
# CMake configuration file for Chrono-Projects
# 
#=============================================================================

cmake_minimum_required(VERSION 2.8)

project(ChronoProjects)

# Set location of additional CMAKE modules
#set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

# Set location of executable
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)

#-----------------------------------------------------------------------------
# Enable CTest
#-----------------------------------------------------------------------------

enable_testing()
include(CTest)

#-----------------------------------------------------------------------------
# Use dependent options
#-----------------------------------------------------------------------------

include(CMakeDependentOption)

#-----------------------------------------------------------------------------
# Always use full RPATH (differentiating between the build and install trees)
#-----------------------------------------------------------------------------

# use, i.e. don't skip the full RPATH for the build tree
set(CMAKE_SKIP_BUILD_RPATH  FALSE)

# when building, don't use the install RPATH already
# (but later on when installing)
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE) 

set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

# add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# the RPATH to be used when installing, but only if it's not a system directory
list(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
if("${isSystemDir}" STREQUAL "-1")
   set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
endif()

#-----------------------------------------------------------------------------
# Force C++11
#-----------------------------------------------------------------------------

#if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
#    set(CH_LINKERFLAG_EXE  "${CH_LINKERFLAG_EXE} -framework IOKit -framework Cocoa -framework OpenGL")
#endif()

#-----------------------------------------------------------------------------
# Fix for VS 2017 15.8 and newer to handle alignment specification with Eigen
#-----------------------------------------------------------------------------

if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
  if(MSVC AND ${MSVC_VERSION} GREATER_EQUAL 1915)
    add_definitions( "-D_ENABLE_EXTENDED_ALIGNED_STORAGE" )
  endif()
endif()

#-----------------------------------------------------------------------------
# Disable some warnings
#-----------------------------------------------------------------------------

if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    if(MSVC)
        add_definitions( "-D_CRT_SECURE_NO_DEPRECATE" )  # avoids deprecation warnings
        add_definitions( "-D_SCL_SECURE_NO_DEPRECATE" )  # avoids deprecation warnings
        add_definitions( "-DNOMINMAX" ) # do not use MSVC's min/max macros
        add_definitions( "-MP" ) # perform parallel builds
    endif(MSVC)
endif()

#-----------------------------------------------------------------------------
# Invoke CMake in subdirectories
#-----------------------------------------------------------------------------

# Options for configuring/building individual sets of programs
option(ENABLE_CONFIG_TESTS "Build the programs for checking Chrono configuration" OFF)
option(ENABLE_METRICS_TESTS "Build the programs for metrics testing" OFF)


# Keep track of all DLLs. Each submodule should append to this list.
list(APPEND ALL_DLLS "")

# Propagate configuration to submodules.
if(ENABLE_CONFIG_TESTS)
  message(STATUS "==== Chrono configuration programs ====\n")
  add_subdirectory(configuration_tests)
endif()
if(ENABLE_METRICS_TESTS)
  message(STATUS "==== Chrono metrics programs ====\n")
  add_subdirectory(metrics_tests)
endif()

add_subdirectory(projects)

# ------------------------------------------------------------------------------
# Copy the data folder to the bin directory
# ------------------------------------------------------------------------------

IF(MSVC)
  FILE( COPY ${CMAKE_SOURCE_DIR}/data/ DESTINATION ${CMAKE_BINARY_DIR}/bin/data/)
ELSEIF(XCODE_VERSION)
  FILE( COPY ${CMAKE_SOURCE_DIR}/data/ DESTINATION ${CMAKE_BINARY_DIR}/bin/data/)
ELSE()
  FILE( COPY ${CMAKE_SOURCE_DIR}/data/ DESTINATION ${CMAKE_BINARY_DIR}/data/)
ENDIF()

SET(CHRONO_DATA_DIR "${CHRONO_DATA_DIR}")
SET(PROJECTS_DATA_DIR "${CMAKE_SOURCE_DIR}/data/" CACHE INTERNAL "")

#-----------------------------------------------------------------------------
# Copy DLLs
#-----------------------------------------------------------------------------

# Add custom command for copying DLLs to the appropriate binary output folder
# (on Windows only).
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
	add_DLL_copy_command("${ALL_DLLS}")
endif()
