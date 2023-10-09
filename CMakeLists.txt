#=============================================================================
# CMake configuration file for Chrono-Projects
# 
#=============================================================================

cmake_minimum_required(VERSION 3.10)

project(ChronoProjects)

# Set location of additional CMAKE modules
#set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

# Set location of executable
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)

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
        add_compile_options(/wd4251) # class needs to have dll-interface
    endif(MSVC)
endif()

#-----------------------------------------------------------------------------
# Invoke CMake in subdirectories
#-----------------------------------------------------------------------------

if(NOT Chrono_DIR)
  set(Chrono_DIR "" CACHE PATH "The directory containing a CMake configuration file for Chrono.")
  return()
endif()

# Options for configuring/building individual sets of programs
option(ENABLE_CONFIG_TESTS "Build the programs for checking Chrono configuration" ON)

# Keep track of all DLLs. Each submodule should append to this list.
list(APPEND ALL_DLL_NAMES "")
list(APPEND ALL_DEPENDENCY_DLLS "")

# Propagate configuration to submodules.
if(ENABLE_CONFIG_TESTS)
  message(STATUS "\n==== Chrono configuration programs ====")
  add_subdirectory(configuration_tests)
endif()

add_subdirectory(projects)

# ------------------------------------------------------------------------------
# Set path to Chrono data directory 
# ------------------------------------------------------------------------------

SET(CHRONO_DATA_DIR "${CHRONO_DATA_DIR}")

#------------------------------------------------------------
# Copy chrono-projects data directory to BUILD tree 
#------------------------------------------------------------

if(MSVC OR XCODE_VERSION)
    file(COPY ${CMAKE_SOURCE_DIR}/data/ DESTINATION ${CMAKE_BINARY_DIR}/bin/data/)
else()
    file(COPY ${CMAKE_SOURCE_DIR}/data/ DESTINATION ${CMAKE_BINARY_DIR}/data/)
endif()


#-----------------------------------------------------------------------------
# Copy DLLs
#-----------------------------------------------------------------------------

# Add custom command for copying DLLs to the appropriate binary output folder
# (on Windows only).
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
  message(STATUS "")
	add_DLL_copy_command("${ALL_DLL_NAMES}" "${ALL_DEPENDENCY_DLLS}")
endif()
