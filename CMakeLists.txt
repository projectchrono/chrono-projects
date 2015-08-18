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

if(${CMAKE_CXX_COMPILER_ID} MATCHES "GNU")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
endif()

if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
	set(CH_LINKERFLAG_EXE  "${CH_LINKERFLAG_EXE} -framework IOKit -framework Cocoa -framework OpenGL")
endif()

#-----------------------------------------------------------------------------
# Invoke CMake in subdirectories
#-----------------------------------------------------------------------------

add_subdirectory(irrlicht)
add_subdirectory(parallel)
#add_subdirectory(parallel_tests)
#add_subdirectory(vehicle)

