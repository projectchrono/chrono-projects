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

#if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
#    set(CH_LINKERFLAG_EXE  "${CH_LINKERFLAG_EXE} -framework IOKit -framework Cocoa -framework OpenGL")
#endif()

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
option(ENABLE_IRRLICHT_PROGRAMS "Build the IRRLICHT programs" OFF)
option(ENABLE_PARALLEL_PROGRAMS "Build the PARALLEL programs" OFF)
option(ENABLE_VEHICLE_PARALLEL_PROGRAMS "Build the VEHICLE + PARALLEL programs" OFF)
option(ENABLE_FEA_PROGRAMS "Build the FEA programs" OFF)
option(ENABLE_FSI_PROGRAMS "Build the FSI programs" OFF)
option(ENABLE_PYPARSER_PROGRAMS "Build the PYPARSER programs" OFF)
option(ENABLE_PARALLEL_TESTS_PROGRAMS "Build the PARALLEL TEST programs" OFF)

# Keep track of all DLLs. Each submodule should append to this list.
list(APPEND ALL_DLLS "")

# Propagate configuration to submodules.
if(ENABLE_IRRLICHT_PROGRAMS)
  message(STATUS "==== Chrono IRRLICHT programs ====")
 add_subdirectory(irrlicht)
endif()
if(ENABLE_PARALLEL_PROGRAMS)
  message(STATUS "==== Chrono PARALLEL programs ====")
  add_subdirectory(parallel)
endif()
if(ENABLE_VEHICLE_PARALLEL_PROGRAMS)
  message(STATUS "==== Chrono VEHICLE_PARALLEL programs ====")
  add_subdirectory(vehicle_parallel)
endif()
if(ENABLE_FEA_PROGRAMS)
  message(STATUS "==== Chrono FEA programs ====")
  add_subdirectory(fea)
endif()
if(ENABLE_FSI_PROGRAMS)
  message(STATUS "==== Chrono FSI programs ====")
  add_subdirectory(fsi)
endif()
if(ENABLE_PYPARSER_PROGRAMS)
  message(STATUS "==== Chrono PYPARSER programs ====")
  add_subdirectory(pyparser)
endif()
#if(ENABLE_PARALLEL_TESTS_PROGRAMS)
#  add_subdirectory(parallel_tests)
#endif()

#--------------------------------------------------------------
# Copy DLLs
#--------------------------------------------------------------

# Add custom command for copying DLLs to the appropriate binary
# output folder (on Windows only).

add_DLL_copy_command("${ALL_DLLS}")


