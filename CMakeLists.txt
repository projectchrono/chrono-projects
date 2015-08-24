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
#if(ENABLE_PARALLEL_TESTS_PROGRAMS)
#  add_subdirectory(parallel_tests)
#endif()

#--------------------------------------------------------------
# Copy DLLs
#--------------------------------------------------------------

# On Windows, traverse the list of DLLs and copy them next to
# the binaries.

if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

  # Remove duplicates in the list of DLLs
  list(REMOVE_DUPLICATES ALL_DLLS)

  message(STATUS "==== Add custom commands for copying DLLs ====")
  message(STATUS "DLL list...")

  # Create custom target for copying DLLs; add it to the default build target
  add_custom_target(COPY_DLLS ALL)

  # Loop over the list of all DLLs and create a custom command
  foreach(DLL ${ALL_DLLS})

    message(STATUS "...${DLL}")

    # Create custom commands, invoked post-build to copy DLLs to the appropriate
    # directory (depending on the configuration selected at build time in VS)
    add_custom_command(
        TARGET COPY_DLLS PRE_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
            "${DLL}"
            "${CMAKE_BINARY_DIR}/bin/$<CONFIGURATION>/"
    )

  endforeach()

endif()
