#=============================================================================
# CMake configuration file for projects requiring the Chrono VEHICLE and
# various other modules
#=============================================================================

#--------------------------------------------------------------
# Find the Chrono package with required components
#--------------------------------------------------------------

# Invoke find_package in CONFIG mode

find_package(Chrono
             COMPONENTS Vehicle 
             OPTIONAL_COMPONENTS Irrlicht MKL Parallel OpenGL
             CONFIG
)

# If Chrono and/or the required component(s) were not found, return now.

if(NOT Chrono_FOUND)
  message("Could not find requirements for VEHICLE projects")
  return()
endif()

#--------------------------------------------------------------
# Include paths and libraries
#--------------------------------------------------------------

# (A) Path to the Chrono include headers
# - If using an installed version of Chrono, this will be the path 
#   to the installed headers (the configuration headers are also
#   available there)
# - If using a build version of Chrono, this will contain both the
#   path to the Chrono sources and the path to the chrono BUILD tree
#   (the latter for the configuration headers)
# 
# (B) Path to the top of the source tree for this project
# - for access to utility headers

include_directories(
    ${CHRONO_INCLUDE_DIRS}
    ${CMAKE_SOURCE_DIR}
)

#--------------------------------------------------------------
# Append to the parent's list of DLLs (and make it visible up)
#--------------------------------------------------------------

list(APPEND ALL_DLLS "${CHRONO_DLLS}")
set(ALL_DLLS "${ALL_DLLS}" PARENT_SCOPE)

#--------------------------------------------------------------
# Include subdirectories, based on their dependencies
#--------------------------------------------------------------

message(STATUS "Programs...")
add_subdirectory(test_steeringControl)
#add_subdirectory(test_pacjekaTire)

if (CHRONO_IRRLICHT_FOUND)
  message(STATUS "Programs using Chrono::Irrlicht...")
  add_subdirectory(test_sprocket)
  add_subdirectory(test_WheeledGenericManeuvers)
  add_subdirectory(test_HMMWV)
endif()

if (CHRONO_PARALLEL_FOUND)
  message(STATUS "Programs using Chrono::Parallel...")
  add_subdirectory(test_M113_parallel)
endif()

message(STATUS "")
