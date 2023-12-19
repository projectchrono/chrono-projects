#=============================================================================
# CMake configuration file for Chrono::Vehicle co-simulation projects.
# Requirements:
#    MPI
#    Chrono::Vehicle and Chrono::Multicore modules
# Optional modules:
#    Chrono::PardisoMKL and Chrono::OpenGL
#=============================================================================

message("\n---- VEHICLE COSIMULATION projects")

#--------------------------------------------------------------
# Find the Chrono package with required components
#--------------------------------------------------------------

# Invoke find_package in CONFIG mode
find_package(Chrono QUIET
             COMPONENTS Vehicle Multicore
             OPTIONAL_COMPONENTS PardisoMKL OpenGL Irrlicht
             CONFIG
)

# If Chrono and/or the required component(s) were not found, return now.
if(NOT Chrono_FOUND)
  message("Could not find requirements for VEHICLE co-simulation projects")
  return()
endif()

if(NOT CHRONO_MPI_FOUND)
  message("The co-simulation projects require MPI which was not found!")
  return()
endif()

# Use same C++ standard as for building Chrono
set(CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_STANDARD ${CHRONO_CXX_STANDARD})

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
# Append to the parent's lists of DLLs (and make them visible)
#--------------------------------------------------------------

list(APPEND ALL_DLL_NAMES       "${CHRONO_DLL_NAMES}")
list(APPEND ALL_DEPENDENCY_DLLS "${CHRONO_DEPENDENCY_DLLS}")
set(ALL_DLL_NAMES       "${ALL_DLL_NAMES}"       PARENT_SCOPE)
set(ALL_DEPENDENCY_DLLS "${ALL_DEPENDENCY_DLLS}" PARENT_SCOPE)

#--------------------------------------------------------------
# Include subdirectories
#--------------------------------------------------------------

message(STATUS "Programs...")

add_subdirectory(test_FEA_cosimulation)

