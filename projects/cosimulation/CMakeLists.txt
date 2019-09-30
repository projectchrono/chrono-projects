#=============================================================================
# CMake configuration file for Chrono::Vehicle co-simulation projects.
# Requirements:
#    MPI
#    Chrono::Vehicle and Chrono::Parallel modules
# Optional modules:
#    Chrono::MKL and Chrono::OpenGL
#=============================================================================

#-----------------------------------------------------------------------------
# MPI support
#-----------------------------------------------------------------------------

message(STATUS "Searching for MPI...")
find_package(MPI)
message(STATUS "  MPI (C++) found:       ${MPI_CXX_FOUND}")
if(MPI_CXX_FOUND)
  message(STATUS "  MPI compiler:          ${MPI_CXX_COMPILER}")
  message(STATUS "  MPI compile flags:     ${MPI_CXX_COMPILE_FLAGS}")
  message(STATUS "  MPI include path:      ${MPI_CXX_INCLUDE_PATH}")
  message(STATUS "  MPI link flags:        ${MPI_CXX_LINK_FLAGS}")
  message(STATUS "  MPI libraries:         ${MPI_CXX_LIBRARIES}")
  message(STATUS "")
  message(STATUS "  MPIEXEC:               ${MPIEXEC}")
  message(STATUS "  MPIEXEC_NUMPROC_FLAG:  ${MPIEXEC_NUMPROC_FLAG}")
  message(STATUS "  MPIEXEC_PREFLAGS:      ${MPIEXEC_PREFLAGS}")
  message(STATUS "  MPIEXEC_POSTFLAGS:     ${MPIEXEC_POSTFLAGS}")
else()
  message(STATUS "  Co-simulation programs disabled...")
endif()

#--------------------------------------------------------------
# Find the Chrono package with required components
#--------------------------------------------------------------

# Invoke find_package in CONFIG mode

find_package(Chrono
             COMPONENTS Vehicle Parallel
             OPTIONAL_COMPONENTS MKL OpenGL Irrlicht
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
# Include subdirectories
#--------------------------------------------------------------

message(STATUS "Programs...")

add_subdirectory(test_FEA_cosimulation)

message(STATUS "")
