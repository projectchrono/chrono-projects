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

#--------------------------------------------------------------
# Append to the parent's lists of DLLs (and make them visible)
#--------------------------------------------------------------

list(APPEND ALL_DLL_NAMES "${CHRONO_DLL_NAMES}")
set(ALL_DLL_NAMES ${ALL_DLL_NAMES} PARENT_SCOPE)

# Make CHRONO_PACKAGE_PREFIX_DIR visible above (needed for add_DLL_copy_command)
set(CHRONO_PACKAGE_PREFIX_DIR ${CHRONO_PACKAGE_PREFIX_DIR} PARENT_SCOPE)

#--------------------------------------------------------------
# Include subdirectories
#--------------------------------------------------------------

message(STATUS "Programs...")

add_subdirectory(test_FEA_cosimulation)

