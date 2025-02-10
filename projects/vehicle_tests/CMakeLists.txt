#=============================================================================
# CMake configuration file for projects requiring the Chrono VEHICLE and
# various other modules
#=============================================================================

message(STATUS "\n---- VEHICLE projects\n")

#--------------------------------------------------------------
# Find the Chrono package with required components
#--------------------------------------------------------------

# Invoke find_package in CONFIG mode
find_package(Chrono QUIET
             COMPONENTS Vehicle 
             OPTIONAL_COMPONENTS Irrlicht PardisoMKL Multicore OpenGL VSG
             CONFIG
)

# If Chrono and/or the required component(s) were not found, return now.
if(NOT Chrono_FOUND)
  message("Could not find requirements for VEHICLE projects")
  return()
endif()

#--------------------------------------------------------------
# Append to the parent's lists of DLLs (and make them visible)
#--------------------------------------------------------------

list(APPEND ALL_DLL_NAMES       "${CHRONO_DLL_NAMES}")
list(APPEND ALL_DEPENDENCY_DLLS "${CHRONO_DEPENDENCY_DLLS}")
set(ALL_DLL_NAMES       "${ALL_DLL_NAMES}"       PARENT_SCOPE)
set(ALL_DEPENDENCY_DLLS "${ALL_DEPENDENCY_DLLS}" PARENT_SCOPE)

# Make CHRONO_PACKAGE_PREFIX_DIR visible above (needed for add_DLL_copy_command)
set(CHRONO_PACKAGE_PREFIX_DIR ${CHRONO_PACKAGE_PREFIX_DIR} PARENT_SCOPE)

#--------------------------------------------------------------
# Include subdirectories, based on their dependencies
#--------------------------------------------------------------

message(STATUS "Programs...")
add_subdirectory(test_steeringControl)
add_subdirectory(test_FEDA)
add_subdirectory(test_BMW_E90)
add_subdirectory(test_Cherokee)

if (CHRONO_IRRLICHT_FOUND)
  message(STATUS "Programs using Chrono::Irrlicht...")
  add_subdirectory(test_sprocket)
  add_subdirectory(test_WheeledGenericManeuvers)
  add_subdirectory(test_HMMWV)
endif()

if (CHRONO_MULTICORE_FOUND)
  message(STATUS "Programs using Chrono::Multicore...")
  add_subdirectory(test_M113_multicore)
endif()
