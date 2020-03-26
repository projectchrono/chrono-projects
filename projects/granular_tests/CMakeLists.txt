#=============================================================================
# CMake configuration file for projects requiring the Chrono Granular module.
# 
# Cannot be used stand-alone (but is mostly self-contained).
#=============================================================================

#--------------------------------------------------------------
# List of all executables
#--------------------------------------------------------------

set(DEMOS
        test_GRAN_DamBreak
        test_GRAN_Bucket
        test_GRAN_rover
        test_GRAN_directshear
        test_GRAN_mixer
        test_GRAN_Coneflow
        test_GRAN_commvessels
        test_GRAN_bulkcompress
        test_GRAN_milsettle
        test_GRAN_testsuite
        test_GRAN_normal
        test_GRAN_slide
        test_GRAN_roll
        test_GRAN_incline
        test_GRAN_repose
)

#--------------------------------------------------------------
# Find the Chrono package with required and optional components
#--------------------------------------------------------------

# Invoke find_package in CONFIG mode

find_package(Chrono
             COMPONENTS Granular
             OPTIONAL_COMPONENTS Vehicle
             CONFIG
)

# If Chrono and/or the required component(s) were not found, return now.
if(NOT Chrono_FOUND)
  message("Could not find requirements for GRANULAR projects")
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
# Append to the parent's list of DLLs
#--------------------------------------------------------------

list(APPEND ALL_DLLS "${CHRONO_DLLS}")
set(ALL_DLLS "${ALL_DLLS}" PARENT_SCOPE)

#--------------------------------------------------------------
# Compilation flags
#--------------------------------------------------------------

set(COMPILE_FLAGS ${CHRONO_CXX_FLAGS})

# Disable some warnings triggered by Irrlicht (Windows only)
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    SET(COMPILE_FLAGS "${COMPILE_FLAGS} /wd4275")
endif()

#--------------------------------------------------------------
# Loop over all demo programs and build them
#--------------------------------------------------------------

message(STATUS "Projects for Chrono::Granular...")

foreach(PROGRAM ${DEMOS})

  message(STATUS "...add ${PROGRAM}")

  add_executable(${PROGRAM}  "${PROGRAM}.cpp")
  source_group(""  FILES "${PROGRAM}.cpp")
  set_target_properties(${PROGRAM} PROPERTIES
    FOLDER demos
    COMPILE_FLAGS "${COMPILE_FLAGS}"
    COMPILE_DEFINITIONS "PROJECTS_DATA_DIR=\"${PROJECTS_DATA_DIR}\";CHRONO_DATA_DIR=\"${CHRONO_DATA_DIR}\""
    LINK_FLAGS "${CHRONO_CXX_FLAGS} ${CHRONO_LINKER_FLAGS}"
  )

  target_link_libraries(${PROGRAM} ${CHRONO_LIBRARIES})

endforeach(PROGRAM)

message(STATUS "")

if (CHRONO_VEHICLE_FOUND)
    message(STATUS "Programs using Chrono::Vehicles...")
    add_subdirectory(HMMWV)
endif()
