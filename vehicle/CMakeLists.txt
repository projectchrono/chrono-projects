#=============================================================================
# CMake configuration file for projects requiring the Chrono unit_VEHICLE
# and unit_PARALLEL
# 
# Cannot be used stand-alone (but is mostly self-contained).
# Depends on Chrono with unit_VEHICLE and unit_PARALLEL enabled and optional
# OpenGL support.
#=============================================================================

message(STATUS "==== Chrono VEHICLE_PARALLEL programs ====")

#--------------------------------------------------------------
# List of all executables
#--------------------------------------------------------------

SET(DEMOS
    hmmwv
    hmmwvDEM
    hmmwvDEM_ditch
)

#--------------------------------------------------------------
# Find the Chrono package with required components
#--------------------------------------------------------------

# Invoke find_package in CONFIG mode

find_package(Chrono
             COMPONENTS Vehicle Parallel
             #OPTIONAL_COMPONENTS OpenGL
             CONFIG
)

# If Chrono and/or the required component(s) were not found, return now.

if(NOT Chrono_FOUND)
  message("Could not find requirements for VEHICLE_PARALLEL projects")
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
# Copy DLLs
#--------------------------------------------------------------

# On Windows, traverse the list of DLLs and copy them next to
# the binaries.

if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

#    TODO

endif()

#--------------------------------------------------------------
# Loop over all demo programs and build them
#--------------------------------------------------------------

MESSAGE(STATUS "Vehicle_Parallel programs...")

FOREACH(PROGRAM ${DEMOS})
    MESSAGE(STATUS "...add ${PROGRAM}")

    ADD_EXECUTABLE(${PROGRAM}  "${PROGRAM}.cpp")
    SOURCE_GROUP(""  FILES "${PROGRAM}.cpp")

    SET_TARGET_PROPERTIES(${PROGRAM} PROPERTIES
        FOLDER demos
        COMPILE_FLAGS "${CHRONO_CFLAGS}"
        COMPILE_DEFINITIONS "CHRONO_DATA_DIR=\"${CHRONO_DATA_DIR}\""
        LINK_FLAGS "${CHRONO_LDFLAGS}"
    )

    TARGET_LINK_LIBRARIES(${PROGRAM} ${CHRONO_LIBRARIES})

    #INSTALL(TARGETS ${PROGRAM} DESTINATION bin)

ENDFOREACH(PROGRAM)

