#=============================================================================
# CMake configuration file for projects requiring the Chrono unit_VEHICLE
#
# Cannot be used stand-alone (but is mostly self-contained).
# Depends on Chrono with unit_VEHICLE and unit_PARALLEL enabled and optional
# OpenGL support.
#=============================================================================

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

# Specify the location of the Chrono project configuration file

# TODO

# Invoke find_package in CONFIG mode

find_package(Chrono
             COMPONENTS Vehicle Parallel
             CONFIG
)

# If Chrono and/or the required component(s) were not found, return now.

if(NOT Chrono_FOUND)
  message("Could not find requirements for VEHICLE projects")
  return()
endif()

# Report what we have found

message(STATUS "Chrono include directories" ${CHRONO_INCLUDE_DIRS})
message(STATUS "Chrono libraries" ${CHRONO_LIBRARIES})
message(STATUS "Chrono DLLs" ${CHRONO_DLLS})

#--------------------------------------------------------------
# Include paths and libraries
#--------------------------------------------------------------

# Path to the Chrono include headers
# - If using an installed version of Chrono, this will be the path 
#   to the installed headers (the configuration headers are also
#   available there)
# - If using a build version of Chrono, this will contain both the
#   path to the Chrono sources and the path to the chrono BUILD tree
#   (the latter for the configuration headers)
# Path to the top of the source tree for this project
# - for access to utility headers

include_directories(
    ${CHRONO_INCLUDE_DIRS}
    ${CMAKE_SOURCE_DIR}
)

# Initialize the list of libraries:
# - main Chrono::Engine library
# - Chrono::Vehicle library
# - Chrono::Parallel library
# - Utility library for Chrono::Parallel

set(LIBRARIES ${CHRONO_LIBRARIES})

# On windows, initialize the list of DLLs to be copied

if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
  set(DLLS ${CHRONO_DLLS})
endif()

# If Chrono was configured with OpenGL support

if(ENABLE_OPENGL)
  # Additional include directories for OpenGL support (GLM, GLEW, GLFW)
  include_directories(${CHRONO_OPENGL_INCLUDE_DIRS})

  # Append additional libraries for OpenGL support
  set(LIBRARIES ${LIBRARIES} ${CHRONO_OPENGL_LIBRARIES})

  # Collect additional DLLs for OpenGL support
  if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    set(DLLS ${DLLS} ${CHRONO_OPENGL_DLLS})
  endif()
endif()

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

MESSAGE(STATUS "Radu demo programs...")

FOREACH(PROGRAM ${DEMOS})
    MESSAGE(STATUS "...add ${PROGRAM}")

    ADD_EXECUTABLE(${PROGRAM}  "${PROGRAM}.cpp")
    SOURCE_GROUP(""  FILES "${PROGRAM}.cpp")

    SET_TARGET_PROPERTIES(${PROGRAM} PROPERTIES
        FOLDER demos
        COMPILE_FLAGS "${CH_BUILDFLAGS}"
        LINK_FLAGS "${CH_LINKERFLAG_EXE}"
    )

    TARGET_LINK_LIBRARIES(${PROGRAM} ${LIBRARIES})
    ADD_DEPENDENCIES(${PROGRAM} ${LIBRARIES})

    INSTALL(TARGETS ${PROGRAM} DESTINATION bin)

ENDFOREACH(PROGRAM)

