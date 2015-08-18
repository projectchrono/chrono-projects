#=============================================================================
# CMake configuration file for projects requiring the Chrono unit_IRRLICHT
# 
# Cannot be used stand-alone (but is mostly self-contained).
# Depends on Chrono with unit_IRRLICHT enabled.
#=============================================================================


#--------------------------------------------------------------
# List of all executables
#--------------------------------------------------------------

set(DEMOS
    demo_collision
    #slider_crank
)

#--------------------------------------------------------------
# Find the Chrono package with required components
#--------------------------------------------------------------

# Invoke find_package in CONFIG mode

find_package(Chrono
             COMPONENTS Irrlicht
             CONFIG
)

# If Chrono and/or the required component(s) were not found, return now.

if(NOT Chrono_FOUND)
  message("Could not find requirements for IRRLICHT projects")
  return()
endif()

# Report what we have found

message(STATUS "Chrono include directories: ${CHRONO_INCLUDE_DIRS}")
message(STATUS "Chrono library names: ${CHRONO_LIB_NAMES}")
message(STATUS "Chrono libraries: ${CHRONO_LIBRARIES}")
message(STATUS "Chrono DLL names: ${CHRONO_DLL_NAMES}")
message(STATUS "Chrono DLLs: ${CHRONO_DLLS}")
message(STATUS "Chrono data directory: ${CHRONO_DATA_DIR}")

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
# - Chrono::Irrlicht library

set(LIBRARIES ${CHRONO_LIBRARIES})

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

message(STATUS "Irrlicht programs...")

foreach(PROGRAM ${DEMOS})

  message(STATUS "...add ${PROGRAM}")

  add_executable(${PROGRAM}  "${PROGRAM}.cpp")
  source_group(""  FILES "${PROGRAM}.cpp")

  set_target_properties(${PROGRAM} PROPERTIES
    FOLDER demos
    COMPILE_FLAGS "${CHRONO_CFLAGS}"
    LINK_FLAGS "${CH_LINKERFLAG_EXE}"
  )

  target_link_libraries(${PROGRAM} ${LIBRARIES})
  
  #install(TARGETS ${PROGRAM} DESTINATION bin)

endforeach(PROGRAM)

