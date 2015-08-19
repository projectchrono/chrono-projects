#=============================================================================
# CMake configuration file for projects requiring the Chrono unit_FEA
# and optionally the unit_IRRLICHT
# 
# Cannot be used stand-alone (but is mostly self-contained).
# Depends on Chrono with unit_IRRLICHT enabled (and optionally unit_IRRLICHT).
#=============================================================================

#--------------------------------------------------------------
# List of all executables
#--------------------------------------------------------------

set(DEMOS_BASIC
    demo_FEAbasic
)

set(DEMOS_IRR
    demo_FEAthermal
)

#--------------------------------------------------------------
# Find the Chrono package with required and optional components
#--------------------------------------------------------------

# Invoke find_package in CONFIG mode

find_package(Chrono
             COMPONENTS FEA
             OPTIONAL_COMPONENTS Irrlicht
             CONFIG
)

# If Chrono and/or the required component(s) were not found, return now.

if(NOT Chrono_FOUND)
  message("Could not find requirements for FEA projects")
  return()
endif()

#--------------------------------------------------------------
# Set the list of executables to be created
#--------------------------------------------------------------

set(DEMOS ${DEMOS_BASIC})

# If the optional component IRRLICHT was found, append to list
if(${CHRONO_IRRLICHT_FOUND})
  set(DEMOS ${DEMOS} ${DEMOS_IRR})
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

set(COMPILE_FLAGS ${CHRONO_CFLAGS})

# Disable some warnings triggered by Irrlicht (Windows only)
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    SET(COMPILE_FLAGS "${COMPILE_FLAGS} /wd4275")
endif()

#--------------------------------------------------------------
# Loop over all demo programs and build them
#--------------------------------------------------------------

message(STATUS "FEA programs...")

foreach(PROGRAM ${DEMOS})

  message(STATUS "...add ${PROGRAM}")

  add_executable(${PROGRAM}  "${PROGRAM}.cpp")
  source_group(""  FILES "${PROGRAM}.cpp")

  set_target_properties(${PROGRAM} PROPERTIES
    FOLDER demos
    COMPILE_FLAGS "${COMPILE_FLAGS}"
    COMPILE_DEFINITIONS "CHRONO_DATA_DIR=\"${CHRONO_DATA_DIR}\""
    LINK_FLAGS "${CHRONO_LDFLAGS}"
  )

  target_link_libraries(${PROGRAM} ${CHRONO_LIBRARIES})
  
  #install(TARGETS ${PROGRAM} DESTINATION bin)

endforeach(PROGRAM)

