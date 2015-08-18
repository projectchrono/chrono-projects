#=============================================================================
# CMake configuration file for projects requiring the Chrono unit_PARALLEL
# 
# Cannot be used stand-alone (but is mostly self-contained).
# Depends on Chrono with unit_PARALLEL enabled and optionally OpenGL support.
#=============================================================================

message(STATUS "==== Chrono PARALLEL programs ====")

#--------------------------------------------------------------
# List of all executables
#--------------------------------------------------------------

set(DEMOS
    demo_crater
    demo_foam
    demo_massflow
    demo_objectdrop
    demo_shear
    demo_soilbin
    demo_suspension
    demo_wheel
    directShear
    mixerDEM
    pressureSinkage
    rad_im_schmutz
    singleWheel
)

#--------------------------------------------------------------
# Find the Chrono package with required components
#--------------------------------------------------------------

# Invoke find_package in CONFIG mode.

find_package(Chrono
             COMPONENTS Parallel
             #OPTIONAL_COMPONENTS OpenGL
             CONFIG
)

# If Chrono and/or the required component(s) were not found, return now.

if(NOT Chrono_FOUND)
  message("Could not find requirements for PARALLEL projects")
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

message(STATUS "Parallel programs...")

foreach(PROGRAM ${DEMOS})

  message(STATUS "...add ${PROGRAM}")

  add_executable(${PROGRAM}  "${PROGRAM}.cpp")
  source_group(""  FILES "${PROGRAM}.cpp")

  set_target_properties(${PROGRAM} PROPERTIES
    FOLDER demos
    COMPILE_FLAGS "${CHRONO_CFLAGS}"
    LINK_FLAGS "${CHRONO_LDFLAGS}"
  )

  target_link_libraries(${PROGRAM} ${CHRONO_LIBRARIES})
  
  #install(TARGETS ${PROGRAM} DESTINATION bin)

endforeach(PROGRAM)

