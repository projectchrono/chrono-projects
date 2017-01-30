#=============================================================================
# CMake configuration file for projects requiring the Chrono PARALLEL module
# 
# Cannot be used stand-alone (but is mostly self-contained).
#=============================================================================

#--------------------------------------------------------------
# List of all executables
#--------------------------------------------------------------

set(DEMOS
    test_PAR_foam
    test_PAR_objectdrop
    test_PAR_soilbin
    test_PAR_suspension
    test_PAR_wheel
    test_PAR_radImSchmutz
    test_PAR_settling
)

set(DEMOS_OPENGL
    test_PAR_active_box
)
#--------------------------------------------------------------
# Find the Chrono package with required components
#--------------------------------------------------------------

# Invoke find_package in CONFIG mode.

find_package(Chrono
             COMPONENTS Parallel
             OPTIONAL_COMPONENTS OpenGL
             CONFIG
)

# If Chrono and/or the required component(s) were not found, return now.

if(NOT Chrono_FOUND)
  message("Could not find requirements for PARALLEL projects")
  return()
endif()

if (CHRONO_OPENGL_FOUND)
    set(DEMOS ${DEMOS} ${DEMOS_OPENGL})
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
# Loop over all demo programs and build them
#--------------------------------------------------------------

message(STATUS "Parallel tests...")

foreach(PROGRAM ${DEMOS})

  message(STATUS "...add ${PROGRAM}")

  add_executable(${PROGRAM}  "${PROGRAM}.cpp")
  source_group(""  FILES "${PROGRAM}.cpp")

  set_target_properties(${PROGRAM} PROPERTIES
    FOLDER demos
    COMPILE_FLAGS "${CHRONO_CXX_FLAGS}"
    COMPILE_DEFINITIONS "CHRONO_DATA_DIR=\"${CHRONO_DATA_DIR}\""
    LINK_FLAGS "${CHRONO_CXX_FLAGS} ${CHRONO_LINKER_FLAGS}"
  )

  target_link_libraries(${PROGRAM} ${CHRONO_LIBRARIES})

endforeach(PROGRAM)

message(STATUS "")
