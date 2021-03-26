#=============================================================================
# CMake configuration file for projects requiring the Chrono Vehicle,
# Synchrono, and Sensor modules
# 
# Cannot be used stand-alone (but is mostly self-contained).
#=============================================================================

#--------------------------------------------------------------
# List of all executables
#--------------------------------------------------------------

set(DEMOS
    demo_SYN_scm
)

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
             COMPONENTS Irrlicht Vehicle Synchrono Sensor
             CONFIG
)

# If Chrono and/or the required component(s) were not found, return now.

if(NOT Chrono_FOUND)
  message("Could not find requirements for Synchrono projects")
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
    ${MPI_CXX_INCLUDE_PATH}
    ${CMAKE_SOURCE_DIR}
)

#--------------------------------------------------------------
# Append to the parent's list of DLLs (and make it visible up)
#--------------------------------------------------------------

list(APPEND ALL_DLLS "${CHRONO_DLLS}")
set(ALL_DLLS "${ALL_DLLS}" PARENT_SCOPE)

#--------------------------------------------------------------
# Compilation flags
#--------------------------------------------------------------

set(COMPILE_FLAGS ${CHRONO_CXX_FLAGS} ${MPI_CXX_COMPILE_FLAGS}) 

# Disable some warnings triggered by Irrlicht (Windows only)
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    SET(COMPILE_FLAGS "${COMPILE_FLAGS} /wd4275")
endif()

#--------------------------------------------------------------
# Loop over all demo programs and build them
#--------------------------------------------------------------

message(STATUS "Configurarion tests for Chrono::Irrlicht...")

foreach(PROGRAM ${DEMOS})

  message(STATUS "...add ${PROGRAM}")

  add_executable(${PROGRAM}  "${PROGRAM}.cpp")
  source_group(""  FILES "${PROGRAM}.cpp")

  set_target_properties(${PROGRAM} PROPERTIES
    FOLDER demos
    COMPILE_FLAGS "${COMPILE_FLAGS}"
    COMPILE_DEFINITIONS "CHRONO_DATA_DIR=\"${CHRONO_DATA_DIR}\";CHRONO_VEHICLE_DATA_DIR=\"${CHRONO_VEHICLE_DATA_DIR}\""
    LINK_FLAGS "${CHRONO_LINKER_FLAGS}"
  )

  target_link_libraries(${PROGRAM} ${CHRONO_LIBRARIES} ${MPI_CXX_LIBRARIES})

endforeach(PROGRAM)

message(STATUS "")
