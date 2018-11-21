#=============================================================================
# CMake configuration file for metrics tests for main Chrono::Engine module
# 
# Cannot be used stand-alone (but is mostly self-contained).
#=============================================================================

#--------------------------------------------------------------
# List of all executables
#--------------------------------------------------------------

set(TESTS
  )

set(TESTS_AVX
    metrics_CH_matMultAVX
)


#--------------------------------------------------------------
# Find the Chrono package with required and optional components
#--------------------------------------------------------------

# Invoke find_package in CONFIG mode

find_package(Chrono
             CONFIG
)

# If Chrono and/or the required component(s) were not found, return now.

if(NOT Chrono_FOUND)
  message("Could not find requirements for CORE metrics")
  return()
endif()

#--------------------------------------------------------------
# Set the list of executables to be created
#--------------------------------------------------------------

# If Chrono has AVX support, append tests
if(CHRONO_AVX_FOUND)
  list(APPEND TESTS ${TESTS_AVX})
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

#--------------------------------------------------------------

if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
  set(WORK_DIR ${PROJECT_BINARY_DIR}/bin/$<CONFIGURATION>)
else()
  set(WORK_DIR ${PROJECT_BINARY_DIR}/bin)
endif()

#--------------------------------------------------------------
# Loop over all demo programs and build them
#--------------------------------------------------------------

message(STATUS "Metrics tests for Chrono::Engine...")

foreach(PROGRAM ${TESTS})

  message(STATUS "...add ${PROGRAM}")

  add_executable(${PROGRAM}  "${PROGRAM}.cpp")
  source_group(""  FILES "${PROGRAM}.cpp")

  set_target_properties(${PROGRAM} PROPERTIES
    FOLDER demos
    COMPILE_FLAGS "${COMPILE_FLAGS}"
    COMPILE_DEFINITIONS "CHRONO_DATA_DIR=\"${CHRONO_DATA_DIR}\""
    LINK_FLAGS "${CHRONO_CXX_FLAGS} ${CHRONO_LINKER_FLAGS}"
  )

  target_link_libraries(${PROGRAM} ${CHRONO_LIBRARIES})

  # Note: this is not intended to work on Windows!
  add_test(NAME ${PROGRAM}
           WORKING_DIRECTORY ${WORK_DIR}
           COMMAND ${WORK_DIR}/${PROGRAM}
           )

endforeach(PROGRAM)

message(STATUS "")
