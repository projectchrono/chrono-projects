#=============================================================================
# CMake configuration file for projects requiring the Chrono MULTICORE module
# 
# Cannot be used stand-alone (but is mostly self-contained).
#=============================================================================

message("\n---- CONSTRAINT FLUID projects")

#--------------------------------------------------------------
# List of all executables
#--------------------------------------------------------------

set(DEMOS
    test_CNFL_settling_test
    test_CNFL_tracked_vehicle_fording
	  test_CNFL_humvee_water_tank
)

#--------------------------------------------------------------
# Find the Chrono package with required components
#--------------------------------------------------------------

# Invoke find_package in CONFIG mode.
set(ZLIB_ROOT "" CACHE PATH "Where is ZLIB located?")
find_package(ZLIB REQUIRED)

if (ZLIB_FOUND)
  message(STATUS "ZLIB include dir: ${ZLIB_INCLUDE_DIRS}")
  message(STATUS "ZLIB library:    ${ZLIB_LIBRARIES}")
else()
  mark_as_advanced(CLEAR ZLIB_ROOT)
endif()

find_package(Chrono QUIET
             COMPONENTS Multicore Vehicle
             OPTIONAL_COMPONENTS OpenGL  
             CONFIG
)

# If Chrono and/or the required component(s) were not found, return now.
if(NOT Chrono_FOUND)
  message("Could not find requirements for CONSTRAINT FLUID projects")
  return()
endif()

# Use same C++ standard as for building Chrono
set(CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_STANDARD ${CHRONO_CXX_STANDARD})

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
    ${ZLIB_INCLUDE_DIRS}
)

#--------------------------------------------------------------
# Append to the parent's lists of DLLs (and make them visible)
#--------------------------------------------------------------

list(APPEND ALL_DLL_NAMES       "${CHRONO_DLL_NAMES}")
list(APPEND ALL_DEPENDENCY_DLLS "${CHRONO_DEPENDENCY_DLLS}")
set(ALL_DLL_NAMES       "${ALL_DLL_NAMES}"       PARENT_SCOPE)
set(ALL_DEPENDENCY_DLLS "${ALL_DEPENDENCY_DLLS}" PARENT_SCOPE)

#--------------------------------------------------------------
# Loop over all demo programs and build them
#--------------------------------------------------------------

message(STATUS "Constraint_fluids tests...")

foreach(PROGRAM ${DEMOS})

  message(STATUS "...add ${PROGRAM}")

  add_executable(${PROGRAM}  "${PROGRAM}.cpp")
  source_group(""  FILES "${PROGRAM}.cpp")

  target_compile_definitions(${PROGRAM} PUBLIC "CHRONO_DATA_DIR=\"${CHRONO_DATA_DIR}\"") 
  target_compile_options(${PROGRAM} PUBLIC ${CHRONO_CXX_FLAGS})
  target_link_options(${PROGRAM} PUBLIC ${CHRONO_LINKER_FLAGS})
  target_link_libraries(${PROGRAM} ${CHRONO_LIBRARIES} ${ZLIB_LIBRARIES})

endforeach()

