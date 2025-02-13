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
  message(STATUS "ZLIB library: ${ZLIB_LIBRARIES}")
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

#--------------------------------------------------------------
# Append to the parent's lists of DLLs (and make them visible)
#--------------------------------------------------------------

list(APPEND ALL_DLL_NAMES "${CHRONO_DLL_NAMES}")
set(ALL_DLL_NAMES ${ALL_DLL_NAMES} PARENT_SCOPE)

# Make CHRONO_PACKAGE_PREFIX_DIR visible above (needed for add_DLL_copy_command)
set(CHRONO_PACKAGE_PREFIX_DIR ${CHRONO_PACKAGE_PREFIX_DIR} PARENT_SCOPE)

#--------------------------------------------------------------
# Loop over all demo programs and build them
#--------------------------------------------------------------

message(STATUS "Constraint_fluids tests...")

foreach(PROGRAM ${DEMOS})
  message(STATUS "...add ${PROGRAM}")

  add_executable(${PROGRAM}  "${PROGRAM}.cpp")
  source_group(""  FILES "${PROGRAM}.cpp")

  set_property(TARGET ${PROGRAM} PROPERTY VS_DEBUGGER_WORKING_DIRECTORY "$<TARGET_FILE_DIR:${PROGRAM}>")

  target_compile_definitions(${PROGRAM} PUBLIC "CHRONO_DATA_DIR=\"${CHRONO_DATA_DIR}\"") 
  target_link_libraries(${PROGRAM} PRIVATE ${CHRONO_TARGETS})
  target_link_libraries(${PROGRAM} ${ZLIB_LIBRARIES})
endforeach()

