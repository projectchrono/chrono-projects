#=============================================================================
# CMake configuration file for projects requiring the Chrono Vehicle,
# Synchrono, and Sensor modules
# 
# Cannot be used stand-alone (but is mostly self-contained).
#=============================================================================

message(STATUS "\n---- SYNCHRONO configuration test\n")

#--------------------------------------------------------------
# List of all executables
#--------------------------------------------------------------

set(DEMOS
    demo_SYN_scm
)

#--------------------------------------------------------------
# Find the Chrono package with required components
#--------------------------------------------------------------

# Invoke find_package in CONFIG mode
find_package(Chrono QUIET
             COMPONENTS Vehicle Synchrono
             OPTIONAL_COMPONENTS Irrlicht Sensor
             CONFIG
)

# If Chrono and/or the required component(s) were not found, return now.
if(NOT Chrono_FOUND)
  message("Could not find requirements for Synchrono projects")
  return()
endif()

#--------------------------------------------------------------
# Append to the parent's lists of DLLs (and make them visible)
#--------------------------------------------------------------

list(APPEND ALL_DLL_NAMES       "${CHRONO_DLL_NAMES}")
list(APPEND ALL_DEPENDENCY_DLLS "${CHRONO_DEPENDENCY_DLLS}")
set(ALL_DLL_NAMES       "${ALL_DLL_NAMES}"       PARENT_SCOPE)
set(ALL_DEPENDENCY_DLLS "${ALL_DEPENDENCY_DLLS}" PARENT_SCOPE)

# Make CHRONO_PACKAGE_PREFIX_DIR visible above (needed for add_DLL_copy_command)
set(CHRONO_PACKAGE_PREFIX_DIR ${CHRONO_PACKAGE_PREFIX_DIR} PARENT_SCOPE)

#--------------------------------------------------------------
# Loop over all demo programs and build them
#--------------------------------------------------------------

message(STATUS "Configuration tests for Chrono::Synchrono...")

foreach(PROGRAM ${DEMOS})
  message(STATUS "...add ${PROGRAM}")

  add_executable(${PROGRAM}  "${PROGRAM}.cpp")
  source_group(""  FILES "${PROGRAM}.cpp")

  set_property(TARGET ${PROGRAM} PROPERTY VS_DEBUGGER_WORKING_DIRECTORY "$<TARGET_FILE_DIR:${PROGRAM}>")

  target_compile_definitions(${PROGRAM} PUBLIC "CHRONO_DATA_DIR=\"${CHRONO_DATA_DIR}\"") 
  target_compile_definitions(${PROGRAM} PUBLIC "CHRONO_VEHICLE_DATA_DIR=\"${CHRONO_VEHICLE_DATA_DIR}\"")
  target_link_libraries(${PROGRAM} PRIVATE ${CHRONO_TARGETS})
endforeach()
