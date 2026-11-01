#=============================================================================
# CMake configuration file for projects requiring the Chrono FEA module
# and optionally the IRRLICHT and VSG modules
# 
# Cannot be used stand-alone (but is mostly self-contained).
#=============================================================================

message(STATUS "\n---- FEA configuration test\n")

#--------------------------------------------------------------
# List of all executables
#--------------------------------------------------------------

set(DEMOS_BASIC
    demo_FEA_basic
)

set(DEMOS_VIS
    demo_FEA_thermal
)

#--------------------------------------------------------------
# Find the Chrono package with required and optional components
#--------------------------------------------------------------

# Invoke find_package in CONFIG mode
find_package(Chrono
             COMPONENTS
             OPTIONAL_COMPONENTS Irrlicht VSG
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

# If one of the optional run-time visualization components was found, append to list
if(CHRONO_IRRLICHT_FOUND OR CHRONO_VSG_FOUND)
  set(DEMOS ${DEMOS} ${DEMOS_VIS})
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

message(STATUS "Configuration tests for FEA...")

foreach(PROGRAM ${DEMOS})
  message(STATUS "...add ${PROGRAM}")

  add_executable(${PROGRAM}  "${PROGRAM}.cpp" FEAvisualization.h)
  source_group(""  FILES "${PROGRAM}.cpp" FEAvisualization.h)

  set_property(TARGET ${PROGRAM} PROPERTY VS_DEBUGGER_WORKING_DIRECTORY "$<TARGET_FILE_DIR:${PROGRAM}>")

  target_compile_definitions(${PROGRAM} PUBLIC "CHRONO_DATA_DIR=\"${CHRONO_DATA_DIR}\"") 
  target_link_libraries(${PROGRAM} PRIVATE ${CHRONO_TARGETS})
endforeach()
