#=============================================================================
# CMake configuration file for projects requiring the Chrono CASCADE and
# IRRLICHT modules.
# 
# Cannot be used stand-alone (but is mostly self-contained).
#=============================================================================

message(STATUS "\n---- CASCADE configuration test\n")

#--------------------------------------------------------------
# List of all executables
#--------------------------------------------------------------

set(DEMOS
    demo_CAS_robot
)

#--------------------------------------------------------------
# Find the Chrono package with required and optional components
#--------------------------------------------------------------

# Invoke find_package in CONFIG mode
find_package(Chrono
             COMPONENTS Cascade Irrlicht VSG
             CONFIG
)

# If Chrono and/or the required component(s) were not found, return now.
if(NOT Chrono_FOUND)
  message("Could not find requirements for CASCADE projects")
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
# Compilation flags
#--------------------------------------------------------------

# Needed for OpenCASCADE
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")	
    add_definitions("/DWNT")
endif()

#--------------------------------------------------------------
# Loop over all demo programs and build them
#--------------------------------------------------------------

message(STATUS "Configuration tests for Chrono::Cascade...")

foreach(PROGRAM ${DEMOS})
  message(STATUS "...add ${PROGRAM}")

  add_executable(${PROGRAM}  "${PROGRAM}.cpp")
  source_group(""  FILES "${PROGRAM}.cpp")

  set_property(TARGET ${PROGRAM} PROPERTY VS_DEBUGGER_WORKING_DIRECTORY "$<TARGET_FILE_DIR:${PROGRAM}>")

  target_compile_definitions(${PROGRAM} PUBLIC "CHRONO_DATA_DIR=\"${CHRONO_DATA_DIR}\"") 
  target_link_libraries(${PROGRAM} PRIVATE ${CHRONO_TARGETS})
endforeach()
