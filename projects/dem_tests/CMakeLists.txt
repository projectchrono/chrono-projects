#=============================================================================
# CMake configuration file for projects requiring the Chrono::Dem module.
# 
# Cannot be used stand-alone (but is mostly self-contained).
#=============================================================================

message(STATUS "\n---- DEM projects\n")

#--------------------------------------------------------------
# List of all executables
#--------------------------------------------------------------

set(TESTS
        test_DEM_DamBreak
        test_DEM_Bucket
        test_DEM_rover
        test_DEM_directshear
        test_DEM_Coneflow
        test_DEM_commvessels
        test_DEM_milsettle
        test_DEM_testsuite
        test_DEM_normal
        test_DEM_slide
        test_DEM_roll
        test_DEM_incline
        #test_DEM_repose
        test_DEM_soilBin
        test_DEM_ballDrop
        test_DEM_compression
)

#--------------------------------------------------------------
# Find the Chrono package with required and optional components
#--------------------------------------------------------------

# Invoke find_package in CONFIG mode
find_package(Chrono QUIET
             COMPONENTS DEM
             OPTIONAL_COMPONENTS Vehicle
             CONFIG
)

# If Chrono and/or the required component(s) were not found, return now.
if(NOT Chrono_FOUND)
  message("Could not find requirements for DEM projects")
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

message(STATUS "Projects for Chrono::Dem...")

foreach(PROGRAM ${TESTS})
  message(STATUS "...add ${PROGRAM}")

  add_executable(${PROGRAM}  "${PROGRAM}.cpp")
  source_group(""  FILES "${PROGRAM}.cpp")

  set_property(TARGET ${PROGRAM} PROPERTY VS_DEBUGGER_WORKING_DIRECTORY "$<TARGET_FILE_DIR:${PROGRAM}>")

  target_compile_definitions(${PROGRAM} PUBLIC "CHRONO_DATA_DIR=\"${CHRONO_DATA_DIR}\"") 
  target_link_libraries(${PROGRAM} PRIVATE ${CHRONO_TARGETS})
endforeach()

