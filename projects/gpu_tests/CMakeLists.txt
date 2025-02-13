#=============================================================================
# CMake configuration file for projects requiring the Chrono::Gpu module.
# 
# Cannot be used stand-alone (but is mostly self-contained).
#=============================================================================

message(STATUS "\n---- GPU projects\n")

#--------------------------------------------------------------
# List of all executables
#--------------------------------------------------------------

set(TESTS
        test_GPU_DamBreak
        test_GPU_Bucket
        test_GPU_rover
        test_GPU_directshear
        test_GPU_Coneflow
        test_GPU_commvessels
        test_GPU_milsettle
        test_GPU_testsuite
        test_GPU_normal
        test_GPU_slide
        test_GPU_roll
        test_GPU_incline
        test_GPU_repose
        test_GPU_soilBin
        test_GPU_ballDrop
        test_GPU_compression
)

#--------------------------------------------------------------
# Find the Chrono package with required and optional components
#--------------------------------------------------------------

# Invoke find_package in CONFIG mode
find_package(Chrono QUIET
             COMPONENTS Gpu OpenGL
             OPTIONAL_COMPONENTS Vehicle
             CONFIG
)

# If Chrono and/or the required component(s) were not found, return now.
if(NOT Chrono_FOUND)
  message("Could not find requirements for GPU projects")
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

message(STATUS "Projects for Chrono::Gpu...")

foreach(PROGRAM ${TESTS})
  message(STATUS "...add ${PROGRAM}")

  add_executable(${PROGRAM}  "${PROGRAM}.cpp")
  source_group(""  FILES "${PROGRAM}.cpp")

  set_property(TARGET ${PROGRAM} PROPERTY VS_DEBUGGER_WORKING_DIRECTORY "$<TARGET_FILE_DIR:${PROGRAM}>")

  target_compile_definitions(${PROGRAM} PUBLIC "CHRONO_DATA_DIR=\"${CHRONO_DATA_DIR}\"") 
  target_link_libraries(${PROGRAM} PRIVATE ${CHRONO_TARGETS})
endforeach()

