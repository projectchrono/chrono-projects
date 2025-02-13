#=============================================================================
# CMake configuration file for projects using the Chrono::FEA module.
# 
# Cannot be used stand-alone (but is mostly self-contained).
#=============================================================================

message(STATUS "\n---- FEA projects\n")

#--------------------------------------------------------------
# List of all executables
#--------------------------------------------------------------

set(TESTS
    test_FEA_shellANCF
    test_FEA_tireCorotational
    test_FEA_tireANCF
    test_FEA_constraints
)

set(TESTS_MKL_IRR
    test_FEA_4HumveeTires
    test_FEA_4HumveeTires_LuGre
    test_FEA_tireANCF_LuGre
    test_FEA_abaqus_wheels
)

set(TESTS_MATLAB_IRR
    test_FEA_miniature
)

#--------------------------------------------------------------
# Find the Chrono package with required components
#--------------------------------------------------------------

# Invoke find_package in CONFIG mode
find_package(Chrono QUIET
             COMPONENTS
             OPTIONAL_COMPONENTS Irrlicht Postprocess PardisoMKL Matlab
             CONFIG
)

# If Chrono and/or the required component(s) were not found, return now.
if(NOT Chrono_FOUND)
  message("Could not find requirements for FEA projects")
  return()
endif()

# If the Chrono::Irrlicht and Chrono::PardisoMKL are available.
if(CHRONO_IRRLICHT_FOUND AND CHRONO_PARDISOMKL_FOUND)
  list(APPEND TESTS ${TESTS_MKL_IRR})
endif()

if (CHRONO_IRRLICHT_FOUND AND CHRONO_MATLAB_FOUND)
	list(APPEND TESTS ${TESTS_MATLAB_IRR})
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

message(STATUS "FEA programs...")

foreach(PROGRAM ${TESTS})
  message(STATUS "...add ${PROGRAM}")

  add_executable(${PROGRAM}  "${PROGRAM}.cpp")
  source_group(""  FILES "${PROGRAM}.cpp")

  set_property(TARGET ${PROGRAM} PROPERTY VS_DEBUGGER_WORKING_DIRECTORY "$<TARGET_FILE_DIR:${PROGRAM}>")

  target_compile_definitions(${PROGRAM} PUBLIC "CHRONO_DATA_DIR=\"${CHRONO_DATA_DIR}\"") 
  target_link_libraries(${PROGRAM} PRIVATE ${CHRONO_TARGETS})
endforeach()
