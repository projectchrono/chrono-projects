#=============================================================================
# CMake configuration file for projects using the Chrono::FEA module.
# 
# Cannot be used stand-alone (but is mostly self-contained).
#=============================================================================

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

find_package(Chrono
             COMPONENTS
             OPTIONAL_COMPONENTS Irrlicht Postprocess MKL MUMPS Matlab
             CONFIG
)

# If Chrono and/or the required component(s) were not found, return now.

if(NOT Chrono_FOUND)
  message("Could not find requirements for FEA projects")
  return()
endif()

# If the Chrono::Irrlicht and Chrono::MKL are available.

if(CHRONO_IRRLICHT_FOUND AND CHRONO_MKL_FOUND)
  list(APPEND TESTS ${TESTS_MKL_IRR})
endif()

if (CHRONO_IRRLICHT_FOUND AND CHRONO_MATLAB_FOUND)
	list(APPEND TESTS ${TESTS_MATLAB_IRR})
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
# Append to the parent's list of DLLs (and make it visible up)
#--------------------------------------------------------------

list(APPEND ALL_DLLS "${CHRONO_DLLS}")
set(ALL_DLLS "${ALL_DLLS}" PARENT_SCOPE)

#--------------------------------------------------------------
# Loop over all demo programs and build them
#--------------------------------------------------------------

MESSAGE(STATUS "FEA programs...")

FOREACH(PROGRAM ${TESTS})
    MESSAGE(STATUS "...add ${PROGRAM}")

    ADD_EXECUTABLE(${PROGRAM}  "${PROGRAM}.cpp")
    SOURCE_GROUP(""  FILES "${PROGRAM}.cpp")

    SET_TARGET_PROPERTIES(${PROGRAM} PROPERTIES
        FOLDER demos
        COMPILE_FLAGS "${CHRONO_CXX_FLAGS}"
        COMPILE_DEFINITIONS "CHRONO_DATA_DIR=\"${CHRONO_DATA_DIR}\""
        LINK_FLAGS "${CHRONO_CXX_FLAGS} ${CHRONO_LINKER_FLAGS}"
    )

    TARGET_LINK_LIBRARIES(${PROGRAM} ${CHRONO_LIBRARIES})

ENDFOREACH(PROGRAM)

message(STATUS "")
