#=============================================================================
# CMake configuration file for projects using the ChronoEngine main module.
# 
# Cannot be used stand-alone (but is mostly self-contained).
#=============================================================================

message("\n---- CORE projects")

#--------------------------------------------------------------
# List of all executables
#--------------------------------------------------------------

set(TESTS
    test_CH_function_recorder
    test_CH_integratorHHT
    test_CH_vector
)

#--------------------------------------------------------------
# Find the Chrono package with required components
#--------------------------------------------------------------

# Invoke find_package in CONFIG mode

find_package(Chrono QUIET
             OPTIONAL_COMPONENTS PardisoMKL
             CONFIG
)

# If Chrono and/or the required component(s) were not found, return now.

if(NOT Chrono_FOUND)
  message("Could not find requirements for CORE projects")
  return()
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

MESSAGE(STATUS "Physics programs...")

FOREACH(PROGRAM ${TESTS})
    MESSAGE(STATUS "...add ${PROGRAM}")

    ADD_EXECUTABLE(${PROGRAM}  "${PROGRAM}.cpp")
    SOURCE_GROUP(""  FILES "${PROGRAM}.cpp")

    target_compile_definitions(${PROGRAM} PUBLIC "CHRONO_DATA_DIR=\"${CHRONO_DATA_DIR}\"") 
    target_compile_options(${PROGRAM} PUBLIC ${CHRONO_CXX_FLAGS})
    target_link_options(${PROGRAM} PUBLIC ${CHRONO_LINKER_FLAGS})
    target_link_libraries(${PROGRAM} ${CHRONO_LIBRARIES})

ENDFOREACH()

