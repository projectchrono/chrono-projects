#=============================================================================
# CMake configuration file for projects using DEM-Engine standalone library.
#=============================================================================

message("\n---- DEME projects")

#--------------------------------------------------------------
# List of all executables
#--------------------------------------------------------------

set(TESTS
    test_DEME_BallDrop_Cosim
)
set(TESTS_VEH_ROBOT
    test_Rover_on_GRC_Incline
    test_Viper_on_SphericalTerrain
)

#--------------------------------------------------------------
# Find the Chrono package with required components
#--------------------------------------------------------------

# Invoke find_package in CONFIG mode

find_package(Chrono QUIET
             OPTIONAL_COMPONENTS Vehicle
             CONFIG
)

find_package(CUDAToolkit REQUIRED)

find_package(ChPF QUIET)
if(NOT ChPF_FOUND)
  message("Could not find ChPF library: It is fine if DEME is built without it; otherwise, find search DEME install dir for it.")
endif()

find_package(DEME REQUIRED)


# If Chrono and/or the required component(s) were not found, return now.
if(NOT Chrono_FOUND)
  message("Could not find Chrono, search Chrono installation dir for it.")
  return()
endif()

if(NOT DEME_FOUND)
  message("Could not find DEM-Engine library, search DEME install dir for it.")
  return()
endif()

if(CHRONO_VEHICLE_FOUND)
  message("Found Chrono::Vehicle, more demos are enabled.")
    list(APPEND TESTS ${TESTS_VEH_ROBOT})
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
    ${DEME_INCLUDE_DIRS}
    ${CMAKE_SOURCE_DIR}
)

#--------------------------------------------------------------
# Append to the parent's list of DLLs (and make it visible up)
#--------------------------------------------------------------

list(APPEND ALL_DLLS "${CHRONO_DLLS}")
set(ALL_DLLS "${ALL_DLLS}" PARENT_SCOPE)

#--------------------------------------------------------------
# DEME requires C++17 to be linked against
#--------------------------------------------------------------

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")

#--------------------------------------------------------------
# Loop over all demo programs and build them
#--------------------------------------------------------------

MESSAGE(STATUS "Physics programs...")

FOREACH(PROGRAM ${TESTS})
    MESSAGE(STATUS "...add ${PROGRAM}")

    ADD_EXECUTABLE(${PROGRAM}  "${PROGRAM}.cpp")
    SOURCE_GROUP(""  FILES "${PROGRAM}.cpp")
    set_property(TARGET ${PROGRAM} PROPERTY C_STANDARD 17)
    set_property(TARGET ${PROGRAM} PROPERTY C_STANDARD_REQUIRED ON)

    target_compile_definitions(${PROGRAM} PUBLIC "CHRONO_DATA_DIR=\"${CHRONO_DATA_DIR}\"") 
    target_compile_definitions(${PROGRAM} PUBLIC "DEME_DATA_DIR=\"${DEME_DATA_DIRS}\"") 
    target_compile_options(${PROGRAM} PUBLIC ${CHRONO_CXX_FLAGS})
    target_link_options(${PROGRAM} PUBLIC ${CHRONO_LINKER_FLAGS})
    target_link_libraries(${PROGRAM} ${CHRONO_LIBRARIES})
    target_link_libraries(${PROGRAM} DEME::simulator_multi_gpu)

ENDFOREACH()

