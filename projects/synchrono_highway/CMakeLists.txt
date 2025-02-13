#=============================================================================
# CMake configuration file for the SynChrono Highway project - an extended
#   demo of a highway environment, interfaced to a human-driven simulator
#=============================================================================

message("\n---- SYNCHRONO HIGHWAY projects")

#--------------------------------------------------------------
# List of all executables
#--------------------------------------------------------------

set(DEMOS
    test_SYN_city
    test_SYN_highway
)

#--------------------------------------------------------------
# Find the Chrono package with required and optional components
#--------------------------------------------------------------

# Invoke find_package in CONFIG mode
find_package(Chrono QUIET
             COMPONENTS SynChrono Vehicle Sensor Irrlicht
             CONFIG
)

# If Chrono and/or the required component(s) were not found, return now.
if(NOT Chrono_FOUND)
  message("Could not find requirements for the SynChrono Highway Project")
  return()
endif()

# If Synchrono was not built with FastDDS support, return now.
if(NOT CHRONO_SYNCHRONO_FASTDDS)
  message("Synchrono built without FastDDS support!")
  return()
endif()

# Get the path to the highway mesh
set(HIGHWAY_DATA_DIR "" CACHE FILEPATH "Path to directory contain highway assets")
add_compile_definitions(HIGHWAY_DATA_DIR=${HIGHWAY_DATA_DIR})

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

# Disable some warnings triggered by Irrlicht (Windows only)
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    add_compile_options(/wd4251) 
endif()

#--------------------------------------------------------------
# Grab the extras
#--------------------------------------------------------------
set(EXTRAS
	# extras/driver/ChCSLDriver.h
	# extras/driver/ChCSLDriver.cpp
	extras/driver/ChLidarWaypointDriver.h
	extras/driver/ChLidarWaypointDriver.cpp
	# extras/driver/joystick.h
	extras/ChCSLSoundEngine.h
	extras/ChCSLSoundEngine.cpp
)

#--------------------------------------------------------------
# Loop over all demo programs and build them
#--------------------------------------------------------------

message(STATUS "Projects for SynChrono Highway...")

foreach(PROGRAM ${DEMOS})
  message(STATUS "...add ${PROGRAM}")

  add_executable(${PROGRAM}  "${PROGRAM}.cpp" ${EXTRAS})
  source_group(""  FILES "${PROGRAM}.cpp")

  set_property(TARGET ${PROGRAM} PROPERTY VS_DEBUGGER_WORKING_DIRECTORY "$<TARGET_FILE_DIR:${PROGRAM}>")

  target_compile_definitions(${PROGRAM} PUBLIC "CHRONO_DATA_DIR=\"${CHRONO_DATA_DIR}\"") 
	target_link_libraries(${PROGRAM} ${CHRONO_TARGETS} "-L/usr/local/cuda/lib64")# -lcudart")
endforeach()
