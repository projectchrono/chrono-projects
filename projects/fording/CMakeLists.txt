#=============================================================================
# CMake configuration file for projects requiring the Chrono PARALLEL module
# 
# Cannot be used stand-alone (but is mostly self-contained).
#=============================================================================

SET(MODEL_FILES
    hmmwv/HMMWV.h
    hmmwv/HMMWV.cpp
    hmmwv/HMMWV_Vehicle.h
    hmmwv/HMMWV_Vehicle.cpp
    hmmwv/HMMWV_VehicleReduced.h
    hmmwv/HMMWV_VehicleReduced.cpp
    hmmwv/HMMWV_DoubleWishbone.h
    hmmwv/HMMWV_DoubleWishbone.cpp
    hmmwv/HMMWV_DoubleWishboneReduced.h
    hmmwv/HMMWV_DoubleWishboneReduced.cpp
    hmmwv/HMMWV_PitmanArm.h
    hmmwv/HMMWV_PitmanArm.cpp
    hmmwv/HMMWV_RackPinion.h
    hmmwv/HMMWV_RackPinion.cpp
    hmmwv/HMMWV_Driveline2WD.h
    hmmwv/HMMWV_Driveline2WD.cpp
    hmmwv/HMMWV_Driveline4WD.h
    hmmwv/HMMWV_Driveline4WD.cpp
    hmmwv/HMMWV_SimplePowertrain.h
    hmmwv/HMMWV_SimplePowertrain.cpp
    hmmwv/HMMWV_Powertrain.h
    hmmwv/HMMWV_Powertrain.cpp
    hmmwv/HMMWV_BrakeSimple.h
    hmmwv/HMMWV_BrakeSimple.cpp
    hmmwv/HMMWV_Wheel.h
    hmmwv/HMMWV_Wheel.cpp
    hmmwv/HMMWV_RigidTire.h
    hmmwv/HMMWV_RigidTire.cpp
)

#--------------------------------------------------------------
# List of all executables
#--------------------------------------------------------------

set(DEMOS
    fording_fluid
)

#--------------------------------------------------------------
# Find the Chrono package with required components
#--------------------------------------------------------------

# Invoke find_package in CONFIG mode.

find_package(Chrono
             COMPONENTS Parallel
             OPTIONAL_COMPONENTS OpenGL 
             CONFIG
)

# If Chrono and/or the required component(s) were not found, return now.

if(NOT Chrono_FOUND)
  message("Could not find requirements for PARALLEL projects")
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

message(STATUS "Parallel tests...")

foreach(PROGRAM ${DEMOS})

  message(STATUS "...add ${PROGRAM}")

  add_executable(${PROGRAM}  "${MODEL_FILES}" "${PROGRAM}.cpp")
  source_group(""  FILES "${PROGRAM}.cpp")

  set_target_properties(${PROGRAM} PROPERTIES
    FOLDER demos
    COMPILE_FLAGS "${CHRONO_CXX_FLAGS}"
    COMPILE_DEFINITIONS "CHRONO_DATA_DIR=\"${CHRONO_DATA_DIR}\""
    LINK_FLAGS "${CHRONO_CXX_FLAGS} ${CHRONO_LINKER_FLAGS}"
  )

  target_link_libraries(${PROGRAM} ${CHRONO_LIBRARIES})

endforeach(PROGRAM)

