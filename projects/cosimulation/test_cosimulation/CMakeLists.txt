#=============================================================================
# CMake configuration file for the vehicle cosimulation tests.
# These example programs requires MPI and the Chrono::Multicore module.
#=============================================================================

IF(NOT MPI_CXX_FOUND OR NOT ENABLE_MODULE_MULTICORE)
    RETURN()
ENDIF()

#--------------------------------------------------------------
# List all tests

SET(TEST_PROGRAMS
  test_VEH_HMMWV_Cosimulation
)

#--------------------------------------------------------------
# Additional include directories

include_directories(${CH_MULTICORE_INCLUDES})
include_directories(${MPI_CXX_INCLUDE_PATH})

#--------------------------------------------------------------
# List of all required libraries

SET(LIBRARIES
    ChronoEngine
    ChronoEngine_multicore
    ChronoEngine_vehicle
    ChronoModels_vehicle
    ${MPI_CXX_LIBRARIES}
)

#--------------------------------------------------------------
# Create the executables

FOREACH(PROGRAM ${TEST_PROGRAMS})
  MESSAGE(STATUS "...add ${PROGRAM}")

  ADD_EXECUTABLE(${PROGRAM} ${PROGRAM}.cpp)
  SOURCE_GROUP(""  FILES  "${PROGRAM}.cpp")

  target_compile_definitions(${PROGRAM} PUBLIC "CHRONO_DATA_DIR=\"${CHRONO_DATA_DIR}\"") 
  target_compile_options(${PROGRAM} PUBLIC ${CHRONO_CXX_FLAGS})
  target_link_options(${PROGRAM} PUBLIC ${CHRONO_LINKER_FLAGS})
                      
  target_link_libraries(${PROGRAM} ${LIBRARIES})

ENDFOREACH()