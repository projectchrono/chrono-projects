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

foreach(PROGRAM ${TEST_PROGRAMS})
  message(STATUS "...add ${PROGRAM}")

  add_executable(${PROGRAM}  "${PROGRAM}.cpp")
  source_group(""  FILES "${PROGRAM}.cpp")

  set_property(TARGET ${PROGRAM} PROPERTY VS_DEBUGGER_WORKING_DIRECTORY "$<TARGET_FILE_DIR:${PROGRAM}>")

  target_compile_definitions(${PROGRAM} PUBLIC "CHRONO_DATA_DIR=\"${CHRONO_DATA_DIR}\"") 
    target_compile_definitions(${PROGRAM} PUBLIC "CHRONO_VEHICLE_DATA_DIR=\"${CHRONO_VEHICLE_DATA_DIR}\"") 
  target_link_libraries(${PROGRAM} PRIVATE ${CHRONO_TARGETS})
endforeach()
