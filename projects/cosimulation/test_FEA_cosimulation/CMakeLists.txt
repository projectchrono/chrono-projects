#=============================================================================
# CMake configuration file for FEA cosimulation tests.
#=============================================================================

if(NOT CHRONO_IRRLICHT_FOUND)
  return()
endif()

#--------------------------------------------------------------
# List all tests

set(TEST_PROGRAMS
  test_FEA_cosimulate_granular
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
