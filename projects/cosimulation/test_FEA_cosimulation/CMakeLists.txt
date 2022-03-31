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

  add_executable(${PROGRAM} ${PROGRAM}.cpp)
  source_group(""  FILES  "${PROGRAM}.cpp")

  target_compile_definitions(${PROGRAM} PUBLIC "CHRONO_DATA_DIR=\"${CHRONO_DATA_DIR}\"") 
  target_compile_options(${PROGRAM} PUBLIC ${CHRONO_CXX_FLAGS})
  target_link_options(${PROGRAM} PUBLIC ${CHRONO_LINKER_FLAGS})
  target_link_libraries(${PROGRAM} ${CHRONO_LIBRARIES})

endforeach()
