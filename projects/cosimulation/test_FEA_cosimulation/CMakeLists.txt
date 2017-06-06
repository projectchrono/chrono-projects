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

  set_target_properties(${PROGRAM} PROPERTIES 
                        COMPILE_FLAGS "${CHRONO_CXX_FLAGS}"
                        COMPILE_DEFINITIONS "CHRONO_DATA_DIR=\"${CHRONO_DATA_DIR}\""
                        LINK_FLAGS "${CHRONO_CXX_FLAGS} ${CHRONO_LINKER_FLAGS}"
                        )
  target_link_libraries(${PROGRAM} ${CHRONO_LIBRARIES})

endforeach()
