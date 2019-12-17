# Test programs for HMMWV simulation
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Additional include paths and libraries
# ------------------------------------------------------------------------------

message(STATUS "test_GRAN_HMMWV temporarily unavailable")
return()

set(PROGRAMS
  test_GRAN_HMMWV
)

# ------------------------------------------------------------------------------
# Add all executables
# ------------------------------------------------------------------------------
message(STATUS "Test programs for GRANULAR HMMWV...")

foreach(PROGRAM ${PROGRAMS})
    message(STATUS "...add ${PROGRAM}")

    add_executable(${PROGRAM}  "${PROGRAM}.cpp" ${MODEL_FILES})
    source_group(""  FILES "${PROGRAM}.cpp" ${MODEL_FILES})

    set_target_properties(${PROGRAM} PROPERTIES
        FOLDER tests
        COMPILE_FLAGS "${CH_CXX_FLAGS} ${CH_PARALLEL_CXX_FLAGS} ${CH_GRANULAR_CXX_FLAGS}"
        COMPILE_DEFINITIONS "CHRONO_DATA_DIR=\"${CHRONO_DATA_DIR}\";CHRONO_VEHICLE_DATA_DIR=\"${CHRONO_VEHICLE_DATA_DIR}\";PROJECTS_DATA_DIR=\"${PROJECTS_DATA_DIR}\""
        LINK_FLAGS "${CH_LINKERFLAG_EXE}"
    )
    target_link_libraries(${PROGRAM} ${LIBRARIES})
	# add_dependencies(${PROGRAM} ${LIBRARIES})
endforeach(PROGRAM)