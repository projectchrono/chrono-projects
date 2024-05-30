#--------------------------------------------------------------
# List of all executables
#--------------------------------------------------------------

SET(DEMOS
    test_VEH_BMW_E90_wall2wall
)

#--------------------------------------------------------------
# Loop over all demo programs and build them
#--------------------------------------------------------------

FOREACH(PROGRAM ${DEMOS})
    MESSAGE(STATUS "...add ${PROGRAM}")

    ADD_EXECUTABLE(${PROGRAM}  "${PROGRAM}.cpp")
    SOURCE_GROUP(""  FILES "${PROGRAM}.cpp")

    target_compile_definitions(${PROGRAM} PUBLIC "CHRONO_DATA_DIR=\"${CHRONO_DATA_DIR}\"") 
    target_compile_definitions(${PROGRAM} PUBLIC "CHRONO_VEHICLE_DATA_DIR=\"${CHRONO_VEHICLE_DATA_DIR}\"") 
    target_compile_options(${PROGRAM} PUBLIC ${CHRONO_CXX_FLAGS})
    target_link_options(${PROGRAM} PUBLIC ${CHRONO_LINKER_FLAGS})
    target_link_libraries(${PROGRAM} ${CHRONO_LIBRARIES})

ENDFOREACH()
