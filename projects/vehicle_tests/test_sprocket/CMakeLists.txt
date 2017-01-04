#--------------------------------------------------------------
# List of all executables
#--------------------------------------------------------------

SET(DEMOS
    test_VEH_sprocketProfile
)

#--------------------------------------------------------------
# Loop over all demo programs and build them
#--------------------------------------------------------------

FOREACH(PROGRAM ${DEMOS})
    MESSAGE(STATUS "...add ${PROGRAM}")

    ADD_EXECUTABLE(${PROGRAM}  "${PROGRAM}.cpp")
    SOURCE_GROUP(""  FILES "${PROGRAM}.cpp")

    SET_TARGET_PROPERTIES(${PROGRAM} PROPERTIES
        FOLDER demos
        COMPILE_FLAGS "${CHRONO_CXX_FLAGS}"
        COMPILE_DEFINITIONS "CHRONO_DATA_DIR=\"${CHRONO_DATA_DIR}\";CHRONO_VEHICLE_DATA_DIR=\"${CHRONO_VEHICLE_DATA_DIR}\""
        LINK_FLAGS "${CHRONO_CXX_FLAGS} ${CHRONO_LINKER_FLAGS}"
    )

    TARGET_LINK_LIBRARIES(${PROGRAM} ${CHRONO_LIBRARIES})

ENDFOREACH(PROGRAM)
