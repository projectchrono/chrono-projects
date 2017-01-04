# Co-simulation tire test rig
# ------------------------------------------------------------------------------

set(PROGRAM test_VEH_tireRig_Cosimulation)

set(TEST_FILES 
    test_VEH_tireRig_Cosimulation.cpp
    BaseNode.h
    BaseNode.cpp
    RigNode.h
    RigNode.cpp
    TerrainNode.h
    TerrainNode.cpp)

source_group("" FILES ${TEST_FILES})

message(STATUS "...add ${PROGRAM}")

include_directories(${MPI_CXX_INCLUDE_PATH})
add_executable(${PROGRAM} ${TEST_FILES})
set_target_properties(${PROGRAM} PROPERTIES 
    COMPILE_FLAGS "${CHRONO_CXX_FLAGS} ${MPI_CXX_COMPILE_FLAGS}"
    COMPILE_DEFINITIONS "CHRONO_DATA_DIR=\"${CHRONO_DATA_DIR}\";CHRONO_VEHICLE_DATA_DIR=\"${CHRONO_VEHICLE_DATA_DIR}\""
    LINK_FLAGS "${CHRONO_LINKER_FLAGS} ${MPI_CXX_LINK_FLAGS}")
target_link_libraries(${PROGRAM} ${CHRONO_LIBRARIES} ${MPI_CXX_LIBRARIES})
