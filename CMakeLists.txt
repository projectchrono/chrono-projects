#--------------------------------------------------------------
# List of all executables

SET(DEMOS
    demo_crater
    demo_massflow
    demo_foam
    demo_soilbin
    demo_wheel
    demo_objectdrop
    rad_im_schmutz
    directShear
    pressureSinkage
    singleWheel
    test_mesh
    test_constraints
    test_utils
    test_collision
    test_scaling
    test_timing
    demo_suspension
    test_fill
    demo_shearDEM
    demo_shearDVI
    test_serial
    test_DVI
)

SET(DEMOS_VEHICLE
    hmmwv
    hmmwvDEM
    hmmwvDEM_ditch
)

#--------------------------------------------------------------
# Additional include paths and libraries

INCLUDE_DIRECTORIES(
    ${CMAKE_SOURCE_DIR}
    ${CMAKE_BINARY_DIR}
)

SET(LIBRARIES
    ${CHRONO_LIB}
    ChronoEngine_Parallel
    ChronoEngine_Utils
)

IF (ENABLE_OPENGL)
    INCLUDE_DIRECTORIES(
        ${GLM_INCLUDE_DIRS}
        ${GLEW_INCLUDE_DIR}
        ${GLFW_INCLUDE_DIR}
    )

    SET(LIBRARIES ${LIBRARIES} ChronoEngine_OpenGL)
ENDIF()

#--------------------------------------------------------------
# Always use full RPATH (differentiating between the build and install trees)

# use, i.e. don't skip the full RPATH for the build tree
SET(CMAKE_SKIP_BUILD_RPATH  FALSE)

# when building, don't use the install RPATH already
# (but later on when installing)
SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE) 

SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

# add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# the RPATH to be used when installing, but only if it's not a system directory
LIST(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
IF("${isSystemDir}" STREQUAL "-1")
   SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
ENDIF("${isSystemDir}" STREQUAL "-1")

#--------------------------------------------------------------

IF("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
	SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
ENDIF()

IF(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
	SET (CH_LINKERFLAG_EXE  "${CH_LINKERFLAG_EXE} -framework IOKit -framework Cocoa -framework OpenGL")
ENDIF()


#--------------------------------------------------------------

MESSAGE(STATUS "Radu demo programs...")

FOREACH(PROGRAM ${DEMOS})
    MESSAGE(STATUS "...add ${PROGRAM}")

    ADD_EXECUTABLE(${PROGRAM}  "${PROGRAM}.cpp")
    SOURCE_GROUP(""  FILES "${PROGRAM}.cpp")

    SET_TARGET_PROPERTIES(${PROGRAM} PROPERTIES
        FOLDER demos
        COMPILE_FLAGS "${CH_BUILDFLAGS}"
        LINK_FLAGS "${CH_LINKERFLAG_EXE}"
    )

    TARGET_LINK_LIBRARIES(${PROGRAM} ${LIBRARIES})
    ADD_DEPENDENCIES(${PROGRAM} ${LIBRARIES})

    INSTALL(TARGETS ${PROGRAM} DESTINATION bin)

ENDFOREACH(PROGRAM)

#--------------------------------------------------------------

IF(ENABLE_VEHICLE)

  MESSAGE(STATUS "Vehicle demo programs...")

  INCLUDE_DIRECTORIES(${CHRONOVEHICLE_INCLUDES})

  SET(LIBRARIES ${LIBRARIES} ${CHRONOVEHICLE_LIBRARIES})

  FOREACH(PROGRAM ${DEMOS_VEHICLE})
    MESSAGE(STATUS "...add ${PROGRAM}")

    ADD_EXECUTABLE(${PROGRAM}  "${PROGRAM}.cpp")
    SOURCE_GROUP(""  FILES "${PROGRAM}.cpp")

    SET_TARGET_PROPERTIES(${PROGRAM} PROPERTIES
        FOLDER demos
        COMPILE_FLAGS "${CH_BUILDFLAGS}"
        LINK_FLAGS "${CH_LINKERFLAG_EXE}"
    )

    TARGET_LINK_LIBRARIES(${PROGRAM} ${LIBRARIES})
    ADD_DEPENDENCIES(${PROGRAM} ${LIBRARIES})

    INSTALL(TARGETS ${PROGRAM} DESTINATION bin)

  ENDFOREACH(PROGRAM)

ENDIF()
