#--------------------------------------------------------------
# List of all executables
#--------------------------------------------------------------

set(TESTS
    test_VEH_FEDA_SSC
    test_VEH_FEDA_RideTool
    test_VEH_FEDA_wall2wall
)

#--------------------------------------------------------------
# Loop over all demo programs and build them
#--------------------------------------------------------------

foreach(PROGRAM ${TESTS})
  message(STATUS "...add ${PROGRAM}")

  add_executable(${PROGRAM}  "${PROGRAM}.cpp")
  source_group(""  FILES "${PROGRAM}.cpp")

  set_property(TARGET ${PROGRAM} PROPERTY VS_DEBUGGER_WORKING_DIRECTORY "$<TARGET_FILE_DIR:${PROGRAM}>")

  target_compile_definitions(${PROGRAM} PUBLIC "CHRONO_DATA_DIR=\"${CHRONO_DATA_DIR}\"") 
  target_link_libraries(${PROGRAM} PRIVATE ${CHRONO_TARGETS})
endforeach()
