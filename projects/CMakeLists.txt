#-----------------------------------------------------------------------------
# Invoke CMake in subdirectories
#-----------------------------------------------------------------------------


option(ENABLE_CONSTRAINT_FLUID_MODELS "Enable the models using constraint fluids, require zlib" OFF)

add_subdirectory(granular_validation)
add_subdirectory(parallel_tests)
add_subdirectory(vehicle_tests)
add_subdirectory(physics_tests)
add_subdirectory(core_tests)

IF (ENABLE_CONSTRAINT_FLUID_MODELS)
   add_subdirectory(constraint_fluids)
ENDIF()

set(ALL_DLLS "${ALL_DLLS}" PARENT_SCOPE)
