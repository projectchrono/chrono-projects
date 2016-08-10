#-----------------------------------------------------------------------------
# Invoke CMake in subdirectories
#-----------------------------------------------------------------------------

add_subdirectory(granular_validation)
add_subdirectory(parallel_tests)
add_subdirectory(vehicle_tests)

set(ALL_DLLS "${ALL_DLLS}" PARENT_SCOPE)
