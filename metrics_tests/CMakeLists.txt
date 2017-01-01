#-----------------------------------------------------------------------------
# Invoke CMake in subdirectories
#-----------------------------------------------------------------------------

add_subdirectory(core)
add_subdirectory(fea)
#add_subdirectory(vehicle)
add_subdirectory(parallel)

set(ALL_DLLS "${ALL_DLLS}" PARENT_SCOPE)
