#-----------------------------------------------------------------------------
# Invoke CMake in subdirectories
#-----------------------------------------------------------------------------

add_subdirectory(fea)
#add_subdirectory(vehicle)
add_subdirectory(multicore)
add_subdirectory(granular)

set(ALL_DLLS "${ALL_DLLS}" PARENT_SCOPE)
