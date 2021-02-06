#-----------------------------------------------------------------------------
# Invoke CMake in subdirectories
#-----------------------------------------------------------------------------

add_subdirectory(fea)
#add_subdirectory(vehicle)
add_subdirectory(multicore)

set(ALL_DLLS "${ALL_DLLS}" PARENT_SCOPE)
