#-----------------------------------------------------------------------------
# Invoke CMake in subdirectories
#-----------------------------------------------------------------------------

add_subdirectory(fea)
add_subdirectory(fsi)
add_subdirectory(irrlicht)
add_subdirectory(parallel)
add_subdirectory(python)
add_subdirectory(cascade)
add_subdirectory(granular)

set(ALL_DLLS "${ALL_DLLS}" PARENT_SCOPE)
