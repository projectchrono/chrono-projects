option(ENABLE_CONFIG_TESTS "Build the programs for checking Chrono configuration" OFF)

if(ENABLE_CONFIG_TESTS)
  message(STATUS "\n==== Chrono configuration programs ====")
else()
  return()
endif()

#-----------------------------------------------------------------------------
# Invoke CMake in subdirectories
#-----------------------------------------------------------------------------

add_subdirectory(fea)
add_subdirectory(fsi)
add_subdirectory(mbs)
add_subdirectory(multicore)
add_subdirectory(parsers)
add_subdirectory(urdf)
add_subdirectory(cascade)
add_subdirectory(gpu)
add_subdirectory(synchrono)
add_subdirectory(vehicle)
add_subdirectory(vsg)

set(ALL_DLL_NAMES       "${ALL_DLL_NAMES}"       PARENT_SCOPE)
set(ALL_DEPENDENCY_DLLS "${ALL_DEPENDENCY_DLLS}" PARENT_SCOPE)
