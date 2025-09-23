if(ENABLE_CONFIG)
  message(STATUS "\n==== Chrono configuration programs ====")
  mark_as_advanced(CLEAR ENABLE_CONFIG_FEA)
  mark_as_advanced(CLEAR ENABLE_CONFIG_FSI)
  mark_as_advanced(CLEAR ENABLE_CONFIG_MBS)
  mark_as_advanced(CLEAR ENABLE_CONFIG_MULTICORE)
  mark_as_advanced(CLEAR ENABLE_CONFIG_PARSERS)
  mark_as_advanced(CLEAR ENABLE_CONFIG_URDF)
  mark_as_advanced(CLEAR ENABLE_CONFIG_CASCADE)
  mark_as_advanced(CLEAR ENABLE_CONFIG_DEM)
  mark_as_advanced(CLEAR ENABLE_CONFIG_SYNCHRONO)
  mark_as_advanced(CLEAR ENABLE_CONFIG_VEHICLE)
  mark_as_advanced(CLEAR ENABLE_CONFIG_VSG)
else()
  mark_as_advanced(FORCE ENABLE_CONFIG_FEA)
  mark_as_advanced(FORCE ENABLE_CONFIG_FSI)
  mark_as_advanced(FORCE ENABLE_CONFIG_MBS)
  mark_as_advanced(FORCE ENABLE_CONFIG_MULTICORE)
  mark_as_advanced(FORCE ENABLE_CONFIG_PARSERS)
  mark_as_advanced(FORCE ENABLE_CONFIG_URDF)
  mark_as_advanced(FORCE ENABLE_CONFIG_CASCADE)
  mark_as_advanced(FORCE ENABLE_CONFIG_DEM)
  mark_as_advanced(FORCE ENABLE_CONFIG_SYNCHRONO)
  mark_as_advanced(FORCE ENABLE_CONFIG_VEHICLE)
  mark_as_advanced(FORCE ENABLE_CONFIG_VSG)
  return()
endif()

#-----------------------------------------------------------------------------
# Invoke CMake in subdirectories
#-----------------------------------------------------------------------------

cmake_dependent_option(ENABLE_CONFIG_FEA
	"Enable configuration tests for FEA module" ON
	"ENABLE_CONFIG" OFF)
cmake_dependent_option(ENABLE_CONFIG_FSI
	"Enable configuration tests for FSI module" ON
	"ENABLE_CONFIG" OFF)
cmake_dependent_option(ENABLE_CONFIG_MBS
	"Enable configuration tests for MBS with Irrlicht module" ON
	"ENABLE_CONFIG" OFF)
cmake_dependent_option(ENABLE_CONFIG_MULTICORE
	"Enable configuration tests for MULTICORE module" ON
	"ENABLE_CONFIG" OFF)
cmake_dependent_option(ENABLE_CONFIG_PARSERS
	"Enable configuration tests for PARSERS module" ON
	"ENABLE_CONFIG" OFF)
cmake_dependent_option(ENABLE_CONFIG_URDF
	"Enable configuration tests for URDF PARSERS module" ON
	"ENABLE_CONFIG" OFF)
cmake_dependent_option(ENABLE_CONFIG_CASCADE
	"Enable configuration tests for CASCADE module" ON
	"ENABLE_CONFIG" OFF)
cmake_dependent_option(ENABLE_CONFIG_DEM
	"Enable configuration tests for DEM module" ON
	"ENABLE_CONFIG" OFF)
cmake_dependent_option(ENABLE_CONFIG_SYNCHRONO
	"Enable configuration tests for SYNCHRONO module" ON
	"ENABLE_CONFIG" OFF)
cmake_dependent_option(ENABLE_CONFIG_VEHICLE
	"Enable configuration tests for VEHICLE module" ON
	"ENABLE_CONFIG" OFF)
cmake_dependent_option(ENABLE_CONFIG_VSG
	"Enable configuration tests for VSG module" ON
	"ENABLE_CONFIG" OFF)

if (ENABLE_CONFIG_FEA)
  add_subdirectory(fea)
endif()
if (ENABLE_CONFIG_FSI)
  add_subdirectory(fsi)
endif()
if (ENABLE_CONFIG_MBS)
  add_subdirectory(mbs)
endif()
if (ENABLE_CONFIG_MULTICORE)
  add_subdirectory(multicore)
endif()
if (ENABLE_CONFIG_PARSERS)
  add_subdirectory(parsers)
endif()
if (ENABLE_CONFIG_URDF)
  add_subdirectory(urdf)
endif()
if (ENABLE_CONFIG_CASCADE)
  add_subdirectory(cascade)
endif()
if (ENABLE_CONFIG_DEM)
  add_subdirectory(dem)
endif()
if (ENABLE_CONFIG_SYNCHRONO)
  add_subdirectory(synchrono)
endif()
if (ENABLE_CONFIG_VEHICLE)
  add_subdirectory(vehicle)
endif()
if (ENABLE_CONFIG_VSG)
  add_subdirectory(vsg)
endif()

set(ALL_DLL_NAMES             "${ALL_DLL_NAMES}"           PARENT_SCOPE)
set(CHRONO_PACKAGE_PREFIX_DIR ${CHRONO_PACKAGE_PREFIX_DIR} PARENT_SCOPE)
