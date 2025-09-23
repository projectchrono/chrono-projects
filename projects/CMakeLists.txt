if(ENABLE_PROJECTS)
  message(STATUS "\n==== Chrono projects ====")
  #mark_as_advanced(CLEAR ENABLE_PROJECTS_DEM_VALIDATION)
  mark_as_advanced(CLEAR ENABLE_PROJECTS_CONSTRAINT_FLUID)
  mark_as_advanced(CLEAR ENABLE_PROJECTS_MULTICORE_TESTS)
  mark_as_advanced(CLEAR ENABLE_PROJECTS_VEHICLE_TESTS)
  mark_as_advanced(CLEAR ENABLE_PROJECTS_CORE_TESTS)
  mark_as_advanced(CLEAR ENABLE_PROJECTS_PHYSICS_TESTS)
  mark_as_advanced(CLEAR ENABLE_PROJECTS_MBS_TESTS)
  mark_as_advanced(CLEAR ENABLE_PROJECTS_FEA_TESTS)
  mark_as_advanced(CLEAR ENABLE_PROJECTS_COSIMULATION)
  mark_as_advanced(CLEAR ENABLE_PROJECTS_DEM_TESTS)
  mark_as_advanced(CLEAR ENABLE_PROJECTS_SYN_HIGHWAY)
else()
  mark_as_advanced(FORCE ENABLE_PROJECTS_DEM_VALIDATION)
  #mark_as_advanced(FORCE ENABLE_PROJECTS_CONSTRAINT_FLUID)
  mark_as_advanced(FORCE ENABLE_PROJECTS_MULTICORE_TESTS)
  mark_as_advanced(FORCE ENABLE_PROJECTS_VEHICLE_TESTS)
  mark_as_advanced(FORCE ENABLE_PROJECTS_CORE_TESTS)
  mark_as_advanced(FORCE ENABLE_PROJECTS_PHYSICS_TESTS)
  mark_as_advanced(FORCE ENABLE_PROJECTS_MBS_TESTS)
  mark_as_advanced(FORCE ENABLE_PROJECTS_FEA_TESTS)
  mark_as_advanced(FORCE ENABLE_PROJECTS_COSIMULATION)
  mark_as_advanced(FORCE ENABLE_PROJECTS_DEM_TESTS)
  mark_as_advanced(FORCE ENABLE_PROJECTS_SYN_HIGHWAY)
  return()
endif()

#-----------------------------------------------------------------------------
# Invoke CMake in subdirectories
#-----------------------------------------------------------------------------

cmake_dependent_option(ENABLE_PROJECTS_DEM_VALIDATION 
	"Enable projects for DEM validation" ON
	"ENABLE_PROJECTS" OFF)
#cmake_dependent_option(ENABLE_PROJECTS_CONSTRAINT_FLUID 
#	"Enable projects using constraint fluids, require zlib" OFF
#	"ENABLE_PROJECTS" OFF)
#cmake_dependent_option(ENABLE_PROJECTS_MULTICORE_TESTS 
#	"Enable projects for Chrono::Multicore tests" ON
#	"ENABLE_PROJECTS" OFF)
cmake_dependent_option(ENABLE_PROJECTS_VEHICLE_TESTS 
	"Enable projects for Chrono::Vehicle tests" ON
	"ENABLE_PROJECTS" OFF)
cmake_dependent_option(ENABLE_PROJECTS_CORE_TESTS 
	"Enable projects for Chrono core tests" ON
	"ENABLE_PROJECTS" OFF)
cmake_dependent_option(ENABLE_PROJECTS_PHYSICS_TESTS 
	"Enable projects for basic physics tests" ON
	"ENABLE_PROJECTS" OFF)
cmake_dependent_option(ENABLE_PROJECTS_MBS_TESTS
	"Enable projects for MBS simulation" ON
	"ENABLE_PROJECTS" OFF)
cmake_dependent_option(ENABLE_PROJECTS_FEA_TESTS 
	"Enable projects for Chrono::FEA tests" ON
	"ENABLE_PROJECTS" OFF)
cmake_dependent_option(ENABLE_PROJECTS_COSIMULATION 
	"Enable projects for Chrono::Vehicle cosimulation" ON
	"ENABLE_PROJECTS" OFF)
cmake_dependent_option(ENABLE_PROJECTS_DEM_TESTS 
	"Enable projects for Chrono::DEM granular tests" ON
	"ENABLE_PROJECTS" OFF)
cmake_dependent_option(ENABLE_PROJECTS_SYN_HIGHWAY 
	"Enable SynChrono Highway test" OFF
	"ENABLE_PROJECTS" OFF)

if (ENABLE_PROJECTS_CORE_TESTS)
  add_subdirectory(core_tests)
endif()

if (ENABLE_PROJECTS_PHYSICS_TESTS)
  add_subdirectory(physics_tests)
endif()

if (ENABLE_PROJECTS_MBS_TESTS)
  add_subdirectory(mbs)
endif()

if (ENABLE_PROJECTS_FEA_TESTS)
  add_subdirectory(fea_tests)
endif()

if (ENABLE_PROJECTS_DEM_VALIDATION)
  add_subdirectory(dem_validation)
endif()

#if (ENABLE_PROJECTS_MULTICORE_TESTS)
#  add_subdirectory(multicore_tests)
#endif()

if (ENABLE_PROJECTS_VEHICLE_TESTS)
  add_subdirectory(vehicle_tests)
endif()

if (ENABLE_PROJECTS_DEM_TESTS)
  add_subdirectory(dem_tests)
endif()

if (ENABLE_PROJECTS_COSIMULATION)
  add_subdirectory(cosimulation)
endif()

#if (ENABLE_PROJECTS_CONSTRAINT_FLUID)
#  add_subdirectory(constraint_fluids)
#endif()

if (ENABLE_PROJECTS_SYN_HIGHWAY)
  add_subdirectory(synchrono_highway)
endif()

set(ALL_DLL_NAMES             "${ALL_DLL_NAMES}"           PARENT_SCOPE)
set(CHRONO_PACKAGE_PREFIX_DIR ${CHRONO_PACKAGE_PREFIX_DIR} PARENT_SCOPE)
