option(ENABLE_PROJECTS "Build the project programs" OFF)

if(ENABLE_PROJECTS)
  message(STATUS "\n==== Chrono projects ====")
else()
  return()
endif()

#-----------------------------------------------------------------------------
# Invoke CMake in subdirectories
#-----------------------------------------------------------------------------

cmake_dependent_option(ENABLE_PRJ_DEM_VALIDATION 
	"Enable projects for DEM validation" ON
	"ENABLE_PROJECTS" OFF)
#cmake_dependent_option(ENABLE_PRJ_CONSTRAINT_FLUID 
#	"Enable projects using constraint fluids, require zlib" OFF
#	"ENABLE_PROJECTS" OFF)
cmake_dependent_option(ENABLE_PRJ_MULTICORE_TESTS 
	"Enable projects for Chrono::Multicore tests" ON
	"ENABLE_PROJECTS" OFF)
cmake_dependent_option(ENABLE_PRJ_VEHICLE_TESTS 
	"Enable projects for Chrono::Vehicle tests" ON
	"ENABLE_PROJECTS" OFF)
cmake_dependent_option(ENABLE_PRJ_CORE_TESTS 
	"Enable projects for Chrono core tests" ON
	"ENABLE_PROJECTS" OFF)
cmake_dependent_option(ENABLE_PRJ_PHYSICS_TESTS 
	"Enable projects for basic physics tests" ON
	"ENABLE_PROJECTS" OFF)
cmake_dependent_option(ENABLE_PRJ_MBS_TESTS
	"Enable projects for MBS simulation" ON
	"ENABLE_PROJECTS" OFF)
cmake_dependent_option(ENABLE_PRJ_FEA_TESTS 
	"Enable projects for Chrono::FEA tests" ON
	"ENABLE_PROJECTS" OFF)
cmake_dependent_option(ENABLE_PRJ_COSIMULATION 
	"Enable projects for Chrono::Vehicle cosimulation" ON
	"ENABLE_PROJECTS" OFF)
cmake_dependent_option(ENABLE_PRJ_GPU_TESTS 
	"Enable projects for Chrono::GPU granular tests" ON
	"ENABLE_PROJECTS" OFF)
cmake_dependent_option(ENABLE_PRJ_SYN_HIGHWAY 
	"Enable SynChrono Highway test" OFF
	"ENABLE_PROJECTS" OFF)

if (ENABLE_PRJ_CORE_TESTS)
  add_subdirectory(core_tests)
endif()

if (ENABLE_PRJ_PHYSICS_TESTS)
  add_subdirectory(physics_tests)
endif()

if (ENABLE_PRJ_MBS_TESTS)
  add_subdirectory(mbs)
endif()

if (ENABLE_PRJ_FEA_TESTS)
  add_subdirectory(fea_tests)
endif()

if (ENABLE_PRJ_DEM_VALIDATION)
  add_subdirectory(dem_validation)
endif()

if (ENABLE_PRJ_MULTICORE_TESTS)
  add_subdirectory(multicore_tests)
endif()

if (ENABLE_PRJ_VEHICLE_TESTS)
  add_subdirectory(vehicle_tests)
endif()

if (ENABLE_PRJ_GPU_TESTS)
  add_subdirectory(gpu_tests)
endif()

if (ENABLE_PRJ_COSIMULATION)
  add_subdirectory(cosimulation)
endif()

#if (ENABLE_PRJ_CONSTRAINT_FLUID)
#  add_subdirectory(constraint_fluids)
#endif()

if (ENABLE_PRJ_SYN_HIGHWAY)
  add_subdirectory(synchrono_highway)
endif()

set(ALL_DLL_NAMES       "${ALL_DLL_NAMES}"       PARENT_SCOPE)
set(ALL_DEPENDENCY_DLLS "${ALL_DEPENDENCY_DLLS}" PARENT_SCOPE)
