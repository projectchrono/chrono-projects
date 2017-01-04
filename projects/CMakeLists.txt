option(ENABLE_PROJECTS "Build the project programs" OFF)

if(ENABLE_PROJECTS)
  message(STATUS "==== Chrono projects ====\n")
endif()

#-----------------------------------------------------------------------------
# Invoke CMake in subdirectories
#-----------------------------------------------------------------------------

cmake_dependent_option(ENABLE_PRJ_GRANULAR_VALIDATION 
	"Enable projects for granular dynamics validation" ON
	"ENABLE_PROJECTS" OFF)
cmake_dependent_option(ENABLE_PRJ_CONSTRAINT_FLUID 
	"Enable projects using constraint fluids, require zlib" OFF
	"ENABLE_PROJECTS" OFF)
cmake_dependent_option(ENABLE_PRJ_PARALLEL_TESTS 
	"Enable projects for Chrono::Parallel tests" ON
	"ENABLE_PROJECTS" OFF)
cmake_dependent_option(ENABLE_PRJ_VEHICLE_TESTS 
	"Enable projects for Chrono::Vehicle tests" ON
	"ENABLE_PROJECTS" OFF)
cmake_dependent_option(ENABLE_PRJ_CORE_TESTS 
	"Enable projects for Chrono::Engine core tests" ON
	"ENABLE_PROJECTS" OFF)
cmake_dependent_option(ENABLE_PRJ_PHYSICS_TESTS 
	"Enable projects for Chrono::Engine physics tests" ON
	"ENABLE_PROJECTS" OFF)
cmake_dependent_option(ENABLE_PRJ_FEA_TESTS 
	"Enable projects for Chrono::FEA tests" ON
	"ENABLE_PROJECTS" OFF)

if (ENABLE_PRJ_CORE_TESTS)
  add_subdirectory(core_tests)
endif()

if (ENABLE_PRJ_PHYSICS_TESTS)
  add_subdirectory(physics_tests)
endif()

if (ENABLE_PRJ_FEA_TESTS)
  add_subdirectory(fea_tests)
endif()

if (ENABLE_PRJ_GRANULAR_VALIDATION)
  add_subdirectory(granular_validation)
endif()

if (ENABLE_PRJ_PARALLEL_TESTS)
  add_subdirectory(parallel_tests)
endif()

if (ENABLE_PRJ_VEHICLE_TESTS)
  add_subdirectory(vehicle_tests)
endif()

if (ENABLE_CONSTRAINT_FLUID_MODELS)
  add_subdirectory(constraint_fluids)
endif()

set(ALL_DLLS "${ALL_DLLS}" PARENT_SCOPE)
