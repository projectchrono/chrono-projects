#include <cstdio>
#include <vector>
#include <cmath>
#include "chrono_parallel/physics/ChSystemParallel.h"
#include "chrono_parallel/solver/ChIterativeSolverParallel.h"

#include "chrono/ChConfig.h"
#include "chrono/utils/ChUtilsCreators.h"
#include "chrono/utils/ChUtilsInputOutput.h"
#include "chrono/utils/ChUtilsGeometry.h"
#include "chrono/utils/ChUtilsGenerators.h"

#include "chrono_parallel/physics/Ch3DOFContainer.h"

#ifdef CHRONO_OPENGL
#include "chrono_opengl/ChOpenGLWindow.h"
#endif

using namespace chrono;
using namespace chrono::collision;

double timestep = 1e-3;
int solver_type = 2;
ChFluidContainer* fluid_container;
std::ofstream ofile;




int setup = 23;
real scale[] = { 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4 ,4, 4 };
real iters[] = { 50, 100, 200, 400, 800, 1600,  50, 100, 200, 400, 800, 1600,  50, 100, 200, 400, 800, 1600,  50, 100, 200, 400, 800, 1600 };
real tolerance =100000;
uint max_iteration = iters[setup];
int mult = scale[setup];
double kernel_radius = .016 * 2;
double thickness = 0.05;

ChVector<> hdim(kernel_radius * 7, kernel_radius * 7, kernel_radius * 7 * mult);

//std::string tolerance_str = tol_string[2];

//test 1
//5k, 10k, 15k, 20k
//constant iteration count, APGD, BB, Jacobi

//test 2
//5k, 10k, 15k, 20k
//constant objective, APGD, BB, Jacobi

//store objective function
//store time taken to solve


void WriteData(ChSystemParallelNSC* msystem, uint i) {
    auto iter_solver = std::static_pointer_cast<ChIterativeSolverParallel>(msystem->GetSolver());
    int iters = iter_solver->GetIterations();
    const std::vector<double>& vhist = iter_solver->GetViolationHistory();
	const std::vector<double>& dhist = iter_solver->GetDeltalambdaHistory();
	real residual = vhist.size() > 0 ? vhist.back() : 0.0;
	real dlambda = dhist.size() > 0 ? dhist.back() : 0.0;

	std::vector<real> density, pressure;
	//std::vector<real3> force;
	//std::vector<chrono::real> force_length;
	//std::vector<chrono::real> marker_height;

	fluid_container->GetFluidDensity(density);
	fluid_container->GetFluidPressure(pressure);
	//fluid_container->GetFluidForce(force);

	//force_length.resize(density.size());
	//marker_height.resize(density.size());

	real avg_density = 0;
	real avg_pressure = 0;
	real KE = 0;
	for (int i = 0; i < density.size(); i++) {
		avg_density += density[i];
		avg_pressure += pressure[i];
		//marker_height[i] = msystem->data_manager->host_data.pos_3dof[i].z;
		//force_length[i] = Length(force[i]);
		KE += .5 * fluid_container->mass * Pow(Length(msystem->data_manager->host_data.vel_3dof[i]), 2);
	}
	avg_density = avg_density / density.size();
	avg_pressure = avg_pressure / pressure.size();
	double SOLVER = msystem->GetTimerSolver();
	printf("%d %d %f %f avg D: %f, avg P:  %f, KE: %f SolverTime: %f\n", i, iters, residual, dlambda, avg_density, avg_pressure,
		KE, SOLVER);

	ofile << i << " " << iters << " " << residual << " " << dlambda << " " << avg_density << " " << avg_pressure
		<< " " << KE <<" "<<SOLVER<<std::endl;

	if (i % 100==0) {
		std::ofstream cfile;
		switch (solver_type) {
		case 0:
			cfile.open("data_fluid_test/fluid_container_data_BB_CONV_"+std::to_string(setup) + "_" + std::to_string(i) +"_.txt");
			break;
		case 1:
			cfile.open("data_fluid_test/fluid_container_data_APGD_CONV_" + std::to_string(setup) + "_" + std::to_string(i) + "_.txt");
			break;
		case 2:
			cfile.open("data_fluid_test/fluid_container_data_JACOBI_CONV_" + std::to_string(setup) + "_" + std::to_string(i) + "_.txt");
			break;
		case 3:
			cfile.open("data_fluid_test/fluid_container_data_GAUSS_SEIDEL_CONV_" + std::to_string(setup) + "_" + std::to_string(i) + "_.txt");
			break;
		default:
			cfile.open("data_fluid_test/fluid_container_data_APGD_CONV_" + std::to_string(setup) + "_" + std::to_string(i) + "_.txt");
			break;
		}

		for (int k = 0; k < vhist.size(); k++) {

			cfile << vhist[k] << " " << dhist[k] << std::endl;
		}
		cfile.close();
	}

}

void AddContainer(ChSystemParallelNSC* sys) {
	// Create a common material
	auto mat = chrono_types::make_shared<ChMaterialSurfaceNSC>();
	mat->SetFriction(1.0);
	utils::CreateBoxContainer(sys, 0, mat, hdim + ChVector<>(0, 0, kernel_radius * 4), thickness, Vector(0, 0, -hdim.z()- kernel_radius), QUNIT, true, false, true,
		true);
}

// -----------------------------------------------------------------------------
// Create the fluid in the shape of a sphere.
// -----------------------------------------------------------------------------
void AddFluid(ChSystemParallelNSC* sys) {
    auto fluid_container = chrono_types::make_shared<ChFluidContainer>();
    sys->Add3DOFContainer(fluid_container);

	fluid_container->tau = timestep * 2;
	fluid_container->epsilon = 1e-8;
	fluid_container->rho = 1000;
	fluid_container->contact_cohesion = 0;

	fluid_container->kernel_radius = kernel_radius;
	fluid_container->mass = .007 * 5.5;
	fluid_container->viscosity = .01;
	fluid_container->enable_viscosity = false;
	fluid_container->alpha = .1;
	fluid_container->contact_compliance = 1e-12;
	fluid_container->contact_mu = 0;

	// msystem.GetSettings()->fluid.max_interactions = 30;
	fluid_container->artificial_pressure = false;
	fluid_container->artificial_pressure_k = .01;
	fluid_container->artificial_pressure_dq = .2 * fluid_container->kernel_radius;
	fluid_container->artificial_pressure_n = 4;
	fluid_container->collision_envelope = 0;  // fluid_container->kernel_ra dius * .05;

	real radius = .2;  //*5
	real dens = 30;
	real3 origin(0, 0, 0);
	real vol;

	std::vector<real3> pos_fluid;
	std::vector<real3> vel_fluid;

	double dist = fluid_container->kernel_radius * .9;
	utils::HCPSampler<> sampler(dist);
	vol = dist * dist * dist * .6;

	utils::Generator::PointVector points = sampler.SampleBox(ChVector<>(0, 0, 0), hdim );

	pos_fluid.resize(points.size());
	vel_fluid.resize(points.size());
	for (int i = 0; i < points.size(); i++) {
		pos_fluid[i] = real3(points[i].x(), points[i].y(), points[i].z()) + origin;
		vel_fluid[i] = real3(0, 0, 0);
	}

	fluid_container->mass = (fluid_container->rho + 400) * vol;
	std::cout << "fluid_mass: " << fluid_container->mass << std::endl;
	fluid_container->UpdatePosition(0);
	fluid_container->AddBodies(pos_fluid, vel_fluid);
}


// -----------------------------------------------------------------------------
// Create the system, specify simulation parameters, and run simulation loop.
// -----------------------------------------------------------------------------
int main(int argc, char* argv[]) {
	double time_end = 5;
	
	
	if (argc == 3) {
		setup = atoi(argv[1]);
		solver_type = atoi(argv[2]);

		 max_iteration = iters[setup];
		 mult = scale[setup];
	}

	// Create system
	// -------------

	ChSystemParallelNSC msystem;
	// Set gravitational acceleration
	msystem.Set_G_acc(ChVector<>(0, 0, -9.81));

	// Set solver parameters
	msystem.GetSettings()->solver.solver_mode = SolverMode::SLIDING;
	msystem.GetSettings()->solver.max_iteration_normal = 0;
	msystem.GetSettings()->solver.max_iteration_sliding = max_iteration;
	msystem.GetSettings()->solver.max_iteration_spinning = 0;
	msystem.GetSettings()->solver.max_iteration_bilateral = 0;
	msystem.GetSettings()->solver.tolerance = tolerance;
	msystem.GetSettings()->solver.test_objective = true;
	msystem.GetSettings()->solver.tolerance_objective = -tolerance;
	msystem.GetSettings()->solver.alpha = 0;
	msystem.GetSettings()->solver.use_full_inertia_tensor = false;
	msystem.GetSettings()->solver.contact_recovery_speed = 100000;
	msystem.GetSettings()->solver.cache_step_length = true;
	msystem.GetSettings()->min_threads = 10;

	switch (solver_type) {
	case 0:
		msystem.ChangeSolverType(SolverType::BB);
		ofile.open("data_fluid_test/fluid_container_data_BB_" + std::to_string(setup) + ".txt");
		break;
	case 1:
		msystem.ChangeSolverType(SolverType::APGD);
		ofile.open("data_fluid_test/fluid_container_data_APGD_" + std::to_string(setup) + ".txt");
		break;
	case 2:
		msystem.ChangeSolverType(SolverType::JACOBI);
		ofile.open("data_fluid_test/fluid_container_data_JACOBI_" + std::to_string(setup) + ".txt");
		break;
	case 3:
		msystem.ChangeSolverType(SolverType::GAUSS_SEIDEL);
		ofile.open("data_fluid_test/fluid_container_data_GAUSS_SEIDEL_" + std::to_string(setup) + ".txt");
		break;
	default:
		msystem.ChangeSolverType(SolverType::APGD);
		ofile.open("data_fluid_test/fluid_container_data_APGD_" + std::to_string(setup) + ".txt");
		break;
	}

	msystem.GetSettings()->collision.narrowphase_algorithm = NarrowPhaseType::NARROWPHASE_HYBRID_MPR;

	AddFluid(&msystem);

	msystem.GetSettings()->collision.collision_envelope = (fluid_container->kernel_radius * .05);
	msystem.GetSettings()->collision.bins_per_axis = vec3(2, 2, 2);
	msystem.SetLoggingLevel(LoggingLevel::LOG_TRACE, true);
	msystem.SetLoggingLevel(LoggingLevel::LOG_INFO, true);
	// Create the fixed and moving bodies
	// ----------------------------------
	AddContainer(&msystem);

	// Perform the simulation
	// ----------------------
	
#if 0
	opengl::ChOpenGLWindow& gl_window = opengl::ChOpenGLWindow::getInstance();
	gl_window.Initialize(1280, 720, "fluid_settling_test", &msystem);
	gl_window.SetCamera(ChVector<>(0, -2, 0), ChVector<>(0, 0, 0), ChVector<>(0, 0, 1), .2);
	gl_window.Pause();

	while (true) {
		if (gl_window.Active()) {
			if (gl_window.DoStepDynamics(timestep)) {

			}
			gl_window.Render();
		}
		else {
			break;
		}
	}
#else
	double exec_time = 0;
	int num_steps = time_end / timestep;

	//=========================================================================================================
	for (int i = 0; i < num_steps; i++) {
		msystem.DoStepDynamics(timestep);
		WriteData(&msystem, i);
		exec_time += msystem.GetTimerStep();
	}
	ofile << exec_time << std::endl;
	printf("Execution Time: %f \n", exec_time);
#endif
	ofile.close();
	return 0;
}
