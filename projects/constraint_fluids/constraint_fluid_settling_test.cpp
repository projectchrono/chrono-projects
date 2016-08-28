#include <stdio.h>
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
ChFluidContainer* fluid_container;
std::ofstream ofile;
real tolerance = 1e-8;
std::string tolerance_str = "1e-8";

int out_frame = 0;
int next_out_frame = 0;

void WriteData(ChSystemParallelDVI* msystem, uint i) {
	double iters = ((ChIterativeSolverParallel*)(msystem->GetSolverSpeed()))->GetTotalIterations();
	const std::vector<double>& vhist = ((ChIterativeSolverParallel*)(msystem->GetSolverSpeed()))->GetViolationHistory();
	const std::vector<double>& dhist = ((ChIterativeSolverParallel*)(msystem->GetSolverSpeed()))->GetDeltalambdaHistory();
	double residual = vhist.size() > 0 ? vhist.back() : 0.0;
	double dlambda = dhist.size() > 0 ? dhist.back() : 0.0;

	std::vector<real> density, pressure;
	std::vector<real3> force;
	std::vector<chrono::real> force_length;
	std::vector<chrono::real> marker_height;

	fluid_container->GetFluidDensity(density);
	fluid_container->GetFluidPressure(pressure);
	fluid_container->GetFluidForce(force);

	force_length.resize(density.size());
	marker_height.resize(density.size());

	real avg_density = 0;
	real avg_pressure = 0;
	real KE = 0;
	for (int i = 0; i < density.size(); i++) {
		avg_density += density[i];
		avg_pressure += pressure[i];
		marker_height[i] = msystem->data_manager->host_data.pos_3dof[i].z;
		force_length[i] = Length(force[i]);
		KE += .5 * fluid_container->mass * Pow(Length(msystem->data_manager->host_data.vel_3dof[i]), 2);
	}
	avg_density = avg_density / density.size();
	avg_pressure = avg_pressure / pressure.size();

	printf("%d %d %f %f avg D, avg P, KE: %f %f %f \n", i, int(iters), residual, dlambda, avg_density, avg_pressure,
		KE);

	ofile << i << " " << int(iters) << " " << residual << " " << dlambda << " " << avg_density << " " << avg_pressure
		<< " " << KE << std::endl;
}

void AddContainer(ChSystemParallelDVI* sys) {
	// Create a common material
	auto mat = std::make_shared<ChMaterialSurface>();
	mat->SetFriction(0.4f);
	ChVector<> hdim(.55, .6, .55);
	utils::CreateBoxContainer(sys, 0, mat, hdim, 0.05, Vector(0, 0, -hdim.z - 0.05 * 4), QUNIT, true, false, true,
		true);
}

// -----------------------------------------------------------------------------
// Create the fluid in the shape of a sphere.
// -----------------------------------------------------------------------------
void AddFluid(ChSystemParallelDVI* sys) {
	fluid_container = new ChFluidContainer(sys);

	fluid_container->tau = timestep * 2;
	fluid_container->epsilon = 1e-8;
	fluid_container->rho = 1000;
	fluid_container->contact_cohesion = 0;

	fluid_container->kernel_radius = .016 * 2;
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
	real3 num_fluid = real3(10, 10, 10);
	real3 origin(0, 0, -.2);
	real vol;

	std::vector<real3> pos_fluid;
	std::vector<real3> vel_fluid;

	double dist = fluid_container->kernel_radius * .9;
	utils::HCPSampler<> sampler(dist);
	vol = dist * dist * dist * .6;
#if 0
	utils::Generator::PointVector points = sampler.SampleBox(ChVector<>(0, 0, 0), radius);
	// vol = 4.0 / 3.0 * CH_C_PI * pow(radius, 3) / real(points.size());

#else
	ChVector<> hdim(.5, .5, .5);
	utils::Generator::PointVector points = sampler.SampleBox(ChVector<>(0, 0, 0), hdim);
	// vol = hdim.x * hdim.y * hdim.z / real(points.size());
#endif

	pos_fluid.resize(points.size());
	vel_fluid.resize(points.size());
	for (int i = 0; i < points.size(); i++) {
		pos_fluid[i] = real3(points[i].x, points[i].y, points[i].z) + origin;
		vel_fluid[i] = real3(0, 0, 0);
	}

	fluid_container->mass = (fluid_container->rho + 400) * vol;
	std::cout << "fluid_mass: " << fluid_container->mass << std::endl;
	fluid_container->UpdatePosition(0);
	fluid_container->AddBodies(pos_fluid, vel_fluid);
}

double tolerances[] = { 100, 10, 1, 1e-1, 1e-3, 1e-5, 1e-8 };
std::string tol_string[] = { "100", "10", "1", "1e-1", "1e-3", "1e-5", "1e-8" };
// -----------------------------------------------------------------------------
// Create the system, specify simulation parameters, and run simulation loop.
// -----------------------------------------------------------------------------
int main(int argc, char* argv[]) {
	double gravity = 9.81;
	double time_end = 10;
	bool useAPGD = false;
	uint max_iteration = 10000;
	if (argc == 3) {
		tolerance = tolerances[atoi(argv[1])];
		tolerance_str = tol_string[atoi(argv[1])];
		useAPGD = atoi(argv[2]);
	}

	// Create system
	// -------------

	ChSystemParallelDVI msystem;
	// Set gravitational acceleration
	msystem.Set_G_acc(ChVector<>(0, 0, -gravity));

	// Set solver parameters
	msystem.GetSettings()->solver.solver_mode = SLIDING;
	msystem.GetSettings()->solver.max_iteration_normal = 0;
	msystem.GetSettings()->solver.max_iteration_sliding = max_iteration;
	msystem.GetSettings()->solver.max_iteration_spinning = 0;
	msystem.GetSettings()->solver.max_iteration_bilateral = 0;
	msystem.GetSettings()->solver.tolerance = tolerance;
	msystem.GetSettings()->solver.alpha = 0;
	msystem.GetSettings()->solver.use_full_inertia_tensor = false;
	msystem.GetSettings()->solver.contact_recovery_speed = 100000;
	msystem.GetSettings()->solver.cache_step_length = true;
	msystem.GetSettings()->min_threads = 8;
	if (useAPGD) {
		msystem.ChangeSolverType(APGD);
	}
	else {
		msystem.ChangeSolverType(BB);
	}

	msystem.GetSettings()->collision.narrowphase_algorithm = NARROWPHASE_HYBRID_MPR;

	AddFluid(&msystem);

	msystem.GetSettings()->collision.collision_envelope = (fluid_container->kernel_radius * .05);
	msystem.GetSettings()->collision.bins_per_axis = vec3(2, 2, 2);
	msystem.SetLoggingLevel(LOG_TRACE, true);
	msystem.SetLoggingLevel(LOG_INFO, true);
	// Create the fixed and moving bodies
	// ----------------------------------
	AddContainer(&msystem);

	// Perform the simulation
	// ----------------------
	ofile.open("data_fluid_test_" + tolerance_str + "/fluid_container_data_" + tolerance_str + ".txt");
#if 0
	opengl::ChOpenGLWindow& gl_window = opengl::ChOpenGLWindow::getInstance();
	gl_window.Initialize(1280, 720, "fluidDVI", &msystem);
	gl_window.SetCamera(ChVector<>(0, -2, 0), ChVector<>(0, 0, 0), ChVector<>(0, 0, 1), .2);
	gl_window.Pause();
	// Uncomment the following two lines for the OpenGL manager to automatically
	// run the simulation in an infinite loop.
	// gl_window.StartDrawLoop(timestep);
	// return 0;
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
	int out_fps = 60;
	double time = 0;
	int sim_frame = 0;

	double exec_time = 0;
	int num_contacts = 0;
	int out_steps = std::ceil((1.0 / timestep) / out_fps);
	int num_steps = time_end / timestep;

	//=========================================================================================================
	int file = 0;
	for (int i = 0; i < num_steps; i++) {
		msystem.DoStepDynamics(timestep);
		std::cout << i << std::endl;
		WriteData(&msystem, i);
		if (i == next_out_frame) {
			std::cout << "write: " << out_frame << std::endl;
			out_frame++;
			next_out_frame += out_steps;
		}

		exec_time += msystem.GetTimerStep();
	}
	ofile << exec_time << std::endl;
	printf("Execution Time: %f \n", exec_time);
#endif

	ofile.close();

	return 0;
}
