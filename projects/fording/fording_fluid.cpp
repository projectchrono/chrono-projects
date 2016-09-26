#include "hmmwv.h"
#include "vehicle.h"
#include "fluid.h"
int threads = 8;
void SetupSystem(ChSystemParallelDVI* system) {
    system->Set_G_acc(ChVector<>(0, 0, -9.81));
    system->SetIntegrationType(ChSystem::INT_ANITESCU);

    system->GetSettings()->min_threads = threads;

    system->GetSettings()->solver.tolerance = tolerance;
    system->GetSettings()->solver.solver_mode = SLIDING;
    system->GetSettings()->solver.max_iteration_normal = max_iteration_normal;
    system->GetSettings()->solver.max_iteration_sliding = max_iteration_sliding;
    system->GetSettings()->solver.max_iteration_spinning = max_iteration_spinning;
    system->GetSettings()->solver.max_iteration_bilateral = 1000;  // make 1000, should be about 220
    system->GetSettings()->solver.max_iteration_fem = 50;
    system->GetSettings()->solver.compute_N = false;
    system->GetSettings()->solver.alpha = 0;
    system->GetSettings()->solver.cache_step_length = true;
    system->GetSettings()->solver.use_full_inertia_tensor = false;
    system->GetSettings()->solver.contact_recovery_speed = contact_recovery_speed;
    system->GetSettings()->solver.bilateral_clamp_speed = 1e8;
    system->ChangeSolverType(BB);
    system->SetLoggingLevel(LOG_INFO);
    system->SetLoggingLevel(LOG_TRACE);

    if (tire_type == FEM) {
        system->GetSettings()->collision.collision_envelope = dof_container->kernel_radius * .01;
    } else {
        system->GetSettings()->collision.collision_envelope = dof_container->kernel_radius * .1;
    }
    system->GetSettings()->collision.bins_per_axis = vec3(100, 20, 25);
    system->GetSettings()->collision.narrowphase_algorithm = NARROWPHASE_HYBRID_MPR;
    system->GetSettings()->collision.fixed_bins = true;
}

// =============================================================================
int main(int argc, char* argv[]) {
    if (argc == 5) {
        simulation_type = static_cast<SIM_TYPE>(atoi(argv[1]));
        tire_type = static_cast<TIRE_TYPE>(atoi(argv[2]));
        stiffness = atof(argv[3]);
        threads = atoi(argv[4]);
    }

    if (tire_type == FEM) {
        container_friction = 1.0;
        tire_mu = .6;
    } else {
        tire_mu = 0.8;
        container_friction = 0.8;
    }

    sim_system = new ChSystemParallelDVI();

    if (simulation_type == SIMFLUID || simulation_type == SIMRIGID) {
        printf("Not using MPM\n");
        rho = 207;
        contact_mu = 0;
    } else if (simulation_type == SIMMPM_FLUID || simulation_type == SIMMPM_RIGID) {
        printf("Using MPM\n");
        rho = 1200;
        contact_mu = .8;
    }

    if (simulation_type == SIMFLUID || simulation_type == SIMMPM_FLUID) {
        printf("Using Fluid\n");
        dof_container = new ChFluidContainer(sim_system);
    } else if (simulation_type == SIMRIGID || simulation_type == SIMMPM_RIGID) {
        printf("Using Rigid\n");
        dof_container = new Ch3DOFRigidContainer(sim_system);
    }

    if (tire_type == FEM) {
        fea_container = new ChFEAContainer(sim_system);

        fea_container->kernel_radius = .015;
        fea_container->material_density = 250;
        fea_container->contact_mu = tire_mu;
        fea_container->contact_cohesion = 0;
        fea_container->alpha = 0.001;
        fea_container->contact_compliance = 1e-7;
        fea_container->beta = 100;
        fea_container->youngs_modulus = 1e7;  // 2e8;
        fea_container->poisson_ratio = .33;
        fea_container->contact_recovery_speed = contact_recovery_speed;
        fea_container->rigid_constraint_recovery_speed = 5;  // not used right now
        fea_container->max_velocity = 10;
        fea_container->SetFamily(6, 4);
    }

    SetupSystem(sim_system);

    CreateContainer(sim_system);
    MyVehicle vehicle(sim_system, "HUMVEE_HULLS.obj");
    CreateFluid(sim_system);
    sim_system->Initialize();

    double time = 0, exec_time = 0;
    int sim_frame = 0, out_frame = 0, next_out_frame = 0;

    forces.resize(10);
    torques.resize(10);
    std::fill(forces.begin(), forces.end(), real3(0));
    std::fill(torques.begin(), torques.end(), real3(0));

#if 0
    opengl::ChOpenGLWindow& gl_window = opengl::ChOpenGLWindow::getInstance();
    gl_window.Initialize(1280, 720, "Fording Simulation", sim_system);
    gl_window.SetCamera(ChVector<>(0, -7.4, 3), ChVector<>(0, -6.4, 2.84), ChVector<>(0, 0, 1), 0.1f);
    gl_window.Pause();

    while (time < time_end) {
        if (gl_window.Active()) {
            if (gl_window.DoStepDynamics(time_step)) {
                vehicle.Update(time);
                time += time_step;
            }
            gl_window.Render();
        } else {
            exit(0);
        }
    }
    return 0;
#else
    std::stringstream sstiff;
    sstiff << stiffness;

    std::string data_output_path = "";
    switch (simulation_type) {
        case (SIMFLUID):
            data_output_path = "data_fluid/";
            break;
        case (SIMMPM_FLUID):
            data_output_path = sstiff.str() + "_data_mpm_fluid/";
            break;
        case (SIMRIGID):
            data_output_path = "data_rigid_no_env/";
            break;
        case (SIMMPM_RIGID):
            data_output_path = sstiff.str() + "_data_mpm_rigid_no_env/";
            break;
    }
    if (tire_type == FEM) {
        data_output_path = "tire_" + data_output_path;
    }
    printf("Writing To: %s \n", data_output_path.c_str());

    while (time < time_end) {
        sim_system->DoStepDynamics(time_step);
        if (sim_frame == next_out_frame) {
            std::cout << "write: " << out_frame << std::endl;
            DumpFluidData(sim_system, data_output_path + "data_" + std::to_string(out_frame) + ".dat", true);
            if (tire_type == FEM) {
                DumpTriangleData(sim_system, data_output_path + "tire_" + std::to_string(out_frame) + ".obj");
            }
            DumpAllObjectsWithGeometryPovray(sim_system,
                                             data_output_path + "vehicle_" + std::to_string(out_frame) + ".dat", true);
            WriteVehicleData(vehicle.my_hmmwv, throttle, braking, forces, torques,
                             data_output_path + "stats_" + std::to_string(out_frame) + ".dat");
            out_frame++;
            next_out_frame += out_steps;
        }
        vehicle.Update(time);

        // Update counters.
        time += time_step;
        sim_frame++;
        exec_time += sim_system->GetTimerStep();
    }
#endif
    cout << "==================================" << endl;
    cout << "Simulation time:   " << exec_time << endl;

    return 0;
}
