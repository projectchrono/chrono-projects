// =============================================================================
// PROJECT CHRONO - http://projectchrono.org
//
// Copyright (c) 2014 projectchrono.org
// All right reserved.
//
// Use of this source code is governed by a BSD-style license that can be found
// in the LICENSE file at the top level of the distribution and at
// http://projectchrono.org/license-chrono.txt.
//
// =============================================================================
// Authors: Daniel Melanz, Radu Serban
// =============================================================================
//
// Chrono::Multicore demo program for mass flow rate studies.
//
// The model simulated here consists of a granular material that flows out of a
// container and the mass of the collected material is measured over time, using
// either penalty or complementarity method for frictional contact.
//
// The global reference frame has Z up.
// All units SI.
// =============================================================================

#include <cstdio>
#include <vector>
#include <cmath>

#include "chrono/ChConfig.h"
#include "chrono/utils/ChUtilsCreators.h"
#include "chrono/utils/ChUtilsGenerators.h"
#include "chrono/utils/ChUtilsInputOutput.h"

#include "chrono_multicore/physics/ChSystemMulticore.h"
#include "chrono_multicore/solver/ChSystemDescriptorMulticore.h"

#include "chrono_thirdparty/filesystem/path.h"

// Control use of OpenGL run-time rendering
// Note: CHRONO_OPENGL is defined in ChConfig.h
////#undef CHRONO_OPENGL

#ifdef CHRONO_OPENGL
#include "chrono_opengl/ChVisualSystemOpenGL.h"
#endif

using namespace chrono;

using std::cout;
using std::flush;
using std::endl;

// -----------------------------------------------------------------------------
// Problem definitions
// -----------------------------------------------------------------------------

// Comment the following line to use NSC contact
#define USE_SMC

enum ProblemType { SETTLING, DROPPING };

ProblemType problem = SETTLING;

// -----------------------------------------------------------------------------
// Global problem definitions
// -----------------------------------------------------------------------------

// Desired number of OpenMP threads (will be clamped to maximum available)
int threads = 8;

// Simulation parameters
double gravity = 9.81;

double time_settling_min = 0.1;
double time_settling_max = 1.0;
double time_dropping_max = 6.0;

#ifdef USE_SMC
double time_step = 1e-5;
int max_iteration = 20;
#else
double time_step = 1e-4;
int max_iteration_normal = 0;
int max_iteration_sliding = 50000;
int max_iteration_spinning = 0;
float contact_recovery_speed = 1.0e30;
#endif

double tolerance = 500.0;

int max_iteration_bilateral = 0;

// Output
#ifdef USE_SMC
const std::string out_dir = "../MASSFLOW_SMC";
#else
const std::string out_dir = "../MASSFLOW_NSC";
#endif

const std::string pov_dir = out_dir + "/POVRAY";
const std::string flow_file = out_dir + "/flow.dat";
const std::string stats_file = out_dir + "/stats.dat";
const std::string checkpoint_file = out_dir + "/settled.dat";

int out_fps_settling = 200;
int out_fps_dropping = 200;

int timing_frame = -1;  // output detailed step timing at this frame

// Parameters for the granular material
double r_g = 0.25e-3;
double rho_g = 2500.0;

float Y_g = 1e7f;
float cr_g = 0.1f;
float mu_g = 0.3f;

// Desired number of particles and X-Y dimensions of the sampling volume for
// granular material.
unsigned int desired_num_particles = 400;

// Parameters for the mechanism material
float Y_c = 2e6f;
float cr_c = 0.1f;
float mu_c = 0.4f;

// Dimensions of mechanism
double height = 6.0e-2;  // height of the cavity
////double width = 0.9525e-2;     // width of the cavity
double width = 0.9398e-2;    // width of the cavity
double thickness = 0.25e-2;  // thickness of walls

double height_insert = sqrt(2 * height * height);
double delta = sqrt(thickness * thickness / 8);

////double speed = 1.5e-3;       // speed of the angled insert
double speed = 1.0e-3;  // speed of the angled insert
double gap = 2e-3;      // size of gap

double time_opening = gap / speed;

// Dimensions of collector
double pos_collector = 4.0e-2;     // position below measuring line
double size_collector = 8.0e-2;    // width and length of collector bin
double height_collector = 1.0e-2;  // height of collector walls

// -----------------------------------------------------------------------------
// Create mechanism
// -----------------------------------------------------------------------------
ChBody* CreateMechanism(ChSystemMulticore* system) {
    // Create the common material
#ifdef USE_SMC
    auto mat_b = chrono_types::make_shared<ChContactMaterialSMC>();
    mat_b->SetYoungModulus(Y_c);
    mat_b->SetFriction(mu_c);
    mat_b->SetRestitution(cr_c);
#else
    auto mat_b = chrono_types::make_shared<ChContactMaterialNSC>();
    mat_b->SetFriction(mu_c);
#endif

    // Angled insert
    auto insert = chrono_types::make_shared<ChBody>();

    insert->SetIdentifier(0);
    insert->SetMass(1);
    insert->SetInertiaXX(ChVector3d(1, 1, 1));
    insert->SetPos(ChVector3d(-0.5 * height - delta, 0, 0.5 * height - delta));
    insert->SetRot(chrono::QuatFromAngleAxis(-CH_C_PI / 4, VECT_Y));
    insert->SetCollide(true);
    insert->SetBodyFixed(true);

    utils::AddBoxGeometry(insert.get(), mat_b, ChVector3d(thickness * 0.5, width * 0.5, height_insert * 0.5));

    system->AddBody(insert);

    // Static slot (back wall)
    auto slot = chrono_types::make_shared<ChBody>();

    slot->SetIdentifier(-1);
    slot->SetMass(1);
    slot->SetInertiaXX(ChVector3d(1, 1, 1));
    slot->SetPos(ChVector3d(0.5 * thickness, 0, 0.5 * height));
    slot->SetRot(ChQuaternion<>(1, 0, 0, 0));
    slot->SetCollide(true);
    slot->SetBodyFixed(true);

    utils::AddBoxGeometry(slot.get(), mat_b, ChVector3d(thickness / 2, width / 2, height / 2), ChVector3d(0, 0, 0));

    system->AddBody(slot);

    // Lateral walls
    auto wall = chrono_types::make_shared<ChBody>();

    wall->SetIdentifier(-2);
    wall->SetMass(1);
    wall->SetInertiaXX(ChVector3d(1, 1, 1));
    wall->SetPos(ChVector3d(0, 0, 0));
    wall->SetRot(ChQuaternion<>(1, 0, 0, 0));
    wall->SetCollide(true);
    wall->SetBodyFixed(true);

    utils::AddBoxGeometry(wall.get(), mat_b, ChVector3d(3 * height / 2, thickness / 2, height),
                          ChVector3d(0, width / 2 + thickness / 2, height / 2));
    utils::AddBoxGeometry(wall.get(), mat_b, ChVector3d(3 * height / 2, thickness / 2, height),
                          ChVector3d(0, -width / 2 - thickness / 2, height / 2));

    system->AddBody(wall);

// Containing bin
#ifdef USE_SMC
    utils::CreateBoxContainer(system, -3, mat_b,
                              ChVector3d(size_collector / 2, size_collector / 2, height_collector / 2), thickness / 2,
                              ChVector3d(0, 0, -pos_collector));
#else
    utils::CreateBoxContainer(system, -3, mat_b,
                              ChVector3d(size_collector / 2, size_collector / 2, height_collector / 2), thickness / 2,
                              ChVector3d(0, 0, -pos_collector));
#endif

    // Return the angled insert body
    return insert.get();
}

// -----------------------------------------------------------------------------
// Create granular material
// -----------------------------------------------------------------------------
void CreateParticles(ChSystemMulticore* system) {
// Create a material for the granular material
#ifdef USE_SMC
    auto mat_g = chrono_types::make_shared<ChContactMaterialSMC>();
    mat_g->SetYoungModulus(Y_g);
    mat_g->SetFriction(mu_g);
    mat_g->SetRestitution(cr_g);
#else
    auto mat_g = chrono_types::make_shared<ChContactMaterialNSC>();
    mat_g->SetFriction(mu_g);
#endif

    // Create a mixture entirely made out of spheres
    double r = 1.01 * r_g;
    utils::PDSampler<double> sampler(2 * r);
    utils::Generator gen(system);

    std::shared_ptr<utils::MixtureIngredient> m1 = gen.AddMixtureIngredient(utils::MixtureType::SPHERE, 1.0);
    m1->setDefaultMaterial(mat_g);
    m1->setDefaultDensity(rho_g);
    m1->setDefaultSize(r_g);

    gen.setBodyIdentifier(1);

    ChVector3d hdims(0.3 * height, 0.3 * width, 0);
    ChVector3d center(-0.4 * height, 0, 0.8 * height);
    ChVector3d vel(0, 0, 0);

    while (gen.getTotalNumBodies() < desired_num_particles) {
        gen.CreateObjectsBox(sampler, center, hdims, vel);
        center.z() += 2 * r;
    }

    std::cout << "Number of particles: " << gen.getTotalNumBodies() << std::endl;
}

// -----------------------------------------------------------------------------
// Find and return the body with specified identifier.
// -----------------------------------------------------------------------------
ChBody* FindBodyById(ChSystemMulticore* sys, int id) {
    for (auto body : sys->Get_bodylist()) {
        if (body->GetIdentifier() == id)
            return body.get();
    }

    return NULL;
}

// -----------------------------------------------------------------------------
// Find the number of particles whose height is below and above, respectively,
// the specified value.
// -----------------------------------------------------------------------------
int GetNumParticlesBelowHeight(ChSystemMulticore* sys, double value) {
    int count = 0;
    for (auto body : sys->Get_bodylist()) {
        if (body->GetIdentifier() > 0 && body->GetPos().z() < value)
            count++;
    }
    return count;
}

int GetNumParticlesAboveHeight(ChSystemMulticore* sys, double value) {
    int count = 0;
    for (auto body : sys->Get_bodylist()) {
        if (body->GetIdentifier() > 0 && body->GetPos().z() > value)
            count++;
    }
    return count;
}

// -----------------------------------------------------------------------------
// Return true if all bodies in the granular mix have a linear velocity whose
// magnitude is below the specified value.
// -----------------------------------------------------------------------------
bool CheckSettled(ChSystemMulticore* sys, double threshold) {
    double t2 = threshold * threshold;

    for (auto body : sys->Get_bodylist()) {
        if (body->GetIdentifier() > 0) {
            double vel2 = body->GetPos_dt().Length2();
            if (vel2 > t2)
                return false;
        }
    }

    return true;
}

// -----------------------------------------------------------------------------
int main(int argc, char* argv[]) {
// Create system
#ifdef USE_SMC
    cout << "Create SMC system" << endl;
    ChSystemMulticoreSMC* sys = new ChSystemMulticoreSMC();
#else
    cout << "Create NSC system" << endl;
    ChSystemMulticoreNSC* sys = new ChSystemMulticoreNSC();
#endif

    sys->SetCollisionSystemType(ChCollisionSystem::Type::MULTICORE);

    // Set number of threads.
    int max_threads = omp_get_num_procs();
    if (threads > max_threads)
        threads = max_threads;
    sys->SetNumThreads(threads);
    cout << "Using " << threads << " threads" << endl;

    // Set gravitational acceleration
    sys->Set_G_acc(ChVector3d(0, 0, -gravity));

    // Edit system settings
    sys->GetSettings()->solver.max_iteration_bilateral = max_iteration_bilateral;
    sys->GetSettings()->solver.tolerance = tolerance;
    sys->GetSettings()->solver.use_full_inertia_tensor = false;

#ifdef USE_SMC
    sys->GetSettings()->collision.narrowphase_algorithm = ChNarrowphase::Algorithm::PRIMS;
#else
    sys->GetSettings()->solver.solver_mode = SolverMode::SLIDING;
    sys->GetSettings()->solver.max_iteration_normal = max_iteration_normal;
    sys->GetSettings()->solver.max_iteration_sliding = max_iteration_sliding;
    sys->GetSettings()->solver.max_iteration_spinning = max_iteration_spinning;
    sys->GetSettings()->solver.alpha = 0;
    sys->GetSettings()->solver.contact_recovery_speed = contact_recovery_speed;
    sys->ChangeSolverType(SolverType::APGDREF);

    sys->GetSettings()->collision.collision_envelope = 0.05 * r_g;
#endif

    sys->GetSettings()->collision.bins_per_axis = vec3(10, 10, 10);

    // Set simulation duration and create bodies (depending on problem type).
    double time_end;
    int out_fps;
    ChBody* insert;

    switch (problem) {
        case SETTLING:
            time_end = time_settling_max;
            out_fps = out_fps_settling;
            insert = CreateMechanism(sys);
            CreateParticles(sys);
            break;

        case DROPPING:
            time_end = time_dropping_max;
            out_fps = out_fps_dropping;
            utils::ReadCheckpoint(sys, checkpoint_file);
            insert = FindBodyById(sys, 0);
            break;
    }

    // Number of steps
    int num_steps = (int)std::ceil(time_end / time_step);
    int out_steps = (int)std::ceil((1.0 / time_step) / out_fps);

    // Zero velocity level for settling check
    double zero_v = 2 * r_g;

    // Create output directories.
    if (!filesystem::create_directory(filesystem::path(out_dir))) {
        cout << "Error creating directory " << out_dir << endl;
        return 1;
    }
    if (!filesystem::create_directory(filesystem::path(pov_dir))) {
        cout << "Error creating directory " << pov_dir << endl;
        return 1;
    }

    // Perform the simulation
    double time = 0;
    int sim_frame = 0;
    int out_frame = 0;
    int next_out_frame = 0;
    double exec_time = 0;
    int num_contacts = 0;
    std::ofstream sfile(stats_file.c_str());
    std::ofstream ffile(flow_file.c_str());

#ifdef CHRONO_OPENGL
    opengl::ChVisualSystemOpenGL vis;
    vis.AttachSystem(sys);
    vis.SetWindowTitle("Foam demo");
    vis.SetWindowSize(1280, 720);
    vis.SetRenderMode(opengl::WIREFRAME);
    vis.Initialize();
    vis.AddCamera(ChVector3d(0, -12 * width, height), ChVector3d(0, 0, height));
    vis.SetCameraVertical(CameraVerticalDir::Z);
#endif

    while (time < time_end) {
        // Output data
        if (sim_frame == next_out_frame) {
            char filename[100];
            sprintf(filename, "%s/data_%03d.dat", pov_dir.c_str(), out_frame + 1);
            utils::WriteVisualizationAssets(sys, filename);

            cout << "------------ Output frame:   " << out_frame << endl;
            cout << "             Sim frame:      " << sim_frame << endl;
            cout << "             Time:           " << time << endl;
            cout << "             Avg. contacts:  " << num_contacts / out_steps << endl;
            cout << "             Execution time: " << exec_time << endl;

            double opening = insert->GetPos().x() + 0.5 * height + delta;
            int count = GetNumParticlesBelowHeight(sys, 0);

            cout << "             Gap:            " << -opening << endl;
            cout << "             Flow:           " << count << endl;

            sfile << time << "  " << exec_time << "  " << num_contacts / out_steps << "\n";

            switch (problem) {
                case SETTLING:
                    // Create a checkpoint from the current state.
                    utils::WriteCheckpoint(sys, checkpoint_file);
                    cout << "             Checkpoint:     " << sys->Get_bodylist().size() << " bodies" << endl;
                    break;
                case DROPPING:
                    // Save current gap opening and number of dropped particles.
                    ffile << time << "  " << -opening << "  " << count << "\n";
                    break;
            }

            out_frame++;
            next_out_frame += out_steps;
            num_contacts = 0;
        }

        // Check for early termination of settling phase.
        if (problem == SETTLING && time > time_settling_min && CheckSettled(sys, zero_v)) {
            cout << "Granular material settled...  time = " << time << endl;
            break;
        }

        // Check for early termination of dropping phase.
        if (problem == DROPPING && time > time_opening &&
            GetNumParticlesAboveHeight(sys, -pos_collector / 2) == 0) {
            cout << "Granular material exhausted... time = " << time << endl;
            break;
        }

// Advance system state by one step.
// Advance simulation by one step
#ifdef CHRONO_OPENGL
        if (vis.Run()) {
            sys->DoStepDynamics(time_step);
            vis.Render();
        } else
            break;
#else
        sys->DoStepDynamics(time_step);
#endif

        // Open the gate until it reaches the specified gap distance.
        if (problem == DROPPING && time <= time_opening) {
            insert->SetPos(ChVector3d(-0.5 * height - delta - time * speed, 0, 0.5 * height - delta));
            insert->SetPos_dt(ChVector3d(-speed, 0, 0));
        }

        time += time_step;
        sim_frame++;
        exec_time += sys->GetTimerStep();
        num_contacts += sys->GetNcontacts();

        // If requested, output detailed timing information for this step
        if (sim_frame == timing_frame)
            sys->PrintStepStats();
    }

    // Create a checkpoint from the last state
    if (problem == SETTLING) {
        cout << "Write checkpoint data to " << checkpoint_file;
        utils::WriteCheckpoint(sys, checkpoint_file);
        cout << "  done.  Wrote " << sys->Get_bodylist().size() << " bodies." << endl;
    }

    // Final stats
    cout << "==================================" << endl;
    cout << "Number of bodies:  " << sys->GetNumBodies() << endl;
    cout << "Simulation time:   " << exec_time << endl;
    cout << "Number of threads: " << threads << endl;

    return 0;
}
