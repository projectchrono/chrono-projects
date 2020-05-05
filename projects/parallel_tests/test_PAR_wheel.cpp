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
// Author: Radu Serban
// =============================================================================
//
// ChronoParallel demo program for testing contact of a wheel shape.
//
// The global reference frame has Z up.
// All units SI.
// =============================================================================

#include <cstdio>
#include <vector>
#include <cmath>

#include "chrono/core/ChStream.h"
#include "chrono/utils/ChUtilsCreators.h"
#include "chrono/utils/ChUtilsGenerators.h"
#include "chrono/utils/ChUtilsInputOutput.h"

#include "chrono_parallel/physics/ChSystemParallel.h"
#include "chrono_parallel/solver/ChSystemDescriptorParallel.h"

#include "chrono_thirdparty/filesystem/path.h"

#include "../utils.h"

using namespace chrono;
using namespace chrono::collision;

using std::cout;
using std::endl;

// =======================================================================
// Note: Run first SETTLING phase, which generates checkpoint file
// SIMULATION phase assumes a checkpoint file exists.

enum ProblemType { SETTLING, SIMULATION };
ProblemType problem = SETTLING;

// =======================================================================
// Global problem definitions

int threads = 8;

// Simulation parameters
double gravity = 9.81;

double time_settling_min = 0.1;
double time_settling_max = 5;
double time_simulation = 10;

double time_step_penalty = 1e-4;
double time_step_complementarity = 1e-3;

int max_iteration = 20;

// Output
const std::string out_dir = "../WHEEL";
const std::string pov_dir = out_dir + "/POVRAY";
const std::string checkpoint_file = out_dir + "/settled.dat";
double out_fps = 60;

// Contact method
ChContactMethod method = ChContactMethod::SMC;

// Parameters for the granular material
int Id_g = 100;
double r_g = 0.04;
double rho_g = 2500;
double vol_g = (4.0 / 3) * CH_C_PI * r_g * r_g * r_g;
double mass_g = rho_g * vol_g;
ChVector<> inertia_g = 0.4 * mass_g * r_g * r_g * ChVector<>(1, 1, 1);

float Y_g = 2e8f;
float mu_g = 0.5f;
float cr_g = 0.1f;
float cohesion_g = 20.0f;

// Parameters for wheel body
int Id_w = 0;
double mass_w = 100;
ChVector<> inertia_w = ChVector<>(2, 2, 4);

float Y_w = 1e8f;
float mu_w = 1.0f;
float cr_w = 0.1f;
float cohesion_w = 0.0f;

// Parameters for the containing bin
int binId = -200;
double hDimX = 1.0;        // length in x direction
double hDimY = 1.0;        // width in y direction
double hDimZ = 0.5;        // height in z direction
double hThickness = 0.04;  // wall thickness

float Y_c = 2e8f;
float mu_c = 1.0f;
float cr_c = 0.1f;
float cohesion_c = 0.0f;

// Height of layer for generator domain
double layerHeight = 0.5;

// =======================================================================

int CreateObjects(ChSystemParallel* system) {
    // Create materials for the granular material and the container
    std::shared_ptr<chrono::ChMaterialSurface> material_g;
    std::shared_ptr<chrono::ChMaterialSurface> material_c;

    switch (method) {
        case ChContactMethod::SMC: {
            auto mat_g = chrono_types::make_shared<ChMaterialSurfaceSMC>();
            mat_g->SetYoungModulus(Y_g);
            mat_g->SetFriction(mu_g);
            mat_g->SetRestitution(cr_g);
            mat_g->SetAdhesion(cohesion_g);

            material_g = mat_g;

            auto mat_c = chrono_types::make_shared<ChMaterialSurfaceSMC>();
            mat_c->SetYoungModulus(Y_c);
            mat_c->SetFriction(mu_c);
            mat_c->SetRestitution(cr_c);
            mat_c->SetAdhesion(cohesion_c);

            material_c = mat_c;

            break;
        }
        case ChContactMethod::NSC: {
            auto mat_g = chrono_types::make_shared<ChMaterialSurfaceNSC>();
            mat_g->SetFriction(mu_g);
            mat_g->SetRestitution(cr_g);
            mat_g->SetCohesion(cohesion_g);

            material_g = mat_g;

            auto mat_c = chrono_types::make_shared<ChMaterialSurfaceNSC>();
            mat_c->SetFriction(mu_c);
            mat_c->SetRestitution(cr_c);
            mat_c->SetCohesion(cohesion_c);

            material_c = mat_c;

            break;
        }
    }

    // Create a mixture entirely made out of spheres
    utils::Generator gen(system);

    std::shared_ptr<utils::MixtureIngredient> m1 = gen.AddMixtureIngredient(utils::MixtureType::SPHERE, 1.0);
    m1->setDefaultMaterial(material_g);
    m1->setDefaultDensity(rho_g);
    m1->setDefaultSize(r_g);

    gen.setBodyIdentifier(Id_g);

    double r = 1.01 * r_g;
    gen.createObjectsBox(utils::SamplingType::POISSON_DISK, 2 * r, ChVector<>(0, 0, r + layerHeight / 2),
                         ChVector<>(hDimX - r, hDimY - r, layerHeight / 2));
    cout << "total granules: " << gen.getTotalNumBodies() << endl;

    // Create the containing bin
    utils::CreateBoxContainer(system, binId, material_c, ChVector<>(hDimX, hDimY, hDimZ), hThickness);

    return gen.getTotalNumBodies();
}

// =======================================================================
// Create the wheel body at the specified height.

std::shared_ptr<ChBody> CreateWheel(ChSystemParallel* system, double z) {
    // Mesh input file
    std::string obj_mesh_file = GetChronoDataFile("wheel_view.obj");
    std::string mesh_name("wheel");

    // Create a material for the wheel
    std::shared_ptr<chrono::ChMaterialSurface> material_w;

    switch (method) {
        case ChContactMethod::SMC: {
            auto mat_w = chrono_types::make_shared<ChMaterialSurfaceSMC>();
            mat_w->SetYoungModulus(Y_w);
            mat_w->SetFriction(mu_w);
            mat_w->SetRestitution(cr_w);
            mat_w->SetAdhesion(cohesion_w);

            material_w = mat_w;

            break;
        }
        case ChContactMethod::NSC: {
            auto mat_w = chrono_types::make_shared<ChMaterialSurfaceNSC>();
            mat_w->SetFriction(mu_w);
            mat_w->SetRestitution(cr_w);
            mat_w->SetCohesion(cohesion_w);

            material_w = mat_w;

            break;
        }
    }

    // Create the wheel body
    auto wheel = std::shared_ptr<ChBody>(system->NewBody());

    wheel->SetIdentifier(Id_w);
    wheel->SetMass(mass_w);
    wheel->SetInertiaXX(inertia_w);
    wheel->SetPos(ChVector<>(0, 0, z));
    wheel->SetRot(ChQuaternion<>(1, 0, 0, 0));
    wheel->SetCollide(true);
    wheel->SetBodyFixed(false);

    wheel->GetCollisionModel()->ClearModel();
    utils::AddTriangleMeshGeometry(wheel.get(), material_w, obj_mesh_file, mesh_name);
    //utils::AddCylinderGeometry(wheel.get(), material_w, 0.3, 0.1);
    wheel->GetCollisionModel()->BuildModel();

    wheel->SetInertiaXX(inertia_w);

    system->AddBody(wheel);

    // Write POV-Ray mesh model.
    utils::WriteMeshPovray(obj_mesh_file, mesh_name, out_dir);

    return wheel;
}

// ========================================================================
// This utility function returns true if all bodies in the granular mix
// have a linear velocity whose magnitude is below the specified value.

bool CheckSettled(ChSystem* sys, double threshold) {
    double t2 = threshold * threshold;
    for (auto body : sys->Get_bodylist()) {
        if (body->GetIdentifier() >= Id_g) {
            double vel2 = body->GetPos_dt().Length2();
            if (vel2 > t2)
                return false;
        }
    }

    return true;
}

// ========================================================================
// These utility functions find the height of the highest or lowest sphere
// in the granular mix, respectively.  We only look at bodies whith
// identifiers larger than Id_g.

double FindHighest(ChSystem* sys) {
    double highest = 0;
    for (auto body : sys->Get_bodylist()) {
        if (body->GetIdentifier() >= Id_g && body->GetPos().z() > highest)
            highest = body->GetPos().z();
    }
    return highest;
}

double FindLowest(ChSystem* sys) {
    double lowest = DBL_MAX;
    for (auto body : sys->Get_bodylist()) {
        if (body->GetIdentifier() >= Id_g && body->GetPos().z() < lowest)
            lowest = body->GetPos().z();
    }
    return lowest;
}

// ========================================================================
int main(int argc, char* argv[]) {
    // Set path to Chrono data
    SetChronoDataPath(CHRONO_DATA_DIR);

    // Create system and set method-specific solver settings
    ChSystemParallel* system;
    switch (method) {
        case ChContactMethod::SMC: {
            ChSystemParallelSMC* sys = new ChSystemParallelSMC;
            sys->GetSettings()->solver.contact_force_model = ChSystemSMC::Hertz;
            sys->GetSettings()->solver.tangential_displ_mode = ChSystemSMC::TangentialDisplacementModel::OneStep;
            sys->GetSettings()->solver.adhesion_force_model = ChSystemSMC::AdhesionForceModel::Constant;
            sys->GetSettings()->solver.use_material_properties = true;

            system = sys;
            break;
        }
        case ChContactMethod::NSC: {
            ChSystemParallelNSC* sys = new ChSystemParallelNSC;
            sys->GetSettings()->solver.solver_type = SolverType::BB;
            sys->GetSettings()->solver.solver_mode = SolverMode::SLIDING;
            sys->GetSettings()->solver.max_iteration_normal = 0;
            sys->GetSettings()->solver.max_iteration_sliding = 200;
            sys->GetSettings()->solver.max_iteration_spinning = 0;
            sys->GetSettings()->solver.alpha = 0;
            sys->GetSettings()->solver.contact_recovery_speed = -1;
            sys->GetSettings()->collision.collision_envelope = r_g / 10;

            system = sys;
            break;
        }
    }

    // Set method-independent solver settings
    system->Set_G_acc(ChVector<>(0, 0, -gravity));
    system->GetSettings()->solver.use_full_inertia_tensor = false;
    system->GetSettings()->solver.tolerance = 1e-3;
    system->GetSettings()->collision.bins_per_axis = vec3(50, 50, 50);
    system->GetSettings()->collision.narrowphase_algorithm = NarrowPhaseType::NARROWPHASE_HYBRID_MPR;

    // Set number of threads.
    int max_threads = omp_get_num_procs();
    if (threads > max_threads)
        threads = max_threads;
    omp_set_num_threads(threads);

    // Create output directories.
    if (!filesystem::create_directory(filesystem::path(out_dir))) {
        cout << "Error creating directory " << out_dir << endl;
        return 1;
    }
    if (!filesystem::create_directory(filesystem::path(pov_dir))) {
        cout << "Error creating directory " << pov_dir << endl;
        return 1;
    }

    // Depending on problem type:
    // - Select end simulation times
    // - Create granular material and container
    // - Create wheel
    double time_end;
    std::shared_ptr<ChBody> wheel;

    switch (problem) {
        case SETTLING:
            time_end = time_settling_max;

            // Create containing bin and the granular material at randomized initial positions
            CreateObjects(system);

            break;

        case SIMULATION:
            time_end = time_simulation;

            // Create the granular material bodies and the container from the checkpoint file.
            cout << "Read checkpoint data from " << checkpoint_file;
            utils::ReadCheckpoint(system, checkpoint_file);
            cout << "  done.  Read " << system->Get_bodylist().size() << " bodies." << endl;

            // Create the wheel.
            double z = FindHighest(system);
            wheel = CreateWheel(system, z + r_g + 0.4);

            break;
    }

    // Set integration step size
    double time_step;
    switch (method) {
        case ChContactMethod::SMC:
            time_step = time_step_penalty;
            break;
        case ChContactMethod::NSC:
            time_step = time_step_complementarity;
            break;
    }
    system->SetStep(time_step);

    // Number of steps
    int num_steps = (int)std::ceil(time_end / time_step);
    int out_steps = (int)std::ceil((1 / time_step) / out_fps);

    // Zero velocity level for settling check (fraction of a grain radius per second)
    double zero_v = 0.9 * r_g;

    // Perform the simulation
    double time = 0;
    int sim_frame = 0;
    int out_frame = 0;
    int next_out_frame = 0;
    double exec_time = 0;

    while (time < time_end) {
        if (sim_frame == next_out_frame) {
            char filename[100];

            sprintf(filename, "%s/data_%03d.dat", pov_dir.c_str(), out_frame);
            utils::WriteShapesPovray(system, filename);

            cout << " --------------------------------- Output frame:   " << out_frame << endl;
            cout << "                                   Sim frame:      " << sim_frame << endl;
            cout << "                                   Time:           " << time << endl;
            cout << "                                   Execution time: " << exec_time << endl;

            // Check if already settled.
            if (problem == SETTLING && time > time_settling_min && CheckSettled(system, zero_v)) {
                cout << "Granular material settled...  time = " << time << endl;
                break;
            }

            // Save checkpoint during settling phase.
            if (problem == SETTLING) {
                utils::WriteCheckpoint(system, checkpoint_file);
            }

            out_frame++;
            next_out_frame += out_steps;
        }

        //TimingOutput(msystem);
        system->DoStepDynamics(time_step);

        time += time_step;
        sim_frame++;
        exec_time += system->GetTimerStep();
    }

    // Create a checkpoint from the last state
    if (problem == SETTLING)
        utils::WriteCheckpoint(system, checkpoint_file);

    // Final stats
    cout << "==================================" << endl;
    cout << "Number of bodies: " << system->Get_bodylist().size() << endl;
    cout << "Simulation time: " << exec_time << endl;
    cout << "Number of threads: " << threads << endl;

    return 0;
}
