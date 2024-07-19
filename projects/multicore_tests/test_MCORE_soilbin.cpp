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
// Author: Dan Melanz, Radu Serban
// =============================================================================
//
// Chrono::Multicore demo program for simulatin of wheel in soilbin.
//
// The global reference frame has Z up.
// All units SI.
// =============================================================================

#include <cstdio>
#include <vector>
#include <cmath>

#include "chrono/utils/ChUtilsGeometry.h"
#include "chrono/utils/ChUtilsCreators.h"
#include "chrono/utils/ChUtilsGenerators.h"
#include "chrono/utils/ChUtilsInputOutput.h"

#include "chrono_multicore/physics/ChSystemMulticore.h"
#include "chrono_multicore/solver/ChSystemDescriptorMulticore.h"

#include "chrono_thirdparty/filesystem/path.h"

using namespace chrono;

using std::cout;
using std::endl;
using std::flush;

// -----------------------------------------------------------------------------
// Problem setup
// -----------------------------------------------------------------------------

// Comment the following line to use NSC contact
#define USE_SMC

// Simulation phase
enum ProblemType {
    SETTLING,
    DROPPING,
};
ProblemType problem = DROPPING;

// -----------------------------------------------------------------------------
// Simulation parameters
// -----------------------------------------------------------------------------

// Desired number of OpenMP threads (will be clamped to maximum available)
int threads = 100;

// Perform dynamic tuning of number of threads?
bool thread_tuning = true;

// Simulation duration.
double time_settling = 5;
double time_dropping = 2;

// Solver parameters
#ifdef USE_SMC
double time_step = 1e-4;
int max_iteration = 20;
#else
double time_step = 1e-4;
int max_iteration_normal = 30;
int max_iteration_sliding = 20;
int max_iteration_spinning = 0;
float contact_recovery_speed = 0.1;
#endif

// Output
#ifdef USE_SMC
const std::string out_dir = "../SOILBIN_SMC";
#else
const std::string out_dir = "../SOILBIN_NSC";
#endif
const std::string pov_dir = out_dir + "/POVRAY";
const std::string checkpoint_file = out_dir + "/settled.dat";
const std::string stats_file = out_dir + "/stats.dat";

int out_fps_settling = 30;
int out_fps_dropping = 60;

// -----------------------------------------------------------------------------
// Parameters for the granular material (identical spheres)
// -----------------------------------------------------------------------------
int tag_particles = 0;
double r_g = 0.1;
double rho_g = 2000;
unsigned int desired_num_particles = 1000;

// -----------------------------------------------------------------------------
// Parameters for the falling object
// -----------------------------------------------------------------------------
// Shape of dropped object
ChCollisionShape::Type shape_o = ChCollisionShape::Type::ROUNDEDCYL;

ChQuaternion<> initRot(1.0, 0.0, 0.0, 0.0);
ChVector3d initLinVel(0.0, 0.0, 0.0);
ChVector3d initAngVel(0.0, 0.0, 0.0);

// -----------------------------------------------------------------------------
// Half-dimensions of the container bin
// -----------------------------------------------------------------------------
double hDimX = 5;
double hDimY = 2;
double hDimZ = 2;

// =============================================================================
// Create container bin.
// =============================================================================
void CreateContainer(ChSystemMulticore* system) {
    double hThickness = 0.1;

#ifdef USE_SMC
    auto mat_c = chrono_types::make_shared<ChContactMaterialSMC>();
    mat_c->SetYoungModulus(2e6f);
    mat_c->SetFriction(0.4f);
    mat_c->SetRestitution(0.1f);

    utils::CreateBoxContainer(system, mat_c, ChVector3d(hDimX, hDimY, hDimZ), hThickness);
#else
    auto mat_c = chrono_types::make_shared<ChContactMaterialNSC>();
    mat_c->SetFriction(0.4f);

    utils::CreateBoxContainer(system, id_c, mat_c, ChVector3d(hDimX, hDimY, hDimZ), hThickness);

#endif
}

// =============================================================================
// Create granular material.
// =============================================================================
void CreateParticles(ChSystemMulticore* system) {
// Create a material for the ball mixture.
#ifdef USE_SMC
    auto mat_g = chrono_types::make_shared<ChContactMaterialSMC>();
    mat_g->SetYoungModulus(1e8f);
    mat_g->SetFriction(0.4f);
    mat_g->SetRestitution(0.1f);
#else
    auto mat_g = chrono_types::make_shared<ChContactMaterialNSC>();
    mat_g->SetFriction(0.4f);
#endif

    // Create a mixture entirely made out of spheres.
    double r = 1.01 * r_g;
    utils::ChPDSampler<double> sampler(2 * r);
    utils::ChGenerator gen(system);

    std::shared_ptr<utils::ChMixtureIngredient> m1 = gen.AddMixtureIngredient(utils::MixtureType::SPHERE, 1.0);
    m1->SetDefaultMaterial(mat_g);
    m1->SetDefaultDensity(rho_g);
    m1->SetDefaultSize(r_g);

    // Create particles, one layer at a time, until the desired number is reached.
    gen.SetStartTag(tag_particles);

    ChVector3d hdims(hDimX - r, hDimY - r, 0);
    ChVector3d center(0, 0, 2 * r);

    while (gen.GetTotalNumBodies() < desired_num_particles) {
        gen.CreateObjectsBox(sampler, center, hdims);
        center.z() += 2 * r;
    }

    cout << "Number of particles: " << gen.GetTotalNumBodies() << endl;
}

// =============================================================================
// Create falling object.
// =============================================================================
void CreateObject(ChSystemMulticore* system, double z) {
    double rho_o = 2000.0;

// -----------------------------------------
// Create a material for the falling object.
// -----------------------------------------

#ifdef USE_SMC
    auto mat_o = chrono_types::make_shared<ChContactMaterialSMC>();
    mat_o->SetYoungModulus(1e8f);
    mat_o->SetFriction(0.4f);
    mat_o->SetRestitution(0.1f);
#else
    auto mat_o = chrono_types::make_shared<ChContactMaterialNSC>();
    mat_o->SetFriction(0.4f);
#endif

    // --------------------------
    // Create the falling object.
    // --------------------------

    auto obj = chrono_types::make_shared<ChBody>();

    obj->EnableCollision(true);
    obj->SetFixed(false);

    // ----------------------------------------------------
    // Depending on the shape of the falling object,
    //    - Calculate bounding radius, volume, and gyration
    //    - Calculate bounding radius, volume, and gyration
    //    - Set contact and visualization shape
    // ----------------------------------------------------

    double rb;
    double vol;
    ChMatrix33<> J;

    switch (shape_o) {
        case ChCollisionShape::Type::SPHERE: {
            double radius = 0.5;
            rb = ChSphere::GetBoundingSphereRadius(radius);
            vol = ChSphere::GetVolume(radius);
            J = ChSphere::GetGyration(radius);
            utils::AddSphereGeometry(obj.get(), mat_o, radius);
        } break;
        case ChCollisionShape::Type::BOX: {
            ChVector3d dims(0.5, 0.75, 1.0);
            rb = ChBox::GetBoundingSphereRadius(dims);
            vol = ChBox::GetVolume(dims);
            J = ChBox::GetGyration(dims);
            utils::AddBoxGeometry(obj.get(), mat_o, dims);
        } break;
        case ChCollisionShape::Type::CAPSULE: {
            double radius = 0.25;
            double len = 0.5;
            rb = ChCapsule::GetBoundingSphereRadius(radius, len);
            vol = ChCapsule::GetVolume(radius, len);
            J = ChCapsule::GetGyration(radius, len);
            utils::AddCapsuleGeometry(obj.get(), mat_o, radius, len);
        } break;
        case ChCollisionShape::Type::CYLINDER: {
            double radius = 0.25;
            double len = 0.5;
            rb = ChCylinder::GetBoundingSphereRadius(radius, len);
            vol = ChCylinder::GetVolume(radius, len);
            J = ChCylinder::GetGyration(radius, len);
            utils::AddCylinderGeometry(obj.get(), mat_o, radius, len);
        } break;
        case ChCollisionShape::Type::ROUNDEDCYL: {
            double radius = 0.25;
            double len = 0.1;
            double srad = 0.1;
            rb = ChRoundedCylinder::GetBoundingSphereRadius(radius, len, srad);
            vol = ChRoundedCylinder::GetVolume(radius, len, srad);
            J = ChRoundedCylinder::GetGyration(radius, len, srad);
            utils::AddRoundedCylinderGeometry(obj.get(), mat_o, radius, len, srad);
        } break;
    }

    // ---------------------
    // Set mass and inertia.
    // ---------------------

    double mass = rho_o * vol;
    obj->SetMass(mass);
    obj->SetInertia(J * mass);

    // ------------------
    // Set initial state.
    // ------------------

    obj->SetPos(ChVector3d(0, 0, z + rb));
    obj->SetRot(initRot);
    obj->SetLinVel(initLinVel);
    obj->SetAngVelLocal(initAngVel);

    // ---------------------
    // Add object to system.
    // ---------------------
    system->AddBody(obj);
}

// =============================================================================
// Find the height of the highest and lowest, respectively, sphere in the
// granular mix, respectively.  We only look at bodies whith stricty positive
// identifiers (to exclude the containing bin).
// =============================================================================
double FindHighest(ChSystem* sys) {
    double highest = 0;
    for (auto body : sys->GetBodies()) {
        if (body->GetTag() >= tag_particles && body->GetPos().z() > highest)
            highest = body->GetPos().z();
    }
    return highest;
}

double FindLowest(ChSystem* sys) {
    double lowest = 1000;
    for (auto body : sys->GetBodies()) {
        if (body->GetTag() >= tag_particles && body->GetPos().z() < lowest)
            lowest = body->GetPos().z();
    }
    return lowest;
}

// =============================================================================
// =============================================================================
int main(int argc, char* argv[]) {
    // --------------------------
    // Create output directories.
    // --------------------------

    if (!filesystem::create_directory(filesystem::path(out_dir))) {
        cout << "Error creating directory " << out_dir << endl;
        return 1;
    }
    if (!filesystem::create_directory(filesystem::path(pov_dir))) {
        cout << "Error creating directory " << pov_dir << endl;
        return 1;
    }

// --------------
// Create system.
// --------------

#ifdef USE_SMC
    cout << "Create SMC system" << endl;
    ChSystemMulticoreSMC* msystem = new ChSystemMulticoreSMC();
#else
    cout << "Create NSC system" << endl;
    ChSystemMulticoreNSC* msystem = new ChSystemMulticoreNSC();
#endif

    msystem->SetCollisionSystemType(ChCollisionSystem::Type::MULTICORE);
    msystem->SetGravitationalAcceleration(ChVector3d(0, 0, -9.81));

    // ----------------------
    // Set number of threads.
    // ----------------------

    int max_threads = omp_get_num_procs();
    if (threads > max_threads)
        threads = max_threads;
    omp_set_num_threads(threads);
    cout << "Using " << threads << " threads" << endl;

    // ---------------------
    // Edit system settings.
    // ---------------------

    msystem->GetSettings()->solver.tolerance = 1e-3;

#ifdef USE_SMC
    msystem->GetSettings()->collision.narrowphase_algorithm = ChNarrowphase::Algorithm::PRIMS;
#else
    msystem->GetSettings()->solver.solver_mode = SolverMode::SLIDING;
    msystem->GetSettings()->solver.max_iteration_normal = max_iteration_normal;
    msystem->GetSettings()->solver.max_iteration_sliding = max_iteration_sliding;
    msystem->GetSettings()->solver.max_iteration_spinning = max_iteration_spinning;
    msystem->GetSettings()->solver.alpha = 0;
    msystem->GetSettings()->solver.contact_recovery_speed = contact_recovery_speed;
    msystem->ChangeSolverType(SolverType::APGDREF);

    msystem->GetSettings()->collision.collision_envelope = 0.05 * r_g;
#endif

    msystem->GetSettings()->collision.bins_per_axis = vec3(10, 10, 10);

    // ----------------------------------------
    // Depending on problem type:
    // - Select end simulation time
    // - Create granular material and container
    // - Create falling object
    // ----------------------------------------

    double time_end;
    int out_fps;

    if (problem == SETTLING) {
        time_end = time_settling;
        out_fps = out_fps_settling;

        cout << "Create granular material" << endl;
        CreateContainer(msystem);
        CreateParticles(msystem);
    } else {
        time_end = time_dropping;
        out_fps = out_fps_dropping;

        // Create the granular material and the container from the checkpoint file.
        cout << "Read checkpoint data from " << checkpoint_file;
        utils::ReadCheckpoint(msystem, checkpoint_file);
        cout << "  done.  Read " << msystem->GetBodies().size() << " bodies." << endl;

        // Create the falling object just above the granular material.
        double z = FindHighest(msystem);
        cout << "Create falling object above height" << z + r_g << endl;
        CreateObject(msystem, z + r_g);
    }

    // Number of steps.
    int num_steps = (int)std::ceil(time_end / time_step);
    int out_steps = (int)std::ceil((1.0 / time_step) / out_fps);

    // -----------------------
    // Perform the simulation.
    // -----------------------
    double time = 0;
    int sim_frame = 0;
    int out_frame = 0;
    int next_out_frame = 0;
    double exec_time = 0;
    int num_contacts = 0;
    std::ofstream sfile(stats_file.c_str());

    while (time < time_end) {
        if (sim_frame == next_out_frame) {
            char filename[100];
            sprintf(filename, "%s/data_%03d.dat", pov_dir.c_str(), out_frame + 1);
            utils::WriteVisualizationAssets(msystem, filename);

            cout << "------------ Output frame:   " << out_frame << endl;
            cout << "             Sim frame:      " << sim_frame << endl;
            cout << "             Time:           " << time << endl;
            cout << "             Lowest point:   " << FindLowest(msystem) << endl;
            cout << "             Avg. contacts:  " << num_contacts / out_steps << endl;
            cout << "             Execution time: " << exec_time << endl;

            sfile << time << "  " << exec_time << "  " << num_contacts / out_steps << "\n";

            // Create a checkpoint from the current state.
            if (problem == SETTLING) {
                cout << "             Write checkpoint data " << flush;
                utils::WriteCheckpoint(msystem, checkpoint_file);
                cout << msystem->GetBodies().size() << " bodies" << endl;
            }

            out_frame++;
            next_out_frame += out_steps;
            num_contacts = 0;
        }

        // Advance dynamics.
        msystem->DoStepDynamics(time_step);

        time += time_step;
        sim_frame++;
        exec_time += msystem->GetTimerStep();
        num_contacts += msystem->GetNumContacts();
    }

    // Create a checkpoint from the last state
    if (problem == SETTLING) {
        cout << "Write checkpoint data to " << checkpoint_file;
        utils::WriteCheckpoint(msystem, checkpoint_file);
        cout << "  done.  Wrote " << msystem->GetBodies().size() << " bodies." << endl;
    }

    // Final stats
    cout << "==================================" << endl;
    cout << "Number of bodies:  " << msystem->GetBodies().size() << endl;
    cout << "Lowest position:   " << FindLowest(msystem) << endl;
    cout << "Simulation time:   " << exec_time << endl;
    cout << "Number of threads: " << threads << endl;

    return 0;
}
