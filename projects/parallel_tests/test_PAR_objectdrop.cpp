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
// ChronoParallel demo program for contact of various shapes.
//
// The global reference frame has Z up.
// All units SI.
// =============================================================================

#include <cstdio>
#include <vector>
#include <cmath>

#include "chrono/ChConfig.h"
#include "chrono/core/ChFileutils.h"
#include "chrono/core/ChStream.h"
#include "chrono/utils/ChUtilsGeometry.h"
#include "chrono/utils/ChUtilsCreators.h"
#include "chrono/utils/ChUtilsInputOutput.h"

#include "chrono_parallel/physics/ChSystemParallel.h"
#include "chrono_parallel/solver/ChSystemDescriptorParallel.h"
#include "chrono_parallel/collision/ChNarrowphaseRUtils.h"

// Note: CHRONO_OPENGL is defined in ChConfig.h
#ifdef CHRONO_OPENGL
#include "chrono_opengl/ChOpenGLWindow.h"
#endif

using namespace chrono;
using namespace chrono::collision;

using std::cout;
using std::endl;

// -----------------------------------------------------------------------------
// Problem setup
// -----------------------------------------------------------------------------

// Comment the following line to use NSC contact
#define USE_SMC

// Parameters for the falling object
collision::ShapeType shape_o = collision::CONE;

ChVector<> initPos(1.0, -1.0, 2.0);
// ChQuaternion<> initRot(1.0, 0.0, 0.0, 0.0);
ChQuaternion<> initRot = Q_from_AngAxis(CH_C_PI / 3, ChVector<>(1, 0, 0));

ChVector<> initLinVel(0.0, 0.0, 0.0);
ChVector<> initAngVel(0.0, 0.0, 0.0);

// Ground contact shapes (SPHERE, CAPSULE, BOX)
collision::ShapeType shape_g = collision::SPHERE;

// -----------------------------------------------------------------------------
// Simulation parameters
// -----------------------------------------------------------------------------

// Desired number of OpenMP threads (will be clamped to maximum available)
int threads = 100;

// Perform dynamic tuning of number of threads?
bool thread_tuning = true;

// Simulation duration.
double time_end = 5;

// Solver parameters
#ifdef USE_SMC
double time_step = 1e-3;
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
const std::string out_dir = "../OBJECTDROP_SMC";
#else
const std::string out_dir = "../OBJECTDROP_NSC";
#endif
const std::string pov_dir = out_dir + "/POVRAY";

int out_fps = 60;

// Continuous loop (only if OpenGL available)
// If true, no output files are generated
bool loop = true;

// =============================================================================
// Create ground body
// =============================================================================
void CreateGround(ChSystemParallel* system) {
// ---------------------------------------
// Create a material and the "ground" body
// ---------------------------------------

#ifdef USE_SMC
    auto mat_g = std::make_shared<ChMaterialSurfaceSMC>();
    mat_g->SetYoungModulus(1e7f);
    mat_g->SetFriction(0.4f);
    mat_g->SetRestitution(0.4f);

    auto ground = std::make_shared<ChBody>(std::make_shared<ChCollisionModelParallel>(), ChMaterialSurface::SMC);
    ground->SetMaterialSurface(mat_g);
#else
    auto mat_g = std::make_shared<ChMaterialSurfaceNSC>();
    mat_g->SetFriction(0.4f);

    auto ground = std::make_shared<ChBody>(std::make_shared<ChCollisionModelParallel>());
    ground->SetMaterialSurface(mat_g);
#endif

    ground->SetIdentifier(-1);
    ground->SetMass(1);
    ground->SetPos(ChVector<>(0, 0, 0));
    ground->SetRot(ChQuaternion<>(1, 0, 0, 0));
    ground->SetBodyFixed(true);
    ground->SetCollide(true);

    // ---------------------------------------------------------
    // Set fixed contact shapes (depending on specified option).
    // ---------------------------------------------------------

    switch (shape_g) {
        case SPHERE:
            // A grid of 5x5 spheres
            {
                double spacing = 1.6;
                double bigR = 2;

                ground->GetCollisionModel()->ClearModel();
                for (int ix = -2; ix < 3; ix++) {
                    for (int iy = -2; iy < 3; iy++) {
                        ChVector<> pos(ix * spacing, iy * spacing, -bigR);
                        utils::AddSphereGeometry(ground.get(), bigR, pos);
                    }
                }
                ground->GetCollisionModel()->BuildModel();
            }
            break;

        case CAPSULE:
            // A set of 7 parallel capsules, rotated by 30 degrees around Z
            {
                double spacing = 1.5;
                double bigR = 1;
                double bigH = 6;

                ChQuaternion<> rot(1, 0, 0, 0);
                rot.Q_from_AngAxis(CH_C_PI / 6, ChVector<>(0, 0, 1));

                ground->GetCollisionModel()->ClearModel();
                for (int ix = -3; ix < 6; ix++) {
                    ChVector<> pos(ix * spacing, 0, -bigR);
                    utils::AddCapsuleGeometry(ground.get(), bigR, bigH, pos, rot);
                }
                ground->GetCollisionModel()->BuildModel();
            }
            break;

        case BOX:
            // A single box
            {
                double bigHx = 6;
                double bigHy = 6;
                double bigHz = 1;

                ground->GetCollisionModel()->ClearModel();
                utils::AddBoxGeometry(ground.get(), ChVector<>(bigHx, bigHy, bigHz), ChVector<>(0, 0, -bigHz));
                ground->GetCollisionModel()->BuildModel();
            }
            break;
    }

    // ------------------------------------
    // Add the "ground" body to the system.
    // ------------------------------------
    system->AddBody(ground);
}

// =============================================================================
// Create falling object
// =============================================================================
void CreateObject(ChSystemParallel* system) {
    double rho_o = 2000.0;

// -----------------------------------------
// Create a material and the falling object.
// -----------------------------------------

#ifdef USE_SMC
    auto mat_o = std::make_shared<ChMaterialSurfaceSMC>();
    mat_o->SetYoungModulus(1e7f);
    mat_o->SetFriction(0.4f);
    mat_o->SetRestitution(0.4f);

    auto obj = std::make_shared<ChBody>(std::make_shared<ChCollisionModelParallel>(), ChMaterialSurface::SMC);
    obj->SetMaterialSurface(mat_o);
#else
    auto mat_o = std::make_shared<ChMaterialSurfaceNSC>();
    mat_o->SetFriction(0.4f);

    auto obj = std::make_shared<ChBody>(std::make_shared<ChCollisionModelParallel>());
    obj->SetMaterialSurface(mat_o);
#endif

    obj->SetIdentifier(1);
    obj->SetCollide(true);
    obj->SetBodyFixed(false);

    // ----------------------------------------------------
    // Depending on the shape of the falling object,
    //    - Calculate bounding radius, volume, and gyration
    //    - Set contact and visualization shape
    // ----------------------------------------------------

    double rb;
    double vol;
    ChMatrix33<> J;

    obj->GetCollisionModel()->ClearModel();

    switch (shape_o) {
        case collision::SPHERE: {
            double radius = 0.3;
            rb = utils::CalcSphereBradius(radius);
            vol = utils::CalcSphereVolume(radius);
            J = utils::CalcSphereGyration(radius);
            utils::AddSphereGeometry(obj.get(), radius);
        } break;
        case collision::BOX: {
            ChVector<> hdims(0.1, 0.2, 0.1);
            rb = utils::CalcBoxBradius(hdims);
            vol = utils::CalcBoxVolume(hdims);
            J = utils::CalcBoxGyration(hdims);
            utils::AddBoxGeometry(obj.get(), hdims);
        } break;
        case collision::CAPSULE: {
            double radius = 0.1;
            double hlen = 0.2;
            rb = utils::CalcCapsuleBradius(radius, hlen);
            vol = utils::CalcCapsuleVolume(radius, hlen);
            J = utils::CalcCapsuleGyration(radius, hlen);
            utils::AddCapsuleGeometry(obj.get(), radius, hlen);
        } break;
        case collision::CYLINDER: {
            double radius = 0.1;
            double hlen = 0.2;
            rb = utils::CalcCylinderBradius(radius, hlen);
            vol = utils::CalcCylinderVolume(radius, hlen);
            J = utils::CalcCylinderGyration(radius, hlen);
            utils::AddCylinderGeometry(obj.get(), radius, hlen);
        } break;
        case collision::ROUNDEDCYL: {
            double radius = 0.1;
            double hlen = 0.2;
            double srad = 0.05;
            rb = utils::CalcRoundedCylinderBradius(radius, hlen, srad);
            vol = utils::CalcRoundedCylinderVolume(radius, hlen, srad);
            J = utils::CalcRoundedCylinderGyration(radius, hlen, srad);
            utils::AddRoundedCylinderGeometry(obj.get(), radius, hlen, srad);
        } break;
        case collision::CONE: {
            double radius = 0.2;
            double height = 0.4;
            rb = utils::CalcConeBradius(radius, height);
            vol = utils::CalcConeVolume(radius, height);
            J = utils::CalcConeGyration(radius, height);
            utils::AddConeGeometry(obj.get(), radius, height);
        } break;
    }

    obj->GetCollisionModel()->BuildModel();

    // ---------------------
    // Set mass and inertia.
    // ---------------------

    double mass = rho_o * vol;
    obj->SetMass(mass);
    obj->SetInertia(J * mass);

    // ------------------
    // Set initial state.
    // ------------------
    assert(initPos.z() > rb);

    obj->SetPos(initPos);
    obj->SetRot(initRot);
    obj->SetPos_dt(initLinVel);
    obj->SetWvel_loc(initAngVel);

    // ---------------------
    // Add object to system.
    // ---------------------
    system->AddBody(obj);
}

// =============================================================================
// =============================================================================
int main(int argc, char* argv[]) {
    // --------------------------
    // Create output directories.
    // --------------------------

    if (ChFileutils::MakeDirectory(out_dir.c_str()) < 0) {
        cout << "Error creating directory " << out_dir << endl;
        return 1;
    }
    if (ChFileutils::MakeDirectory(pov_dir.c_str()) < 0) {
        cout << "Error creating directory " << pov_dir << endl;
        return 1;
    }

// --------------
// Create system.
// --------------
    char title[100];
#ifdef USE_SMC
    sprintf(title, "Object Drop >> SMC");
    cout << "Create SMC system" << endl;
    ChSystemParallelSMC* msystem = new ChSystemParallelSMC();
#else
    sprintf(title, "Object Drop >> NSC");
    cout << "Create NSC system" << endl;
    ChSystemParallelNSC* msystem = new ChSystemParallelNSC();
#endif

    msystem->Set_G_acc(ChVector<>(0, 0, -9.81));

    // ----------------------
    // Set number of threads.
    // ----------------------

    int max_threads = omp_get_num_procs();
    if (threads > max_threads)
        threads = max_threads;
    msystem->SetParallelThreadNumber(threads);
    omp_set_num_threads(threads);
    cout << "Using " << threads << " threads" << endl;

    // ---------------------
    // Edit system settings.
    // ---------------------

    msystem->GetSettings()->solver.tolerance = 1e-3;

#ifdef USE_SMC
    msystem->GetSettings()->collision.narrowphase_algorithm = NarrowPhaseType::NARROWPHASE_HYBRID_MPR;
#else
    msystem->GetSettings()->solver.solver_mode = SolverMode::SLIDING;
    msystem->GetSettings()->solver.max_iteration_normal = max_iteration_normal;
    msystem->GetSettings()->solver.max_iteration_sliding = max_iteration_sliding;
    msystem->GetSettings()->solver.max_iteration_spinning = max_iteration_spinning;
    msystem->GetSettings()->solver.alpha = 0;
    msystem->GetSettings()->solver.contact_recovery_speed = contact_recovery_speed;
    msystem->ChangeSolverType(SolverType::APGDREF);

    msystem->GetSettings()->solver.contact_recovery_speed = 1;
#endif

    msystem->GetSettings()->collision.bins_per_axis = vec3(10, 10, 10);

    // --------------
    // Create bodies.
    // --------------
    CreateGround(msystem);
    CreateObject(msystem);

// -----------------------
// Perform the simulation.
// -----------------------

#ifdef CHRONO_OPENGL
    // Initialize OpenGL
    opengl::ChOpenGLWindow& gl_window = opengl::ChOpenGLWindow::getInstance();
    gl_window.Initialize(1280, 720, title, msystem);
    gl_window.SetCamera(ChVector<>(0, -10, 0), ChVector<>(0, 0, 0), ChVector<>(0, 0, 1));

    // Let the OpenGL manager run the simulation until interrupted.
    if (loop) {
        gl_window.StartDrawLoop(time_step);
        return 0;
    }
#endif

    // Run simulation for specified time.
    int out_steps = (int)std::ceil((1.0 / time_step) / out_fps);

    double time = 0;
    int sim_frame = 0;
    int out_frame = 0;
    int next_out_frame = 0;
    double exec_time = 0;
    int num_contacts = 0;

    while (time < time_end) {
        if (sim_frame == next_out_frame) {
            char filename[100];
            sprintf(filename, "%s/data_%03d.dat", pov_dir.c_str(), out_frame + 1);
            utils::WriteShapesPovray(msystem, filename);

            cout << "------------ Output frame:   " << out_frame << endl;
            cout << "             Sim frame:      " << sim_frame << endl;
            cout << "             Time:           " << time << endl;
            cout << "             Avg. contacts:  " << num_contacts / out_steps << endl;
            cout << "             Execution time: " << exec_time << endl;

            out_frame++;
            next_out_frame += out_steps;
            num_contacts = 0;
        }

// Advance dynamics.
#ifdef CHRONO_OPENGL
        if (gl_window.Active()) {
            gl_window.DoStepDynamics(time_step);
            gl_window.Render();
        } else
            break;
#else
        msystem->DoStepDynamics(time_step);
#endif

        // Update counters.
        time += time_step;
        sim_frame++;
        exec_time += msystem->GetTimerStep();
        num_contacts += msystem->GetNcontacts();
    }

    // Final stats
    cout << "==================================" << endl;
    cout << "Simulation time:   " << exec_time << endl;
    cout << "Number of threads: " << threads << endl;

    return 0;
}
